#include "ram/algorithm.hpp"

#include <deque>

#include "biosoup/nucleic_acid.hpp"

namespace ram {

std::vector<Kmer> Minimize(
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    MinimizeConfig config) {
  if (sequence->inflated_len < config.kmer_length) {
    return std::vector<Kmer>{};
  }

  std::uint64_t mask = (1ULL << (config.kmer_length * 2)) - 1;

  auto hash = [&](std::uint64_t key) -> std::uint64_t {
    key = ((~key) + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;
    return key;
  };

  std::deque<Kmer> window;
  auto window_add = [&](std::uint64_t value, std::uint64_t location) -> void {
    while (!window.empty() && window.back().value > value) {
      window.pop_back();
    }
    window.emplace_back(value, location);
  };
  auto window_update = [&](std::uint32_t position) -> void {
    while (!window.empty() && (window.front().position()) < position) {
      window.pop_front();
    }
  };

  std::uint64_t shift = (config.kmer_length - 1) * 2;
  std::uint64_t minimizer = 0;
  std::uint64_t reverse_minimizer = 0;
  std::uint64_t id = static_cast<std::uint64_t>(sequence->id) << 32;
  std::uint64_t is_stored = 1ULL << 63;

  std::vector<Kmer> dst;

  for (std::uint32_t i = 0; i < sequence->inflated_len; ++i) {
    std::uint64_t c = sequence->Code(i);
    minimizer = ((minimizer << 2) | c) & mask;
    reverse_minimizer = (reverse_minimizer >> 2) | ((c ^ 3) << shift);
    if (i >= config.kmer_length - 1U) {
      if (minimizer < reverse_minimizer) {
        window_add(hash(minimizer), (i - (config.kmer_length - 1U)) << 1 | 0);
      } else if (minimizer > reverse_minimizer) {
        window_add(hash(reverse_minimizer),
                   (i - (config.kmer_length - 1U)) << 1 | 1);
      }
    }
    if (i >= (config.kmer_length - 1U) + (config.window_length - 1U)) {
      for (auto it = window.begin(); it != window.end(); ++it) {
        if (it->value != window.front().value) {
          break;
        }
        if (it->origin & is_stored) {
          continue;
        }
        dst.emplace_back(it->value, id | it->origin);
        it->origin |= is_stored;
      }
      window_update(i - (config.kmer_length - 1U) -
                    (config.window_length - 1U) + 1);
    }
  }

  if (config.minhash) {
    RadixSort(std::span<Kmer>(dst), config.kmer_length * 2,
              KmerValueProjection);
    dst.resize(sequence->inflated_len / config.kmer_length);
    RadixSort(std::span<Kmer>(dst), 64, KmerOriginProjection);
  }

  return dst;
}

std::vector<biosoup::Overlap> Map(
    const std::unique_ptr<biosoup::NucleicAcid>& lhs,
    const std::unique_ptr<biosoup::NucleicAcid>& rhs,
    MinimizeConfig minimize_config, ChainConfig chain_config) {
  auto lhs_sketch = ::ram::Minimize(lhs, minimize_config);
  if (lhs_sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  auto rhs_sketch = ::ram::Minimize(
      rhs, MinimizeConfig{.kmer_length = minimize_config.kmer_length,
                          .window_length = minimize_config.window_length});
  if (rhs_sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  RadixSort(std::span<Kmer>(lhs_sketch), minimize_config.kmer_length * 2,
            KmerValueProjection);
  RadixSort(std::span<Kmer>(rhs_sketch), minimize_config.kmer_length * 2,
            KmerValueProjection);

  std::uint64_t rhs_id = rhs->id;

  std::vector<Match> matches;
  for (std::uint32_t i = 0, j = 0; i < lhs_sketch.size(); ++i) {
    while (j < rhs_sketch.size()) {
      if (lhs_sketch[i].value < rhs_sketch[j].value) {
        break;
      } else if (lhs_sketch[i].value == rhs_sketch[j].value) {
        for (std::uint32_t k = j; k < rhs_sketch.size(); ++k) {
          if (lhs_sketch[i].value != rhs_sketch[k].value) {
            break;
          }

          std::uint64_t strand =
              (lhs_sketch[i].strand() & 1) == (rhs_sketch[k].strand() & 1);
          std::uint64_t lhs_pos = lhs_sketch[i].position();
          std::uint64_t rhs_pos = rhs_sketch[k].position();
          std::uint64_t diagonal =
              !strand ? rhs_pos + lhs_pos : rhs_pos - lhs_pos + (3ULL << 30);

          matches.emplace_back((((rhs_id << 1) | strand) << 32) | diagonal,
                               (lhs_pos << 32) | rhs_pos);
        }
        break;
      } else {
        ++j;
      }
    }
  }

  return Chain(lhs->id, std::move(matches), chain_config);
}

std::vector<biosoup::Overlap> Chain(std::uint64_t lhs_id,
                                    std::vector<Match>&& matches,
                                    ChainConfig config) {
  RadixSort(std::span<Match>(matches), 64, MatchGroupProjection);
  matches.emplace_back(-1, -1);  // stop dummy

  std::vector<std::pair<std::uint64_t, std::uint64_t>> intervals;
  for (std::uint64_t i = 1, j = 0; i < matches.size(); ++i) {  // NOLINT
    if (matches[i].group - matches[j].group > config.bandwidth) {
      if (i - j >= 4) {
        if (!intervals.empty() && intervals.back().second > j) {  // extend
          intervals.back().second = i;
        } else {  // new
          intervals.emplace_back(j, i);
        }
      }
      ++j;
      while (j < i && matches[i].group - matches[j].group > config.bandwidth) {
        ++j;
      }
    }
  }

  std::vector<biosoup::Overlap> dst;
  for (const auto& it : intervals) {
    std::uint64_t j = it.first;
    std::uint64_t i = it.second;

    if (i - j < config.chain) {
      continue;
    }

    RadixSort(std::span(matches.begin() + j, matches.begin() + i), 64,
              MatchPositionProjection);

    std::uint64_t strand = matches[j].strand();

    std::vector<std::uint64_t> indices;
    if (strand) {                         // same strand
      indices = LongestMatchSubsequence(  // increasing
          std::span(matches.cbegin() + j, matches.cbegin() + i),
          [](std::uint32_t lhs, std::uint32_t rhs) noexcept -> bool {
            return lhs < rhs;
          });
    } else {                              // different strand
      indices = LongestMatchSubsequence(  // decreasing
          std::span(matches.cbegin() + j, matches.cbegin() + i),
          [](std::uint32_t lhs, std::uint32_t rhs) noexcept -> bool {
            return lhs > rhs;
          });
    }

    if (indices.size() < config.chain) {
      continue;
    }

    indices.emplace_back(matches.size() - 1 - j);  // stop dummy from above
    for (std::uint64_t k = 1, l = 0; k < indices.size(); ++k) {
      if (matches[j + indices[k]].lhs_position() -
              matches[j + indices[k - 1]].lhs_position() >
          config.gap) {
        if (k - l < config.chain) {
          l = k;
          continue;
        }

        std::uint32_t lhs_matches = 0;
        std::uint32_t lhs_begin = 0;
        std::uint32_t lhs_end = 0;
        std::uint32_t rhs_matches = 0;
        std::uint32_t rhs_begin = 0;
        std::uint32_t rhs_end = 0;

        for (std::uint64_t m = l; m < k; ++m) {
          std::uint32_t lhs_pos = matches[j + indices[m]].lhs_position();
          if (lhs_pos > lhs_end) {
            lhs_matches += lhs_end - lhs_begin;
            lhs_begin = lhs_pos;
          }
          lhs_end = lhs_pos + config.kmer_length;

          std::uint32_t rhs_pos = matches[j + indices[m]].rhs_position();
          rhs_pos = strand ? rhs_pos
                           : (1U << 31) - (rhs_pos + config.kmer_length - 1);
          if (rhs_pos > rhs_end) {
            rhs_matches += rhs_end - rhs_begin;
            rhs_begin = rhs_pos;
          }
          rhs_end = rhs_pos + config.kmer_length;
        }
        lhs_matches += lhs_end - lhs_begin;
        rhs_matches += rhs_end - rhs_begin;
        if (std::min(lhs_matches, rhs_matches) < config.min_matches) {
          l = k;
          continue;
        }

        dst.emplace_back(
            lhs_id, matches[j + indices[l]].lhs_position(),
            config.kmer_length + matches[j + indices[k - 1]].lhs_position(),
            matches[j].rhs_id(),
            strand ? matches[j + indices[l]].rhs_position()
                   : matches[j + indices[k - 1]].rhs_position(),
            config.kmer_length +
                (strand ? matches[j + indices[k - 1]].rhs_position()
                        : matches[j + indices[l]].rhs_position()),
            std::min(lhs_matches, rhs_matches), strand);

        l = k;
      }
    }
  }
  return dst;
}

}  // namespace ram
