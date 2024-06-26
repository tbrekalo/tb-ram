#include "ram/algorithm.hpp"

#include <cmath>
#include <deque>
#include <format>

#include "biosoup/nucleic_acid.hpp"
#include "glog/logging.h"

namespace ram {

namespace {

// Requirements for the comparison operator passed to a
// function for finding match subsequences
template <class T>
concept BinaryMatchComparator = requires(const T& t, std::uint32_t arg) {
  { t(arg, arg) } noexcept -> std::same_as<bool>;
};

struct ChainEntry {
  std::int64_t match_idx;
  double score;
};

struct ScoredMatchSequences {
  std::vector<std::uint64_t> indices;
  double score;
};

auto CreateMatchSequenceComparator(std::span<Match> matches,
                                   std::uint64_t strand) {
  auto comp =
      // same strand; increasing
      strand ? +[](std::uint32_t lhs, std::uint32_t rhs) noexcept -> bool {
    return lhs < rhs;
  }
  // different strand; decreasing
  : +[](std::uint32_t lhs, std::uint32_t rhs) noexcept -> bool {
      return lhs > rhs;
    };

  return [matches, comp](ChainEntry head, const Match& cur_match) -> bool {
    return matches[head.match_idx].lhs_position() < cur_match.lhs_position() &&
           comp(matches[head.match_idx].rhs_position(),
                cur_match.rhs_position());
  };
}

auto CreateScoreCalculator(ChainConfig config, std::span<Match> matches,
                           std::uint64_t strand) {
  return [matches, strand, w = config.kmer_length, g = config.gap,
          bandwidth = config.bandwidth](const std::uint32_t prev_idx,
                                        const std::uint32_t cur_idx) -> double {
    DCHECK(prev_idx < cur_idx);
    auto beta = [&] -> double {
      if (matches[cur_idx].lhs_position() <= matches[prev_idx].lhs_position()) {
        return -(1e9 + 11);
      }

      auto calc_pos_dif = [](std::int32_t prev_pos,
                             std::int32_t cur_pos) -> std::int32_t {
        return std::max(cur_pos - prev_pos, 0);
      };

      auto query_gap = calc_pos_dif(matches[prev_idx].lhs_position() + w,
                                    matches[cur_idx].lhs_position());
      auto target_gap = strand
                            ? calc_pos_dif(matches[prev_idx].rhs_position() + w,
                                           matches[cur_idx].rhs_position())
                            : calc_pos_dif(matches[cur_idx].rhs_position() + w,
                                           matches[prev_idx].rhs_position());

      auto abs_diff =
          static_cast<std::uint32_t>(std::abs(query_gap - target_gap));
      if (std::max(query_gap, target_gap) >= g || abs_diff > bandwidth) {
        return -(1e9 + 11);
      }

      return abs_diff;
    }();

    if (beta < 0) {
      return -(1e9 + 11);
    }

    auto alpha = [&] {
      DCHECK(matches[cur_idx].lhs_position() >=
             matches[prev_idx].lhs_position());
      if (strand) {
        DCHECK(matches[cur_idx].rhs_position() >=
               matches[prev_idx].rhs_position())
            << std::format("strand={} matches[{}]={} matches[{}]={}", strand,
                           cur_idx, matches[cur_idx].rhs_position(), prev_idx,
                           matches[prev_idx].rhs_position());
      } else {
        DCHECK(matches[cur_idx].rhs_position() <=
               matches[prev_idx].rhs_position())
            << std::format("strand={} matches[{}]={} matches[{}]={} ", strand,
                           cur_idx, matches[cur_idx].rhs_position(), prev_idx,
                           matches[prev_idx].rhs_position());
      }

      // THIS IS WRONG!!!! TODO: FIX IT!!!
      return std::min(std::min(matches[cur_idx].lhs_position() -
                                   matches[prev_idx].lhs_position(),
                               strand ? matches[cur_idx].rhs_position() -
                                            matches[prev_idx].rhs_position()
                                      : matches[prev_idx].rhs_position() -
                                            matches[cur_idx].rhs_position()),
                      w);
    }();

    auto gama = (beta > 0 ? 0.01 * w * beta + 0.5 * std::log2(beta) : 0.);

    return alpha - gama;
  };
}

}  // namespace

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

std::vector<Index> ConstructIndices(
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> sequences,
    MinimizeConfig minimize_config,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool) {
  std::vector<Index> indices(
      1uz << std::min(14uz, 2uz * minimize_config.kmer_length));
  if (sequences.empty()) {
    return indices;
  }

  std::vector<std::vector<Kmer>> minimizers(indices.size());
  {
    std::uint64_t mask = indices.size() - 1;
    std::vector<std::vector<Kmer>> buffer;
    for (auto idx = 0uz; idx != sequences.size();) {
      const auto batch_first_idx = idx;
      idx = [sequences](std::size_t start_idx) {
        auto curr_idx = start_idx;
        for (auto batch_size = 0uz;
             curr_idx < sequences.size() && batch_size < 5'0000'000uz;
             ++curr_idx) {
          batch_size += sequences[curr_idx]->inflated_len;
        }
        return curr_idx;
      }(idx);

      buffer.resize(idx - batch_first_idx);
      std::vector<std::future<void>> futures;
      for (auto i = batch_first_idx; i < idx; ++i) {
        futures.push_back(thread_pool->Submit(
            [&sequences, &minimize_config, &buffer,
             batch_first_idx](std::size_t minimize_idx) -> void {
              buffer[minimize_idx - batch_first_idx] =
                  Minimize(sequences[minimize_idx], minimize_config);
            },
            i));
      }

      for (auto& future : futures) {
        future.wait();
      }

      for (auto& kmers : buffer) {
        for (auto kmer : kmers) {
          auto& m = minimizers[kmer.value & mask];
          if (m.capacity() == m.size()) {
            m.reserve(m.capacity() * 1.5);
          }
          m.push_back(kmer);
        }

        std::vector<Kmer>{}.swap(kmers);
      }
    }
  }

  {
    std::vector<std::future<std::pair<std::size_t, std::size_t>>> futures;
    for (std::uint32_t i = 0; i < minimizers.size(); ++i) {
      futures.emplace_back(thread_pool->Submit(
          [&](std::uint32_t i) -> std::pair<std::size_t, std::size_t> {
            if (minimizers[i].empty()) {
              return std::make_pair(0, 0);
            }

            RadixSort<Kmer>(std::span<Kmer>(minimizers[i]),
                            minimize_config.kmer_length * 2,
                            KmerValueProjection);

            minimizers[i].emplace_back(-1, -1);  // stop dummy

            std::size_t num_origins = 0;
            std::size_t num_keys = 0;

            for (std::uint64_t j = 1, c = 1; j < minimizers[i].size();
                 ++j, ++c) {
              if (minimizers[i][j - 1].value != minimizers[i][j].value) {
                if (c > 1) {
                  num_origins += c;
                }
                ++num_keys;
                c = 0;
              }
            }

            return std::make_pair(num_origins, num_keys);
          },
          i));
    }
    for (std::uint32_t i = 0; i < minimizers.size(); ++i) {
      if (minimizers[i].empty()) {
        continue;
      }

      auto [num_origins, num_keys] = futures[i].get();
      indices[i].origins.reserve(num_origins);
      indices[i].locator.reserve(num_keys);

      for (std::uint64_t j = 1, c = 1; j < minimizers[i].size(); ++j, ++c) {
        if (minimizers[i][j - 1].value != minimizers[i][j].value) {
          if (c == 1) {
            indices[i].locator.emplace(minimizers[i][j - 1].value << 1 | 1,
                                       minimizers[i][j - 1].origin);
          } else {
            indices[i].locator.emplace(minimizers[i][j - 1].value << 1,
                                       indices[i].origins.size() << 32 | c);
            for (std::uint64_t k = j - c; k < j; ++k) {
              indices[i].origins.emplace_back(minimizers[i][k].origin);
            }
          }
          c = 0;
        }
      }

      std::vector<Kmer>().swap(minimizers[i]);
    }
  }

  return indices;
}

std::uint32_t CalculateKmerThreshold(std::vector<Index> indices,
                                     double frequency) {
  if (!(0 <= frequency && frequency <= 1)) {
    throw std::invalid_argument(
        "[ram::MinimizerEngine::Filter] error: invalid frequency");
  }

  if (frequency == 0) {
    return -1;
  }

  std::vector<std::uint32_t> occurrences;
  for (const auto& it : indices) {
    for (const auto& jt : it.locator) {
      if (jt.first & 1) {
        occurrences.emplace_back(1);
      } else {
        occurrences.emplace_back(static_cast<std::uint32_t>(jt.second));
      }
    }
  }

  if (occurrences.empty()) {
    return -1;
  }

  std::nth_element(occurrences.begin(),
                   occurrences.begin() + (1 - frequency) * occurrences.size(),
                   occurrences.end());
  return occurrences[(1 - frequency) * occurrences.size()] + 1;
}

std::vector<Match> MatchPairs(const std::unique_ptr<biosoup::NucleicAcid>& lhs,
                              const std::unique_ptr<biosoup::NucleicAcid>& rhs,
                              MinimizeConfig minimize_config) {
  auto lhs_sketch = Minimize(lhs, minimize_config);
  if (lhs_sketch.empty()) {
    return std::vector<Match>{};
  }

  auto rhs_sketch = Minimize(
      rhs, MinimizeConfig{.kmer_length = minimize_config.kmer_length,
                          .window_length = minimize_config.window_length});
  if (rhs_sketch.empty()) {
    return std::vector<Match>{};
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

  return matches;
}

std::vector<biosoup::Overlap> MapPairs(
    const std::unique_ptr<biosoup::NucleicAcid>& lhs,
    const std::unique_ptr<biosoup::NucleicAcid>& rhs,
    MinimizeConfig minimize_config, ChainConfig chain_config) {
  return Chain(lhs->id, MatchPairs(lhs, rhs, minimize_config), chain_config);
}

static std::vector<std::pair<std::uint64_t, std::uint64_t>>
FindMatchGroupIntervals(ChainConfig config, std::span<const Match> matches) {
  std::vector<std::pair<std::uint64_t, std::uint64_t>> intervals;
  for (std::uint64_t i = 1, j = 0; i < matches.size(); ++i) {
    if (matches[i].group - matches[j].group > config.bandwidth) {
      if (auto len = i - j; len >= config.chain) {
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

  DCHECK(std::ranges::all_of(
      intervals,
      [matches](std::pair<std::int64_t, std::int64_t> interval) -> bool {
        return std::ranges::all_of(
            matches.subspan(interval.first, interval.second - interval.first),
            [rhs_id = matches[interval.first].rhs_id()](
                std::uint32_t id) -> bool { return rhs_id == id; },
            &Match::rhs_id);
      }));

  return intervals;
}

static std::vector<std::pair<std::uint64_t, std::uint64_t>>
FindMatchReadIntervals(ChainConfig config, std::span<const Match> matches) {
  std::vector<std::pair<std::uint64_t, std::uint64_t>> intervals;
  for (std::uint64_t i = 1, j = 0; i < matches.size(); ++i) {
    if (matches[i].rhs_id() != matches[j].rhs_id() ||
        matches[i].strand() != matches[j].strand()) {
      if (i - j > config.chain) {
        intervals.emplace_back(j, i);
      }
      j = i;
    }
  }

  return intervals;
}

static std::vector<ScoredMatchSequences> CreateMatchSequences(
    ChainConfig config, std::span<Match> matches, std::uint64_t strand) {
  DCHECK(matches.empty());
  DCHECK(std::ranges::all_of(
      matches,
      [rhs_id = matches.front().rhs_id()](std::uint64_t id) -> bool {
        return rhs_id == id;
      },
      &Match::rhs_id));
  std::vector<ScoredMatchSequences> dst;

  std::vector<std::int64_t> predecessor(matches.size(), -1);
  std::vector<ChainEntry> heads(matches.size());
  std::vector<ChainEntry> sorted_matches;

  auto score_fn = CreateScoreCalculator(config, matches, strand);
  auto cmp_fn = CreateMatchSequenceComparator(matches, strand);
  auto assert_prev_idx = [strand, matches, &sorted_matches](
                             std::int64_t sorted_matches_idx,
                             std::int64_t cur_match_idx) -> void {
    if (sorted_matches_idx > 0 && strand) {
      DCHECK(
          matches[cur_match_idx].rhs_position() >=
          matches[sorted_matches[sorted_matches_idx].match_idx].rhs_position())
          << std::format("strand={} matches[{}]={} matches[{}]={} ", strand,
                         cur_match_idx, matches[cur_match_idx].rhs_position(),
                         sorted_matches[sorted_matches_idx].match_idx,
                         matches[sorted_matches[sorted_matches_idx].match_idx]
                             .rhs_position());
    } else if (sorted_matches_idx > 0 && !strand) {
      DCHECK(
          matches[cur_match_idx].rhs_position() <=
          matches[sorted_matches[sorted_matches_idx].match_idx].rhs_position())
          << std::format("strand={} matches[{}]={} matches[{}]={} ", strand,
                         cur_match_idx, matches[cur_match_idx].rhs_position(),
                         sorted_matches[sorted_matches_idx].match_idx,
                         matches[sorted_matches[sorted_matches_idx].match_idx]
                             .rhs_position());
    }
  };

  sorted_matches.push_back({.match_idx = -1, .score = 0.});
  for (std::int64_t match_idx = 0;
       match_idx < static_cast<std::int64_t>(matches.size()); ++match_idx) {
    std::int64_t cur_sorted_matches_idx =
        std::lower_bound(sorted_matches.begin() + 1, sorted_matches.end(),
                         matches[match_idx], cmp_fn) -
        sorted_matches.begin();

    auto prev_sorted_matches_idx = cur_sorted_matches_idx - 1;
    if (static_cast<std::uint32_t>(cur_sorted_matches_idx) ==
        sorted_matches.size()) {
      sorted_matches.resize(sorted_matches.size() + 1);
    }

    assert_prev_idx(prev_sorted_matches_idx, match_idx);
    auto score =
        prev_sorted_matches_idx > 0
            ? score_fn(sorted_matches[prev_sorted_matches_idx].match_idx,
                       match_idx) +
                  sorted_matches[prev_sorted_matches_idx].score
            : config.kmer_length;

    for (std::int64_t candidate_prev = prev_sorted_matches_idx - 1, k = 50;
         candidate_prev > 0 && k > 0; --candidate_prev, --k) {
      if (strand and
          matches[sorted_matches[candidate_prev].match_idx].rhs_position() >=
              matches[match_idx].rhs_position()) {
        continue;
      }

      if (!strand and
          matches[sorted_matches[candidate_prev].match_idx].rhs_position() <=
              matches[match_idx].rhs_position()) {
        continue;
      }

      assert_prev_idx(candidate_prev, match_idx);
      if (auto new_score =
              score_fn(sorted_matches[candidate_prev].match_idx, match_idx) +
              sorted_matches[candidate_prev].score;
          new_score > score) {
        score = new_score;
        prev_sorted_matches_idx = candidate_prev;
        break;
      }
    }

    sorted_matches[cur_sorted_matches_idx] = heads[match_idx] = {
        .match_idx = match_idx, .score = std::max(score, 0.)};
    predecessor[match_idx] =
        score > 0 ? sorted_matches[prev_sorted_matches_idx].match_idx : -1;
  }

  std::ranges::sort(heads, std::ranges::greater{}, &ChainEntry::score);
  std::vector<std::uint8_t> visited(matches.size(), 0);
  for (auto [head_idx, head_score] : heads) {
    if (head_score < config.min_matches) {
      break;
    }

    auto unvisit = [&] {
      for (auto chain_node = head_idx; chain_node != -1 && !visited[chain_node];
           chain_node = predecessor[chain_node]) {
        visited[chain_node] = 0;
      }
    };

    std::vector<std::uint64_t> chain;
    for (auto chain_node = head_idx; chain_node != -1 && !visited[chain_node];
         chain_node = predecessor[chain_node]) {
      chain.push_back(chain_node);
      visited[chain_node] = 1;
    }

    if (chain.size() < config.chain) {
      unvisit();
      continue;
    }

    double score = 0.;
    std::ranges::reverse(chain);

    for (auto chain_idx = 1uz; chain_idx < chain.size(); ++chain_idx) {
      DCHECK(matches[chain[chain_idx - 1]].lhs_position() <=
             matches[chain[chain_idx]].lhs_position());
      if (strand) {
        DCHECK(matches[chain[chain_idx - 1]].rhs_position() <=
               matches[chain[chain_idx]].rhs_position());
      } else {
        DCHECK(matches[chain[chain_idx - 1]].rhs_position() >=
               matches[chain[chain_idx]].rhs_position());
      }

      score += score_fn(chain[chain_idx - 1], chain[chain_idx]);
    }

    if (score < config.min_matches) {
      unvisit();
      continue;
    }

    dst.push_back({.indices = std::move(chain), .score = score});
  }

  return dst;
}

static std::vector<ScoredMatchSequences> FindChainIndices(
    ChainConfig config, std::span<Match> matches, std::uint64_t lhs_idx,
    std::uint64_t rhs_idx, std::uint64_t strand) {
  auto match_sequence = CreateMatchSequences(
      config, std::span(matches.begin() + lhs_idx, matches.begin() + rhs_idx),
      strand);

  {
    auto [first, last] = std::ranges::remove_if(
        match_sequence, [config](const ScoredMatchSequences& arg) -> bool {
          return arg.indices.size() < config.chain;
        });
    match_sequence.erase(first, last);
  }

  return match_sequence;
}

std::vector<MatchChain> FindChainMatches([[maybe_unused]] std::uint32_t lhs_id,
                                         std::vector<Match>&& matches,
                                         ChainConfig config) {
  std::ranges::sort(matches, [](const Match& a, const Match& b) -> bool {
    return a.rhs_id() != b.rhs_id() ? a.rhs_id() < b.rhs_id()
                                    : a.strand() < b.strand();
  });
  matches.emplace_back(-1, -1);  // stop dummy

  auto intervals = FindMatchReadIntervals(config, matches);
  std::vector<MatchChain> match_chains;

  for (const auto& it : intervals) {
    std::uint64_t lhs_idx = it.first;
    std::uint64_t rhs_idx = it.second;
    RadixSort(std::span(matches.begin() + lhs_idx, matches.begin() + rhs_idx),
              64, MatchPositionProjection);
    std::uint64_t strand = matches[lhs_idx].strand();
    auto indices = FindChainIndices(config, matches, lhs_idx, rhs_idx, strand);

    for (const auto& [index_vec, score] : indices) {
      DCHECK(not index_vec.empty());
      for (auto ik = 1uz; ik < index_vec.size(); ++ik) {
        DCHECK(matches[lhs_idx + index_vec[ik - 1]].lhs_position() <=
               matches[lhs_idx + index_vec[ik]].lhs_position());
        if (strand) {
          DCHECK(matches[lhs_idx + index_vec[ik - 1]].rhs_position() <=
                 matches[lhs_idx + index_vec[ik]].rhs_position());
          continue;
        }

        DCHECK(matches[lhs_idx + index_vec[ik - 1]].rhs_position() >=
               matches[lhs_idx + index_vec[ik]].rhs_position());
      }

      auto local_matches = std::vector<Match>();
      std::uint32_t lhs_matches = 0;
      std::uint32_t lhs_begin = 0;
      std::uint32_t lhs_end = 0;
      std::uint32_t rhs_matches = 0;
      std::uint32_t rhs_begin = 0;
      std::uint32_t rhs_end = 0;

      for (auto idx : index_vec) {
        std::uint32_t lhs_pos = matches[lhs_idx + idx].lhs_position();
        if (lhs_pos > lhs_end) {
          lhs_matches += lhs_end - lhs_begin;
          lhs_begin = lhs_pos;
        }
        lhs_end = lhs_pos + config.kmer_length;

        std::uint32_t rhs_pos = matches[lhs_idx + idx].rhs_position();
        rhs_pos =
            strand ? rhs_pos : (1U << 31) - (rhs_pos + config.kmer_length - 1);
        if (rhs_pos > rhs_end) {
          rhs_matches += rhs_end - rhs_begin;
          rhs_begin = rhs_pos;
        }
        rhs_end = rhs_pos + config.kmer_length;
        local_matches.push_back(matches[lhs_idx + idx]);
      }
      lhs_matches += lhs_end - lhs_begin;
      rhs_matches += rhs_end - rhs_begin;

      auto lhs_overlap_len = [&] {
        DCHECK(matches[lhs_idx + index_vec.back()].lhs_position() >=
               matches[lhs_idx + index_vec.front()].lhs_position())
            << std::format("lhs_end={} lhs_start={}",
                           matches[lhs_idx + index_vec.back()].lhs_position(),
                           matches[lhs_idx + index_vec.front()].lhs_position());

        return matches[lhs_idx + index_vec.back()].lhs_position() -
               matches[lhs_idx + index_vec.front()].lhs_position();
      }();

      auto rhs_overlap_len = [&] {
        if (strand) {
          DCHECK(matches[lhs_idx + index_vec.back()].rhs_position() >=
                 matches[lhs_idx + index_vec.front()].rhs_position())
              << std::format("rhs_end={} rhs_start={}",
                             matches[index_vec.back()].rhs_position(),
                             matches[index_vec.front()].rhs_position());

          return matches[lhs_idx + index_vec.back()].rhs_position() -
                 matches[lhs_idx + index_vec.front()].rhs_position();
        }

        DCHECK(matches[lhs_idx + index_vec.front()].rhs_position() >=
               matches[lhs_idx + index_vec.back()].rhs_position())
            << std::format("rhs_end={} rhs_start={}",
                           matches[lhs_idx + index_vec.front()].rhs_position(),
                           matches[lhs_idx + index_vec.back()].rhs_position());

        return matches[lhs_idx + index_vec.front()].rhs_position() -
               matches[lhs_idx + index_vec.back()].rhs_position();
      }();

      const auto lhs_matches_normed = 1. * lhs_matches / lhs_overlap_len;
      const auto rhs_matches_normed = 1. * rhs_matches / rhs_overlap_len;

      if (std::fabs(rhs_matches_normed - lhs_matches_normed) >= 0.005) {
        continue;
      }

      match_chains.push_back(MatchChain{.matches = std::move(local_matches),
                                        .lhs_matches = lhs_matches,
                                        .rhs_matches = rhs_matches,
                                        .score = score});
    }
  }

  return match_chains;
}

std::vector<OverlapAI> ChainAI(std::uint32_t lhs_id,
                               std::vector<Match>&& matches,
                               ChainConfig config) {
  std::vector<OverlapAI> dst;
  std::vector<std::uint32_t> diffs;
  for (const auto& [chain_matches, lhs_matches, rhs_matches, score] :
       FindChainMatches(lhs_id, std::move(matches), config)) {
    if (chain_matches.empty()) {
      continue;
    }

    diffs.resize(chain_matches.size() - 1);
    for (auto i = 0uz; i + 1uz < chain_matches.size(); ++i) {
      diffs[i] = chain_matches[i + 1uz].lhs_position() -
                 chain_matches[i].lhs_position();
    }
    std::ranges::sort(diffs);

    const auto strand = chain_matches.front().strand();
    dst.push_back(OverlapAI{
        .lhs_id = lhs_id,
        .lhs_begin = chain_matches.front().lhs_position(),
        .lhs_end = config.kmer_length + chain_matches.back().lhs_position(),
        .lhs_matches = lhs_matches,

        .strand = strand,

        .rhs_id = chain_matches.front().rhs_id(),
        .rhs_begin = strand ? chain_matches.front().rhs_position()
                            : chain_matches.back().rhs_position(),
        .rhs_end = config.kmer_length +
                   (strand ? chain_matches.back().rhs_position()
                           : chain_matches.front().rhs_position()),
        .rhs_matches = rhs_matches,

        .score = score,
        .diff_mean =
            [&diffs] {
              return std::accumulate(diffs.begin(), diffs.end(), 0.,
                                     std::plus<double>{}) /
                     diffs.size();
            }(),

        .q75 = diffs[diffs.size() * 0.75],
        .q90 = diffs[diffs.size() * 0.90],
        .q95 = diffs[diffs.size() * 0.95],
        .q98 = diffs[diffs.size() * 0.98],
    });
  }

  return dst;
}

std::vector<biosoup::Overlap> Chain(std::uint32_t lhs_id,
                                    std::vector<Match>&& matches,
                                    ChainConfig config) {
  std::vector<biosoup::Overlap> dst;
  for (const auto& jt : FindChainMatches(lhs_id, std::move(matches), config)) {
    if (jt.matches.empty()) {
      continue;
    }

    const auto strand = jt.matches.front().strand();
    dst.emplace_back(
        lhs_id, jt.matches.front().lhs_position(),
        config.kmer_length + jt.matches.back().lhs_position(),
        jt.matches.front().rhs_id(),
        strand ? jt.matches.front().rhs_position()
               : jt.matches.back().rhs_position(),
        config.kmer_length + (strand ? jt.matches.back().rhs_position()
                                     : jt.matches.front().rhs_position()),
        std::min(jt.lhs_matches, jt.rhs_matches), strand);
  }

  return dst;
}

std::vector<Match> MatchToIndex(
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    std::span<const Index> indices, MapToIndexConfig map_config,
    MinimizeConfig minimize_config, ChainConfig chain_config,
    std::vector<std::uint32_t>* filtered) {
  auto sketch = Minimize(sequence, minimize_config);
  if (sketch.empty()) {
    return std::vector<Match>{};
  }

  std::vector<Match> matches;
  auto add_match = [&](const Kmer& kmer, uint64_t origin) -> void {
    auto id = [](std::uint64_t origin) -> std::uint32_t {
      return static_cast<std::uint32_t>(origin >> 32);
    };
    auto position = [](std::uint64_t origin) -> std::uint32_t {
      return static_cast<std::uint32_t>(origin) >> 1;
    };
    auto strand = [](std::uint64_t origin) -> bool { return origin & 1; };

    if (map_config.avoid_equal && sequence->id == id(origin)) {
      return;
    }
    if (map_config.avoid_symmetric && sequence->id > id(origin)) {
      return;
    }

    std::uint64_t rhs_id = id(origin);
    std::uint64_t strand_ = kmer.strand() == strand(origin);
    std::uint64_t lhs_pos = kmer.position();
    std::uint64_t rhs_pos = position(origin);
    std::uint64_t diagonal =
        !strand_ ? rhs_pos + lhs_pos : rhs_pos - lhs_pos + (3ULL << 30);

    matches.emplace_back((((rhs_id << 1) | strand_) << 32) | diagonal,
                         (lhs_pos << 32) | rhs_pos);
  };

  struct Hit {
    const Kmer* kmer;
    std::uint32_t n;
    const uint64_t* origins;

    Hit(const Kmer* kmer, std::uint32_t n, const uint64_t* origins)
        : kmer(kmer), n(n), origins(origins) {}

    bool operator<(const Hit& other) const { return n < other.n; }
  };
  std::vector<Hit> filtered_hits;

  std::uint64_t mask = indices.size() - 1;
  std::uint32_t prev = 0;

  sketch.emplace_back(-1, sequence->inflated_len << 1);  // stop dummy

  for (std::uint32_t i = 0; i < sketch.size(); ++i) {
    const auto& kmer = sketch[i];
    std::uint32_t j = kmer.value & mask;
    const uint64_t* origins = nullptr;
    auto n = indices[j].Find(kmer.value, &origins);
    if (n > map_config.occurrence && i + 1 != sketch.size()) {
      filtered_hits.emplace_back(&kmer, n, origins);
      if (filtered) {
        filtered->emplace_back(kmer.position());
      }
      continue;
    }

    std::size_t rescuees =
        std::min(static_cast<std::size_t>(kmer.position() - prev) /
                     chain_config.bandwidth,
                 filtered_hits.size());
    if (rescuees) {
      std::partial_sort(filtered_hits.begin(), filtered_hits.begin() + rescuees,
                        filtered_hits.end());
      for (auto it = filtered_hits.begin(); rescuees; rescuees--, ++it) {
        for (; it->n; it->n--, ++it->origins) {
          add_match(*it->kmer, *it->origins);
        }
      }
    }
    filtered_hits.clear();
    prev = kmer.position();

    for (; n; n--, ++origins) {
      add_match(kmer, *origins);
    }
  }

  return matches;
}

std::vector<MatchChain> ChainOnIndex(
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    std::span<const Index> indices, MapToIndexConfig map_config,
    MinimizeConfig minimize_config, ChainConfig chain_config,
    std::vector<std::uint32_t>* filtered) {
  return FindChainMatches(sequence->id,
                          MatchToIndex(sequence, indices, map_config,
                                       minimize_config, chain_config, filtered),
                          chain_config);
}

std::vector<biosoup::Overlap> MapToIndex(
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    std::span<const Index> indices, MapToIndexConfig map_config,
    MinimizeConfig minimize_config, ChainConfig chain_config,
    std::vector<std::uint32_t>* filtered) {
  return Chain(sequence->id,
               MatchToIndex(sequence, indices, map_config, minimize_config,
                            chain_config, filtered),
               chain_config);
}

std::vector<OverlapAI> MapToIndexAI(
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    std::span<const Index> indices, MapToIndexConfig map_config,
    MinimizeConfig minimize_config, ChainConfig chain_config,
    std::vector<std::uint32_t>* filtered) {
  return ChainAI(sequence->id,
                 MatchToIndex(sequence, indices, map_config, minimize_config,
                              chain_config, filtered),
                 chain_config);
}

}  // namespace ram
