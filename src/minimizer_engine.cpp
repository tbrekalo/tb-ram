// Copyright (c) 2020 Robert Vaser

#include "ram/minimizer_engine.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <span>
#include <stdexcept>
#include <vector>

namespace ram {

namespace {

template <class T, class Compare>
[[gnu::noinline]] void RadixSort(std::span<T> source, std::uint8_t max_bits,
                                 Compare comp) {  //  unary comparison function
  if (source.empty()) {
    return;
  }

  std::vector<T> buffer(source.size());
  auto to_sort = std::span(source);
  auto sorted = std::span(buffer);

  std::uint8_t shift = 0;
  std::uint64_t counts[0x100]{};
  std::uint64_t buckets[0x100]{};
  for (; shift < max_bits; shift += 8) {
    std::memset(counts, 0, sizeof(counts));
    for (auto it : to_sort) {
      ++counts[std::invoke(comp, it) >> shift & 0xFF];
    }
    if (counts[std::invoke(comp, to_sort.front()) >> shift & 0xFF] ==
        source.size()) {
      continue;
    }

    for (std::uint64_t i = 0, j = 0; i < 0x100; j += counts[i++]) {
      buckets[i] = j;
    }
    for (auto it : to_sort) {
      sorted[buckets[std::invoke(comp, it) >> shift & 0xFF]++] = it;
    }

    std::swap(sorted, to_sort);
  }

  if (to_sort.data() != source.data()) {  // copy the sorted array for odd cases
    std::ranges::copy(to_sort, source.begin());
  }
}

}  // namespace

MinimizerEngine::MinimizerEngine(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, std::uint32_t k,
    std::uint32_t w, std::uint32_t bandwidth, std::uint32_t chain,
    std::uint32_t matches, std::uint32_t gap)
    : k_(std::min(std::max(k, 1U), 31U)),
      w_(w),
      bandwidth_(bandwidth),
      chain_(chain),
      matches_(matches),
      gap_(gap),
      occurrence_(-1),
      index_(1U << std::min(14U, 2 * k_)),
      thread_pool_(thread_pool ? thread_pool
                               : std::make_shared<thread_pool::ThreadPool>(1)) {
}

std::uint32_t MinimizerEngine::Index::Find(std::uint64_t key,
                                           const std::uint64_t** dst) const {
  auto it = locator.find(key << 1);
  if (it == locator.end()) {
    return 0;
  }
  if (it->first & 1) {
    *dst = &(it->second);
    return 1;
  }
  *dst = &(origins[it->second >> 32]);
  return static_cast<std::uint32_t>(it->second);
}

void MinimizerEngine::Minimize(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator first,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator last,
    bool minhash) {
  for (auto& it : index_) {
    it.origins.clear();
    it.locator.clear();
  }

  if (first >= last) {
    return;
  }

  std::vector<std::vector<Kmer>> minimizers(index_.size());
  {
    std::uint64_t mask = index_.size() - 1;

    while (first != last) {
      std::size_t batch_size = 0;
      std::vector<std::future<std::vector<Kmer>>> futures;
      for (; first != last && batch_size < 50000000; ++first) {
        batch_size += (*first)->inflated_len;
        futures.emplace_back(thread_pool_->Submit(
            [&](decltype(first) it) -> std::vector<Kmer> {
              return Minimize(*it, minhash);
            },
            first));
      }
      for (auto& it : futures) {
        for (const auto& jt : it.get()) {
          auto& m = minimizers[jt.value & mask];
          if (m.capacity() == m.size()) {
            m.reserve(m.capacity() * 1.5);
          }
          m.emplace_back(jt);
        }
      }
    }
  }

  {
    std::vector<std::future<std::pair<std::size_t, std::size_t>>> futures;
    for (std::uint32_t i = 0; i < minimizers.size(); ++i) {
      futures.emplace_back(thread_pool_->Submit(
          [&](std::uint32_t i) -> std::pair<std::size_t, std::size_t> {
            if (minimizers[i].empty()) {
              return std::make_pair(0, 0);
            }

            RadixSort(std::span(minimizers[i]), k_ * 2, &Kmer::value);
            minimizers[i].emplace_back(-1, -1);  // stop dummy

            std::size_t num_origins = 0;
            std::size_t num_keys = 0;

            for (std::uint64_t j = 1, c = 1; j < minimizers[i].size();
                 ++j, ++c) {  // NOLINT
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
      auto num_entries = futures[i].get();
      if (minimizers[i].empty()) {
        continue;
      }

      index_[i].origins.reserve(num_entries.first);
      index_[i].locator.reserve(num_entries.second);

      for (std::uint64_t j = 1, c = 1; j < minimizers[i].size(); ++j, ++c) {
        if (minimizers[i][j - 1].value != minimizers[i][j].value) {
          if (c == 1) {
            index_[i].locator.emplace(minimizers[i][j - 1].value << 1 | 1,
                                      minimizers[i][j - 1].origin);
          } else {
            index_[i].locator.emplace(minimizers[i][j - 1].value << 1,
                                      index_[i].origins.size() << 32 | c);
            for (std::uint64_t k = j - c; k < j; ++k) {
              index_[i].origins.emplace_back(minimizers[i][k].origin);
            }
          }
          c = 0;
        }
      }

      std::vector<Kmer>().swap(minimizers[i]);
    }
  }
}

void MinimizerEngine::Filter(double frequency) {
  if (!(0 <= frequency && frequency <= 1)) {
    throw std::invalid_argument(
        "[ram::MinimizerEngine::Filter] error: invalid frequency");
  }

  if (frequency == 0) {
    occurrence_ = -1;
    return;
  }

  std::vector<std::uint32_t> occurrences;
  for (const auto& it : index_) {
    for (const auto& jt : it.locator) {
      if (jt.first & 1) {
        occurrences.emplace_back(1);
      } else {
        occurrences.emplace_back(static_cast<std::uint32_t>(jt.second));
      }
    }
  }

  if (occurrences.empty()) {
    occurrence_ = -1;
    return;
  }

  std::nth_element(occurrences.begin(),
                   occurrences.begin() + (1 - frequency) * occurrences.size(),
                   occurrences.end());
  occurrence_ = occurrences[(1 - frequency) * occurrences.size()] + 1;
}

std::vector<biosoup::Overlap> MinimizerEngine::Map(
    const std::unique_ptr<biosoup::NucleicAcid>& sequence, bool avoid_equal,
    bool avoid_symmetric, bool minhash,
    std::vector<std::uint32_t>* filtered) const {
  auto sketch = Minimize(sequence, minhash);
  if (sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
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

    if (avoid_equal && sequence->id == id(origin)) {
      return;
    }
    if (avoid_symmetric && sequence->id > id(origin)) {
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
    Kmer kmer;
    std::uint32_t n;
    const uint64_t* origins;

    Hit() = default;
    Hit(Kmer kmer, std::uint32_t n, const uint64_t* origins)
        : kmer(kmer), n(n), origins(origins) {}

    bool operator<(const Hit& other) const { return n < other.n; }
  };
  std::vector<Hit> filtered_hits;

  std::uint64_t mask = index_.size() - 1;
  std::uint32_t prev = 0;

  RadixSort(std::span(sketch), k_ * 2, &Kmer::value);
  sketch.emplace_back(-1, sequence->inflated_len << 1);  // stop dummy
  std::vector<Hit> to_add, hits;

  for (const auto kmer : sketch) {
    std::uint32_t group = kmer.value & mask;
    const uint64_t* origins = nullptr;
    auto n = index_[group].Find(kmer.value, &origins);
    hits.emplace_back(kmer, n, origins);
  }

  RadixSort(
      std::span(hits.begin(), hits.end() - 1), 32,
      [](Hit const& hit) -> std::uint32_t { return hit.kmer.position(); });

  for (auto [kmer, n, origins] : hits) {
    if (n > occurrence_) {
      filtered_hits.emplace_back(kmer, n, origins);
      if (filtered) {
        filtered->emplace_back(kmer.position());
      }
      continue;
    }

    std::size_t rescuees =
        std::min(static_cast<std::size_t>(kmer.position() - prev) / bandwidth_,
                 filtered_hits.size());
    if (rescuees) {
      std::partial_sort(filtered_hits.begin(), filtered_hits.begin() + rescuees,
                        filtered_hits.end());
      for (auto it = filtered_hits.begin(); rescuees; rescuees--, ++it) {
        to_add.emplace_back(it->kmer, it->n, it->origins);
      }
    }
    filtered_hits.clear();
    prev = kmer.position();
    to_add.emplace_back(kmer, n, origins);
  }

  RadixSort(std::span(to_add), 2 * k_,
            [](Hit const& hit) -> bool { return hit.kmer.value; });
  for (auto [kmer, n, origins] : to_add) {
    for (; n > 0; --n, ++origins) {
      add_match(kmer, *origins);
    }
  }

  return Chain(sequence->id, std::move(matches));
}

std::vector<biosoup::Overlap> MinimizerEngine::Map(
    const std::unique_ptr<biosoup::NucleicAcid>& lhs,
    const std::unique_ptr<biosoup::NucleicAcid>& rhs, bool minhash) const {
  auto lhs_sketch = Minimize(lhs, minhash);
  if (lhs_sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  auto rhs_sketch = Minimize(rhs);
  if (rhs_sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  RadixSort(std::span(lhs_sketch), k_ * 2, &Kmer::value);
  RadixSort(std::span(rhs_sketch), k_ * 2, &Kmer::value);

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

  return Chain(lhs->id, std::move(matches));
}

std::vector<biosoup::Overlap> MinimizerEngine::Chain(
    std::uint64_t lhs_id, std::vector<Match>&& matches) const {
  RadixSort(std::span(matches), 64, &Match::group);
  matches.emplace_back(-1, -1);  // stop dummy

  std::vector<std::pair<std::uint64_t, std::uint64_t>> intervals;
  for (std::uint64_t i = 1, j = 0; i < matches.size(); ++i) {  // NOLINT
    if (matches[i].group - matches[j].group > bandwidth_) {
      if (i - j >= 4) {
        if (!intervals.empty() && intervals.back().second > j) {  // extend
          intervals.back().second = i;
        } else {  // new
          intervals.emplace_back(j, i);
        }
      }
      ++j;
      while (j < i && matches[i].group - matches[j].group > bandwidth_) {
        ++j;
      }
    }
  }

  std::vector<biosoup::Overlap> dst;
  for (const auto& it : intervals) {
    std::uint64_t j = it.first;
    std::uint64_t i = it.second;

    if (i - j < chain_) {
      continue;
    }

    RadixSort(std::span(matches.begin() + j, matches.begin() + i), 64,
              &Match::positions);
    std::uint64_t strand = matches[j].strand();

    std::vector<std::uint64_t> indices;
    if (strand) {                    // same strand
      indices = LongestSubsequence(  // increasing
          matches.begin() + j, matches.begin() + i, std::less<std::uint64_t>());
    } else {                         // different strand
      indices = LongestSubsequence(  // decreasing
          matches.begin() + j, matches.begin() + i,
          std::greater<std::uint64_t>());
    }

    if (indices.size() < chain_) {
      continue;
    }

    indices.emplace_back(matches.size() - 1 - j);  // stop dummy from above
    for (std::uint64_t k = 1, l = 0; k < indices.size(); ++k) {
      if (matches[j + indices[k]].lhs_position() -
              matches[j + indices[k - 1]].lhs_position() >
          gap_) {
        if (k - l < chain_) {
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
          lhs_end = lhs_pos + k_;

          std::uint32_t rhs_pos = matches[j + indices[m]].rhs_position();
          rhs_pos = strand ? rhs_pos : (1U << 31) - (rhs_pos + k_ - 1);
          if (rhs_pos > rhs_end) {
            rhs_matches += rhs_end - rhs_begin;
            rhs_begin = rhs_pos;
          }
          rhs_end = rhs_pos + k_;
        }
        lhs_matches += lhs_end - lhs_begin;
        rhs_matches += rhs_end - rhs_begin;
        if (std::min(lhs_matches, rhs_matches) < matches_) {
          l = k;
          continue;
        }

        dst.emplace_back(
            lhs_id, matches[j + indices[l]].lhs_position(),
            k_ + matches[j + indices[k - 1]].lhs_position(),
            matches[j].rhs_id(),
            strand ? matches[j + indices[l]].rhs_position()
                   : matches[j + indices[k - 1]].rhs_position(),
            k_ + (strand ? matches[j + indices[k - 1]].rhs_position()
                         : matches[j + indices[l]].rhs_position()),
            std::min(lhs_matches, rhs_matches), strand);

        l = k;
      }
    }
  }
  return dst;
}

std::vector<MinimizerEngine::Kmer> MinimizerEngine::Minimize(
    const std::unique_ptr<biosoup::NucleicAcid>& sequence, bool minhash) const {
  if (sequence->inflated_len < k_ + w_ - 1U) {
    return std::vector<Kmer>{};
  }

  std::uint64_t mask = (1ULL << (k_ * 2)) - 1;

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

  std::int32_t w_begin = 0, w_end = 0;
  std::vector<Kmer> kmers(sequence->inflated_len - k_ + 1U);
  auto window_add = [&](std::uint64_t value, std::uint64_t location) -> void {
    while (w_begin < w_end && kmers[w_end - 1].value > value) {
      --w_end;
    }
    kmers[w_end++] = Kmer(value, location);
  };
  auto window_update = [&](std::uint32_t position) -> void {
    while (w_begin < w_end && kmers[w_begin].position() < position) {
      ++w_begin;
    }
  };

  std::uint64_t shift = (k_ - 1) * 2;
  std::uint64_t minimizer = 0;
  std::uint64_t reverse_minimizer = 0;
  std::uint64_t id = static_cast<std::uint64_t>(sequence->id) << 32;
  std::uint64_t is_stored = 1ULL << 63;

  for (std::uint32_t i = 0, j = 0; i < sequence->inflated_len; ++i) {
    std::uint64_t c = sequence->Code(i);
    minimizer = ((minimizer << 2) | c) & mask;
    reverse_minimizer = (reverse_minimizer >> 2) | ((c ^ 3) << shift);
    if (i + 1U < k_) [[unlikely]] {
      continue;
    }

    std::int32_t flag = minimizer < reverse_minimizer;
    kmers[j++] = Kmer(hash(minimizer * flag + reverse_minimizer * (1 - flag)),
                      id | (i - (k_ - 1U)) << 1 | (1 - flag));
  }

  std::int32_t n = 0;
  std::vector<Kmer> dst(kmers.size() - w_ + 1U);
  for (std::uint32_t i = 0; i < kmers.size(); ++i) {
    window_add(kmers[i].value, kmers[i].origin);
    if (i + 1U < w_) [[unlikely]] {
      continue;
    }
    for (std::int32_t idx = w_begin; idx < w_end; ++idx) {
      if (kmers[idx].value != kmers[w_begin].value) {
        break;
      }
      if (kmers[idx].origin & is_stored) {
        continue;
      }
      dst[n++] = kmers[idx];
      kmers[idx].origin |= is_stored;
    }
    window_update(i - (w_ - 1U) + 1);
  }

  dst.resize(n);
  if (minhash) {
    RadixSort(std::span(dst), k_ * 2, &Kmer::value);
    dst.resize(sequence->inflated_len / k_);
    RadixSort(std::span(dst), 64, &Kmer::origin);
  }

  return dst;
}

template <typename Compare>
std::vector<std::uint64_t> MinimizerEngine::LongestSubsequence(
    std::vector<Match>::const_iterator first,
    std::vector<Match>::const_iterator last,
    Compare comp) {  // binary comparison function
  if (first >= last) {
    return std::vector<std::uint64_t>{};
  }

  std::vector<std::uint64_t> minimal(last - first + 1, 0);
  std::vector<std::uint64_t> predecessor(last - first, 0);

  std::uint64_t longest = 0;
  for (auto it = first; it != last; ++it) {
    std::uint64_t lo = 1, hi = longest;
    while (lo <= hi) {
      std::uint64_t mid = lo + (hi - lo) / 2;
      if ((first + minimal[mid])->lhs_position() < it->lhs_position() &&
          comp((first + minimal[mid])->rhs_position(), it->rhs_position())) {
        lo = mid + 1;
      } else {
        hi = mid - 1;
      }
    }

    predecessor[it - first] = minimal[lo - 1];
    minimal[lo] = it - first;
    longest = std::max(longest, lo);
  }

  std::vector<std::uint64_t> dst;
  for (std::uint64_t i = 0, j = minimal[longest]; i < longest; ++i) {
    dst.emplace_back(j);
    j = predecessor[j];
  }
  std::reverse(dst.begin(), dst.end());

  return dst;
}

}  // namespace ram
