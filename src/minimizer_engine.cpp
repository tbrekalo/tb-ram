// Copyright (c) 2020 Robert Vaser

#include "ram/minimizer_engine.hpp"

#include <stdexcept>

#include "ram/algorithm.hpp"

namespace ram {

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
              return ::ram::Minimize(*it, MinimizeConfig{
                                              .kmer_length = k_,
                                              .window_length = w_,
                                              .minhash = minhash,
                                          });
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

            RadixSort<Kmer>(std::span<Kmer>(minimizers[i]), k_ * 2,
                            KmerValueProjection);

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
  auto sketch = ::ram::Minimize(sequence, MinimizeConfig{.kmer_length = k_,
                                                         .window_length = w_,
                                                         .minhash = minhash});
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
    const Kmer* kmer;
    std::uint32_t n;
    const uint64_t* origins;

    Hit(const Kmer* kmer, std::uint32_t n, const uint64_t* origins)
        : kmer(kmer), n(n), origins(origins) {}

    bool operator<(const Hit& other) const { return n < other.n; }
  };
  std::vector<Hit> filtered_hits;

  std::uint64_t mask = index_.size() - 1;
  std::uint32_t prev = 0;

  sketch.emplace_back(-1, sequence->inflated_len << 1);  // stop dummy

  for (const auto& kmer : sketch) {
    std::uint32_t i = kmer.value & mask;
    const uint64_t* origins = nullptr;
    auto n = index_[i].Find(kmer.value, &origins);
    if (n > occurrence_) {
      filtered_hits.emplace_back(&kmer, n, origins);
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

  return Chain(sequence->id, std::move(matches),
               ChainConfig{.kmer_length = k_,
                           .bandwidth = bandwidth_,
                           .chain = chain_,
                           .min_matches = matches_,
                           .gap = gap_});
}

}  // namespace ram
