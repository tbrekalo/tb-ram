#ifndef RAM_ALGORITHM_HPP_
#define RAM_ALGORITHM_HPP_

#include <algorithm>
#include <concepts>
#include <memory>
#include <span>
#include <vector>

#include "biosoup/overlap.hpp"
#include "ram/types.hpp"

namespace biosoup {
class NucleicAcid;
}

namespace ram {

// Projects a type T into an integral type which can be used in a radix sort.
template <class T, class P>
concept RadixProjection = requires(const T& t, const P& p) {
  { p(t) } noexcept -> std::same_as<std::uint64_t>;
};

// Performns an in-place radix sort.
template <class T, class Proj>
  requires(RadixProjection<T, Proj>)
inline void RadixSort(std::span<T> values, std::uint8_t max_bits, Proj proj) {
  if (values.empty()) {
    return;
  }

  std::vector<T> buffer(values.size());
  std::span sorted_values(buffer);

  std::uint64_t bucket_indices[0x100]{};  // 256 b
  std::uint8_t shift = 0;
  for (; shift < max_bits; shift += 8) {
    std::uint64_t counts[0x100]{};  // 256 b
    for (const auto& it : values) {
      ++counts[(proj(it) >> shift) & 0xFF];
    }

    // exclusive scan
    for (std::uint64_t bucket_idx = 0, count = 0; bucket_idx < 0x100;
         count += counts[bucket_idx++]) {
      bucket_indices[bucket_idx] = count;
    }

    // copies values from old array into sorted buckets
    for (const auto& it : values) {
      sorted_values[bucket_indices[(proj(it) >> shift) & 0xFF]++] =
          std::move(it);
    }

    // std::span<T> performs a shallow swap.
    std::swap(values, sorted_values);
  }

  if ((shift / 8) & 1) {
    std::ranges::copy(values, sorted_values.begin());
  }
}

// Requirements for the comparison operator passed to the
// LongestMatchSubsequence function.
template <class T>
concept BinaryMatchComparator = requires(const T& t, std::uint32_t arg) {
  { t(arg, arg) } noexcept -> std::same_as<bool>;
};

// Expects the input matches to be sorted by the lhs position.
template <BinaryMatchComparator C>
inline std::vector<std::uint64_t> LongestMatchSubsequence(
    std::span<const Match> matches, C comp) {
  if (matches.empty()) {
    return std::vector<std::uint64_t>{};
  }

  std::vector<std::uint64_t> minimal(matches.size() + 1, 0);
  std::vector<std::uint64_t> predecessor(matches.size(), 0);

  std::uint64_t longest = 0;
  for (auto idx = 0uz; idx < matches.size(); ++idx) {
    std::uint64_t lo = 1, hi = longest;
    while (lo <= hi) {
      std::uint64_t mid = lo + (hi - lo) / 2;
      if (matches[minimal[mid]].lhs_position() < matches[idx].lhs_position() &&
          comp(matches[minimal[mid]].rhs_position(),
               matches[idx].rhs_position())) {
        lo = mid + 1;
      } else {
        hi = mid - 1;
      }
    }

    predecessor[idx] = minimal[lo - 1];
    minimal[lo] = idx;
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

struct MinimizeConfig {
  std::uint32_t kmer_length = 15;
  std::uint32_t window_length = 5;
  bool minhash = false;
};

std::vector<Kmer> Minimize(
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    MinimizeConfig config);

struct ChainConfig {
  std::uint32_t kmer_length = 15;
  std::uint32_t bandwidth = 500;
  std::uint32_t chain = 4;
  std::uint32_t min_matches = 100;
  std::uint64_t gap = 10'000;
};

std::vector<biosoup::Overlap> Chain(std::uint64_t lhs_id,
                                    std::vector<Match>&& matches,
                                    ChainConfig config);

// Find overlaps between a pair of sequences.
// Minhash argument from configuration will only be applied on the lhs sequence.
std::vector<biosoup::Overlap> Map(
    const std::unique_ptr<biosoup::NucleicAcid>& lhs,
    const std::unique_ptr<biosoup::NucleicAcid>& rhs,
    MinimizeConfig minimize_config, ChainConfig chain_config);

}  // namespace ram

#endif  // RAM_ALGORITHM_HPP_