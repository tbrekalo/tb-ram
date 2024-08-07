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

  constexpr auto kNBuckets = 0x100;
  std::uint64_t bucket_indices[kNBuckets]{};  // 256 b
  std::uint8_t shift = 0;
  for (; shift < max_bits; shift += 8) {
    std::uint64_t counts[kNBuckets]{};  // 256 b
    for (const auto& it : values) {
      ++counts[(proj(it) >> shift) & 0xFF];
    }

    // exclusive scan
    [&bucket_indices, &counts]<std::size_t... Is>(std::index_sequence<Is...>) {
      std::uint64_t count = 0;
      (..., [&](std::size_t bucket_idx) {
        bucket_indices[bucket_idx] = count;
        count += counts[bucket_idx];
      }(Is));
    }(std::make_index_sequence<kNBuckets>{});

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

std::vector<Kmer> Minimize(
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    MinimizeConfig config);

std::vector<Index> ConstructIndices(
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> sequences,
    MinimizeConfig minimize_config,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool);

std::uint32_t CalculateKmerThreshold(std::vector<Index> indices,
                                     double frequency);

// Find matches between a pair of sequences.
// Minhash argument from configuration will only be applied on the lhs sequence.
std::vector<Match> MatchPairs(const std::unique_ptr<biosoup::NucleicAcid>& lhs,
                              const std::unique_ptr<biosoup::NucleicAcid>& rhs,
                              MinimizeConfig minimize_config);

std::vector<Match> MatchToIndex(
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    std::span<const Index> indices, MapToIndexConfig map_config,
    MinimizeConfig minimize_config, ChainConfig chain_config,
    std::vector<std::uint32_t>* filtered);

std::vector<MatchChain> FindChainMatches(
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    std::vector<Match>&& matches, ChainConfig config);

std::vector<biosoup::Overlap> ChainLCSKpp(
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    std::vector<Match>&& matches, ChainConfig config);

// Chain matches into overlaps.
std::vector<biosoup::Overlap> Chain(
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    std::vector<Match>&& matches, ChainConfig config);

// Chain matches into overlaps for ai worlkoad.
std::vector<OverlapAI> ChainAI(
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    std::vector<Match>&& matches, ChainConfig config);

std::vector<MatchChain> ChainOnIndex(
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    std::span<const Index> indices, MapToIndexConfig map_config,
    MinimizeConfig minimize_config, ChainConfig chain_config,
    std::vector<std::uint32_t>* filtered);

std::vector<biosoup::Overlap> MapToIndex(
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    std::span<const Index> indices, MapToIndexConfig map_config,
    MinimizeConfig minimize_config, ChainConfig chain_config,
    std::vector<std::uint32_t>* filtered);

std::vector<biosoup::Overlap> LCSKppToIndex(
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    std::span<const Index> indices, MapToIndexConfig map_config,
    MinimizeConfig minimize_config, ChainConfig chain_config,
    std::vector<std::uint32_t>* filtered);

std::vector<OverlapAI> MapToIndexAI(
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    std::span<const Index> indices, MapToIndexConfig map_config,
    MinimizeConfig minimize_config, ChainConfig chain_config,
    std::vector<std::uint32_t>* filtered);

}  // namespace ram

#endif  // RAM_ALGORITHM_HPP_
