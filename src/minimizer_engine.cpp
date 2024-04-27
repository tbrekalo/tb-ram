#include "ram/minimizer_engine.hpp"

#include "ram/algorithm.hpp"
#include "ram/types.hpp"
#include "thread_pool/thread_pool.hpp"

namespace ram {

struct ram::MinimizerEngine::Impl {
  AlgoConfig algo_config;
  std::vector<Index> indices;
  std::uint32_t occurences;
};

MinimizerEngine::MinimizerEngine(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, std::uint32_t k,
    std::uint32_t w, std::uint32_t bandwidth, std::uint32_t chain,
    std::uint32_t matches, std::uint32_t gap)
    : pimpl_(std::unique_ptr<Impl>(new Impl{
          .algo_config =
              AlgoConfig{.thread_pool = std::move(thread_pool),
                         .minimize_config =
                             MinimizeConfig{
                                 .kmer_length = k,
                                 .window_length = w,
                             },
                         .chain_config = ChainConfig{.kmer_length = k,
                                                     .bandwidth = bandwidth,
                                                     .chain = chain,
                                                     .min_matches = matches,
                                                     .gap = gap}},
          .indices = {},
          .occurences = 0})) {}

void MinimizerEngine::Minimize(
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    bool minhash) {
  auto minimize_cfg = MinimizeConfig{
      .kmer_length = pimpl_->algo_config.minimize_config.kmer_length,
      .window_length = pimpl_->algo_config.minimize_config.window_length,
      .minhash = minhash,
  };
  pimpl_->indices =
      ConstructIndices(targets, minimize_cfg, pimpl_->algo_config.thread_pool);
}

void MinimizerEngine::Filter(double frequency) {
  pimpl_->occurences = CalculateKmerThreshold(pimpl_->indices, frequency);
}

std::vector<biosoup::Overlap> MinimizerEngine::Map(
    const std::unique_ptr<biosoup::NucleicAcid>& sequence, bool avoid_equal,
    bool avoid_symmetric, bool minhash,
    std::vector<std::uint32_t>* filtered) const {
  auto map_to_index_cfg = MapToIndexConfig{.avoid_equal = avoid_equal,
                                           .avoid_symmetric = avoid_symmetric,
                                           .occurrence = pimpl_->occurences};

  auto minimize_cfg = MinimizeConfig{
      .kmer_length = pimpl_->algo_config.minimize_config.kmer_length,
      .window_length = pimpl_->algo_config.minimize_config.window_length,
      .minhash = minhash};

  return MapToIndex(sequence, pimpl_->indices, map_to_index_cfg, minimize_cfg,
                    pimpl_->algo_config.chain_config, filtered);
}

std::vector<biosoup::Overlap> MinimizerEngine::Map(
    const std::unique_ptr<biosoup::NucleicAcid>& lhs,
    const std::unique_ptr<biosoup::NucleicAcid>& rhs, bool minhash) const {
  return MapPairs(
      lhs, rhs,
      MinimizeConfig{
          .kmer_length = pimpl_->algo_config.minimize_config.kmer_length,
          .window_length = pimpl_->algo_config.minimize_config.window_length,
          .minhash = minhash},
      pimpl_->algo_config.chain_config);
}

}  // namespace ram
