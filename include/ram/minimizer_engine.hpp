#ifndef RAM_MINIMIZER_ENGINE_HPP_
#define RAM_MINIMIZER_ENGINE_HPP_

#include <span>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "thread_pool/thread_pool.hpp"

namespace ram {

class MinimizerEngine {
 public:
  MinimizerEngine(
      std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr,
      std::uint32_t k = 15,  // element of [1, 31]
      std::uint32_t w = 5, std::uint32_t bandwidth = 500, std::uint32_t chain = 4,
      std::int32_t matches = 100, std::int32_t gap = 10000);

  MinimizerEngine(const MinimizerEngine&) = delete;
  MinimizerEngine& operator=(const MinimizerEngine&) = delete;

  MinimizerEngine(MinimizerEngine&&) = default;
  MinimizerEngine& operator=(MinimizerEngine&&) = default;

  ~MinimizerEngine();

  // transform set of sequences to minimizer index
  // minhash = pick only the smallest sequence->data.size() / k minimizers
  void Minimize(std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
                bool minhash = false);

  // set occurrence frequency threshold
  void Filter(double frequency);

  // find overlaps in preconstructed minimizer index
  std::vector<biosoup::Overlap> Map(
      const std::unique_ptr<biosoup::NucleicAcid>& sequence,
      bool avoid_equal,      // ignore overlaps in which lhs_id == rhs_id
      bool avoid_symmetric,  // ignore overlaps in which lhs_id > rhs_id
      bool minhash = false,  // only lhs
      std::vector<std::uint32_t>* filtered = nullptr) const;

  // find overlaps between a pair of sequences
  std::vector<biosoup::Overlap> Map(
      const std::unique_ptr<biosoup::NucleicAcid>& lhs,
      const std::unique_ptr<biosoup::NucleicAcid>& rhs,
      bool minhash = false) const;  // only lhs

 private:
  struct Impl;
  std::unique_ptr<Impl> pimpl_;
};

}  // namespace ram

#endif  // RAM_MINIMIZER_ENGINE_HPP_
