// Copyright (c) 2020 Robet Vaser

#ifndef RAM_MINIMIZER_ENGINE_HPP_
#define RAM_MINIMIZER_ENGINE_HPP_

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "ram/types.hpp"
#include "thread_pool/thread_pool.hpp"

namespace ram {

class MinimizerEngine {
 public:
  MinimizerEngine(
      std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr,
      std::uint32_t k = 15,  // element of [1, 31]
      std::uint32_t w = 5, std::uint32_t bandwidth = 500,
      std::uint32_t chain = 4, std::uint32_t matches = 100,
      std::uint32_t gap = 10000);

  MinimizerEngine(const MinimizerEngine&) = delete;
  MinimizerEngine& operator=(const MinimizerEngine&) = delete;

  MinimizerEngine(MinimizerEngine&&) = default;
  MinimizerEngine& operator=(MinimizerEngine&&) = default;

  ~MinimizerEngine() = default;

  // transform set of sequences to minimizer index
  // minhash = pick only the smallest sequence->data.size() / k minimizers
  void Minimize(
      std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator first,
      std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator last,
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
  class Index {
   public:
    Index() = default;

    std::uint32_t Find(std::uint64_t key, const std::uint64_t** dst) const;

    struct Hash {
      std::size_t operator()(std::uint64_t key) const {
        return std::hash<std::uint64_t>()(key >> 1);
      }
    };
    struct KeyEqual {
      bool operator()(std::uint64_t lhs, std::uint64_t rhs) const {
        return (lhs >> 1) == (rhs >> 1);
      }
    };

    std::vector<std::uint64_t> origins;
    std::unordered_map<std::uint64_t, std::uint64_t, Hash, KeyEqual> locator;
  };

  std::uint32_t k_;
  std::uint32_t w_;
  std::uint32_t bandwidth_;
  std::uint32_t chain_;
  std::uint32_t matches_;
  std::uint64_t gap_;
  std::uint32_t occurrence_;
  std::vector<Index> index_;
  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
};

}  // namespace ram

#endif  // RAM_MINIMIZER_ENGINE_HPP_
