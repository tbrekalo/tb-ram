#ifndef RAM_TYPES_HPP_
#define RAM_TYPES_HPP_

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "thread_pool/thread_pool.hpp"

namespace ram {

struct MinimizeConfig {
  std::uint32_t kmer_length = 15;
  std::uint32_t window_length = 5;
  bool minhash = false;
};

struct ChainConfig {
  std::uint32_t kmer_length = 15;
  std::uint32_t bandwidth = 500;
  std::uint32_t chain = 4;
  std::int32_t min_matches = 100;
  std::int64_t gap = 10'000;
  std::uint32_t n_chains = 5;
};

struct MapToIndexConfig {
  bool avoid_equal = true;
  bool avoid_symmetric = true;
  std::uint32_t occurrence = -1;
};

struct AlgoConfig {
  std::shared_ptr<thread_pool::ThreadPool> thread_pool;
  MinimizeConfig minimize_config;
  ChainConfig chain_config;
};

struct Kmer {
 public:
  constexpr Kmer() {}

  constexpr Kmer(std::uint64_t value, std::uint64_t origin)
      : value(value), origin(origin) {}

  constexpr std::uint32_t id() const noexcept {
    return static_cast<std::uint32_t>(origin >> 32);
  }

  constexpr std::uint32_t position() const noexcept {
    return static_cast<std::uint32_t>(origin) >> 1;
  }

  constexpr bool strand() const noexcept { return origin & 1; }

  std::uint64_t value;
  std::uint64_t origin;
};

inline constexpr std::uint64_t KmerValueProjection(const Kmer& kmer) noexcept {
  return kmer.value;
}

inline constexpr std::uint64_t KmerOriginProjection(const Kmer& kmer) noexcept {
  return kmer.origin;
}

struct Match {
 public:
  constexpr Match() {}

  constexpr Match(std::uint64_t group, std::uint64_t positions)
      : group(group), positions(positions) {}

  constexpr std::uint32_t rhs_id() const noexcept {
    return static_cast<std::uint32_t>(group >> 33);
  }

  constexpr bool strand() const noexcept { return (group >> 32) & 1; }

  constexpr std::uint32_t diagonal() const noexcept {
    return static_cast<std::uint32_t>(group);
  }

  [[gnu::always_inline]] constexpr std::uint32_t lhs_position() const noexcept {
    return static_cast<std::uint32_t>(positions >> 32);
  }

  [[gnu::always_inline]] constexpr std::uint32_t rhs_position() const noexcept {
    return static_cast<std::uint32_t>(positions);
  }

  std::uint64_t group;
  std::uint64_t positions;
};

inline constexpr std::uint64_t MatchGroupProjection(
    const Match& match) noexcept {
  return match.group;
}

inline constexpr std::uint64_t MatchPositionProjection(
    const Match& match) noexcept {
  return match.positions;
}

class Index {
 public:
  Index() = default;

  std::uint32_t Find(std::uint64_t key, const std::uint64_t** dst) const {
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

  struct Hash {
    std::size_t operator()(std::uint64_t key) const noexcept {
      return std::hash<std::uint64_t>()(key >> 1);
    }
  };
  struct KeyEqual {
    bool operator()(std::uint64_t lhs, std::uint64_t rhs) const noexcept {
      return (lhs >> 1) == (rhs >> 1);
    }
  };

  std::vector<std::uint64_t> origins;
  std::unordered_map<std::uint64_t, std::uint64_t, Hash, KeyEqual> locator;
};

struct MatchChain {
  std::vector<Match> matches;
  std::uint32_t lhs_matches;
  std::uint32_t rhs_matches;
  double score;
};

struct OverlapAI {
  std::uint32_t lhs_id;
  std::uint32_t lhs_begin;
  std::uint32_t lhs_end;
  std::uint32_t lhs_matches;

  bool strand;  // Watcon-Crick strand

  std::uint32_t rhs_id;
  std::uint32_t rhs_begin;
  std::uint32_t rhs_end;
  std::uint32_t rhs_matches;

  double score;
  double diff_mean;
  std::uint32_t q75;
  std::uint32_t q90;
  std::uint32_t q95;
  std::uint32_t q98;
};

}  // namespace ram

#endif  // RAM_TYPES_HPP_
