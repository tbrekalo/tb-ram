#ifndef RAM_TYPES_HPP_
#define RAM_TYPES_HPP_

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace ram {

struct Kmer {
 public:
  Kmer();
  Kmer(std::uint64_t value, std::uint64_t origin);

  std::uint32_t id() const noexcept;

  std::uint32_t position() const noexcept;

  bool strand() const noexcept;

  std::uint64_t value;
  std::uint64_t origin;
};

[[using gnu: always_inline, const, hot]] std::uint64_t KmerValueProjection(
    const Kmer& kmer) noexcept;

[[using gnu: always_inline, const, hot]] std::uint64_t KmerOriginProjection(
    const Kmer& kmer) noexcept;

struct Match {
 public:
  Match();
  Match(std::uint64_t group, std::uint64_t positions);

  [[gnu::always_inline]] std::uint32_t rhs_id() const noexcept;

  [[gnu::always_inline]] bool strand() const noexcept;

  [[gnu::always_inline]] std::uint32_t diagonal() const noexcept;

  [[gnu::always_inline]] std::uint32_t lhs_position() const noexcept;

  [[gnu::always_inline]] std::uint32_t rhs_position() const noexcept;

  std::uint64_t group;
  std::uint64_t positions;
};

[[using gnu: always_inline, const, hot]] std::uint64_t MatchGroupProjection(
    const Match& match) noexcept;

[[using gnu: always_inline, const, hot]] std::uint64_t MatchPositionProjection(
    const Match& match) noexcept;

class Index {
 public:
  Index();

  std::uint32_t Find(std::uint64_t key, const std::uint64_t** dst) const;

  struct Hash {
    std::size_t operator()(std::uint64_t key) const noexcept;
  };
  struct KeyEqual {
    bool operator()(std::uint64_t lhs, std::uint64_t rhs) const noexcept;
  };

  std::vector<std::uint64_t> origins;
  std::unordered_map<std::uint64_t, std::uint64_t, Hash, KeyEqual> locator;
};

}  // namespace ram

#endif  // RAM_TYPES_HPP_