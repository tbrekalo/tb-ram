#include "ram/types.hpp"

namespace ram {

Kmer::Kmer() {}

std::uint32_t Kmer::id() const noexcept {
  return static_cast<std::uint32_t>(origin >> 32);
}

std::uint32_t Kmer::position() const noexcept {
  return static_cast<std::uint32_t>(origin) >> 1;
}

bool Kmer::strand() const noexcept { return origin & 1; }

Kmer::Kmer(std::uint64_t value, std::uint64_t origin)
    : value(value), origin(origin) {}

std::uint64_t KmerValueProjection(const Kmer& kmer) noexcept {
  return kmer.value;
}

std::uint64_t KmerOriginProjection(const Kmer& kmer) noexcept {
  return kmer.origin;
}

Match::Match() {}

Match::Match(std::uint64_t group, std::uint64_t positions)
    : group(group), positions(positions) {}

std::uint32_t Match::rhs_id() const noexcept {
  return static_cast<std::uint32_t>(group >> 33);
}

bool Match::strand() const noexcept { return (group >> 32) & 1; }

std::uint32_t Match::diagonal() const noexcept {
  return static_cast<std::uint32_t>(group);
}

std::uint32_t Match::lhs_position() const noexcept {
  return static_cast<std::uint32_t>(positions >> 32);
}

std::uint32_t Match::rhs_position() const noexcept {
  return static_cast<std::uint32_t>(positions);
}

std::uint64_t MatchGroupProjection(const Match& match) noexcept {
  return match.group;
}

std::uint64_t MatchPositionProjection(const Match& match) noexcept {
  return match.positions;
}

Index::Index() {}

std::size_t Index::Hash::operator()(std::uint64_t key) const noexcept {
  return std::hash<std::uint64_t>()(key >> 1);
}

bool Index::KeyEqual::operator()(std::uint64_t lhs,
                                 std::uint64_t rhs) const noexcept {
  return (lhs >> 1) == (rhs >> 1);
}

std::uint32_t Index::Find(std::uint64_t key, const std::uint64_t** dst) const {
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

}  // namespace ram
