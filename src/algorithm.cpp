#include "ram/algorithm.hpp"

#include <deque>

#include "biosoup/nucleic_acid.hpp"

namespace ram {

std::vector<Kmer> Minimize(
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    MinimizeConfig config) {
  if (sequence->inflated_len < config.kmer_length) {
    return std::vector<Kmer>{};
  }

  std::uint64_t mask = (1ULL << (config.kmer_length * 2)) - 1;

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

  std::deque<Kmer> window;
  auto window_add = [&](std::uint64_t value, std::uint64_t location) -> void {
    while (!window.empty() && window.back().value > value) {
      window.pop_back();
    }
    window.emplace_back(value, location);
  };
  auto window_update = [&](std::uint32_t position) -> void {
    while (!window.empty() && (window.front().position()) < position) {
      window.pop_front();
    }
  };

  std::uint64_t shift = (config.kmer_length - 1) * 2;
  std::uint64_t minimizer = 0;
  std::uint64_t reverse_minimizer = 0;
  std::uint64_t id = static_cast<std::uint64_t>(sequence->id) << 32;
  std::uint64_t is_stored = 1ULL << 63;

  std::vector<Kmer> dst;

  for (std::uint32_t i = 0; i < sequence->inflated_len; ++i) {
    std::uint64_t c = sequence->Code(i);
    minimizer = ((minimizer << 2) | c) & mask;
    reverse_minimizer = (reverse_minimizer >> 2) | ((c ^ 3) << shift);
    if (i >= config.kmer_length - 1U) {
      if (minimizer < reverse_minimizer) {
        window_add(hash(minimizer), (i - (config.kmer_length - 1U)) << 1 | 0);
      } else if (minimizer > reverse_minimizer) {
        window_add(hash(reverse_minimizer),
                   (i - (config.kmer_length - 1U)) << 1 | 1);
      }
    }
    if (i >= (config.kmer_length - 1U) + (config.window_length - 1U)) {
      for (auto it = window.begin(); it != window.end(); ++it) {
        if (it->value != window.front().value) {
          break;
        }
        if (it->origin & is_stored) {
          continue;
        }
        dst.emplace_back(it->value, id | it->origin);
        it->origin |= is_stored;
      }
      window_update(i - (config.kmer_length - 1U) -
                    (config.window_length - 1U) + 1);
    }
  }

  if (config.minhash) {
    RadixSort(std::span<Kmer>(dst), config.kmer_length * 2,
              KmerValueProjection);
    dst.resize(sequence->inflated_len / config.kmer_length);
    RadixSort(std::span<Kmer>(dst), 64, KmerOriginProjection);
  }

  return dst;
}

}  // namespace ram
