// Copyright (c) 2020 Robert Vaser

#include <cstdlib>
#include <iostream>
#include <mutex>
#include <string_view>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/progress_bar.hpp"
#include "biosoup/timer.hpp"
#include "cxxopts.hpp"
#include "ram/algorithm.hpp"
#include "ram/io.hpp"
#include "tbb/tbb.h"

enum class Mode {
  kMatch,
  kOverlap,
};

std::ostream& operator<<(std::ostream& ostrm, Mode mode) {
  switch (mode) {
    case Mode::kMatch:
      return ostrm << "match";
    case Mode::kOverlap:
      return ostrm << "overlap";
  }
}

std::istream& operator>>(std::istream& istrm, Mode& mode) {
  using namespace std::literals;
  std::string repr;
  istrm >> repr;

  if (repr == "match"sv) {
    mode = Mode::kMatch;
    return istrm;
  }

  if (repr == "overlap"sv) {
    mode = Mode::kOverlap;
    return istrm;
  }

  istrm.setstate(std::ios::failbit);
  return istrm;
}

struct BatchContext {
  ram::AlgoConfig algo_cfg;
  ram::MapToIndexConfig map2index_cfg;
  std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets;
  std::span<const std::unique_ptr<biosoup::NucleicAcid>> sequences;
  std::span<ram::Index> indices;
  std::function<void()> update_progress;
};

template <Mode mode>
static const auto ExecuteBatchImpl = [](BatchContext ctx) -> void {
  auto [operation, prototype, print] = [] {
    if constexpr (mode == Mode::kMatch) {
      return std::tuple{&ram::MatchToIndex, std::vector<ram::Match>{},
                        &ram::PrintMatchBatch};
    } else {
      return std::tuple{&ram::MapToIndex, std::vector<biosoup::Overlap>{},
                        &ram::PrintOverlapBatch};
    }
  }();

  std::vector<decltype(prototype)> values(ctx.sequences.size());
  tbb::parallel_for(0uz, ctx.sequences.size(), [&](std::size_t idx) -> void {
    values[idx] = operation(ctx.sequences[idx], ctx.indices, ctx.map2index_cfg,
                            ctx.algo_cfg.minimize, ctx.algo_cfg.chain, nullptr);
    ctx.update_progress();
  });

  print(std::cout, ctx.targets, ctx.sequences, values);
};

static const auto ExecuteBatch = [](Mode mode, BatchContext ctx) -> void {
  switch (mode) {
    case Mode::kMatch:
      ExecuteBatchImpl<Mode::kMatch>(ctx);
      break;
    case Mode::kOverlap:
      ExecuteBatchImpl<Mode::kOverlap>(ctx);
      break;
  }
};

int main(int argc, char** argv) {
  cxxopts::Options options("ram", "sequence mapping tool");

  /* clang-format off */
  options.add_options()
    ("inputs", "target and query seuqnces",
      cxxopts::value<std::vector<std::string>>());
  options.add_options("mode")
    ("mode",
      "ram operating mode (match|overlap)",
      cxxopts::value<Mode>()->default_value("overlap"));
  options.add_options("algorithm")
    ("k,kmer-length",
      "length of minimizers",
      cxxopts::value<std::uint32_t>()->default_value("15"))
    ("w,window-length",
      "length of sliding window from which minimizers are sampled",
      cxxopts::value<std::uint32_t>()->default_value("5"))
    ("f,frequency-threshold",
      "threshold for ignoring most frequent minimizers",
      cxxopts::value<double>()->default_value("0.001"))
    ("bandwidth",
      "size of bandwidth in which minimizer hits can be chained",
      cxxopts::value<std::uint32_t>()->default_value("500"))
    ("chain",
      "minimal number of chained minimizer hits in overlap",
      cxxopts::value<std::uint32_t>()->default_value("4"))
    ("matches",
      "minimal number of matching bases in overlap",
      cxxopts::value<std::uint32_t>()->default_value("100"))
    ("gap",
      "maximal gap between minimizer hits in a chain",
      cxxopts::value<std::uint64_t>()->default_value("10000"))
    ("minhash",
      "use only a portion of all minimizers")
    ("t,threads",
      "number of threads",
      cxxopts::value<std::uint32_t>()->default_value("1"));
  options.add_options("info")
    ("v,version", "print version and exit early")
    ("h,help", "print help and exit early");
  options.positional_help("<target> [<query>]");
  /* clang-format on */

  try {
    auto early_quit = false;
    options.parse_positional({"inputs"});
    auto parsed_options = options.parse(argc, argv);
    if (parsed_options.count("version")) {
      std::cerr << VERSION << std::endl;
      early_quit = true;
    }

    if (parsed_options.count("help")) {
      std::cerr << options.help() << std::endl;
      early_quit = true;
    }

    if (early_quit) {
      return 0;
    }

    auto input_paths = parsed_options["inputs"].as<std::vector<std::string>>();
    auto tparser = ram::CreateParser(input_paths[0]);

    auto is_ava = input_paths.size() == 1uz;
    std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> sparser = nullptr;
    if (input_paths.size() > 1) {
      sparser = ram::CreateParser(input_paths[1]);
      if (sparser == nullptr) {
        return 1;
      }
      is_ava = input_paths[0] == input_paths[1];
    } else {
      sparser = ram::CreateParser(input_paths[0]);
      is_ava = true;
    }

    auto num_threads = parsed_options["threads"].as<std::uint32_t>();
    auto frequency = parsed_options["frequency-threshold"].as<double>();

    auto cfg = ram::AlgoConfig{
        .minimize =
            ram::MinimizeConfig{
                .kmer_length =
                    parsed_options["kmer-length"].as<std::uint32_t>(),
                .window_length =
                    parsed_options["window-length"].as<std::uint32_t>(),
                .minhash = parsed_options["minhash"].as<bool>(),
            },
        .chain = ram::ChainConfig{
            .kmer_length = parsed_options["kmer-length"].as<std::uint32_t>(),
            .bandwidth = parsed_options["bandwidth"].as<std::uint32_t>(),
            .chain = parsed_options["chain"].as<std::uint32_t>(),
            .min_matches = parsed_options["matches"].as<std::uint32_t>(),
            .gap = parsed_options["gap"].as<std::uint64_t>(),
        }};

    auto arena = tbb::task_arena(num_threads);
    biosoup::Timer timer{};

    arena.execute([&] {
      while (true) {
        timer.Start();
        std::vector<std::unique_ptr<biosoup::NucleicAcid>> targets;
        targets = tparser->Parse(1ULL << 32);
        if (targets.empty()) {
          break;
        }

        std::cerr << "[ram::] parsed " << targets.size() << " targets "
                  << std::fixed << timer.Stop() << "s" << std::endl;

        timer.Start();
        auto indices = ram::ConstructIndices(targets, cfg.minimize);
        auto map_to_index_cfg = ram::MapToIndexConfig{
            .avoid_equal = is_ava,
            .avoid_symmetric = is_ava,
            .occurrence = ram::CalculateKmerThreshold(indices, frequency),
        };

        std::cerr << "[ram::] minimized targets " << std::fixed << timer.Stop()
                  << "s" << std::endl;

        std::uint64_t num_targets = biosoup::NucleicAcid::num_objects;
        biosoup::NucleicAcid::num_objects = 0;
        while (true) {
          timer.Start();
          std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
          sequences = sparser->Parse(1U << 30);
          if (sequences.empty()) {
            break;
          }

          std::mutex bar_mtx;
          biosoup::ProgressBar bar{static_cast<std::uint32_t>(sequences.size()),
                                   16};
          ExecuteBatchImpl<Mode::kOverlap>(
              {cfg, map_to_index_cfg, targets, sequences, indices,
               [&timer, &bar, &bar_mtx, n = sequences.size()] {
                 std::lock_guard lk{bar_mtx};
                 if (++bar) {
                   std::cerr << "[ram::] batch progress "
                             << 100. * bar.event_counter() / n << "% [" << bar
                             << "] " << std::fixed << timer.Lap() << "s"
                             << "\r";
                 }
               }});
          std::cerr << std::endl;
          timer.Stop();

          if (is_ava && biosoup::NucleicAcid::num_objects >= num_targets) {
            break;
          }
        }

        sparser->Reset();
        biosoup::NucleicAcid::num_objects = num_targets;
      }
    });

    std::cerr << "[ram::] " << timer.elapsed_time() << "s" << std::endl;
  } catch (const std::exception& exception) {
    std::cerr << exception.what() << std::endl;
    return 1;
  }

  return 0;
}
