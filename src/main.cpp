// Copyright (c) 2020 Robert Vaser

#include <cstdlib>
#include <iostream>
#include <mutex>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/progress_bar.hpp"
#include "biosoup/timer.hpp"
#include "cxxopts.hpp"
#include "ram/algorithm.hpp"
#include "ram/io.hpp"
#include "tbb/tbb.h"

int main(int argc, char** argv) {
  cxxopts::Options options("ram", "sequence mapping tool");

  /* clang-format off */
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
  options.add_options("input")
    ("inputs", "target and query seuqnces",
      cxxopts::value<std::vector<std::string>>());
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

    auto minimize_cfg = ram::MinimizeConfig{
        .kmer_length = parsed_options["kmer-length"].as<std::uint32_t>(),
        .window_length = parsed_options["window-length"].as<std::uint32_t>(),
        .minhash = parsed_options["minhash"].as<bool>(),
    };

    auto chain_cfg = ram::ChainConfig{
        .kmer_length = parsed_options["kmer-length"].as<std::uint32_t>(),
        .bandwidth = parsed_options["bandwidth"].as<std::uint32_t>(),
        .chain = parsed_options["chain"].as<std::uint32_t>(),
        .min_matches = parsed_options["matches"].as<std::uint32_t>(),
        .gap = parsed_options["gap"].as<std::uint64_t>(),
    };

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
        auto indices = ram::ConstructIndices(targets, minimize_cfg);
        auto occurrence = ram::CalculateKmerThreshold(indices, frequency);

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

          biosoup::ProgressBar bar{static_cast<std::uint32_t>(sequences.size()),
                                   16};
          auto update_progress = [&timer, &bar, mtx = std::mutex{}] mutable {
            std::lock_guard lk{mtx};
            if (++bar) {
              std::cerr << "[ram::] mapped " << bar.event_counter()
                        << " sequences "
                        << "[" << bar << "] " << std::fixed << timer.Lap()
                        << "s"
                        << "\r";
            }
          };

          std::vector<std::vector<biosoup::Overlap>> overlaps(sequences.size());
          tbb::parallel_for(
              0uz, sequences.size(), [&](std::size_t idx) -> void {
                overlaps[idx] =
                    ram::MapSeqToIndex(sequences[idx], indices,
                                       ram::MapToIndexConfig{
                                           .avoid_equal = is_ava,
                                           .avoid_symmetric = is_ava,
                                           .occurrence = occurrence,
                                       },
                                       minimize_cfg, chain_cfg, nullptr);
                update_progress();
              });

          std::uint64_t rhs_offset = targets.front()->id;
          std::uint64_t lhs_offset = sequences.front()->id;
          for (auto& it : overlaps) {
            for (const auto& jt : it) {
              std::cout << sequences[jt.lhs_id - lhs_offset]->name << "\t"
                        << sequences[jt.lhs_id - lhs_offset]->inflated_len
                        << "\t" << jt.lhs_begin << "\t" << jt.lhs_end << "\t"
                        << (jt.strand ? "+" : "-") << "\t"
                        << targets[jt.rhs_id - rhs_offset]->name << "\t"
                        << targets[jt.rhs_id - rhs_offset]->inflated_len << "\t"
                        << jt.rhs_begin << "\t" << jt.rhs_end << "\t"
                        << jt.score << "\t"
                        << std::max(jt.lhs_end - jt.lhs_begin,
                                    jt.rhs_end - jt.rhs_begin)
                        << "\t" << 255 << std::endl;
            }
          }
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
