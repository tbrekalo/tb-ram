// Copyright (c) 2020 Robert Vaser

#include <cstdlib>
#include <format>
#include <iostream>
#include <mutex>
#include <ranges>
#include <string_view>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/progress_bar.hpp"
#include "cxxopts.hpp"
#include "glog/logging.h"
#include "ram/algorithm.hpp"
#include "ram/io.hpp"
#include "thread_pool/thread_pool.hpp"

enum class Mode {
  kChain,
  kMatch,
  kOverlap,
  kOverlapAI,
};

std::ostream& operator<<(std::ostream& ostrm, Mode mode) {
  switch (mode) {
    case Mode::kChain:
      return ostrm << "chain";
    case Mode::kMatch:
      return ostrm << "match";
    case Mode::kOverlap:
      return ostrm << "overlap";
    case Mode::kOverlapAI:
      return ostrm << "overlap-ai";
  }
}

std::istream& operator>>(std::istream& istrm, Mode& mode) {
  using namespace std::literals;
  std::string repr;
  istrm >> repr;

  if (repr == "chain"sv) {
    mode = Mode::kChain;
    return istrm;
  }

  if (repr == "match"sv) {
    mode = Mode::kMatch;
    return istrm;
  }

  if (repr == "overlap"sv) {
    mode = Mode::kOverlap;
    return istrm;
  }

  if (repr == "overlap-ai") {
    mode = Mode::kOverlapAI;
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

static auto PrintChainBatch =
    [chain_idx = 0uz](
        std::ostream& ostrm,
        std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
        std::span<const std::unique_ptr<biosoup::NucleicAcid>> queries,
        std::span<const std::vector<ram::MatchChain>> match_chains) mutable
    -> void {
  ram::PrintMatchChainBatch(ostrm, targets, queries, match_chains);
  auto foo = std::ranges::fold_left(
      match_chains |
          std::ranges::views::transform(
              [](const std::vector<ram::MatchChain>& vec) -> std::size_t {
                return vec.size();
              }),
      chain_idx, std::plus<std::size_t>{});
};

template <Mode mode>
static const auto ExecuteBatchImpl = [](BatchContext ctx) -> void {
  auto [operation, prototype, print] = [] {
    if constexpr (mode == Mode::kChain) {
      return std::tuple{&ram::ChainOnIndex, std::vector<ram::MatchChain>{},
                        PrintChainBatch};
    } else if constexpr (mode == Mode::kMatch) {
      return std::tuple{&ram::MatchToIndex, std::vector<ram::Match>{},
                        &ram::PrintMatchBatch};
    } else if constexpr (mode == Mode::kOverlap) {
      return std::tuple{&ram::MapToIndex, std::vector<biosoup::Overlap>{},
                        &ram::PrintOverlapBatch};
    } else {
      return std::tuple{&ram::MapToIndexAI, std::vector<ram::OverlapAI>{},
                        &ram::PrintOverlapAIBatch};
    }
  }();

  std::vector<decltype(prototype)> values(ctx.sequences.size());

  std::vector<std::future<void>> futures;
  for (auto idx = 0uz; idx < ctx.sequences.size(); ++idx) {
    futures.push_back(ctx.algo_cfg.thread_pool->Submit(
        [&](std::size_t jdx) -> void {
          values[jdx] =
              operation(ctx.sequences[jdx], ctx.indices, ctx.map2index_cfg,
                        ctx.algo_cfg.minimize, ctx.algo_cfg.chain, nullptr);
          ctx.update_progress();
        },
        idx));
  }

  for (auto& future : futures) {
    future.wait();
  }

  print(std::cout, ctx.targets, ctx.sequences, values);
};

static const auto ExecuteBatch = [](Mode mode, BatchContext ctx) -> void {
  switch (mode) {
    case Mode::kChain:
      ExecuteBatchImpl<Mode::kChain>(ctx);
      break;
    case Mode::kMatch:
      ExecuteBatchImpl<Mode::kMatch>(ctx);
      break;
    case Mode::kOverlap:
      ExecuteBatchImpl<Mode::kOverlap>(ctx);
      break;
    case Mode::kOverlapAI:
      ExecuteBatchImpl<Mode::kOverlapAI>(ctx);
      break;
  }
};

int main(int argc, char** argv) {
  /* clang-format on */
  FLAGS_logtostderr = true;
  google::InitGoogleLogging(argv[0]);
  google::InstallPrefixFormatter(
      +[](std::ostream& ostrm, const google::LogMessage& m, void*) -> void {
        ostrm << std::format(
            "datetime={:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}.{:06d} "
            "level={} "
            "threadid={} "
            "loc={}:{}",

            m.time().year(), m.time().month(), m.time().day(), m.time().hour(),
            m.time().min(), m.time().sec(), m.time().usec(),

            google::GetLogSeverityName(m.severity())[0],

            m.thread_id(),

            m.basename(), m.line());
      });
  /* clang-format on */
  cxxopts::Options options("ram", "sequence mapping tool");

  /* clang-format off */
  options.add_options()
    ("inputs", "target and query seuqnces",
      cxxopts::value<std::vector<std::string>>());
  options.add_options("io")
    ("mode",
      "ram operating mode (chain|match|overlap)",
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
    ("query-batch-size",
     "maximum query batch size in mebibites",
     cxxopts::value<std::uint64_t>()->default_value("1024"))
    ("target-batch-size",
      "maximum target batch size in mebibytes",
      cxxopts::value<std::uint64_t>()->default_value("4096"))
    ("t,threads",
      "number of threads",
      cxxopts::value<std::uint64_t>()->default_value("1"));
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
      google::ShutdownGoogleLogging();
      return 0;
    }

    auto input_paths = parsed_options["inputs"].as<std::vector<std::string>>();
    auto tparser = ram::CreateParser(input_paths[0]);

    auto is_ava = input_paths.size() == 1uz;
    std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> sparser = nullptr;
    if (input_paths.size() > 1) {
      sparser = ram::CreateParser(input_paths[1]);
      if (sparser == nullptr) {
        google::ShutdownGoogleLogging();
        return 1;
      }
      is_ava = input_paths[0] == input_paths[1];
    } else {
      sparser = ram::CreateParser(input_paths[0]);
      is_ava = true;
    }

    auto num_threads = parsed_options["threads"].as<std::uint64_t>();
    auto frequency = parsed_options["frequency-threshold"].as<double>();

    auto cfg = ram::AlgoConfig{
        .thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads),
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

    const auto mode = parsed_options["mode"].as<Mode>();
    const auto query_batch_size =
        parsed_options["query-batch-size"].as<std::uint64_t>() *
        (1ULL << 20ULL);
    const auto target_batch_size =
        parsed_options["target-batch-size"].as<std::uint64_t>() *
        (1ULL << 20ULL);

    std::ios_base::sync_with_stdio(false);
    while (true) {
      std::vector<std::unique_ptr<biosoup::NucleicAcid>> targets;
      targets = tparser->Parse(target_batch_size);
      if (targets.empty()) {
        break;
      }

      LOG(INFO) << std::format("event=parsed-targets value={}", targets.size());
      auto indices =
          ram::ConstructIndices(targets, cfg.minimize, cfg.thread_pool);
      auto map_to_index_cfg = ram::MapToIndexConfig{
          .avoid_equal = is_ava,
          .avoid_symmetric = is_ava,
          .occurrence = ram::CalculateKmerThreshold(indices, frequency),
      };

      LOG(INFO) << std::format("event=minimized-targets value={}",
                               targets.size());

      std::uint64_t num_targets = biosoup::NucleicAcid::num_objects;
      biosoup::NucleicAcid::num_objects = 0;
      while (true) {
        std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
        sequences = sparser->Parse(query_batch_size);
        if (sequences.empty()) {
          break;
        }

        std::mutex bar_mtx;
        biosoup::ProgressBar bar{static_cast<std::uint32_t>(sequences.size()),
                                 16};
        ExecuteBatch(mode, {cfg, map_to_index_cfg, targets, sequences, indices,
                            [&bar, &bar_mtx, n = sequences.size()] {
                              std::lock_guard lk{bar_mtx};
                              if (++bar) {
                                LOG(INFO) << std::format(
                                    "event=batch-progress value={}/{}",
                                    bar.event_counter(), n);
                              }
                            }});

        if (is_ava && biosoup::NucleicAcid::num_objects >= num_targets) {
          break;
        }
      }

      sparser->Reset();
      biosoup::NucleicAcid::num_objects = num_targets;
    }
  } catch (const std::exception& exception) {
    std::cerr << exception.what() << std::endl;
    google::ShutdownGoogleLogging();
    return 1;
  }

  google::ShutdownGoogleLogging();
  return 0;
}
