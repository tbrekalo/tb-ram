#include "ram/io.hpp"

#include <exception>
#include <format>
#include <sstream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

static constexpr auto kFastaSuffxies =
    std::array<char const*, 4>{".fasta", "fasta.gz", ".fa", ".fa.gz"};

static constexpr auto kFastqSuffixes =
    std::array<char const*, 4>{".fastq", ".fastq.gz", ".fq", ".fq.gz"};

static auto IsSuffixFor(std::string_view const suffix,
                        std::string_view const query) -> bool {
  return suffix.length() <= query.length()
             ? suffix == query.substr(query.length() - suffix.length())
             : false;
}

namespace ram {

std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> CreateParser(
    const std::string& path) {
  using namespace std::placeholders;
  if (std::any_of(kFastaSuffxies.cbegin(), kFastaSuffxies.cend(),
                  std::bind(IsSuffixFor, _1, path.c_str()))) {
    return bioparser::Parser<biosoup::NucleicAcid>::Create<
        bioparser::FastaParser>(path);
  }
  if (std::any_of(kFastqSuffixes.cbegin(), kFastqSuffixes.cend(),
                  std::bind(IsSuffixFor, _1, path.c_str()))) {
    return bioparser::Parser<biosoup::NucleicAcid>::Create<
        bioparser::FastqParser>(path);
  }

  throw std::runtime_error([path] {
    auto ostrm = std::ostringstream{};
    ostrm << "[ram::CreateParser] error: file " << path
          << " has unsupported format extension (valid extensions: .fasta, "
          << ".fasta.gz, .fna, .fna.gz, .fa, .fa.gz, .fastq, .fastq.gz, "
          << ".fq, .fq.gz)";
    return ostrm.str();
  }());
}

void PrintOverlapBatch(
    std::ostream& ostrm,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> queries,
    std::span<const std::vector<biosoup::Overlap>> overlaps) {
  std::uint64_t rhs_offset = targets.front()->id;
  std::uint64_t lhs_offset = queries.front()->id;
  for (const auto& it : overlaps) {
    for (const auto& jt : it) {
      /* clang-format off */
      ostrm << std::format(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        queries[jt.lhs_id - lhs_offset]->name,
        queries[jt.lhs_id - lhs_offset]->inflated_len,
        jt.lhs_begin,
        jt.lhs_end,

        (jt.strand ? '+' : '-'),

        targets[jt.rhs_id - rhs_offset]->name,
        targets[jt.rhs_id - rhs_offset]->inflated_len,
        jt.rhs_begin,
        jt.rhs_end,

        jt.score,
        std::max(jt.lhs_end - jt.lhs_begin, jt.rhs_end - jt.rhs_begin),
        255
      );
      /* clang-format on */
    }
  }

  std::flush(ostrm);
}

void PrintOverlapAIBatch(
    std::ostream& ostrm,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> queries,
    std::span<const std::vector<OverlapAI>> overlaps) {
  std::uint64_t rhs_offset = targets.front()->id;
  std::uint64_t lhs_offset = queries.front()->id;
  for (const auto& it : overlaps) {
    for (const auto& jt : it) {
      /* clang-format off */
      ostrm << std::format(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        queries[jt.lhs_id - lhs_offset]->name,
        queries[jt.lhs_id - lhs_offset]->inflated_len,
        jt.lhs_begin,
        jt.lhs_end,
        jt.lhs_matches,

        (jt.strand ? '+' : '-'),

        targets[jt.rhs_id - rhs_offset]->name,
        targets[jt.rhs_id - rhs_offset]->inflated_len,
        jt.rhs_begin,
        jt.rhs_end,
        jt.rhs_matches,

        jt.score,
        jt.diff_mean,
        jt.q75,
        jt.q90,
        jt.q95,
        jt.q98
      );
      /* clang-format on */
    }
  }

  std::flush(ostrm);
}

void PrintMatchBatch(
    std::ostream& ostrm,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> queries,
    std::span<const std::vector<Match>> matches) {
  std::uint64_t rhs_offset = targets.front()->id;
  for (auto lhs_idx = 0uz; lhs_idx < matches.size(); ++lhs_idx) {
    /* clang-format off */
    for (const auto& match : matches[lhs_idx]) {
      ostrm << std::format(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        queries[lhs_idx]->name,
        queries[lhs_idx]->inflated_len,
        match.lhs_position(),

        (match.strand() ? '+' : '-'),

        targets[match.rhs_id() - rhs_offset]->name,
        targets[match.rhs_id() - rhs_offset]->inflated_len,
        match.rhs_position()
      );
    }
    /* clang-format on */
  }

  std::flush(ostrm);
}

void PrintMatchChainBatch(
    std::ostream& ostrm,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> queries,
    std::span<const std::vector<MatchChain>> match_chains) {
  std::uint64_t chain_id = 0;
  std::uint64_t rhs_offset = targets.front()->id;
  for (auto lhs_idx = 0uz; lhs_idx < match_chains.size(); ++lhs_idx) {
    for (const auto& match_chain : match_chains[lhs_idx]) {
      for (const auto& match : match_chain.matches) {
        /* clang-format off */
        ostrm << std::format(
          "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
          queries[lhs_idx]->name,
          queries[lhs_idx]->inflated_len,
          match.lhs_position(),

          (match.strand() ? '+' : '-'),

          targets[match.rhs_id() - rhs_offset]->name,
          targets[match.rhs_id() - rhs_offset]->inflated_len,
          match.rhs_position(),

          match_chain.lhs_matches,
          match_chain.rhs_matches,
          chain_id
        );
        /* clang-format on */
      }
      ++chain_id;
    }
  }

  std::flush(ostrm);
}

}  // namespace ram
