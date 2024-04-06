#include "ram/io.hpp"

#include <exception>
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
      ostrm << queries[jt.lhs_id - lhs_offset]->name
            << '\t' << queries[jt.lhs_id - lhs_offset]->inflated_len
            << '\t' << jt.lhs_begin
            << '\t' << jt.lhs_end

            << '\t' << (jt.strand ? '+' : '-')

            << "\t" << targets[jt.rhs_id - rhs_offset]->name
            << '\t' << targets[jt.rhs_id - rhs_offset]->inflated_len
            << '\t' << jt.rhs_begin
            << '\t' << jt.rhs_end

            << '\t' << jt.score
            << '\t' << std::max(jt.lhs_end - jt.lhs_begin, jt.rhs_end - jt.rhs_begin)
            << '\t' << 255
            << '\n';
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
      ostrm << queries[jt.lhs_id - lhs_offset]->name
            << '\t' << queries[jt.lhs_id - lhs_offset]->inflated_len
            << '\t' << jt.lhs_begin
            << '\t' << jt.lhs_end
            << '\t' << jt.lhs_matches

            << '\t' << (jt.strand ? '+' : '-')

            << "\t" << targets[jt.rhs_id - rhs_offset]->name
            << '\t' << targets[jt.rhs_id - rhs_offset]->inflated_len
            << '\t' << jt.rhs_begin
            << '\t' << jt.rhs_end
            << '\t' << jt.rhs_matches

            << '\t' << jt.diff_mean
            << '\t' << jt.q75
            << '\t' << jt.q90
            << '\t' << jt.q95
            << '\t' << jt.q98
            << '\n';
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
      ostrm << queries[lhs_idx]->name
            << '\t' << queries[lhs_idx]->inflated_len
            << '\t' << match.lhs_position()

            << '\t' << (match.strand() ? '+' : '-')

            << '\t' << targets[match.rhs_id() - rhs_offset]->name
            << '\t' << targets[match.rhs_id() - rhs_offset]->inflated_len
            << '\t' << match.rhs_position()
            << '\n';
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
        ostrm << queries[lhs_idx]->name
              << '\t' << queries[lhs_idx]->inflated_len
              << '\t' << match.lhs_position()

              << '\t' << (match.strand() ? '+' : '-')

              << '\t' << targets[match.rhs_id() - rhs_offset]->name
              << '\t' << targets[match.rhs_id() - rhs_offset]->inflated_len
              << '\t' << match.rhs_position()

              << '\t' << match_chain.lhs_matches
              << '\t' << match_chain.rhs_matches
              << '\t' << chain_id
              << '\n';
        /* clang-format on */
      }
      ++chain_id;
    }
  }

  std::flush(ostrm);
}

}  // namespace ram
