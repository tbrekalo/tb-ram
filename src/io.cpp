#include "ram/io.hpp"

#include <exception>
#include <sstream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace ram {
std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> CreateParser(
    const std::string& path) {
  auto is_suffix = [](const std::string& s, const std::string& suff) {
    return s.size() < suff.size()
               ? false
               : s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta") || is_suffix(path, ".fasta.gz") ||
      is_suffix(path, ".fna") || is_suffix(path, ".fna.gz") ||
      is_suffix(path, ".fa") || is_suffix(path, ".fa.gz")) {
    return bioparser::Parser<biosoup::NucleicAcid>::Create<
        bioparser::FastaParser>(path);  // NOLINT
  }
  if (is_suffix(path, ".fastq") || is_suffix(path, ".fastq.gz") ||
      is_suffix(path, ".fq") || is_suffix(path, ".fq.gz")) {
    return bioparser::Parser<biosoup::NucleicAcid>::Create<
        bioparser::FastqParser>(path);  // NOLINT
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
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> sequences,
    std::span<const std::vector<biosoup::Overlap>> overlaps) {
  std::uint64_t rhs_offset = targets.front()->id;
  std::uint64_t lhs_offset = sequences.front()->id;
  for (const auto& it : overlaps) {
    for (const auto& jt : it) {
      ostrm << sequences[jt.lhs_id - lhs_offset]->name << '\t'
            << sequences[jt.lhs_id - lhs_offset]->inflated_len << '\t'
            << jt.lhs_begin << '\t' << jt.lhs_end << '\t'
            << (jt.strand ? '+' : '-') << "\t"
            << targets[jt.rhs_id - rhs_offset]->name << '\t'
            << targets[jt.rhs_id - rhs_offset]->inflated_len << '\t'
            << jt.rhs_begin << '\t' << jt.rhs_end << '\t' << jt.score << '\t'
            << std::max(jt.lhs_end - jt.lhs_begin, jt.rhs_end - jt.rhs_begin)
            << '\t' << 255 << '\n';
    }
  }

  std::flush(ostrm);
}

void PrintMatchBatch(
    std::ostream& ostrm,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> sequences,
    std::span<const std::vector<Match>> matches) {
  std::uint64_t rhs_offset = targets.front()->id;
  for (auto lhs_idx = 0uz; lhs_idx < matches.size(); ++lhs_idx) {
    for (const auto& match : matches[lhs_idx]) {
      ostrm << sequences[lhs_idx]->name << '\t'
            << sequences[lhs_idx]->inflated_len << '\t' << match.lhs_position()
            << '\t' << (match.strand() ? '+' : '-') << '\t'
            << targets[match.rhs_id() - rhs_offset]->name << '\t'
            << targets[match.rhs_id() - rhs_offset]->inflated_len << '\t'
            << match.rhs_position() << '\n';
    }
  }

  std::flush(ostrm);
}

}  // namespace ram
