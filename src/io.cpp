#include "ram/io.hpp"

#include <exception>
#include <sstream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"

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
}  // namespace ram
