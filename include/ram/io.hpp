#ifndef RAM_IO_HPP_
#define RAM_IO_HPP_

#include <memory>
#include <ostream>
#include <span>

#include "bioparser/parser.hpp"
#include "ram/types.hpp"

namespace biosoup {

class NucleicAcid;
struct Overlap;

}  // namespace biosoup

namespace ram {

std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> CreateParser(
    const std::string& path);

void PrintOverlapBatch(
    std::ostream& ostrm,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> queries,
    std::span<const std::vector<biosoup::Overlap>> overlaps);

void PrintOverlapAIBatch(
    std::ostream& ostrm,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> queries,
    std::span<const std::vector<OverlapAI>> overlaps);

void PrintMatchBatch(
    std::ostream& ostrm,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> queries,
    std::span<const std::vector<Match>> matches);

void PrintMatchChainBatch(
    std::ostream& ostrm,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> targets,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> queries,
    std::span<const std::vector<MatchChain>> match_chains);

}  // namespace ram

#endif  // RAM_IO_HPP_
