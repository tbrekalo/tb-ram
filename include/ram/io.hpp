#ifndef RAM_IO_HPP_
#define RAM_IO_HPP_

#include <memory>
#include <ostream>
#include <span>

#include "bioparser/parser.hpp"

namespace biosoup {

class NucleicAcid;
struct Overlap;

}  // namespace biosoup

namespace ram {

std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> CreateParser(
    const std::string& path);

void PrintOverlapBatch(
    std::ostream& ostrm,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> batch_targets,
    std::span<const std::unique_ptr<biosoup::NucleicAcid>> batch_queries,
    std::span<const std::vector<biosoup::Overlap>> overlaps);

}  // namespace ram

#endif  // RAM_IO_HPP_
