#ifndef RAM_IO_HPP_
#define RAM_IO_HPP_

#include <memory>

#include "bioparser/parser.hpp"

namespace biosoup {

class NucleicAcid;

}

namespace ram {

std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> CreateParser(
    const std::string& path);

}

#endif  // RAM_IO_HPP_
