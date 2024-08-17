#ifndef LCSKPP_H_
#define LCSKPP_H_

#include <string_view>
#include <vector>

namespace biosoup {
class NucleicAcid;
}

namespace ram {

struct MatchInterval {
  std::int32_t query_first;  // inclusive
  std::int32_t query_last;   // exclusive

  std::int32_t target_first;  // inclusive
  std::int32_t target_last;   // exclusive

  friend constexpr auto operator<=>(MatchInterval, MatchInterval) = default;
};

struct LCSKppResult {
  std::int32_t score;
  std::vector<MatchInterval> match_intervals;
};

struct ArgNucleicAcid {
  biosoup::NucleicAcid* ptr;
  std::int32_t start;
  std::int32_t end;
  bool is_rc;
};

LCSKppResult LCSKpp(ArgNucleicAcid lhs, ArgNucleicAcid rhs, std::int32_t k);

LCSKppResult LCSKpp(const std::string_view rows, const std::string_view cols,
                    std::int32_t k);

}  // namespace ram

#endif  // LCSKPP_H_
