#ifndef LCSKPP_H_
#define LCSKPP_H_

#include <string>
#include <vector>

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

// Returns length of LCSk++ of given strings a and b.
// LCSK++ indices will be reconstructed to vector pointed to by
// `reconstruction` arg pointer. If it's set to NULL reconstruction is skipped.
LCSKppResult lcskpp(const std::string& rows, const std::string& cols, int k);

}  // namespace ram

#endif  // LCSKPP_H_
