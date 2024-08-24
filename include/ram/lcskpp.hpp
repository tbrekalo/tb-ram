#ifndef LCSKPP_H_
#define LCSKPP_H_

#include <string_view>
#include <vector>

namespace biosoup {
class NucleicAcid;
}

namespace ram {

struct LCSKppMatch {
  int row_idx;  // row end position (inclusive)
  int col_idx;  // col end position (inclusive)

  friend auto operator<=>(LCSKppMatch lhs, LCSKppMatch rhs) = default;
};

struct MatchInterval {
  std::int32_t query_first;  // inclusive
  std::int32_t query_last;   // exclusive

  std::int32_t target_first;  // inclusive
  std::int32_t target_last;   // exclusive

  friend constexpr auto operator<=>(MatchInterval, MatchInterval) = default;

  [[gnu::always_inline]] std::int32_t lhs_position() const noexcept {
    return query_first;
  }
  [[gnu::always_inline]] std::int32_t rhs_position() const noexcept {
    return target_first;
  }
};

struct LCSKppResult {
  std::int32_t score;
  std::vector<MatchInterval> match_intervals;
};

struct ArgNucleicAcid {
  biosoup::NucleicAcid* ptr;
  std::uint32_t first;
  std::uint32_t last;
  bool is_rc;
};

LCSKppResult LCSKpp(ArgNucleicAcid lhs, ArgNucleicAcid rhs, std::int32_t k);

LCSKppResult LCSKpp(const std::string_view rows, const std::string_view cols,
                    std::int32_t k);

LCSKppResult LCSKpp(const std::vector<LCSKppMatch>& matches, std::int32_t k,
                    std::int32_t n_cols);

}  // namespace ram

#endif  // LCSKPP_H_
