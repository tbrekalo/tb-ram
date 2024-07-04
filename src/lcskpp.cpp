// Copyright 2017: I. Katanic, G. Matula
#include "ram/lcskpp.hpp"

#include <optional>

#include "biosoup/nucleic_acid.hpp"
#include "glog/logging.h"

namespace ram {

namespace {

// If there was no overflow, same as:
//   *c = a*b
//   return *c < 2^63
constexpr std::optional<uint64_t> product_fits_in_63bits(uint64_t a,
                                                         uint64_t b) {
  if (a <= LLONG_MAX / b) {
    return a * b;
  }
  return std::nullopt;
}

// If there was no overflow, same as:
//   *c = pow(a, b)
//   return *c < 2^63
constexpr std::optional<uint64_t> power_fits_in_63bits(uint64_t a, uint64_t b) {
  uint64_t ret = 1;
  while (b > 0) {
    if (auto opt_prod = product_fits_in_63bits(ret, a); opt_prod) {
      ret = *opt_prod;
      --b;

      continue;
    }

    return std::nullopt;
  }
  return ret;
}

constexpr uint64_t next_pw2(uint64_t x) {
  uint64_t power_2 = 1;
  while (power_2 < x) power_2 *= 2;
  return power_2;
}

constexpr int log2(uint64_t x) { return x <= 1 ? 0 : 1 + log2(x / 2); }

// Sorts vector of pairs v, by v.first using radix sort.
// Argument maks should be equal to largest v.first.
void radix_sort(std::vector<std::pair<uint64_t, int>>& v, uint64_t maks) {
  const int len = 9;
  const int sz = 1 << len;
  const int mask = sz - 1;
  int pos[sz + 1];

  std::vector<std::pair<uint64_t, int>> w(v.size());

  int n_blocks = 0;
  while (maks > 1) {
    n_blocks++;
    maks >>= len;
  }

  for (int block = 0; block < n_blocks; ++block) {
    memset(pos, 0, sizeof pos);

    for (auto& x : v) {
      ++pos[((x.first >> (block * len)) & mask) + 1];
    }

    for (int i = 0; i < sz; ++i) {
      pos[i + 1] += pos[i];
    }

    for (auto& x : v) {
      w[pos[(x.first >> (block * len)) & mask]++] = x;
    }

    v.swap(w);
  }
}

struct PrefixDb {
  int column_idx;
  int match_idx;

  friend auto operator<=>(PrefixDb lhs, PrefixDb rhs) = default;
};

struct Match {
  int row_idx;  // row end position (inclusive)
  int col_idx;  // col end position (inclusive)

  friend auto operator<=>(Match lhs, Match rhs) = default;
};

// Finds all k-matches in given pair of strings.
std::vector<Match> calc_matches(const std::string& a, const std::string& b,
                                int k) {
  std::vector<Match> matches;
  if (a.empty() || b.empty()) {
    return {};
  }
  // First, remap characters to interval [0, sigma>,
  // where sigma is alphabet size.
  int n = a.size();
  int m = b.size();
  constexpr int sigma = 4;

  // sigma_to_k is number of possible k length strings
  uint64_t sigma_to_k;
  if (auto opt_sigma_to_k = power_fits_in_63bits(sigma, k);
      not opt_sigma_to_k) {
    fprintf(stderr, "This implementation works only when sigma^k < 2^63.");
    return {};
  } else {
    sigma_to_k = *opt_sigma_to_k;
  }
  // increase sigma_to_k to the next power of two, for cheaper mod operation
  sigma_to_k = next_pw2(sigma_to_k);
  uint64_t mod_mask = sigma_to_k - 1;

  //  a) throw k-length substrings hashes from both string, into 16 buckets,
  //      based on first 4 bits of the hash.
  //  b) radix sort each bucket
  //  c) find consecutive hashes of the same value and save pairs of hashes
  //      from different strings.

  constexpr int nBuckets = 16;
  // std::pair<uint64_t, int> -> (hash, end_position)
  std::vector<std::vector<std::pair<uint64_t, int>>> hashes(nBuckets);
  for (int p = 0; p < nBuckets; ++p) {
    hashes[p].reserve((n + m) / nBuckets);
  }

  uint64_t current_hash = 0;
  for (int i = 0; i < n; ++i) {
    current_hash = (current_hash * sigma +
                    biosoup::kNucleotideCoder[static_cast<uint8_t>(a[i])]) &
                   mod_mask;
    if (i >= k - 1) {
      hashes[current_hash % nBuckets].push_back({current_hash / nBuckets, i});
    }
  }

  current_hash = 0;

  for (int i = 0; i < m; ++i) {
    current_hash = (current_hash * sigma +
                    biosoup::kNucleotideCoder[static_cast<uint8_t>(b[i])]) &
                   mod_mask;
    if (i >= k - 1)
      hashes[current_hash % nBuckets].push_back(
          {current_hash / nBuckets, (n + i)});
  }

  for (auto& bucket : hashes) {
    radix_sort(bucket, sigma_to_k / nBuckets);

    for (std::size_t i = 0, j = 0; i < bucket.size(); i = j) {
      for (j = i + 1; j < bucket.size() && bucket[j].first == bucket[i].first;
           ++j);

      if (j - i > 1) {
        std::size_t s = i;
        while (s < j && bucket[s].second < n) ++s;

        for (std::size_t k1 = i; k1 < s; ++k1) {
          for (std::size_t k2 = s; k2 < j; ++k2) {
            matches.push_back({bucket[k1].second, bucket[k2].second - n});
          }
        }
      }
    }
  }

  // Now we have to sort matching pairs. If there is small number of pairs
  // sort them using std::sort, otherwise use the pigeonhole sort:
  if (matches.size() * log2(matches.size()) < 3u * n) {
    std::ranges::sort(matches);
  } else {
    std::vector<int> pos(n + 1, 0);
    std::vector<Match> tmp(matches.size());
    for (auto& match : matches) {
      pos[match.row_idx + 1]++;
    }
    for (int i = 0; i < n; ++i) {
      pos[i + 1] += pos[i];
    }
    for (auto& p : matches) {
      tmp[pos[p.row_idx]++] = p;
    }
    matches.swap(tmp);
  }
  return matches;
}

}  // namespace

LCSKppResult lcskpp(const std::string& rows, const std::string& cols, int k) {
  if (rows.empty() || cols.empty()) {
    return {};
  }

  int n_cols = cols.size();
  auto matches = calc_matches(rows, cols, k);
  if (matches.empty()) {
    return {};
  }

  int n_matches = matches.size();
  std::vector<PrefixDb> MinYPrefix(n_cols + 1, {n_cols + 1, -1});
  MinYPrefix[0].column_idx = -1;

  int r = 0;      // current lcsk++ length
  int log_r = 1;  // log(r+1) + 1

  std::vector<int> match_dp(matches.size());
  std::vector<int> predecessor(matches.size());

  int match_idx = 0;
  int query_match_idx = 0;
  int cont_idx = 0;
  int last_idx = -1;

  // We process matches row by row (while loop).
  // When at row i, we also handle queries on MinYPrefix for matches
  // in row i+k-1 (query_idx).
  // cont_idx points to ptr's possible continuation (applicable only to k > 1).
  while (match_idx < n_matches) {
    int update_row_index = matches[match_idx].row_idx;  // current row index
    int matches_cur_row_idx = match_idx;

    // Processing queries for future rows untill row row_index+k-1
    while (query_match_idx < n_matches &&
           matches[query_match_idx].row_idx < update_row_index + k) {
      int query_row_idx = matches[query_match_idx].row_idx;

      // count matches in row query_row_i
      int query_row_end = query_match_idx;
      while (query_row_end < n_matches &&
             matches[query_row_end].row_idx == query_row_idx)
        query_row_end++;

      // decide between binary search for each match (Hunt and Szymanski) and
      // sweep through MinYPrefix (Kuo and Cross).
      if ((query_row_end - query_match_idx) * log_r * 6 < r) {
        int last_l = 0;
        while (query_match_idx < query_row_end) {
          int j = matches[query_match_idx].col_idx;

          if (MinYPrefix[last_l].column_idx < j - k + 1) {
            int lo = last_l + 1, hi = r + 1;
            while (lo < hi) {
              int mid = (lo + hi) / 2;
              if (MinYPrefix[mid].column_idx < j - k + 1)
                lo = mid + 1;
              else
                hi = mid;
            }
            last_l = lo;
          }

          predecessor[query_match_idx] = MinYPrefix[last_l - 1].match_idx;
          match_dp[query_match_idx++] = last_l - 1 + k;
        }
      } else {
        int len = 0;
        while (query_match_idx < query_row_end) {
          int col_idx = matches[query_match_idx].col_idx;
          while (MinYPrefix[len].column_idx < col_idx - k + 1) len++;
          predecessor[query_match_idx] = MinYPrefix[len - 1].match_idx;
          match_dp[query_match_idx++] = len - 1 + k;
        }
      }
    }

    if (k > 1) {
      while (cont_idx < matches_cur_row_idx &&
             matches[cont_idx].row_idx < update_row_index - 1)
        cont_idx++;
    }

    // now the main loop through matches in row i, finishing off their
    // dp calculation and doing updates to MinYPrefix
    while (match_idx < n_matches &&
           matches[match_idx].row_idx == update_row_index) {
      int col_idx = matches[match_idx].col_idx;
      int& cur_dp = match_dp[match_idx];

      // update MinYPrefix
      for (int s = cur_dp; s > cur_dp - k && MinYPrefix[s].column_idx > col_idx;
           --s) {
        MinYPrefix[s] = {col_idx, match_idx};
      }

      // Try continuation
      if (k > 1) {
        while (cont_idx < matches_cur_row_idx &&
               matches[cont_idx].col_idx < col_idx - 1)
          cont_idx++;
        if (cont_idx < matches_cur_row_idx &&
            matches[cont_idx].col_idx == col_idx - 1 &&
            match_dp[cont_idx] + 1 > cur_dp) {
          cur_dp = match_dp[cont_idx] + 1;
          predecessor[match_idx] = cont_idx;
          MinYPrefix[cur_dp] =
              std::min(MinYPrefix[cur_dp], PrefixDb{col_idx, match_idx});
        }
      }

      if (cur_dp > r) {
        r = cur_dp;
        last_idx = match_idx;
        while ((1 << log_r) < r + 1) log_r++;
      }

      match_idx++;
    }
  }

  std::vector<MatchInterval> match_intervals;
  MatchInterval cur_interval{
      .query_first = matches[last_idx].row_idx - k + 1,
      .query_last = matches[last_idx].row_idx + 1,

      .target_first = matches[last_idx].col_idx - k + 1,
      .target_last = matches[last_idx].col_idx + 1,
  };

  int i = last_idx;
  while (i != -1) {
    DCHECK(0 <= i && i < matches.size());
    if ((predecessor[i] != -1 &&
         matches[predecessor[i]].row_idx + 1 == matches[i].row_idx &&
         matches[predecessor[i]].col_idx + 1 == matches[i].col_idx)) {
      i = predecessor[i];
      continue;
    }

    // non-overlapping match
    DCHECK(predecessor[i] == -1 ||
           (matches[predecessor[i]].row_idx + k <= matches[i].row_idx &&
            matches[predecessor[i]].col_idx + k <= matches[i].col_idx));

    cur_interval.query_first = matches[i].row_idx - k + 1;
    cur_interval.target_first = matches[i].col_idx - k + 1;
    match_intervals.push_back(cur_interval);

    if (predecessor[i] != -1) {
      cur_interval = {
          .query_first = matches[predecessor[i]].row_idx - k + 1,
          .query_last = matches[predecessor[i]].row_idx + 1,

          .target_first = matches[predecessor[i]].col_idx - k + 1,
          .target_last = matches[predecessor[i]].col_idx + 1,
      };
    }

    i = predecessor[i];
  }

  std::ranges::reverse(match_intervals);
  return {.score = r, .match_intervals = std::move(match_intervals)};
}

}  // namespace ram
