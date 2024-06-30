// Copyright 2017: I. Katanic, G. Matula
#ifndef LCSKPP_H_
#define LCSKPP_H_

#include <cassert>
#include <format>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace ram {

namespace detail {

// If there was no overflow, same as:
//   *c = a*b
//   return *c < 2^63
inline bool product_fits_in_63bits(uint64_t a, uint64_t b, uint64_t* c) {
  if (a <= LLONG_MAX / b) {
    *c = a * b;
    return true;
  }
  return false;
}

// If there was no overflow, same as:
//   *c = pow(a, b)
//   return *c < 2^63
inline bool power_fits_in_63bits(uint64_t a, uint64_t b, uint64_t* c) {
  *c = 1;
  while (b > 0) {
    if (!product_fits_in_63bits(*c, a, c)) return false;
    b--;
  }
  return true;
}

inline uint64_t next_pw2(uint64_t x) {
  uint64_t power_2 = 1;
  while (power_2 < x) power_2 *= 2;
  return power_2;
}

inline int log2(uint64_t x) { return x <= 1 ? 0 : 1 + log2(x / 2); }

// Sorts vector of pairs v, by v.first using radix sort.
// Argument maks should be equal to largest v.first.
inline void radix_sort(std::vector<std::pair<uint64_t, int>>& v,
                       uint64_t maks) {
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

    for (auto& x : v) ++pos[((x.first >> (block * len)) & mask) + 1];

    for (int i = 0; i < sz; ++i) pos[i + 1] += pos[i];

    for (auto& x : v) w[pos[(x.first >> (block * len)) & mask]++] = x;

    v.swap(w);
  }
}

// Sorts vector of integers using radix sort.
// Compares only first log2(max_value-1) bits of integers.
inline void radix_sort(std::vector<uint64_t>& v, uint64_t max_value) {
  const int max_len = 9;

  std::vector<uint64_t> w(v.size());
  int n_bits = log2(max_value);

  int shift = 0;
  while (n_bits > 0) {
    std::vector<std::size_t> pos((1 << max_len) + 1, 0);

    int len = std::min(max_len, n_bits);
    int sz = 1 << len;
    int mask = sz - 1;

    for (auto& x : v) {
      ++pos[((x >> shift) & mask) + 1];
    }

    for (int i = 0; i < sz; ++i) {
      pos[i + 1] += pos[i];
    }

    for (auto& x : v) {
      w[pos[(x >> shift) & mask]++] = x;
    }

    shift += len;
    n_bits -= len;
    v.swap(w);
  }
}

// Finds all k-matches in given pair of strings.
inline void calc_matches(const std::string& a, const std::string& b, int k,
                         std::vector<std::pair<int, int>>* matches) {
  if (a.empty() || b.empty()) {
    return;
  }
  // First, remap characters to interval [0, sigma>,
  // where sigma is alphabet size.
  std::vector<int> map(256, -1);
  int sigma = 0;
  int n = a.size();
  int m = b.size();
  for (int i = 0; i < n; ++i) {
    if (map[(uint8_t)a[i]] == -1) {
      map[(uint8_t)a[i]] = sigma++;
    }
  }
  for (int i = 0; i < m; ++i) {
    if (map[(uint8_t)b[i]] == -1) {
      map[(uint8_t)b[i]] = sigma++;
    }
  }
  // sigma_to_k is number of possible k length strings
  uint64_t sigma_to_k;
  if (!power_fits_in_63bits(sigma, k, &sigma_to_k)) {
    fprintf(stderr, "This implementation works only when sigma^k < 2^63.");
    return;
  }
  // increase sigma_to_k to the next power of two, for cheaper mod operation
  sigma_to_k = next_pw2(sigma_to_k);
  uint64_t mod_mask = sigma_to_k - 1;

  // Now to find actual pairs, we use perfect hash-table if sigma_to_k is small
  // enough, and otherwise:
  //   a) throw k-length substrings hashes from both string, into 16 buckets,
  //      based on first 4 bits of the hash.
  //   b) radix sort each bucket
  //   c) find consecutive hashes of the same value and save pairs of hashes
  //      from different strings.

  if (sigma_to_k <= (1 << 20) && sigma_to_k <= 4 * (n + m)) {
    // for each value of hash build a linked list containing such substrings
    // from the second string
    std::vector<int> last(sigma_to_k, -1);
    std::vector<int> prev(m);
    uint64_t current_hash = 0;
    for (int i = 0; i < m; ++i) {
      current_hash = (current_hash * sigma + map[(uint8_t)b[i]]) & mod_mask;
      if (i >= k - 1) {
        prev[i] = last[current_hash];
        last[current_hash] = i;
      }
    }

    // now for each substring of the first string, fetch equal substrings from
    // the linked list
    current_hash = 0;
    for (int i = 0; i < n; ++i) {
      current_hash = (current_hash * sigma + map[(uint8_t)a[i]]) & mod_mask;
      if (i >= k - 1) {
        int sz = matches->size();
        for (int j = last[current_hash]; j != -1; j = prev[j]) {
          matches->push_back({i, j});
        }
        reverse(matches->begin() + sz, matches->end());
      }
    }
  } else {
    // Now we need to store pairs (hash, position in string), if both things
    // can fit into single 64-bit integer, we do so. Otherwise we use std::pair.

    uint64_t _tmp;
    if (product_fits_in_63bits(sigma_to_k, n + m, &_tmp)) {
      const int P = 16;  // number of buckets
      std::vector<uint64_t> hashes[P];
      for (int p = 0; p < P; ++p) {
        hashes[p].reserve((n + m) / P);
      }

      // we use bitwise shifting instead of div/mul by sigma_to_k
      int shift = log2(sigma_to_k);

      uint64_t current_hash = 0;
      for (int i = 0; i < n; ++i) {
        current_hash = (current_hash * sigma + map[(uint8_t)a[i]]) & mod_mask;
        if (i >= k - 1)
          hashes[current_hash % P].push_back((i * 1LLU << shift) +
                                             current_hash / P);
      }

      current_hash = 0;
      for (int i = 0; i < m; ++i) {
        current_hash = (current_hash * sigma + map[(uint8_t)b[i]]) & mod_mask;
        if (i >= k - 1)
          hashes[current_hash % P].push_back(((n + i) * 1LLU << shift) +
                                             current_hash / P);
      }

      for (int p = 0; p < P; ++p) {
        radix_sort(hashes[p], sigma_to_k / P);

        int sz = hashes[p].size();
        for (int i = 0, j = 0; i < sz; i = j) {
          for (j = i + 1;
               j < sz && (hashes[p][j] & mod_mask) == (hashes[p][i] & mod_mask);
               ++j);

          if (j - i > 1) {
            int s = i;
            while (s < j && (hashes[p][s] >> shift) < n) ++s;

            for (int k1 = i; k1 < s; ++k1) {
              for (int k2 = s; k2 < j; ++k2) {
                matches->push_back(
                    {hashes[p][k1] >> shift, (hashes[p][k2] >> shift) - n});
              }
            }
          }
        }
      }
    } else {
      const int P = 16;
      std::vector<std::pair<uint64_t, int>> hashes[P];
      for (int p = 0; p < P; ++p) {
        hashes[p].reserve((n + m) / P);
      }

      uint64_t current_hash = 0;
      for (int i = 0; i < n; ++i) {
        current_hash = (current_hash * sigma + map[(uint8_t)a[i]]) & mod_mask;
        if (i >= k - 1)
          hashes[current_hash % P].push_back({current_hash / P, i});
      }

      current_hash = 0;

      for (int i = 0; i < m; ++i) {
        current_hash = (current_hash * sigma + map[(uint8_t)b[i]]) & mod_mask;
        if (i >= k - 1)
          hashes[current_hash % P].push_back({current_hash / P, (n + i)});
      }

      for (int p = 0; p < P; ++p) {
        radix_sort(hashes[p], sigma_to_k / P);

        int sz = hashes[p].size();
        for (int i = 0, j = 0; i < sz; i = j) {
          int s = 0;
          for (j = i + 1; j < sz && hashes[p][j].first == hashes[p][i].first;
               ++j);

          if (j - i > 1) {
            int s = i;
            while (s < j && hashes[p][s].second < n) ++s;

            for (int k1 = i; k1 < s; ++k1) {
              for (int k2 = s; k2 < j; ++k2) {
                matches->push_back(
                    {hashes[p][k1].second, hashes[p][k2].second - n});
              }
            }
          }
        }
      }
    }

    // Now we have to sort matching pairs. If there is small number of pairs
    // sort them using std::sort, otherwise use the pigeonhole sort:
    if (matches->size() * log2(matches->size()) < 3 * n) {
      sort(matches->begin(), matches->end());
    } else {
      std::vector<int> pos(n + 1, 0);
      std::vector<std::pair<int, int>> tmp(matches->size());
      for (auto& p : *matches) pos[p.first + 1]++;
      for (int i = 0; i < n; ++i) pos[i + 1] += pos[i];
      for (auto& p : *matches) tmp[pos[p.first]++] = p;
      matches->swap(tmp);
    }
  }
}

}  // namespace detail

// Returns length of LCSk++ of given strings a and b.
// LCSK++ indices will be reconstructed to vector pointed to by
// `reconstruction` arg pointer. If it's set to NULL reconstruction is skipped.
inline int lcskpp(const std::string& a, const std::string& b, int k,
                  std::vector<std::pair<int, int>>* reconstruction) {
  if (a.empty() || b.empty()) {
    return 0;
  }

  int m = b.size();
  std::vector<std::pair<int, int>> matches;
  detail::calc_matches(a, b, k, &matches);

  int n_matches = matches.size();
  std::vector<std::pair<int, int>> MinYPrefix(m + 1, {m + 1, -1});
  MinYPrefix[0].first = -1;

  int r = 0;      // current lcsk++ length
  int log_r = 1;  // log(r+1) + 1

  std::vector<int> match_dp(matches.size());
  std::vector<int> recon(matches.size());

  int ptr = 0;
  int query_ptr = 0;
  int cont_ptr = 0;
  int last_idx = -1;

  // We process matches row by row (while loop).
  // When at row i, we also handle queries on MinYPrefix for matches
  // in row i+k-1 (query_ptr).
  // cont_ptr points to ptr's possible continuation (applicable only to k > 1).
  while (ptr < n_matches) {
    int i = matches[ptr].first;  // current row index
    int i_row_ptr = ptr;

    while (query_ptr < n_matches && matches[query_ptr].first < i + k) {
      int query_row_i = matches[query_ptr].first;

      // count matches in row query_row_i
      int query_row_end = query_ptr;
      while (query_row_end < n_matches &&
             matches[query_row_end].first == query_row_i)
        query_row_end++;
      int row_count = query_row_end - query_ptr;

      // decide between binary search for each match (Hunt and Szymanski) and
      // sweep through MinYPrefix (Kuo and Cross).
      if (row_count * log_r * 6 < r) {
        int last_l = 0;
        while (query_ptr < query_row_end) {
          int j = matches[query_ptr].second;

          if (MinYPrefix[last_l].first < j - k + 1) {
            int lo = last_l + 1, hi = r + 1;
            while (lo < hi) {
              int mid = (lo + hi) / 2;
              if (MinYPrefix[mid].first < j - k + 1)
                lo = mid + 1;
              else
                hi = mid;
            }
            last_l = lo;
          }

          recon[query_ptr] = MinYPrefix[last_l - 1].second;
          match_dp[query_ptr++] = last_l - 1 + k;
        }
      } else {
        int l = 0;
        while (query_ptr < query_row_end) {
          int j = matches[query_ptr].second;
          while (MinYPrefix[l].first < j - k + 1) l++;
          recon[query_ptr] = MinYPrefix[l - 1].second;
          match_dp[query_ptr++] = l - 1 + k;
        }
      }
    }

    if (k > 1) {
      while (cont_ptr < i_row_ptr && matches[cont_ptr].first < i - 1)
        cont_ptr++;
    }

    // now the main loop through matches in row i, finishing off their
    // dp calculation and doing updates to MinYPrefix
    while (ptr < n_matches && matches[ptr].first == i) {
      int j = matches[ptr].second;
      int& ptr_dp = match_dp[ptr];

      // update MinYPrefix
      for (int s = ptr_dp; s > ptr_dp - k && MinYPrefix[s].first > j; --s) {
        MinYPrefix[s] = {j, ptr};
      }

      // Try continuation
      if (k > 1) {
        while (cont_ptr < i_row_ptr && matches[cont_ptr].second < j - 1)
          cont_ptr++;
        if (cont_ptr < i_row_ptr && matches[cont_ptr].second == j - 1 &&
            match_dp[cont_ptr] + 1 > ptr_dp) {
          ptr_dp = match_dp[cont_ptr] + 1;
          recon[ptr] = cont_ptr;
          MinYPrefix[ptr_dp] = min(MinYPrefix[ptr_dp], {j, ptr});
        }
      }

      if (ptr_dp > r) {
        r = ptr_dp;
        last_idx = ptr;
        while ((1 << log_r) < r + 1) log_r++;
      }

      ptr++;
    }
  }

  if (reconstruction) {
    reconstruction->clear();

    int i = last_idx;
    while (i != -1) {
      assert(0 <= i && i < matches.size());

      if (recon[i] != -1 && matches[recon[i]].first + 1 == matches[i].first &&
          matches[recon[i]].second + 1 == matches[i].second) {
        // continuation
        reconstruction->push_back({matches[i].first, matches[i].second});
      } else {
        // non-overlapping match
        assert(recon[i] == -1 ||
               (matches[recon[i]].first + k <= matches[i].first &&
                matches[recon[i]].second + k <= matches[i].second));
        for (int j = 0; j < k; ++j) {
          reconstruction->push_back(
              {matches[i].first - j, matches[i].second - j});
        }
      }

      i = recon[i];
    }
    reverse(reconstruction->begin(), reconstruction->end());
  }

  return r;
}

}  // namespace ram

#endif  // LCSKPP_H_
