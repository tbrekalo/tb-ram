// Copyright (c) 2020 Robert Vaser

#include "bioparser/fasta_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "gtest/gtest.h"
#include "ram/algorithm.hpp"
#include "thread_pool/thread_pool.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace ram {
namespace test {

class RamMinimizerEngineTest : public ::testing::Test {
 public:
  void SetUp() override {
    thread_pool = std::make_shared<thread_pool::ThreadPool>(1);
    biosoup::NucleicAcid::num_objects = 0;
    auto p =
        bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(
            TEST_DATA);
    s = p->Parse(-1);
    EXPECT_EQ(2, s.size());
  }

  std::shared_ptr<thread_pool::ThreadPool> thread_pool;
  std::vector<std::unique_ptr<biosoup::NucleicAcid>> s;
};

TEST_F(RamMinimizerEngineTest, Map) {
  auto indices = ConstructIndices(thread_pool, s, MinimizeConfig{});
  auto occurrence = CalculateKmerThreshold(indices, 0.001);
  std::vector<std::uint32_t> filtered;

  auto o = MapSeqToIndex(s.front(), indices,
                         MapToIndexConfig{.occurrence = occurrence},
                         MinimizeConfig{}, ChainConfig{}, &filtered);
  EXPECT_EQ(1, o.size());
  EXPECT_EQ(0, o.front().lhs_id);
  EXPECT_EQ(30, o.front().lhs_begin);
  EXPECT_EQ(1869, o.front().lhs_end);
  EXPECT_EQ(1, o.front().rhs_id);
  EXPECT_EQ(0, o.front().rhs_begin);
  EXPECT_EQ(1893, o.front().rhs_end);
  EXPECT_EQ(585, o.front().score);
  EXPECT_TRUE(o.front().strand);

  o = MapSeqToIndex(s.back(), indices,
                    MapToIndexConfig{.occurrence = occurrence},
                    MinimizeConfig{}, ChainConfig{}, &filtered);
  EXPECT_TRUE(o.empty());

  o = MapSeqToIndex(
      s.back(), indices,
      MapToIndexConfig{.avoid_symmetric = false, .occurrence = occurrence},
      MinimizeConfig{}, ChainConfig{}, &filtered);

  EXPECT_EQ(1, o.size());
  EXPECT_EQ(1, o.front().lhs_id);
  EXPECT_EQ(0, o.front().lhs_begin);
  EXPECT_EQ(1893, o.front().lhs_end);
  EXPECT_EQ(0, o.front().rhs_id);
  EXPECT_EQ(30, o.front().rhs_begin);
  EXPECT_EQ(1869, o.front().rhs_end);
  EXPECT_EQ(585, o.front().score);
  EXPECT_TRUE(o.front().strand);

  o = MapSeqToIndex(s.front(), indices,
                    MapToIndexConfig{.avoid_equal = false,
                                     .avoid_symmetric = true,
                                     .occurrence = occurrence},
                    MinimizeConfig{}, ChainConfig{}, &filtered);
  EXPECT_EQ(2, o.size());
  EXPECT_EQ(0, o.front().lhs_id);
  EXPECT_EQ(2, o.front().lhs_begin);
  EXPECT_EQ(1897, o.front().lhs_end);
  EXPECT_EQ(0, o.front().rhs_id);
  EXPECT_EQ(2, o.front().rhs_begin);
  EXPECT_EQ(1897, o.front().rhs_end);
  EXPECT_EQ(1895, o.front().score);
  EXPECT_TRUE(o.front().strand);
}

TEST_F(RamMinimizerEngineTest, Pair) {
  auto o = MapPairs(s.front(), s.back(), MinimizeConfig{}, ChainConfig{});
  EXPECT_EQ(1, o.size());
  EXPECT_EQ(0, o.front().lhs_id);
  EXPECT_EQ(30, o.front().lhs_begin);
  EXPECT_EQ(1869, o.front().lhs_end);
  EXPECT_EQ(1, o.front().rhs_id);
  EXPECT_EQ(0, o.front().rhs_begin);
  EXPECT_EQ(1893, o.front().rhs_end);
  EXPECT_EQ(585, o.front().score);
  EXPECT_TRUE(o.front().strand);

  o = MapPairs(s.back(), s.front(), MinimizeConfig{}, ChainConfig{});
  EXPECT_EQ(1, o.size());
  EXPECT_EQ(1, o.front().lhs_id);
  EXPECT_EQ(0, o.front().lhs_begin);
  EXPECT_EQ(1893, o.front().lhs_end);
  EXPECT_EQ(0, o.front().rhs_id);
  EXPECT_EQ(30, o.front().rhs_begin);
  EXPECT_EQ(1869, o.front().rhs_end);
  EXPECT_EQ(585, o.front().score);
  EXPECT_TRUE(o.front().strand);
}

TEST_F(RamMinimizerEngineTest, Filter) {
  auto minimize_config = MinimizeConfig{.kmer_length = 9, .window_length = 3};
  auto indices = ConstructIndices(thread_pool, s, minimize_config);
  auto occurrence_0_001 = CalculateKmerThreshold(indices, 0.001);
  std::vector<std::uint32_t> filtered;

  auto o = MapSeqToIndex(
      s.front(), indices, MapToIndexConfig{.occurrence = occurrence_0_001},
      minimize_config, ChainConfig{.kmer_length = 9}, &filtered);

  EXPECT_EQ(1, o.size());
  EXPECT_EQ(0, o.front().lhs_id);
  EXPECT_EQ(31, o.front().lhs_begin);
  EXPECT_EQ(1888, o.front().lhs_end);
  EXPECT_EQ(1, o.front().rhs_id);
  EXPECT_EQ(1, o.front().rhs_begin);
  EXPECT_EQ(1914, o.front().rhs_end);
  EXPECT_EQ(994, o.front().score);
  EXPECT_TRUE(o.front().strand);

  auto occurrence_0_1 = CalculateKmerThreshold(indices, 0.1);

  o = MapSeqToIndex(s.front(), indices,
                    MapToIndexConfig{.occurrence = occurrence_0_1},
                    minimize_config, ChainConfig{.kmer_length = 9}, &filtered);
  EXPECT_EQ(1, o.size());
  EXPECT_EQ(0, o.front().lhs_id);
  EXPECT_EQ(31, o.front().lhs_begin);
  EXPECT_EQ(1888, o.front().lhs_end);
  EXPECT_EQ(1, o.front().rhs_id);
  EXPECT_EQ(1, o.front().rhs_begin);
  EXPECT_EQ(1914, o.front().rhs_end);
  EXPECT_EQ(980, o.front().score);
  EXPECT_TRUE(o.front().strand);
}

TEST_F(RamMinimizerEngineTest, Micromize) {
  auto indices = ConstructIndices(thread_pool, s, MinimizeConfig{});
  std::vector<std::uint32_t> filtered;

  auto o =
      MapSeqToIndex(s.front(), indices, MapToIndexConfig{},
                    MinimizeConfig{.minhash = true}, ChainConfig{}, &filtered);

  EXPECT_EQ(1, o.size());
  EXPECT_EQ(0, o.front().lhs_id);
  EXPECT_EQ(80, o.front().lhs_begin);
  EXPECT_EQ(1857, o.front().lhs_end);
  EXPECT_EQ(1, o.front().rhs_id);
  EXPECT_EQ(55, o.front().rhs_begin);
  EXPECT_EQ(1881, o.front().rhs_end);
  EXPECT_EQ(242, o.front().score);
  EXPECT_TRUE(o.front().strand);

  o = MapPairs(s.front(), s.back(), MinimizeConfig{.minhash = true},
               ChainConfig{});
  EXPECT_EQ(1, o.size());
  EXPECT_EQ(0, o.front().lhs_id);
  EXPECT_EQ(80, o.front().lhs_begin);
  EXPECT_EQ(1857, o.front().lhs_end);
  EXPECT_EQ(1, o.front().rhs_id);
  EXPECT_EQ(55, o.front().rhs_begin);
  EXPECT_EQ(1881, o.front().rhs_end);
  EXPECT_EQ(242, o.front().score);
  EXPECT_TRUE(o.front().strand);
}

}  // namespace test
}  // namespace ram
