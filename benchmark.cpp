
#include "omp/Random.h"
#include "omp/Hand.h"
#include "omp/HandEvaluator.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <tbb/tbb.h>

#if OMP_BENCHMARK_3RD_PARTY
    // Include last because of all the badly named macros. I hate C programmers.
    #include "evaluators/SKPokerEval/SevenEval.h"
    #include "evaluators/twoplustwo/twoplustwo.h"
    #include "evaluators/ACE_eval/ace_eval.h"
    #include "evaluators/HoldemShowdown/HandEval.h"
    #include "evaluators/poker-eval/include/poker_defs.h"
    #include "evaluators/poker-eval/include/inlines/eval.h"
#endif

using namespace std;
using namespace tbb;

// Base class for evaluator adaptors utilizing CRTP.
template<class TEval, class THand = nullptr_t>
class AdaptorBase
{
public:
    typedef THand Hand;
    // Card index format:
    static const unsigned CARD_OFFSET = 0;
    static const unsigned RANK_MAJOR = true;
    static const unsigned ASCENDING_RANKS = true;

    void initHand(THand& h) const
    {
    }

    void addCard(THand& h, unsigned cardIdx) const
    {
    }

    static unsigned cardIdxToCanonical(unsigned idx)
    {
        idx -= TEval::CARD_OFFSET;
        unsigned rank = TEval::RANK_MAJOR ? idx >> 2 : idx % 13;
        unsigned suit = TEval::RANK_MAJOR ? idx & 3 : idx / 13;
        rank = TEval::ASCENDING_RANKS ? rank : 12 - rank;
        return 4 * rank + suit;
    }

    static unsigned cardIdxFromCanonical(unsigned idx)
    {
        unsigned rank = idx >> 2;
        unsigned suit = idx & 3;
        rank = TEval::ASCENDING_RANKS ? rank : 12 - rank;
        return (TEval::RANK_MAJOR ? 4 * rank + suit : 13 * suit + rank) + TEval::CARD_OFFSET;
    }
};

// OMPEval
class Omp : public omp::HandEvaluator, public AdaptorBase<Omp, omp::Hand>
{
public:
    void initHand(Hand& h) const
    {
        h = Hand::empty();
    }

    void addCard(Hand& h, unsigned cardIdx) const
    {
        h += cardIdx;
    }

    unsigned evaluate(const Hand& h, unsigned c1, unsigned c2, unsigned c3, unsigned c4,
            unsigned c5, unsigned c6, unsigned c7) const
    {
        return omp::HandEvaluator::evaluate(h);
    }
};

#if OMP_BENCHMARK_3RD_PARTY

// SKPokerEval
class Skpe : public SevenEval, public AdaptorBase<Skpe>
{
public:
    static const unsigned ASCENDING_RANKS = false;

    unsigned evaluate(const Hand& h, unsigned c1, unsigned c2, unsigned c3, unsigned c4,
            unsigned c5, unsigned c6, unsigned c7) const
    {
        return GetRank(c1, c2, c3, c4, c5, c6, c7);
    }
};

// 2+2 Evaluator
// This requires some modifications to 2+2 code so that we don't have to read the table from a file.
class Tpt: public AdaptorBase<Tpt>
{
public:
    typedef unsigned Hand;
    static const unsigned CARD_OFFSET = 1; // Card have indexes 1-52.

    Tpt()
    {
        static bool initVar = (generateArrays(), initVar);
    }

    void initHand(Hand& h) const
    {
        h = 53;
    }

    void addCard(Hand& h, unsigned cardIdx) const
    {
        h = HR[h + cardIdx];
    }

    unsigned evaluate(const Hand& h, unsigned c1, unsigned c2, unsigned c3, unsigned c4,
            unsigned c5, unsigned c6, unsigned c7) const
    {
        return h;
    }
};

// ACE Evaluator
class Ace: public AdaptorBase<Ace, array<Card,ACEHAND>>
{
public:
    static const unsigned RANK_MAJOR = false;

    void initHand(Hand& h) const
    {
        h = Hand{};
    }

    void addCard(Hand& h, unsigned cardIdx) const
    {
        static uint32_t ACE_CARDS[52] = {
            0x41, 0x101, 0x401, 0x1001, 0x4001, 0x10001, 0x40001, 0x100001,
                0x400001, 0x1000001, 0x4000001, 0x10000001, 0x40000001,
            0x42, 0x102, 0x402, 0x1002, 0x4002, 0x10002, 0x40002, 0x100002,
                0x400002, 0x1000002, 0x4000002, 0x10000002, 0x40000002,
            0x44, 0x104, 0x404, 0x1004, 0x4004, 0x10004, 0x40004, 0x100004,
                0x400004, 0x1000004, 0x4000004, 0x10000004, 0x40000004,
            0x48, 0x108, 0x408, 0x1008, 0x4008, 0x10008, 0x40008, 0x100008,
                0x400008, 0x1000008, 0x4000008, 0x10000008, 0x40000008
        };
        ACE_addcard(h, ACE_CARDS[cardIdx]);
    }

    unsigned evaluate(const Hand& h, unsigned c1, unsigned c2, unsigned c3, unsigned c4,
            unsigned c5, unsigned c6, unsigned c7) const
    {
        Card E(Card h[]); // Works if ACE compiled as c++, otherwise we need extern C block.
        return ACE_evaluate(const_cast<Card*>(&h[0]));
    }
};

// Steve Brecher's Holdem Showdown
class Sbhs: public AdaptorBase<Sbhs, Hand_T>
{
public:
    static const unsigned RANK_MAJOR = false;

    Sbhs()
    {
        static bool initVar = (Init_Hand_Eval(), initVar);
    }

    void initHand(Hand& h) const
    {
        ZeroHand(h);
    }

    void addCard(Hand& h, unsigned cardIdx) const
    {
        Hand c;
        c.as64Bits = IndexToMask(cardIdx);
        AddHandTo(h, c);
    }

    unsigned evaluate(const Hand& h, unsigned c1, unsigned c2, unsigned c3, unsigned c4,
            unsigned c5, unsigned c6, unsigned c7) const
    {
        return Hand_7_Eval(h);
    }
};

// Pokersource poker-eval
class Pse: public AdaptorBase<Sbhs, StdDeck_CardMask>
{
public:
    void initHand(Hand& h) const
    {
        h = Hand{};
    }

    void addCard(Hand& h, unsigned cardIdx) const
    {
        StdDeck_CardMask_OR(h, h, StdDeck_MASK(cardIdx));
    }

    unsigned evaluate(const Hand& h, unsigned c1, unsigned c2, unsigned c3, unsigned c4,
            unsigned c5, unsigned c6, unsigned c7) const
    {
        return StdDeck_StdRules_EVAL_N(h, 7);
    }
};

#endif //OMP_BENCHMARK_3RD_PARTY

// Benchmark template that allows testing any evaluator with a proper adaptor.
template<typename TEval>
class Benchmark
{
public:
    typedef typename TEval::Hand Hand;

    // Run test and benchmarks.
    void run()
    {
        cout << endl;
        if (!is_same<TEval,Omp>::value)
            test(Omp());
        sequential<false>();
        random1();
        if (!is_same<Hand,nullptr_t>::value)
            random2();
        sequential<true>();
    }

private:
    // Benchmark sequential evaluation.
    template<bool tSingleSuit>
    void sequential()
    {
        cout << "Sequential evaluation" << (tSingleSuit ? " (flush hands):" : ":") <<  endl;

        static const unsigned END = (tSingleSuit ? 13 : 52) + TEval::CARD_OFFSET;
        
        unsigned sum = 0;
        uint64_t count = 0;

        auto t1 = chrono::high_resolution_clock::now();
        
        /*
        for (unsigned i = 0; i < (tSingleSuit ? 200000 : 5); ++i) {
            for (unsigned c1 = TEval::CARD_OFFSET; c1 < END; c1++) {
                for (unsigned c2 = c1 + 1; c2 < END; c2++) {
                    for (unsigned c3 = c2 + 1; c3 < END; c3++) {
                        for (unsigned c4 = c3 + 1; c4 < END; c4++) {
                            Hand h4;
                            mEval.initHand(h4);
                            mEval.addCard(h4, c1);
                            mEval.addCard(h4, c2);
                            mEval.addCard(h4, c3);
                            mEval.addCard(h4, c4);
                            for (unsigned c5 = c4 + 1; c5 < END; c5++) {
                                Hand h5 = h4;
                                mEval.addCard(h5, c5);
                                for (unsigned c6 = c5 + 1; c6 < END; c6++) {
                                    Hand h6 = h5;
                                    mEval.addCard(h6, c6);
                                    for (unsigned c7 = c6 + 1; c7 < END; c7++) {
                                        Hand h7 = h6;
                                        mEval.addCard(h7, c7);
                                        sum += mEval.evaluate(h7, c1, c2, c3, c4, c5, c6, c7);
                                        ++count;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        */

        unsigned suit_count = (tSingleSuit ? 200000 : 5);
        unsigned* sums = new unsigned[suit_count*END];
        uint64_t* counts = new uint64_t[suit_count*END];
        //init 0
        std::fill(sums,sums+suit_count*END,0);
        std::fill(counts,counts+suit_count*END,0);

        tbb::parallel_for(tbb::blocked_range2d<unsigned>(0, suit_count, TEval::CARD_OFFSET, END),
        [&](const tbb::blocked_range2d<unsigned>& range) {
            for (unsigned i = range.rows().begin(); i < range.rows().end(); ++i) {
                    for (unsigned c1 = range.cols().begin(); c1 < range.cols().end(); ++c1) {
                        //do parall
                        for (unsigned c2 = c1 + 1; c2 < END; c2++) {
                            for (unsigned c3 = c2 + 1; c3 < END; c3++) {
                                for (unsigned c4 = c3 + 1; c4 < END; c4++) {
                                    Hand h4;
                                    mEval.initHand(h4);
                                    mEval.addCard(h4, c1);
                                    mEval.addCard(h4, c2);
                                    mEval.addCard(h4, c3);
                                    mEval.addCard(h4, c4);
                                    for (unsigned c5 = c4 + 1; c5 < END; c5++) {
                                        Hand h5 = h4;
                                        mEval.addCard(h5, c5);
                                        for (unsigned c6 = c5 + 1; c6 < END; c6++) {
                                            Hand h6 = h5;
                                            mEval.addCard(h6, c6);
                                            for (unsigned c7 = c6 + 1; c7 < END; c7++) {
                                                Hand h7 = h6;
                                                mEval.addCard(h7, c7);
                                                sums[i*END+c1] += mEval.evaluate(h7, c1, c2, c3, c4, c5, c6, c7);
                                                ++counts[i*END+c1];
                                            }
                                        }
                                    }
                                }
                            }
                        }//end of tbb parallel       
                    } //for 2
            }//for 1
        });

        
        sum = tbb::parallel_reduce(
            tbb::blocked_range<unsigned>(0, suit_count*END), 0, 
                    [&](const tbb::blocked_range<unsigned>& rg, unsigned init) {
                for (unsigned r = rg.begin(); r != rg.end(); ++r) {
                    init += sums[r];
                }
                return init; 
            },
            [](unsigned x, unsigned y) {
            return x + y; 
            }
        );

        count = tbb::parallel_reduce(
            tbb::blocked_range<unsigned>(0, suit_count*END), 0, 
                    [&](const tbb::blocked_range<unsigned>& rg, uint64_t init) {
                for (unsigned r = rg.begin(); r != rg.end(); ++r) {
                    init += counts[r];
                }
                return init; 
            },
            [](uint64_t x, uint64_t y) {
            return x + y; 
            }
        );

        auto t2 = chrono::high_resolution_clock::now();
        double t = 1e-9 * chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count();
        cout << "   " << count << " evals  " << (1e-6 * count / t) << "M/s  " << t << "s  " << sum << endl;

        delete []sums;
        delete []counts;

    }

    // Benchmark random order evaluation using card arrays.
    void random1()
    {
        cout << "Random order evaluation (card arrays):" << endl;
        omp::XoroShiro128Plus rng(0);
        omp::FastUniformIntDistribution<unsigned> rnd(0, 51);
        
        uint64_t count = 0;
        unsigned sum = 0;
        const uint64_t ITERS = 50;
        const uint64_t RHANDS = 10000000;
        
        unsigned* sums = new unsigned[ITERS*RHANDS];
        unsigned* counts = new unsigned[ITERS*RHANDS];
        //unsigned tests[1][1];
        //all_sums[0][0]=0;
        //all_counts[0][0]=0;
        //all_sums[49][9999999]=0;
        //all_counts[49][9999999]=0;
        //std::fill(all_sums,all_sums+hands_count,0);
        //std::fill(all_counts,all_counts+hands_count,0);

        vector<array<uint8_t,7>> table = generateRandomHands(RHANDS);

        auto t1 = chrono::high_resolution_clock::now();
        /*
        for (int i = 0; i < 50; ++i) {
            for (auto& cards : table) {
                Hand h;
                mEval.initHand(h);
                mEval.addCard(h, cards[0]);
                mEval.addCard(h, cards[1]);
                mEval.addCard(h, cards[2]);
                mEval.addCard(h, cards[3]);
                mEval.addCard(h, cards[4]);
                mEval.addCard(h, cards[5]);
                mEval.addCard(h, cards[6]);
                sum += mEval.evaluate(h, cards[0], cards[1], cards[2], cards[3], cards[4], cards[5], cards[6]);
                ++count;
            }
        }
        */

        tbb::parallel_for(tbb::blocked_range2d<uint64_t>(0, ITERS, 0, RHANDS),
        [&](const tbb::blocked_range2d<uint64_t>& range) {
            for (uint64_t i = range.rows().begin(); i < range.rows().end(); ++i) {
                    for (uint64_t j = range.cols().begin(); j < range.cols().end(); ++j) {
                        
                        Hand h;
                        mEval.initHand(h);
                        mEval.addCard(h, table[j][0]);
                        mEval.addCard(h, table[j][1]);
                        mEval.addCard(h, table[j][2]);
                        mEval.addCard(h, table[j][3]);
                        mEval.addCard(h, table[j][4]);
                        mEval.addCard(h, table[j][5]);
                        mEval.addCard(h, table[j][6]);
                        sums[i*RHANDS+j] = mEval.evaluate(h, table[j][0], table[j][1], table[j][2], table[j][3], table[j][4], table[j][5], table[j][6]);
                        counts[i*RHANDS+j]=1;
                        
                    }
            }
        });

        
        sum = tbb::parallel_reduce(
            tbb::blocked_range<unsigned>(0, ITERS*RHANDS), 0, 
                    [&](const tbb::blocked_range<unsigned>& rg, unsigned init) {
                for (unsigned r = rg.begin(); r != rg.end(); ++r) {
                    init += sums[r];
                }
                return init; 
            },
            [](unsigned x, unsigned y) {
            return x + y; 
            }
        );

        count = tbb::parallel_reduce(
            tbb::blocked_range<unsigned>(0, ITERS*RHANDS), 0, 
                    [&](const tbb::blocked_range<unsigned>& rg, uint64_t init) {
                for (unsigned r = rg.begin(); r != rg.end(); ++r) {
                    init += counts[r];
                }
                return init; 
            },
            [](uint64_t x, uint64_t y) {
            return x + y; 
            }
        );
        

        auto t2 = chrono::high_resolution_clock::now();
        double t = 1e-9 * chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count();
        cout << "   " << count << " evals  " << (1e-6 * count / t) << "M/s  " << t << "s  " << sum << endl;
        
        delete []sums;
        delete []counts;
    }

    // Benchmark random order evaluation using Hand objects.
    void random2()
    {
        cout << "Random order evaluation (precalculated Hand objects):" << endl;
        omp::XoroShiro128Plus rng(0);
        omp::FastUniformIntDistribution<unsigned> rnd(0, 51);
        uint64_t count = 0;
        unsigned sum = 0;

        vector<Hand> table;
        for (unsigned i = 0; i < 10000000; ++i) {
            uint64_t usedCardsMask = 0;
            table.emplace_back();
            mEval.initHand(table.back());
            for(unsigned j = 0; j < 7; ++j) {
                unsigned card;
                uint64_t cardMask;
                do {
                    card = rnd(rng) + TEval::CARD_OFFSET;
                    cardMask = 1ull << card;
                } while (usedCardsMask & cardMask);
                usedCardsMask |= cardMask;
                mEval.addCard(table.back(), card);
            }
        }

        uint64_t hands_size = table.size();
        const uint64_t ITERS = 50;
        unsigned* sums = new unsigned[ITERS*hands_size];
        unsigned* counts = new unsigned[ITERS*hands_size];

        auto t1 = chrono::high_resolution_clock::now();

        /*
        for (int i = 0; i < ITERS; ++i) {
            for (auto& hand: table) {
                sum += mEval.evaluate(hand, 0, 0, 0, 0, 0, 0, 0);
                ++count;
            }
        }
        */
        
        tbb::parallel_for(tbb::blocked_range2d<uint64_t>(0, ITERS, 0, hands_size),
        [&](const tbb::blocked_range2d<uint64_t>& range) {
            for (uint64_t i = range.rows().begin(); i < range.rows().end(); ++i) {
                    for (uint64_t j = range.cols().begin(); j < range.cols().end(); ++j) {
                        sums[i*hands_size+j] = mEval.evaluate(table[j], 0, 0, 0, 0, 0, 0, 0);
                        counts[i*hands_size+j]=1;
                    }
            }
        });
        
        sum = tbb::parallel_reduce(
            tbb::blocked_range<unsigned>(0, ITERS*hands_size), 0, 
                    [&](const tbb::blocked_range<unsigned>& rg, unsigned init) {
                for (unsigned r = rg.begin(); r != rg.end(); ++r) {
                    init += sums[r];
                }
                return init; 
            },
            [](unsigned x, unsigned y) {
            return x + y; 
            }
        );

        count = tbb::parallel_reduce(
            tbb::blocked_range<unsigned>(0, ITERS*hands_size), 0, 
                    [&](const tbb::blocked_range<unsigned>& rg, uint64_t init) {
                for (unsigned r = rg.begin(); r != rg.end(); ++r) {
                    init += counts[r];
                }
                return init; 
            },
            [](uint64_t x, uint64_t y) {
            return x + y; 
            }
        );

        auto t2 = chrono::high_resolution_clock::now();
        double t = 1e-9 * chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count();
        cout << "   " << count << " evals  " << (1e-6 * count / t) << "M/s  " << t << "s  " << sum << endl;

        delete []sums;
        delete []counts;
    }

    // Generate a vector of random hands. The random seed is deterministic on purpose.
    static vector<array<uint8_t,7>> generateRandomHands(size_t count)
    {
        omp::XoroShiro128Plus rng(0);
        omp::FastUniformIntDistribution<unsigned> rnd(0, 51);

        vector<array<uint8_t,7>> table;
        for (unsigned i = 0; i < count; ++i) {
            uint64_t usedCardsMask = 0;
            table.emplace_back();
            for (auto& card : table.back()) {
                uint64_t cardMask;
                do {
                    card = rnd(rng);
                    cardMask = 1ull << card;
                } while (usedCardsMask & cardMask);
                usedCardsMask |= cardMask;
            }
        }

        return table;
    }

    // Test this evaluator against another one by comparing the order in which they rank random hands.
    template<class TEval2>
    void test(const TEval2& eval2)
    {
        auto table = generateRandomHands(100000);
        auto sorted1 = sortHands(mEval, table);
        auto sorted2 = sortHands(eval2, table);
        for (size_t i = 0; i < sorted1.size(); ++i) {
            if (sorted1[i].first != sorted2[i].first) {
                cout << "Incorrectly ranked hand: " << i << endl;
                return;
            }
        }
        cout << "Test OK." << endl;
    }

    // Sort a list of hands based on their rank using another evaluator.
    template<class TEval2>
    static vector<pair<size_t,unsigned>> sortHands(const TEval2& eval, const vector<array<uint8_t,7>>& table)
    {
        vector<pair<size_t,unsigned>> handValues;
        handValues.reserve(table.size());
        size_t i = 0;
        for (auto& cards: table) {
            typename TEval2::Hand h;
            eval.initHand(h);
            for (uint8_t c : cards)
                eval.addCard(h, TEval2::cardIdxFromCanonical(c));
            unsigned hv = eval.evaluate(h,
                    TEval2::cardIdxFromCanonical(cards[0]),
                    TEval2::cardIdxFromCanonical(cards[1]),
                    TEval2::cardIdxFromCanonical(cards[2]),
                    TEval2::cardIdxFromCanonical(cards[3]),
                    TEval2::cardIdxFromCanonical(cards[4]),
                    TEval2::cardIdxFromCanonical(cards[5]),
                    TEval2::cardIdxFromCanonical(cards[6])
                    );
            handValues.emplace_back(i++, hv);
        }
        sort(handValues.begin(), handValues.end(), [](const pair<size_t,unsigned>& lhs,
                const pair<size_t,unsigned>& rhs){
            if (lhs.second != rhs.second)
                return lhs.second < rhs.second;
            return lhs.first < rhs.first;
        });
        return handValues;
    }

    TEval mEval;
};

void benchmark()
{
    // Benchmark only one at a time because there's some weird performance interference.
    Benchmark<Omp>().run();
    //Benchmark<Skpe>().run();
    //Benchmark<Tpt>().run();
    //Benchmark<Ace>().run();
    //Benchmark<Sbhs>().run();
    //Benchmark<Pse>().run();
}
