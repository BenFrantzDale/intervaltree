#ifndef __INTERVAL_TREE_H
#define __INTERVAL_TREE_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>
#include <cassert>

template <class Scalar, typename Value>
class Interval {
public:
    Scalar start;
    Scalar stop;
    Value value;
    Interval(const Scalar& s, const Scalar& e, const Value& v)
    : start(std::min(s, e))
    , stop(std::max(s, e))
    , value(v) 
    {}
};

template <class Scalar, typename Value>
Value intervalStart(const Interval<Scalar,Value>& i) {
    return i.start;
}

template <class Scalar, typename Value>
Value intervalStop(const Interval<Scalar, Value>& i) {
    return i.stop;
}

template <class Scalar, typename Value>
std::ostream& operator<<(std::ostream& out, const Interval<Scalar, Value>& i) {
    out << "Interval(" << i.start << ", " << i.stop << "): " << i.value;
    return out;
}

template <typename IntervalType>
struct MinMaxFunction {
    using Scalar = decltype(IntervalType::start);
    auto start(const IntervalType& x) const -> decltype(x.start) const& { return x.start; }
    auto stop(const IntervalType& x) const -> decltype(x.stop) const& { return x.stop; }
};

template <class IntervalType, class StartStopFn = MinMaxFunction<IntervalType>>
class IntervalTreeBase {
public:
    using interval = IntervalType;
    using Scalar = typename StartStopFn::Scalar;
    typedef std::vector<interval> interval_vector;

    struct IntervalStartCmp {
        const StartStopFn& ssf;
        bool operator()(const interval& a, const interval& b) {
            return ssf.start(a) < ssf.start(b);
        }
    };

    struct IntervalStopCmp {
        const StartStopFn& ssf;
        bool operator()(const interval& a, const interval& b) {
            return ssf.stop(a) < ssf.stop(b);
        }
    };

    explicit IntervalTreeBase(StartStopFn ssf = StartStopFn())
        : startStopFn(ssf)
    {}

    ~IntervalTreeBase() = default;

    std::unique_ptr<IntervalTreeBase> clone() const {
        return std::unique_ptr<IntervalTreeBase>(new IntervalTreeBase(*this));
    }

    IntervalTreeBase(const IntervalTreeBase& other)
    :   startStopFn(other.startStopFn),
        tree(other.tree)
    {}

    IntervalTreeBase& operator=(IntervalTreeBase&&) = default;
    IntervalTreeBase(IntervalTreeBase&&) = default;

    IntervalTreeBase& operator=(const IntervalTreeBase& other) {
        startStopFn = other.startStopFn;
        tree = other.tree;
        return *this;
    }

    IntervalTreeBase(
            interval_vector&& ivals,
            StartStopFn ssf = StartStopFn(),
            std::size_t depth = 16,
            std::size_t minbucket = 64,
            std::size_t maxbucket = 512, 
            Scalar leftextent = 0,
            Scalar rightextent = 0)
      : startStopFn(ssf),
        tree(ivals.begin(), ivals.end(),
             ssf,
             depth,
             minbucket, maxbucket,
             leftextent,
             rightextent)
    {}

    Scalar min() const {
        return tree.min(startStopFn);
    }
    Scalar max() const {
        return tree.max(startStopFn);
    }
    // Call f on all intervals near the range [start, stop]:
    template <class UnaryFunction>
    void visit_near(const Scalar& start, const Scalar& stop, UnaryFunction f) const {
        tree.visit_near(start, stop, f, startStopFn);
    }

    // Call f on all intervals crossing pos
    template <class UnaryFunction>
    void visit_overlapping(const Scalar& pos, UnaryFunction f) const {
        visit_overlapping(pos, pos, f);
    }

    // Call f on all intervals overlapping [start, stop]
    template <class UnaryFunction>
    void visit_overlapping(const Scalar& start, const Scalar& stop, UnaryFunction f) const {
        auto filterF = [&](const interval& intrvl) {
            if (startStopFn.stop(intrvl) >= start && startStopFn.start(intrvl) <= stop) {
                // Only apply f if overlapping
                f(intrvl);
            }
        };
        visit_near(start, stop, filterF);
    }

    // Call f on all intervals contained within [start, stop]
    template <class UnaryFunction>
    void visit_contained(const Scalar& start, const Scalar& stop, UnaryFunction f) const {
        auto filterF = [&](const interval& intrvl) {
            if (start <= startStopFn.start(intrvl) && startStopFn.stop(intrvl) <= stop) {
                f(intrvl);
            }
        };
        visit_near(start, stop, filterF);
    }

    interval_vector findOverlapping(const Scalar& start, const Scalar& stop) const {
        interval_vector result;
        visit_overlapping(start, stop,
                          [&](const interval& intrvl) {
                            result.emplace_back(intrvl);
                          });
        return result;
    }

    interval_vector findContained(const Scalar& start, const Scalar& stop) const {
        interval_vector result;
        visit_contained(start, stop,
                        [&](const interval& intrvl) {
                          result.push_back(intrvl);
                        });
        return result;
    }
    bool empty() const {
        return tree.empty();
    }

    template <class UnaryFunction>
    void visit_all(UnaryFunction f) const {
        tree.visit_all(f);
    }

    std::pair<Scalar, Scalar> extentBruitForce() const {
        struct Extent {
            StartStopFn& startStopFn;
            std::pair<Scalar, Scalar> x = {std::numeric_limits<Scalar>::max(),
                                                       std::numeric_limits<Scalar>::min() };
            void operator()(const interval & intrvl) {
                x.first  = std::min(x.first,  startStopFn.start(intrvl));
                x.second = std::max(x.second, startStopFn.stop(intrvl));
            }
                                                                };
        auto extent = Extent{startStopFn};

        visit_all([&](const interval & intrvl) { extent(intrvl); });
        return extent.x;
                                            }

    // Check all constraints.
    // If first is false, second is invalid.
    std::pair<bool, std::pair<Scalar, Scalar>> is_valid() const {
        return tree.is_valid(startStopFn);
    }

    friend std::ostream& operator<<(std::ostream& os, const IntervalTreeBase& itree) {
        return writeOut(os, itree);
    }

    friend std::ostream& writeOut(std::ostream& os, const IntervalTreeBase& itree,
                                  std::size_t depth = 0) {
        auto pad = [&]() { for (std::size_t i = 0; i != depth; ++i) { os << ' '; } };
        pad(); os << "center: " << itree.center << '\n';
        for (const interval & inter : itree.intervals) {
            pad(); os << inter << '\n';
        }
        if (itree.left) {
            pad(); os << "left:\n";
            writeOut(os, *itree.left, depth + 1);
        } else {
            pad(); os << "left: nullptr\n";
        }
        if (itree.right) {
            pad(); os << "right:\n";
            writeOut(os, *itree.right, depth + 1);
        } else {
            pad(); os << "right: nullptr\n";
        }
        return os;
    }

private:
    StartStopFn startStopFn;
    struct Node {
        Scalar center;
        interval_vector intervals;
        std::unique_ptr<Node> left;
        std::unique_ptr<Node> right;

        template <typename IntervalIter>
        Node(const IntervalIter begIntervals,
             const IntervalIter endIntervals,
             StartStopFn& startStopFn = StartStopFn(),
             std::size_t depth = 16,
             std::size_t minbucket = 64,
             std::size_t maxbucket = 512,
             Scalar leftextent = 0,
             Scalar rightextent = 0)
        {
            assert(begIntervals <= endIntervals); // Shouldn't be here with empty range.
            const auto nIntervals = static_cast<std::size_t>(std::distance(begIntervals, endIntervals));
            --depth;
            const auto minmaxStop = std::minmax_element(begIntervals, endIntervals,
                                                        IntervalStopCmp{startStopFn});
            const auto minmaxStart = std::minmax_element(begIntervals, endIntervals,
                                                         IntervalStartCmp{startStopFn});
            if (begIntervals != endIntervals) {
                center = (startStopFn.start(*minmaxStart.first) + startStopFn.stop(*minmaxStop.second)) / 2;
            }
            if (leftextent == 0 && rightextent == 0) {
                // sort intervals by start
                std::sort(begIntervals, endIntervals, IntervalStartCmp{startStopFn});
            } else {
                assert(std::is_sorted(begIntervals, endIntervals, IntervalStartCmp{startStopFn}));
            }
            if (depth == 0 || (nIntervals < minbucket && nIntervals < maxbucket)) {
                std::sort(begIntervals, endIntervals, IntervalStartCmp{startStopFn});
                intervals = interval_vector(begIntervals, endIntervals);
                assert(is_valid(startStopFn).first);
                return;
            } else {
                Scalar leftp = 0;
                Scalar rightp = 0;

                if (leftextent || rightextent) {
                    leftp = leftextent;
                    rightp = rightextent;
                } else {
                    leftp = startStopFn.start(*begIntervals);
                    rightp = startStopFn.stop(*std::max_element(begIntervals, endIntervals,
                                                                IntervalStopCmp{startStopFn}));
                }

                //auto leftEnds = std::partition(begIntervals, endIntervals, [&](const interval& intrvl) { return startStopFn.stop(intrvl) < center; });
                // Should be sorted?
                auto leftEnds = std::find_if(begIntervals, endIntervals, [&](const interval& intrvl) { return !(startStopFn.stop(intrvl) < center); });
                assert(begIntervals <= leftEnds);
                assert(leftEnds <= endIntervals);
                // Now partition the stuff that belongs to this tree node to the left of the stuff that goes in the right child:
                //auto rightBeg = std::partition(leftEnds, endIntervals, [&](const interval& intrvl) { return !(startStopFn.start(intrvl) > center); });
                auto rightBeg = std::find_if(leftEnds, endIntervals, [&](const interval& intrvl) { return startStopFn.start(intrvl) > center; });
                assert(begIntervals <= rightBeg);
                assert(leftEnds <= rightBeg);
                assert(rightBeg <= endIntervals);
                assert(rightBeg <= endIntervals);
                intervals = interval_vector(leftEnds, rightBeg);
                if (begIntervals != leftEnds) {
                    left.reset(new Node(begIntervals, leftEnds,
                                        startStopFn,
                                        depth, minbucket, maxbucket,
                                        leftp, center));
                }
                if (rightBeg != endIntervals) {
                    right.reset(new Node(rightBeg, endIntervals,
                                         startStopFn,
                                         depth, minbucket, maxbucket,
                                         center, rightp));
                }
            }
            assert(is_valid(startStopFn).first);
        }
        bool empty() const {
            if (left && !left->empty()) {
                return false;
            }
            if (!intervals.empty()) {
                return false;
            }
            if (right && !right->empty()) {
                return false;
            }
            return true;
        }

        Scalar min(const StartStopFn& startStopFn) const {
            assert(!empty());
            if (left) { return left->min(startStopFn); }
            return startStopFn.start(*std::min_element(intervals.begin(), intervals.end(),
                                                       IntervalStartCmp{startStopFn}));
        }
        Scalar max(const StartStopFn& startStopFn) const {
            assert(!empty());
            if (right) { return right->max(startStopFn); }
            return startStopFn.stop(*std::max_element(intervals.begin(), intervals.end(),
                                                      IntervalStopCmp{startStopFn}));
        }
        template <class UnaryFunction>
        void visit_near(const Scalar& start, const Scalar& stop, UnaryFunction f, const StartStopFn& startStopFn) const {
            if (!intervals.empty() && ! (stop < startStopFn.start(intervals.front()))) {
                for (auto & i : intervals) {
                    f(i);
                }
            }
            if (left && start <= center) {
                left->visit_near(start, stop, f, startStopFn);
            }
            if (right && stop >= center) {
                right->visit_near(start, stop, f, startStopFn);
            }
        }
        template <class UnaryFunction>
        void visit_all(UnaryFunction f) const {
            if (left) {
                left->visit_all(f);
            }
            std::for_each(intervals.begin(), intervals.end(), f);
            if (right) {
                right->visit_all(f);
            }
        }
        // Check all constraints.
        // If first is false, second is invalid.
        std::pair<bool, std::pair<Scalar, Scalar>> is_valid(const StartStopFn& startStopFn) const {
            const auto minmaxStop = std::minmax_element(intervals.begin(), intervals.end(),
                                                        IntervalStopCmp{startStopFn});
            const auto minmaxStart = std::minmax_element(intervals.begin(), intervals.end(),
                                                         IntervalStartCmp{startStopFn});

            std::pair<bool, std::pair<Scalar, Scalar>> result = {true, { std::numeric_limits<Scalar>::max(),
                                                                         std::numeric_limits<Scalar>::min() }};
            if (!intervals.empty()) {
                result.second.first   = std::min(result.second.first,  startStopFn.start(*minmaxStart.first));
                result.second.second  = std::min(result.second.second, startStopFn.stop(*minmaxStop.second));
            }
            if (left) {
                auto valid = left->is_valid(startStopFn);
                result.first &= valid.first;
                result.second.first   = std::min(result.second.first,  valid.second.first);
                result.second.second  = std::min(result.second.second, valid.second.second);
                if (!result.first) { return result; }
                if (valid.second.second >= center) {
                    result.first = false;
                    return result;
                }
            }
            if (right) {
                auto valid = right->is_valid(startStopFn);
                result.first &= valid.first;
                result.second.first   = std::min(result.second.first,  valid.second.first);
                result.second.second  = std::min(result.second.second, valid.second.second);
                if (!result.first) { return result; }
                if (valid.second.first <= center) {
                    result.first = false;
                    return result;
                }
            }
            if (!std::is_sorted(intervals.begin(), intervals.end(), IntervalStartCmp{startStopFn})) {
                result.first = false;
            }
            return result;
        }
    };
    Node tree;
};


template <typename Scalar, typename Value>
using IntervalTree = IntervalTreeBase<Interval<Scalar, Value>>;

#endif
