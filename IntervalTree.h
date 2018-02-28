#ifndef __INTERVAL_TREE_H
#define __INTERVAL_TREE_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>

template <class Scalar, typename Value = std::size_t>
class Interval {
public:
    Value start;
    Value stop;
    Scalar value;
    Interval(Value s, Value e, const Scalar& v)
        : start(s)
        , stop(e)
        , value(v)
    { }
};

template <class Scalar, typename Value>
Value intervalStart(const Interval<Scalar,Value>& i) {
    return i.start;
}

template <class Scalar, typename Value>
Value intervalStop(const Interval<Scalar,Value>& i) {
    return i.stop;
}

template <class Scalar, typename Value>
  std::ostream& operator<<(std::ostream& out, Interval<Scalar,Value>& i) {
    out << "Interval(" << i.start << ", " << i.stop << "): " << i.value;
    return out;
}

template <class Scalar, typename Value = std::size_t>
class IntervalStartSorter {
public:
    bool operator()(const Interval<Scalar,Value>& a, const Interval<Scalar,Value>& b) {
        return a.start < b.start;
    }
};

template <class Scalar, typename Value = std::size_t>
class IntervalTree {

public:
    typedef Interval<Scalar,Value> interval;
    typedef std::vector<interval> interval_vector;

    interval_vector intervals;
    std::unique_ptr<IntervalTree> left;
    std::unique_ptr<IntervalTree> right;
    Value center;

    IntervalTree<Scalar,Value>(void)
        : left(nullptr)
        , right(nullptr)
        , center(0)
    { }

    IntervalTree(IntervalTree&&) = default;
    IntervalTree& operator=(IntervalTree&&) = default;

    std::size_t bytes_used() const {
        return (sizeof(*this) +
                sizeof(interval_vector::value_type) * intervals.capacity() +
                left  ? left->bytes_used()  : 0 +
                right ? right->bytes_used() : 0);
    }
private:
    std::unique_ptr<IntervalTree> clone() const {
        return std::unique_ptr<IntervalTree>(new IntervalTree(*this));
    }
public:

    IntervalTree<Scalar,Value>(const IntervalTree& other)
    :   intervals(other.intervals),
        left(other.left ? other.left->clone() : nullptr),
        right(other.right ? other.right->clone() : nullptr),
        center(other.center)
    {
    }

public:

    IntervalTree<Scalar,Value>& operator=(const IntervalTree& other) {
        center = other.center;
        intervals = other.intervals;
        left = other.left ? other.left->clone() : nullptr;
        right = other.right ? other.right->clone() : nullptr;
        return *this;
    }

    // Note: changes the order of ivals
    IntervalTree<Scalar,Value>(
            interval_vector& ivals,
            std::size_t depth = 16,
            std::size_t minbucket = 64,
            Value leftextent = 0,
            Value rightextent = 0,
            std::size_t maxbucket = 512
            )
        : left(nullptr)
        , right(nullptr)
    {

        --depth;
        IntervalStartSorter<Scalar,Value> intervalStartSorter;
        if (depth == 0 || (ivals.size() < minbucket && ivals.size() < maxbucket)) {
            std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
            intervals = ivals;
        } else {
            if (leftextent == 0 && rightextent == 0) {
                // sort intervals by start
              std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
            }

            Value leftp = 0;
            Value rightp = 0;
            Value centerp = 0;

            if (leftextent || rightextent) {
                leftp = leftextent;
                rightp = rightextent;
            } else {
                leftp = ivals.front().start;
                std::vector<Value> stops;
                stops.resize(ivals.size());
                transform(ivals.begin(), ivals.end(), stops.begin(), intervalStop<Scalar,Value>);
                rightp = *max_element(stops.begin(), stops.end());
            }

            centerp = ivals.at(ivals.size() / 2).start;
            center = centerp;

            interval_vector lefts;
            interval_vector rights;

            for (const interval& interval_i : ivals) {
                if (interval_i.stop < center) {
                    lefts.push_back(interval_i);
                } else if (interval_i.start > center) {
                    rights.push_back(interval_i);
                } else {
                    intervals.push_back(interval_i);
                }
            }

            if (!lefts.empty()) {
                left = std::unique_ptr<IntervalTree>(new IntervalTree(lefts, depth, minbucket, leftp, centerp));
            }
            if (!rights.empty()) {
                right = std::unique_ptr<IntervalTree>(new IntervalTree(rights, depth, minbucket, centerp, rightp));
            }
        }
    }

    interval_vector findOverlapping(const Scalar & pos) const {
      return findOverlapping(pos, pos);
    }

    interval_vector findOverlapping(const Scalar & start, const Scalar & stop) const {
      interval_vector ov;
      findOverlapping(start, stop, [&ov](const interval& i) { ov.push_back(i); });
      return ov;
    }

    template <typename F>
    F findOverlapping(const Scalar & pos, F overlapping) const {
        return findOverlapping(pos, pos, std::move(overlapping));
    }

    template <typename F>
    F findOverlapping(const Scalar & start, const Scalar & stop, F overlapping) const {
        if (left && start <= center) {
            overlapping = left->findOverlapping(start, stop, overlapping);
        }
        if (!intervals.empty() && stop >= intervals.front().start) {
            for (const interval& interval : intervals) {
                if (interval.stop >= start && interval.start <= stop) {
                    overlapping = overlapping(interval);
                }
            }
        }
        if (right && stop >= center) {
            overlapping = right->findOverlapping(start, stop, overlapping);
        }
        return std::move(overlapping);
    }

    interval_vector findContained(Value start, Value stop) const {
        interval_vector contained;
        findContained(start, stop, contained);
        return contained;
    }

    template <typename F>
    F findContained(Value start, Value stop, F contained) const {
        if (!intervals.empty() && !(stop < intervals.front().start)) {
            for (const interval & interv : intervals) {
                if (interv.start >= start && interv.stop <= stop) {
                    contained(interv);
                }
            }
        }

        if (left && start <= center) {
            contained = left->findContained(start, stop, std::move(contained));
        }

        if (right && stop >= center) {
            contained = right->findContained(start, stop, std::move(contained));
        }
        return std::move(contained); // Return std::move(f) just like std::for_each does.
    }

    ~IntervalTree(void) = default;

};

#endif
