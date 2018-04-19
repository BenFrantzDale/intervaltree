# intervaltree

## Overview

An interval tree can be used to efficiently find a set of numeric intervals overlapping or containing another interval.

This library provides a basic implementation of an interval tree using C++ templates, allowing the insertion of arbitrary types into the tree.

## Usage

Add `#include "IntervalTree.h"` to the source files in which you will use the interval tree.

To make an IntervalTree where interval endpoints are `int`s containing `char`s:

```c++
std::vector<Interval<int, char>> intervals;
intervals.push_back(Interval<int, char>(2, 10, 'a'));
intervals.push_back(Interval<int, char>(3, 4, 'b'));
intervals.push_back(Interval<int, char>(20, 100, 'c'));
IntervalTree<int, char> tree;
assert(tree.empty()); // So far it's empty.
tree = IntervalTree<int, char>(std::move(intervals));
```

Member functions `visit_contained` and `visit_overlapping` let you call an arbitrary
function on each interval that is either fully contained by the given range or overalps it.
For example:

```c++
tree.visit_contained(2, 7, [](const Interval<int, char>& i) { std::cout << i.value << std::endl; });
```

will print out "`b`" because only the `Interval<int, char>(3, 4, 'b')` interval is fully contained in [2, 7].
Likewise, 

```c++
tree.visit_overlapping(8, 50, [](const Interval<int, char>& i) { std::cout << i.value << std::endl; });
```

will print out both "`a`" and "`c`" because [2, 10] and [20, 100] both overlap [8, 50] while [3, 4] does not.
Of course, `visit_overlapping` can be used with repeated arguments to find intervals that cross a point and so

```c++
tree.visit_overlapping(4, 4, [](const Interval<int, char>& i) { std::cout << i.value << std::endl; });
```

will output both "`a`" and "`b`".

Methods `findContained(start, stop)` and `findOverlapping(start, stop)` are provided for convenience.
Rather than taking a function, they return a `std::vector<Interval<Scalar, Value>>` containing copies
of the intervals they find:

```c++
auto results = tree.findContained(start, stop);
std::cout << "found " << results.size() << " overlapping intervals" << std::endl;
```

Likewise `IntervalTree::findOverlapping` finds the intervals which are contained or partially overlap the interval [start, stop].

### Author: Erik Garrison <erik.garrison@gmail.com>

### License: MIT
