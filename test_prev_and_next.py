#!/usr/bin/env python

import itertools as it

def prev_and_next(iterable):
    prevs, items, nexts = it.tee(iterable, 3)
    prevs = it.chain([None], prevs)
    nexts = it.chain(it.islice(nexts, 1, None), [None])
    return zip(prevs, items, nexts)

def prev_and_next_circular(iterable):
    # Number of elements
    n = sum(1 for _ in iterable)
    prevs, items, nexts = it.tee(iterable, 3)
    prevs = it.islice(it.cycle(prevs), n-1, None)
    nexts = it.islice(it.cycle(nexts), 1, None)
    return zip(prevs, items, nexts)

a = range(0,10)

print('# prev_and_next')
for p, t, n in prev_and_next(a):
    print(p, t, n)
print()

print('# prev_and_next_circular')
for p, t, n in prev_and_next_circular(a):
    print(p, t, n)
