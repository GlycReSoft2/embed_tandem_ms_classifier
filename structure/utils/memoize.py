from functools import wraps
from collections import OrderedDict


def memoize(maxsize=1000):
    """Make a memoization decorator. A negative value of `maxsize` means
    no size limit."""
    def deco(f):
        """Memoization decorator. Items of `kwargs` must be hashable."""
        memo = {}

        @wraps(f)
        def func(*args, **kwargs):
            key = (args, frozenset(kwargs.items()))
            if key not in memo:
                if len(memo) == maxsize:
                    memo.popitem()
                memo[key] = f(*args, **kwargs)
            return memo[key]
        return func
    return deco


def memoize_partial_sequence(maxsize=1000, slice_backs=None):
    if slice_backs is None:
        slice_backs = [1]

    def deco(f):
        memo = OrderedDict()

        @wraps(f)
        def func(sequence, **kwargs):
            kwds = frozenset(kwargs.items())
            key = (sequence, kwds)
            if key not in memo:
                for s in slice_backs:
                    if (sequence[:-s], kwds) in memo:
                        stem = memo[(sequence[:-s], kwds)]
                        term = f(sequence[-s:], **kwargs)
                        memo[key] = stem + term
                        break
                else:
                    memo[key] = f(sequence, **kwargs)
            res = memo[key]
            if len(memo) >= maxsize:
                memo.popitem(False)
            return res
        return func
    return deco
