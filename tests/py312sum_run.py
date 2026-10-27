#!/usr/bin/env python3
"""Run a script with the built-in sum() of Python >= 3.12 (compensated float summation, Neumaier) even on an
older interpreter, so `lqcore theta --py-sum neumaier` can be checked here. Usage: py312sum_run.py SCRIPT ARGS..."""
import builtins, math, runpy, sys

_orig = builtins.sum


def sum312(iterable, start=0):
    if sys.version_info >= (3, 12):
        return _orig(iterable, start)
    it = iter(iterable)
    result = start
    if type(result) is int:
        for item in it:
            if type(item) in (int, bool):
                result += item
                continue
            result = result + item
            break
        else:
            return result
    if type(result) is float:
        f, c = result, 0.0
        for item in it:
            if type(item) is float:
                t = f + item
                if abs(f) >= abs(item):
                    c += (f - t) + item
                else:
                    c += (item - t) + f
                f = t
                continue
            if type(item) is int:
                f += float(item)
                continue
            r = f + c if (c and math.isfinite(c)) else f
            r = r + item
            for x in it:
                r = r + x
            return r
        if c and math.isfinite(c):
            f += c
        return f
    for item in it:
        result = result + item
    return result


builtins.sum = sum312
script = sys.argv[1]
sys.argv = sys.argv[1:]
runpy.run_path(script, run_name="__main__")
