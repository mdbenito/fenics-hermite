# coding: utf-8

from functools import reduce
from decorator import FunctionMaker
from inspect import getargspec


def fnop(op, doc, *fs):
    """ Combine an arbitrary number of functions with op
        in a new function.

        Arguments:
            op: operator name as a string"""
    n = len(getargspec(fs[0]).args)
    correct = True
    for f in fs[1:]:
        correct = correct and n == len(getargspec(f).args)
    if not correct:
        raise ValueError('All functions must have the same arity.')
    arg_names = ['x'+str(i) for i in range(n)]
    args = '(' + ', '.join(arg_names) + ')'  # '(x0, x1, ..., x(n-1))'
    body = str('return reduce(lambda x,y: x ' + op + ' y, [f' +
               args + ' for f in fs])')
    res = FunctionMaker.create('res' + args, body,
                               {'fs': fs, 'reduce': reduce},
                               addsource=True, doc=doc)
    return res


def fnand(*fs):
    """ Combine an arbitrary number of functions with `and`
        in a new function"""
    def func_name(f):
        arg_names = getargspec(f).args
        args = '(' + ', '.join(arg_names) + ')'
        return f.__name__ + args
    funcs = ', '.join(list(map(func_name, fs)))
    doc = "Returns the 'and' concatenation of %s" % funcs
    return fnop('and', doc, *fs)


def fnnot(f):
    """ Negation of a function."""
    n = len(getargspec(f).args)
    arg_names = ['x'+str(i) for i in range(n)]
    args = '(' + ', '.join(arg_names) + ')'  # '(x0, x1, ..., x(n-1))'
    doc = "Returns the negation of %s%s" % (f.__name__, args)
    res = FunctionMaker.create('res' + args, 'return not f' + args, {'f': f},
                               addsource=True, doc=doc)
    return res


def color_reset(color_func):
    def reset_string(s):
        return color_func(s) + "\x1b[0m"
    return reset_string


@color_reset
def red(s):
    return "\x1b[31m" + s


@color_reset
def green(s):
    return "\x1b[32m" + s


@color_reset
def yellow(s):
    return "\x1b[33m" + s


@color_reset
def blue(s):
    return "\x1b[34m" + s
