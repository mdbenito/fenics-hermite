# coding: utf-8

from functools import reduce
from decorator import FunctionMaker
from inspect import getargspec
from time import time

__all__ = ['fnand', 'fnor', 'fnnot', 'red', 'green', 'yellow', 'blue']

##############################################################################
# Boolean operations on functions

def fnop(op, doc, *fs):
    """ Combine an arbitrary number of functions using binary operator op
        into a new function of the same arity.

        Arguments:
        ----------
            op: operator name as a string ("and", "or", "+", ...)
            doc: docstring for the returned callable.

        Returns:
        --------
            A callable of the same arity as the arguments.
        """
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
        into a new function.

    Arguments:
    ----------
        *fs: one or more Callables of equal arity and boolean
             return value.
        
    Returns:
    --------
        A Callable which returns True iff all of the *fs return
        True on the arguments.
    """
    def func_name(f):
        arg_names = getargspec(f).args
        args = '(' + ', '.join(arg_names) + ')'
        return f.__name__ + args
    funcs = ', '.join(list(map(func_name, fs)))
    doc = "Returns the 'and' concatenation of %s" % funcs
    return fnop('and', doc, *fs)


def fnor(*fs):
    """ Combine an arbitrary number of functions with `or`
        into a new function.
    
    Arguments:
    ----------
        *fs: one or more Callables of equal arity and boolean
             return value.
        
    Returns:
    --------
        A Callable which returns True iff any of the *fs return
        True on the arguments.
    """
    def func_name(f):
        arg_names = getargspec(f).args
        args = '(' + ', '.join(arg_names) + ')'
        return f.__name__ + args
    funcs = ', '.join(list(map(func_name, fs)))
    doc = "Returns the 'or' concatenation of %s" % funcs
    return fnop('or', doc, *fs)


def fnnot(f):
    """ Negation of a function.
    
    Arguments:
    ----------
        *fs: one Callable of the same arity as f and boolean
             return value.
        
    Returns:
    --------
        A Callable which returns True iff f(...) returns False.
    """
    n = len(getargspec(f).args)
    arg_names = ['x'+str(i) for i in range(n)]
    args = '(' + ', '.join(arg_names) + ')'  # '(x0, x1, ..., x(n-1))'
    doc = "Returns the negation of %s%s" % (f.__name__, args)
    res = FunctionMaker.create('res' + args, 'return not f' + args, {'f': f},
                               addsource=True, doc=doc)
    return res


##############################################################################
# Taking derivatives of functions

import autograd as ad
#import autograd.numpy as np

def make_derivatives(fun):
    """ Decorates fun with an additional argument to compute derivatives.
    The decorated function accepts two arguments: an np.ndarray with the
    point of evaluation and a tuple "derivatives" with a multiindex.
    Example:
        @make_derivatives
        def fun2d(x, y):
            return x*y**2
        
        fun2d([1., 2.])  --> 4.
        fun2d([0., 3.], derivatives=(1,0))  --> 9.
    """
    def diff_fun(x, derivatives=()):
        ret = fun
        for i, n in enumerate(derivatives):
            for j in range(n):
                ret = ad.grad(ret, i)
        return ret(*x)

    diff_fun.__doc__ = "Evaluates %s and its derivatives "\
                       "(use multiindex notation)\n\n"\
                       "Example: %s(x, derivatives=(2,0))\n"\
                       "         Returns the 2nd partial derivative of %s"\
                       " wrt. the first variable,\n"\
                       "         evaluated at x "\
                       " (which must be a numpy.ndarray)." % (fun.__name__,
                                                              fun.__name__,
                                                              fun.__name__)
    return diff_fun


##############################################################################
# Pretty printing, pretty useless.

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


##############################################################################

class Msg(object):
    """ A simple context manager to time chunks of code and 
    print status messages.
    
    Usage:
        with Msg("Doing amazing stuff"):
            do_amazing_stuff()
            do_more_of_it()
    
    Output will be:
        Doing amazing stuff... done in XX.XXs.
    
    Configure with either a second parameter or a global variable:
        Msg(message, level:int): will print messages only if level < DEBUG
                                 set level to a huge number to disable output 
        global DEBUG: set to 0 to disable output for *all* Msg()

    TODO: capture stdout/stderr and paste after the message,
          log it or whatnot
          
    """
    def __init__(self, msg:str, level:int=100):
        global DEBUG
        try:
            self.output = level < DEBUG
        except (NameError, TypeError): # TypeError because of unorderable types
            self.output = False
        self.msg = msg
        self.begin = time()

    def __enter__(self):
        if self.output:
            print(self.msg + "... ", end='')

    def __exit__(self, type, value, traceback):
        if self.output:
            print("done in %.3fs." % (time()-self.begin))