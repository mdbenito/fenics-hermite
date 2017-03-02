# coding: utf-8

from functools import reduce

def fnand(*fs):
    def res(*args):
        return reduce(lambda x,y: x and y, [f(*args) for f in fs])
    return res

def fnor(*fs):
    def res(*args):
        return reduce(lambda x,y: x or y, [f(*args) for f in fs])
    return res

def fnnot(f):
    def res(*args):
        return not f(*args)
    return res

def one_arg(f):
    def res(x):
        return f(x)
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