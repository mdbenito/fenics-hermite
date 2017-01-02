
# coding: utf-8

# In[ ]:

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

