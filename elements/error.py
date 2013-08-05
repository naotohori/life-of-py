#!/usr/bin/env python
'''
@author: Naoto Hori
'''

class MyError(Exception) :
    def __init__(self, classname='unknown', funcname='unknown', title='unknown'):
        self._class = classname
        self._func = funcname
        self._title = title
        
    def show(self):
        print 'class:', self._class
        print 'function:', self._func
        print 'matter:', self._title
