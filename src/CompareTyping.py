#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import os

def compare_typing( typing_files ):
    files = sort_files( typing_files )
    print files

def sort_files( typing_files ):
    files = {}
    for filename in typing_files:
        filepath = os.path.abspath( filename )
        parts = filepath.split('/')
        if len(parts) > 1:
            name = parts[-2]
            files[name] = filepath
    return files

if __name__ == '__main__':
    import sys

    typing_files = sys.argv[1:]

    compare_typing( typing_files )