import argparse
import os
import sys
version = '0.3.0'

class abspathAction(argparse.Action):
    def __call__(self, parser, namespace, value, option_string):
        '''Check if file exists, then convert to absolute path'''
        if not os.path.exists(value):
            print('Error! File "' + value + '" not found. Cannot continue', file=sys.stderr)
            sys.exit(1)

        setattr(namespace, self.dest, os.path.abspath(value))

