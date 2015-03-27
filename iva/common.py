import argparse
import os
import sys
import subprocess
version = '0.11.5'

class abspathAction(argparse.Action):
    def __call__(self, parser, namespace, value, option_string):
        '''Check if file exists, then convert to absolute path'''
        if not os.path.exists(value):
            print('Error! File "' + value + '" not found. Cannot continue', file=sys.stderr)
            sys.exit(1)

        setattr(namespace, self.dest, os.path.abspath(value))


def syscall(cmd, allow_fail=False, verbose=False):
    if verbose:
        print('syscall:', cmd, flush=True)
    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        if allow_fail:
            return False
        else:
            print('The following command failed with exit code', error.returncode, file=sys.stderr)
            print(cmd, file=sys.stderr)
            print('\nThe output was:\n', file=sys.stderr)
            print(error.output.decode(), file=sys.stderr)
            sys.exit(1)

    return True


def decode(x):
    try:
        s = x.decode()
    except:
        return x
    return s
