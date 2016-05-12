# Copyright (c) 2014-2016 Genome Research Ltd.
#
# This file is part of IVA.
#
# IVA is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.
import argparse
import os
import sys
import subprocess
version = '1.0.6'

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
