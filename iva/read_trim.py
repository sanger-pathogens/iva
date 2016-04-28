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
import os
from iva import common

def run_trimmomatic(reads1, reads2, outprefix, trimmo_jar, adapters, minlen=50, verbose=0, threads=1, qual_trim=''):
    cmd = ' '.join([
        'java -Xmx1000m -jar',
        trimmo_jar,
        'PE',
        '-threads', str(threads),
        reads1,
        reads2,
        outprefix + '_1.fq',
        outprefix + '.unpaired_1.fq',
        outprefix + '_2.fq',
        outprefix + '.unpaired_2.fq',
        'ILLUMINACLIP:' + os.path.abspath(adapters) + ':2:10:7:1',
        qual_trim,
        'MINLEN:' + str(minlen)
    ])

    if verbose:
        print('Run trimmomatic:', cmd)
    common.syscall(cmd)
    os.unlink(outprefix + '.unpaired_1.fq')
    os.unlink(outprefix + '.unpaired_2.fq')
