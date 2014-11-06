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
