import shutil
import subprocess
import re
import sys
import pyfastaq
from iva import common

class Error (Exception): pass


def is_in_path(prog):
    return shutil.which(prog) is not None


prog_to_version_cmd = {
    'blastn': ('blastn -version', re.compile('^blastn: (.*)$')),
    'makeblastdb': ('makeblastdb -version', re.compile('makeblastdb: (.*)$')),
    'kmc': ('kmc', re.compile('^K-Mer Counter \(KMC\) ver\. (.*) \(.*\)$')),
    'kmc_dump': ('kmc_dump', re.compile('^KMC dump ver. (.*) \(.*\)$')),
    'kraken': ('kraken --version', re.compile('^Kraken version (.*)$')),
    'kraken-build': ('kraken-build --version', re.compile('^Kraken version (.*)$')),
    'nucmer': ('nucmer --version', re.compile('^NUCmer \(NUCleotide MUMmer\) version (.*)$')),
    'R': ('R --version', re.compile('^R version (.*) \(.*\) --')),
    'smalt': ('smalt version', re.compile('^Version: (.*)$')),
    'samtools': ('samtools', re.compile('^Version: (.*)$')),
}


assembly_progs = [
    'kmc',
    'kmc_dump',
    'nucmer',
    'smalt',
    'samtools',
]


qc_progs = [
    'nucmer',
    'R',
    'smalt',
    'samtools',
]


qc_progs_optional = [
    'blastn',
    'makeblastdb',
    'kraken',
    'kraken-build',
]


qc_make_db_progs = [
    'kraken',
    'kraken-build',
]


def get_version(prog, must_be_in_path=True):
    assert prog in prog_to_version_cmd
    if not is_in_path(prog):
        if must_be_in_path:
            raise Error('Error getting version of ' + prog + ' - not found in path.')
        else:
            return 'UNKNOWN - not in path'

    cmd, regex = prog_to_version_cmd[prog]
    cmd_output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    cmd_output = common.decode(cmd_output[0]).split('\n')[:-1] + common.decode(cmd_output[1]).split('\n')[:-1]
    for line in cmd_output:
        hits = regex.search(line)
        if hits:
            return hits.group(1)
    return 'UNKNOWN ...\n    I tried running this to get the version: "' + cmd + '"\n    and the output didn\'t match this regular expression: "' + regex.pattern + '"'


def get_all_versions(progs, must_be_in_path=True):
    info = []
    for prog in sorted(progs):
        version = get_version(prog, must_be_in_path=must_be_in_path)
        info.append(' '.join(['Using', prog, 'version', version]))
    return info


def write_prog_info(script, filename):
    if script == 'iva':
        required = assembly_progs
        optional = None
    elif script == 'iva_qc':
        required = qc_progs
        optional = qc_progs_optional
    elif script == 'iva_qc_make_db':
        required = qc_make_db_progs
        optional = None
    else:
        raise Error('Script ' + script + ' not recognised')

    f = pyfastaq.utils.open_file_write(filename)
    print(' '.join(sys.argv), file=f)
    print('IVA version', common.version, file=f)
    print('\n'.join(get_all_versions(required)), file=f)
    if optional is not None:
        print('\n'.join(get_all_versions(optional, must_be_in_path=False)), file=f)
    pyfastaq.utils.close(f)

