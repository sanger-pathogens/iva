import shutil
import subprocess
import re

class Error (Exception): pass


def is_in_path(prog):
    return shutil.which(prog) is not None


prog_to_version_cmd = {
    'kmc': ('kmc', re.compile('^K-Mer Counter \(KMC\) ver\. (.*) \(.*\)$')),
    'kmc_dump': ('kmc_dump', re.compile('^KMC dump ver. (.*) \(.*\)$')),
    'kraken': ('kraken --version', re.compile('^Kraken version (.*)$')),
    'kraken-build': ('kraken-build --version', re.compile('^Kraken version (.*)$')),
    'nucmer': ('nucmer --version', re.compile('^NUCmer \(NUCleotide MUMmer\) version (.*)$')),
    'R': ('R --version', re.compile('^R version (.*) \(.*\) --')),
    'reapr': ('reapr', re.compile('^REAPR version: (.*)$')),
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
    'kraken',
    'nucmer',
    'R',
    'smalt',
    'samtools',
]


qc_progs_optional = [
    'kraken',
    'kraken-build',
    'reapr',
]


qc_make_db_progs = [
    'kraken',
    'kraken-build',
]


def get_version(prog, must_be_in_path=True):
    assert prog in prog_to_version_cmd
    if not is_in_path(prog):
        if must_be_in_path:
            raise Error('Error getting version of', prog, '- not found in path.')
        else:
            return 'UNKNOWN - not in path'

    cmd, regex = prog_to_version_cmd[prog]
    cmd_output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    cmd_output = cmd_output[0].decode().split('\n')[:-1] + cmd_output[1].decode().split('\n')[:-1]
    for line in cmd_output:
        hits = regex.search(line)
        if hits:
            return hits.group(1)
    return 'UNKNOWN ...\n    I tried running this to get the version: "' + cmd + '"\n    and the output didn\'t match this regular expression: "' + regex.pattern + '"'


def print_all_versions(progs, must_be_in_path=True):
    for prog in sorted(progs):
        version = get_version(prog, must_be_in_path=must_be_in_path)
        print('Using', prog, 'version', version)

