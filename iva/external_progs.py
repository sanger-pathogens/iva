import shutil
import subprocess
import re

class Error (Exception): pass


def is_in_path(prog):
    return shutil.which(prog) is not None


prog_to_version_cmd = {
    'kmc': ('kmc', re.compile('^K-Mer Counter \(KMC\) ver\. (.*) \(.*\)$')),
    'kmc_dump': ('kmc_dump', re.compile('^KMC dump ver. (.*) \(.*\)$')),
    'nucmer': ('nucmer --version', re.compile('^NUCmer \(NUCleotide MUMmer\) version (.*)$')),
    'smalt': ('smalt version', re.compile('^Version: (.*)$')),
    'samtools': ('samtools', re.compile('^Version: (.*)$')),
}


def get_version(prog):
    assert prog in prog_to_version_cmd
    assert is_in_path(prog)
    cmd, regex = prog_to_version_cmd[prog]
    cmd_output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    cmd_output = cmd_output[0].decode().split('\n')[:-1] + cmd_output[1].decode().split('\n')[:-1]
    for line in cmd_output:
        hits = regex.search(line)
        if hits:
            return hits.group(1)
    return 'UNKNOWN ...\n    I tried running this to get the version: "' + cmd + '"\n    and the output didn\'t match this regular expression: "' + regex.pattern + '"'


def print_all_versions():
    for prog in sorted(prog_to_version_cmd):
        version = get_version(prog)
        print('Using', prog, 'version', version)

