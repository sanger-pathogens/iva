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
import zipfile
import tempfile
import shutil
import os

class Error (Exception): pass


class Extractor:
    def __init__(self, egg):
        self.egg = os.path.abspath(egg)
        if not os.path.exists(self.egg):
            raise Error('Egg not found: ' + self.egg)

        if os.path.isdir(self.egg):
            self.zip_file = None
        elif os.path.isfile(self.egg):
            try:
                self.zip_file = zipfile.ZipFile(self.egg)
            except:
                raise Error('Error opening zip egg file: ' + self.egg)
            self.zip_filenames = set(self.zip_file.namelist())
        else:
            raise Error('Egg not recognised as dir or file: ' + self.egg)


    def _copy_file_unzipped(self, infile, outfile):
        to_copy = os.path.join(self.egg, infile)
        try:
            shutil.copyfile(to_copy, outfile)
        except:
            raise Error('Error copying file:' + to_copy + ' -> ' + outfile)


    def _copy_file_zipped(self, infile, outfile):
        if infile not in self.zip_filenames:
            raise Error('Error copying file. Not found. ' + infile)
        tmpdir = tempfile.mkdtemp(prefix='tmp.copy_file_zipped.', dir=os.getcwd())
        self.zip_file.extract(infile, path=tmpdir) 
        os.rename(os.path.join(tmpdir, infile), outfile)
        shutil.rmtree(tmpdir)


    def copy_file(self, infile, outfile):
        if self.zip_file is None:
            self._copy_file_unzipped(infile, outfile)
        else:
            self._copy_file_zipped(infile, outfile)


    def _copy_dir_unzipped(self, indir, outdir):
        to_copy = os.path.join(self.egg, indir)
        try:
            shutil.copytree(to_copy, outdir)
        except:
            raise Error('Error copying directory ', indir)


    def _copy_dir_zipped(self, indir, outdir):
        tmpdir = tempfile.mkdtemp(prefix='tmp.copy_dir_zipped', dir=os.getcwd())
        if not indir.endswith('/'):
            indir += '/'

        files = [x for x in self.zip_filenames if x.startswith(indir)]

        for f in files:
            self.zip_file.extract(f, path=tmpdir)

        os.rename(os.path.join(tmpdir, indir), outdir)
        shutil.rmtree(tmpdir)


    def copy_dir(self, indir, outdir):
        if self.zip_file is None:
            self._copy_dir_unzipped(indir, outdir)
        else:
            self._copy_dir_zipped(indir, outdir)
