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
