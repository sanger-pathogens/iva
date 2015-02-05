import unittest
import sys
import shutil
import os
import filecmp
from iva import egg_extract

modules_dir = os.path.dirname(os.path.abspath(egg_extract.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestEggExtract(unittest.TestCase):
    def setUp(self):
        self.zipped = os.path.join(data_dir, 'egg_extract_test_dir.zip')
        self.not_zipped = os.path.join(data_dir, 'egg_extract_test_dir')
        self.egg_zipped = egg_extract.Extractor(self.zipped)
        self.egg_not_zipped = egg_extract.Extractor(self.not_zipped)


    def test_init_and_write(self):
        '''test __init__'''
        with self.assertRaises(egg_extract.Error):
            egg_extract.Extractor('notafilethrowerror')

        self.assertTrue(self.egg_zipped.zip_file is not None)
        self.assertTrue(self.egg_not_zipped.zip_file is None)


    def test_copy_file_unzipped(self):
        '''test _copy_file_unzipped'''
        to_copy = os.path.join('iva', 'gage', 'gage1')
        outfile = 'tmp.copy_file_unzipped'
        if os.path.exists(outfile):
            os.unlink(outfile)
        self.egg_not_zipped._copy_file_unzipped(to_copy, outfile)
        self.assertTrue(os.path.exists(outfile))
        os.unlink(outfile)


    def test_copy_file_zipped(self):
        '''test _copy_file_zipped'''
        outfile = 'tmp.copy_file_zipped'
        if os.path.exists(outfile):
            os.unlink(outfile)
        to_copy = os.path.join('iva', 'gage', 'gage1')
        self.egg_zipped._copy_file_zipped(to_copy, outfile)
        self.assertTrue(os.path.exists(outfile))
        os.unlink(outfile)


    def test_copy_file(self):
        '''test copy_file'''
        outfile = 'tmp.copy_file'
        for egg in [self.egg_zipped, self.egg_not_zipped]:
            if os.path.exists(outfile):
                os.unlink(outfile)
            to_copy = os.path.join('iva', 'gage', 'gage1')
            egg.copy_file(to_copy, outfile)
            self.assertTrue(os.path.exists(outfile))
            os.unlink(outfile)


    def test_copy_dir_unzipped(self):
        '''test copy_dir_unzipped'''
        outdir = 'tmp.copy_dir_unzipped'
        to_copy = os.path.join('iva', 'gage')
        if os.path.exists(outdir):
            shutil.rmtree(outdir)

        self.egg_not_zipped._copy_dir_unzipped(to_copy, outdir)
        self.assertTrue(os.path.exists(outdir))
        self.assertTrue(os.path.isdir(outdir))
        self.assertTrue(os.path.exists(os.path.join(outdir, 'gage1')))
        self.assertTrue(os.path.exists(os.path.join(outdir, 'gage2')))
        shutil.rmtree(outdir)


    def test_copy_dir_zipped(self):
        '''test copy_dir_zipped'''
        outdir = 'tmp.copy_dir_zipped'
        to_copy = os.path.join('iva', 'gage')
        if os.path.exists(outdir):
            shutil.rmtree(outdir)

        self.egg_zipped._copy_dir_zipped(to_copy, outdir)
        self.assertTrue(os.path.exists(outdir))
        self.assertTrue(os.path.isdir(outdir))
        self.assertTrue(os.path.exists(os.path.join(outdir, 'gage1')))
        self.assertTrue(os.path.exists(os.path.join(outdir, 'gage2')))
        shutil.rmtree(outdir)


    def test_copy_dir(self):
        '''test copy_dir'''
        outdir = 'tmp.copy_dir'
        for egg in [self.egg_zipped, self.egg_not_zipped]:
            if os.path.exists(outdir):
                shutil.rmtree(outdir)
            to_copy = os.path.join('iva', 'gage')
            egg.copy_dir(to_copy, outdir)
            self.assertTrue(os.path.exists(outdir))
            self.assertTrue(os.path.isdir(outdir))
            self.assertTrue(os.path.exists(os.path.join(outdir, 'gage1')))
            self.assertTrue(os.path.exists(os.path.join(outdir, 'gage2')))
            shutil.rmtree(outdir)

