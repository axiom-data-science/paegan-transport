import unittest
import os
import warnings
import paegan.transport.export as ex

try:
    import tables
except:
    warnings.warn("PyTables(HDF5) export tests will fail without Pytables installed (`pip install tables`)")


class ExportTest(unittest.TestCase):

    def setUp(self):
        self.output_path = "/data/lm/tests/output"

    def test_h5_trackline(self):
        ex.H5Trackline.export(folder=self.output_path, h5_file=os.path.normpath(os.path.join(os.path.dirname(__file__), "./resources/files/hdf5_output.h5")))

    def test_h5_particle_trackline(self):
        ex.H5ParticleTracklines.export(folder=self.output_path, h5_file=os.path.normpath(os.path.join(os.path.dirname(__file__), "./resources/files/hdf5_output.h5")))

    def test_h5_shapefile(self):
        ex.H5GDALShapefile.export(folder=self.output_path, h5_file=os.path.normpath(os.path.join(os.path.dirname(__file__), "./resources/files/hdf5_output.h5")))
