import os
import unittest
from pytest import raises
from datetime import datetime, timedelta
from paegan.transport.models.transport import Transport
from paegan.transport.exceptions import ModelError, BaseDataControllerError
from paegan.transport.model_controller import CachingModelController, BaseModelController
from shapely.geometry import Point
import logging
from paegan.logger.easy_logger import EasyLogger


class CachingModelControllerTest(unittest.TestCase):

    def setUp(self):
        self.log = EasyLogger('testlog.txt', level=logging.PROGRESS)
        self.log.logger.info("**************************************")
        self.log.logger.info(self.id())

        self.start_lat = 60.75
        self.start_lon = -147
        self.start_depth = 0
        self.num_particles = 4
        self.time_step = 3600
        self.num_steps = 10
        self.start_time = datetime(2014, 1, 2, 00)
        self.transport = Transport(horizDisp=0.05, vertDisp=0.0003)

        output_dir = "/data/lm/tests/output"
        self.output_path = os.path.join(output_dir, self.id())

        cache_dir = "/data/lm/tests/cache"
        self.cache_path = os.path.join(cache_dir, self.id())

        self.bathy_file = "/data/lm/bathy/ETOPO1_Bed_g_gmt4.grd"

        self.shoreline_path = "/data/lm/shore"

    def tearDown(self):
        self.log.logger.info("**************************************")
        self.log.close()

    def test_run_from_multiple_files_with_cache(self):
        models = [self.transport]

        p = Point(self.start_lon, self.start_lat)

        model = CachingModelController(geometry=p,
                                       depth=self.start_depth,
                                       start=self.start_time,
                                       step=self.time_step,
                                       nstep=self.num_steps,
                                       npart=self.num_particles,
                                       models=models,
                                       use_bathymetry=False,
                                       se_shoreline=False,
                                       time_chunk=10,
                                       horiz_chunk=4,
                                       bathy_path=self.bathy_file)

        particles = model.run("/data/lm/tests/pws_das_2014*.nc", output_formats = ['NetCDF'], output_path=self.output_path, cache_path=self.cache_path, remove_cache=False)
        self.assertEquals(len(particles), self.num_particles)
        self.assertTrue(os.path.exists(self.cache_path))
        os.remove(self.cache_path)
        self.assertTrue(os.path.exists(os.path.join(self.output_path, "trajectories.nc")))

    def test_run_from_multiple_files_without_cache(self):
        models = [self.transport]

        p = Point(self.start_lon, self.start_lat)

        model = BaseModelController(geometry=p,
                                    depth=self.start_depth,
                                    start=self.start_time,
                                    step=self.time_step,
                                    nstep=self.num_steps,
                                    npart=self.num_particles,
                                    models=models,
                                    use_bathymetry=False,
                                    use_shoreline=False)

        particles = model.run("/data/lm/tests/pws_das_2014*.nc", output_formats = ['NetCDF', 'trackline'], output_path=self.output_path)
        self.assertEquals(len(particles), self.num_particles)
        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))
        self.assertTrue(os.path.exists(os.path.join(self.output_path, "trajectories.nc")))

    def test_run_from_dap_with_cache(self):
        models = [self.transport]

        p = Point(self.start_lon, self.start_lat)

        model = CachingModelController(geometry=p,
                                       depth=self.start_depth,
                                       start=self.start_time,
                                       step=self.time_step,
                                       nstep=200,
                                       npart=self.num_particles,
                                       models=models,
                                       use_bathymetry=False,
                                       use_shoreline=False,
                                       time_chunk=24,
                                       horiz_chunk=4)

        particles = model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", output_formats = ['NetCDF'], output_path=self.output_path, cache_path=self.cache_path, remove_cache=False)
        self.assertEquals(len(particles), self.num_particles)
        self.assertTrue(os.path.exists(self.cache_path))
        os.remove(self.cache_path)
        self.assertTrue(os.path.exists(os.path.join(self.output_path, "trajectories.nc")))

    def test_run_from_dap_without_cache(self):
        models = [self.transport]

        p = Point(self.start_lon, self.start_lat)

        model = CachingModelController(geometry=p,
                                       depth=self.start_depth,
                                       start=self.start_time,
                                       step=self.time_step,
                                       nstep=self.num_steps,
                                       npart=self.num_particles,
                                       models=models,
                                       use_bathymetry=False,
                                       use_shoreline=False)

        particles = model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", output_formats = ['NetCDF'], output_path=self.output_path)
        self.assertEquals(len(particles), self.num_particles)
        # We didn't pass remove_cache=False, so it should have been removed by the CachingModelController.
        self.assertFalse(os.path.exists(self.cache_path))
        self.assertTrue(os.path.exists(os.path.join(self.output_path, "trajectories.nc")))

    def test_run_from_polygon(self):
        models = [self.transport]

        p = Point(self.start_lon, self.start_lat).buffer(0.001)

        model = BaseModelController(geometry=p,
                                    depth=self.start_depth,
                                    start=datetime(2014, 1, 2, 0),
                                    step=self.time_step,
                                    nstep=self.num_steps,
                                    npart=self.num_particles,
                                    models=models,
                                    use_bathymetry=False,
                                    use_shoreline=False,
                                    time_chunk=10,
                                    horiz_chunk=4)

        particles = model.run("/data/lm/tests/pws_das_2014*.nc", output_formats=['NetCDF'], output_path=self.output_path)
        self.assertEquals(len(particles), self.num_particles)
        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))
        self.assertTrue(os.path.exists(os.path.join(self.output_path, "trajectories.nc")))

    def test_run_from_point(self):
        models = [self.transport]

        p = Point(self.start_lon, self.start_lat)

        model = BaseModelController(geometry=p,
                                    depth=self.start_depth,
                                    start=self.start_time,
                                    step=self.time_step,
                                    nstep=self.num_steps,
                                    npart=self.num_particles,
                                    models=models,
                                    use_bathymetry=False,
                                    use_shoreline=False)

        particles = model.run("/data/lm/tests/pws_das_2014*.nc", output_formats = ['NetCDF'], output_path=self.output_path)
        self.assertEquals(len(particles), self.num_particles)
        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))
        self.assertTrue(os.path.exists(os.path.join(self.output_path, "trajectories.nc")))

    def test_run_from_point_with_wfs_shoreline(self):
        models = [self.transport]

        p = Point(self.start_lon, self.start_lat)

        model = BaseModelController(geometry=p,
                                    depth=self.start_depth,
                                    start=self.start_time,
                                    step=self.time_step,
                                    nstep=self.num_steps,
                                    npart=self.num_particles,
                                    models=models,
                                    use_bathymetry=False,
                                    use_shoreline=True,
                                    shoreline_path='http://geo.asascience.com/geoserver/shorelines/ows',
                                    shoreline_feature='shorelines:10m_land_polygons')

        particles = model.run("/data/lm/tests/pws_das_2014*.nc", output_formats = ['NetCDF'], output_path=self.output_path)
        self.assertEquals(len(particles), self.num_particles)
        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))
        self.assertTrue(os.path.exists(os.path.join(self.output_path, "trajectories.nc")))

    def test_time_method_interp(self):
        models = [self.transport]

        p = Point(self.start_lon, self.start_lat)

        model = BaseModelController(geometry=p,
                                    depth=self.start_depth,
                                    start=self.start_time,
                                    step=self.time_step,
                                    nstep=self.num_steps,
                                    npart=self.num_particles,
                                    models=models,
                                    use_bathymetry=False,
                                    use_shoreline=False,
                                    time_method="interp")

        particles = model.run("/data/lm/tests/pws_das_2014*.nc", output_formats = ['NetCDF'], output_path=self.output_path)
        self.assertEquals(len(particles), self.num_particles)
        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))
        self.assertTrue(os.path.exists(os.path.join(self.output_path, "trajectories.nc")))

    def test_time_method_nearest(self):
        models = [self.transport]

        p = Point(self.start_lon, self.start_lat)

        model = BaseModelController(geometry=p,
                                    depth=self.start_depth,
                                    start=self.start_time,
                                    step=self.time_step,
                                    nstep=self.num_steps,
                                    npart=self.num_particles,
                                    models=models,
                                    use_bathymetry=False,
                                    use_shoreline=False,
                                    time_method="nearest")

        particles = model.run("/data/lm/tests/pws_das_2014*.nc", output_formats = ['NetCDF'], output_path=self.output_path)
        self.assertEquals(len(particles), self.num_particles)
        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))
        self.assertTrue(os.path.exists(os.path.join(self.output_path, "trajectories.nc")))

    def test_time_method_bad(self):
        models = [self.transport]

        p = Point(self.start_lon, self.start_lat)

        with self.assertRaises(TypeError):
            BaseModelController(geometry=p,
                                depth=self.start_depth,
                                start=self.start_time,
                                step=self.time_step,
                                nstep=self.num_steps,
                                npart=self.num_particles,
                                models=models,
                                use_bathymetry=False,
                                use_shoreline=False,
                                time_method="umm_what_am_i")

    def test_start_on_land_from_lat_lon(self):
        # Set the start position and time for the models
        start_lat = 60.15551950079041
        start_lon = -148.1999130249019

        models = [self.transport]

        model = BaseModelController(latitude=start_lat,
                                    longitude=start_lon,
                                    depth=self.start_depth,
                                    start=self.start_time,
                                    step=self.time_step,
                                    nstep=self.num_steps,
                                    npart=self.num_particles,
                                    models=models,
                                    use_bathymetry=False,
                                    use_shoreline=True)

        with raises(ModelError):
            model.run("/data/lm/tests/pws_das_2014*.nc")

    def test_start_on_land_from_point_no_depth(self):
        # Set the start position and time for the models
        start_lat = 60.15551950079041
        start_lon = -148.1999130249019

        p = Point(start_lon, start_lat)

        models = [self.transport]

        model = BaseModelController(geometry=p,
                                    start=self.start_time,
                                    step=self.time_step,
                                    nstep=self.num_steps,
                                    npart=self.num_particles,
                                    models=models,
                                    use_bathymetry=False,
                                    use_shoreline=True)

        with raises(ModelError):
            model.run("/data/lm/tests/pws_das_2014*.nc")

    def test_start_on_land_from_point_with_depth(self):
        # Set the start position and time for the models
        start_lat = 60.15551950079041
        start_lon = -148.1999130249019
        depth     = -10

        p = Point(start_lon, start_lat)

        models = [self.transport]

        model = BaseModelController(geometry=p,
                                    depth=depth,
                                    start=self.start_time,
                                    step=self.time_step,
                                    nstep=self.num_steps,
                                    npart=self.num_particles,
                                    models=models,
                                    use_bathymetry=False,
                                    use_shoreline=True)

        with raises(ModelError):
            model.run("/data/lm/tests/pws_das_2014*.nc")

    def test_bad_dataset(self):
        models = [self.transport]

        model = BaseModelController(latitude=self.start_lat,
                                    longitude=self.start_lon,
                                    depth=self.start_depth,
                                    start=self.start_time,
                                    step=self.time_step,
                                    nstep=self.num_steps,
                                    npart=self.num_particles,
                                    models=models,
                                    use_bathymetry=False,
                                    use_shoreline=False)

        with raises(BaseDataControllerError):
            model.run("http://example.com/thisisnotadataset.nc")

    def test_timechunk_greater_than_timestep(self):
        models = [self.transport]

        model = CachingModelController(latitude=self.start_lat,
                                       longitude=self.start_lon,
                                       depth=self.start_depth,
                                       start=self.start_time,
                                       step=self.time_step,
                                       nstep=self.num_steps,
                                       npart=self.num_particles,
                                       models=models,
                                       use_bathymetry=False,
                                       use_shoreline=False,
                                       time_chunk=48,
                                       horiz_chunk=2)

        particles = model.run("/data/lm/tests/pws_das_2014*.nc", output_formats = ['NetCDF'], output_path=self.output_path, cache_path=self.cache_path, remove_cache=False)
        self.assertEquals(len(particles), self.num_particles)
        self.assertTrue(os.path.exists(self.cache_path))
        os.remove(self.cache_path)
        self.assertTrue(os.path.exists(os.path.join(self.output_path, "trajectories.nc")))

    def test_no_local_data_for_requested_run(self):
        models = [self.transport]

        # Start is after available time
        model = BaseModelController(latitude=self.start_lat,
                                    longitude=self.start_lon,
                                    depth=self.start_depth,
                                    start=datetime.utcnow() + timedelta(days=30),
                                    step=self.time_step,
                                    nstep=self.num_steps,
                                    npart=self.num_particles,
                                    models=models,
                                    use_bathymetry=False,
                                    use_shoreline=False)

        with self.assertRaises(BaseDataControllerError):
            model.run("/data/lm/tests/pws_das_2014*.nc", output_formats = ['NetCDF'], output_path=self.output_path, cache_path=self.cache_path, remove_cache=False)

        # Start is OK but Ending is after available time
        model = BaseModelController(latitude=self.start_lat,
                                    longitude=self.start_lon,
                                    depth=self.start_depth,
                                    start=datetime(2014, 1, 1, 0),
                                    step=self.time_step,
                                    nstep=500,
                                    npart=self.num_particles,
                                    models=models,
                                    use_bathymetry=False,
                                    use_shoreline=False)

        with self.assertRaises(BaseDataControllerError):
            model.run("/data/lm/tests/pws_das_2014*.nc", output_formats = ['NetCDF'], output_path=self.output_path, cache_path=self.cache_path, remove_cache=False)

    def test_no_dap_data_for_requested_run(self):
        models = [self.transport]

        # Start is after available time
        model = CachingModelController(latitude=self.start_lat,
                                       longitude=self.start_lon,
                                       depth=self.start_depth,
                                       start=datetime.utcnow() + timedelta(days=30),
                                       step=self.time_step,
                                       nstep=self.num_steps,
                                       npart=self.num_particles,
                                       models=models,
                                       use_bathymetry=False,
                                       use_shoreline=False)

        with self.assertRaises(BaseDataControllerError):
            model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", output_formats = ['NetCDF'], output_path=self.output_path, cache_path=self.cache_path, remove_cache=False)

        # Start is OK but Ending is after available time
        model = CachingModelController(latitude=self.start_lat,
                                       longitude=self.start_lon,
                                       depth=self.start_depth,
                                       start=datetime.utcnow() - timedelta(days=2),
                                       step=self.time_step,
                                       nstep=500,
                                       npart=self.num_particles,
                                       models=models,
                                       use_bathymetry=False,
                                       use_shoreline=False)

        with self.assertRaises(BaseDataControllerError):
            model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", output_formats = ['NetCDF'], output_path=self.output_path, cache_path=self.cache_path, remove_cache=False)

    def test_run_10m_shoreline(self):
        models = [self.transport]

        p = Point(self.start_lon, self.start_lat)

        model = BaseModelController(geometry=p,
                                    depth=self.start_depth,
                                    start=self.start_time,
                                    step=self.time_step,
                                    nstep=300,
                                    npart=2,
                                    models=models,
                                    use_bathymetry=False,
                                    use_shoreline=True,
                                    shoreline_index_buffer=0.05)

        particles = model.run("/data/lm/tests/pws_das_2014*.nc")
        self.assertEquals(len(particles), 2)
        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))

    def test_run_west_coast_shoreline(self):

        models = [self.transport]

        p = Point(self.start_lon, self.start_lat)

        model = BaseModelController(geometry=p,
                                    depth=self.start_depth,
                                    start=self.start_time,
                                    step=self.time_step,
                                    nstep=300,
                                    npart=2,
                                    models=models,
                                    use_bathymetry=False,
                                    use_shoreline=True,
                                    shoreline_path='/data/lm/shore/westcoast/New_Land_Clean.shp',
                                    shoreline_index_buffer=0.05)

        particles = model.run("/data/lm/tests/pws_das_2014*.nc")
        self.assertEquals(len(particles), 2)
        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))
