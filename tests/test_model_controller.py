import os
import unittest
from pytest import raises
from datetime import datetime, timedelta
from paegan.transport.models.transport import Transport
from paegan.transport.exceptions import ModelError, BaseDataControllerError
from paegan.transport.controllers import CachingModelController, BaseModelController, IPythonClusterModelController, DistributedModelController
from shapely.geometry import Point
import logging
from paegan.logger.easy_logger import EasyLogger
import paegan.transport.export as ex

import fiona
import gj2ascii
from shapely.geometry import shape


class ModelControllerTest(unittest.TestCase):

    def setUp(self):
        self.log = EasyLogger('testlog.txt', level=logging.PROGRESS)
        self.log.logger.handlers.append(logging.StreamHandler())
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
        test_dir = self.id().split('.')[-1]
        self.output_path = os.path.join(output_dir, test_dir)
        if not os.path.isdir(self.output_path):
            os.makedirs(self.output_path)
        self.output_formats = [ex.H5Trackline]

        cache_dir = "/data/lm/tests/cache"
        self.cache_path = os.path.join(cache_dir, test_dir, 'cache.nc')

        self.bathy_file = "/data/lm/bathy/global/ETOPO1_Bed_g_gmt4.grd"

        self.shoreline_path = "/data/lm/shore"

    def tearDown(self):
        self.log.close()

    def draw_trackline(self, geojson_file):
        with fiona.open(geojson_file) as src:
            self.log.logger.info((gj2ascii.render(src, 100, fill='.', char='o', bbox=shape(src.next()['geometry']).buffer(0.1).envelope.bounds)))

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
                                       horiz_chunk=4)

        model.setup_run("/data/lm/tests/pws_das_2014*.nc", cache_path=self.cache_path, remove_cache=False)
        model.run(output_formats=self.output_formats, output_path=self.output_path)

        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))
        self.assertTrue(os.path.exists(self.cache_path))
        os.remove(self.cache_path)

    def test_run_from_multiple_files_without_cache_on_ipy_cluster(self):
        try:
            from IPython.parallel import Client
            client = Client()

            pool = client.load_balanced_view()
        except:
            raise unittest.SkipTest("Cluster connection failed")

        models = [self.transport]
        p = Point(self.start_lon, self.start_lat)
        model = IPythonClusterModelController(geometry=p,
                                              depth=self.start_depth,
                                              start=self.start_time,
                                              step=self.time_step,
                                              nstep=self.num_steps,
                                              npart=self.num_particles,
                                              models=models,
                                              use_bathymetry=False,
                                              use_shoreline=False,
                                              pool=pool)

        model.setup_run("/data/lm/tests/pws_das_2014*.nc")
        model.run(output_formats=self.output_formats, output_path=self.output_path)

        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))
        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))

    def test_run_from_multiple_files_without_cache(self, pool=None):
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
                                    pool=pool)

        model.setup_run("/data/lm/tests/pws_das_2014*.nc")
        model.run(output_formats=self.output_formats, output_path=self.output_path)

        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))

        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))

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

        model.setup_run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", cache_path=self.cache_path, remove_cache=False)
        model.run(output_formats=self.output_formats, output_path=self.output_path)

        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))

        self.assertTrue(os.path.exists(self.cache_path))
        os.remove(self.cache_path)

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

        model.setup_run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc")
        model.run(output_formats=self.output_formats, output_path=self.output_path)

        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))

        # We didn't pass remove_cache=False, so it should have been removed by the CachingModelController.
        self.assertFalse(os.path.exists(self.cache_path))

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

        model.setup_run("/data/lm/tests/pws_das_2014*.nc")
        model.run(output_formats=self.output_formats, output_path=self.output_path)

        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))

        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))

    def test_distributed_from_polygon(self):
        models = [self.transport]
        p = Point(self.start_lon, self.start_lat).buffer(0.001)
        model = DistributedModelController(geometry=p,
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

        model.setup_run("/data/lm/tests/pws_das_2014*.nc", redis_url='redis://127.0.0.1:6379/150', redis_results_channel='test_distributed_from_polygon:results')
        model.run(output_formats=self.output_formats, output_path=self.output_path)

        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))

        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))

    @unittest.skip("Need to find a new WFS server to test against")
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

        model.setup_run("/data/lm/tests/pws_das_2014*.nc")
        model.run(output_formats=self.output_formats, output_path=self.output_path)

        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))

        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))

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

        model.setup_run("/data/lm/tests/pws_das_2014*.nc")
        model.run(output_formats=self.output_formats, output_path=self.output_path)

        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))

        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))

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

        model.setup_run("/data/lm/tests/pws_das_2014*.nc")
        model.run(output_formats=self.output_formats, output_path=self.output_path)

        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))

        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))

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
            model.setup_run("/data/lm/tests/pws_das_2014*.nc")
            model.run()

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
            model.setup_run("/data/lm/tests/pws_das_2014*.nc")
            model.run()

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
            model.setup_run("/data/lm/tests/pws_das_2014*.nc")
            model.run()

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
            model.setup_run("http://example.com/thisisnotadataset.nc")
            model.run()

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

        model.setup_run("/data/lm/tests/pws_das_2014*.nc", cache_path=self.cache_path, remove_cache=False)
        model.run(output_formats=self.output_formats, output_path=self.output_path)

        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))

        self.assertTrue(os.path.exists(self.cache_path))
        os.remove(self.cache_path)

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
            model.setup_run("/data/lm/tests/pws_das_2014*.nc", cache_path=self.cache_path, remove_cache=False)
            model.run(output_formats=self.output_formats, output_path=self.output_path)

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
            model.setup_run("/data/lm/tests/pws_das_2014*.nc", cache_path=self.cache_path, remove_cache=False)
            model.run(output_formats=self.output_formats, output_path=self.output_path)

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
            model.setup_run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", cache_path=self.cache_path, remove_cache=False)
            model.run(output_formats=self.output_formats, output_path=self.output_path)

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
            model.setup_run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", cache_path=self.cache_path, remove_cache=False)
            model.run(output_formats=self.output_formats, output_path=self.output_path)

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

        model.setup_run("/data/lm/tests/pws_das_2014*.nc")
        model.run()

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

        model.setup_run("/data/lm/tests/pws_das_2014*.nc")
        model.run()

        # Not a caching controller, no cache path should exist
        self.assertFalse(os.path.exists(self.cache_path))
