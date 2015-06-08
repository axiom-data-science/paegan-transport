import unittest
from paegan.transport.models.transport import Transport
import os
import json
import urllib
from datetime import datetime
import logging
from paegan.logger.easy_logger import EasyLogger
import paegan.transport.export as ex
from paegan.transport.models.behavior import LarvaBehavior
from paegan.transport.controllers import CachingModelController
import pytz

import fiona
import gj2ascii
from shapely.geometry import shape


class LarvalBehaviourTest(unittest.TestCase):

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
        self.output_formats = [ex.H5GDALShapefile, ex.H5Trackline]

        cache_dir = "/data/lm/tests/cache"
        self.cache_path = os.path.join(cache_dir, test_dir, 'cache.nc')

        self.bathy_file = "/data/lm/bathy/global/ETOPO1_Bed_g_gmt4.grd"

        self.shoreline_path = "/data/lm/shore"

    def tearDown(self):
        self.log.close()

    def draw_trackline(self, geojson_file):
        with fiona.open(geojson_file) as src:
            self.log.logger.info((gj2ascii.render(src, 100, fill='.', char='o', bbox=shape(src.next()['geometry']).buffer(0.1).envelope.bounds)))

    def test_behavior_growth_and_settlement(self):
        # 6 days
        num_steps = 144
        num_particles = 2

        # Behavior
        behavior_config = open(os.path.normpath(os.path.join(os.path.dirname(__file__), "./resources/files/behavior_for_run_testing.json"))).read()
        lb = LarvaBehavior(json=behavior_config)

        models = [self.transport]
        models.append(lb)

        model = CachingModelController(
            latitude=self.start_lat,
            longitude=self.start_lon,
            depth=self.start_depth,
            start=self.start_time,
            step=self.time_step,
            nstep=num_steps,
            npart=num_particles,
            models=models,
            use_bathymetry=True,
            use_shoreline=True,
            time_chunk=24,
            horiz_chunk=4,
            time_method='nearest',
            bathy_path=self.bathy_file)

        model.setup_run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", cache=self.cache_path)
        model.run(output_formats=self.output_formats, output_path=self.output_path)
        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))

    def test_quick_settlement(self):
        num_steps = 24
        num_particles = 4

        # Behavior
        behavior_config = open(os.path.normpath(os.path.join(os.path.dirname(__file__), "./resources/files/behavior_quick_settle.json"))).read()
        lb = LarvaBehavior(json=behavior_config)

        models = [self.transport]
        models.append(lb)

        model = CachingModelController(
            latitude=self.start_lat,
            longitude=self.start_lon,
            depth=self.start_depth,
            start=self.start_time,
            step=self.time_step,
            nstep=num_steps,
            npart=num_particles,
            models=models,
            use_bathymetry=True,
            use_shoreline=True,
            time_chunk=12,
            horiz_chunk=2,
            time_method='nearest',
            bathy_path=self.bathy_file)

        model.setup_run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_DAS.nc", cache=self.cache_path)
        model.run(output_path=self.output_path, output_formats=self.output_formats)
        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))

    @unittest.skip("Lost behavior file")
    def test_kayak_island(self):
        # 6 days
        num_steps = 1632
        num_particles = 100
        time_step = 3600

        behavior_config = json.loads(urllib.urlopen("http://behaviors.larvamap.asascience.com/library/50ef1bb1cc7b61000700001d.json").read())
        lb = LarvaBehavior(data=behavior_config[u'results'][0])

        models = [Transport(horizDisp=0.01, vertDisp=0.001)]
        models.append(lb)

        start_time = datetime(2011, 5, 2, 00, tzinfo=pytz.utc)
        start_lat = 59.93517413488866
        start_lon = -144.496213677788
        depth = -1

        shoreline_path = os.path.join(self.shoreline_path, "westcoast", "New_Land_Clean.shp")

        model = CachingModelController(
            latitude=start_lat,
            longitude=start_lon,
            depth=depth,
            start=start_time,
            step=time_step,
            nstep=num_steps,
            npart=num_particles,
            models=models,
            use_bathymetry=True,
            use_shoreline=True,
            time_chunk=24,
            horiz_chunk=5,
            time_method='interp',
            shoreline_path=shoreline_path,
            bathy_path=self.bathy_file)

        model.setup_run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L1_FCST.nc", cache=self.cache_path)
        model.run(output_path=self.output_path, output_formats=self.output_formats)
        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))

    @unittest.skip("Lost behavior file")
    def test_sheep_bay(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: Sheep Bay")

        # 6 days
        num_steps = 1632

        num_particles = 100

        time_step = 3600

        behavior_config = json.loads(urllib.urlopen("http://behaviors.larvamap.asascience.com/library/50ef1bb1cc7b61000700001d.json").read())
        lb = LarvaBehavior(data=behavior_config[u'results'][0])

        models = [Transport(horizDisp=0.01, vertDisp=0.001)]
        models.append(lb)

        start_time = datetime(2011, 5, 2, 00, tzinfo=pytz.utc)

        start_lat = 60.60899655733162
        start_lon = -145.97402533055956
        depth = -1

        shoreline_path = os.path.join(self.shoreline_path, "westcoast", "New_Land_Clean.shp")

        model = CachingModelController(
            latitude=start_lat,
            longitude=start_lon,
            depth=depth,
            start=start_time,
            step=time_step,
            nstep=num_steps,
            npart=num_particles,
            models=models,
            use_bathymetry=True,
            use_shoreline=True,
            time_chunk=24,
            horiz_chunk=5,
            time_method='interp',
            shoreline_path=shoreline_path,
            bathy_path=self.bathy_file)

        model.setup_run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", cache=self.cache_path)
        model.run(output_path=self.output_path, output_formats=self.output_formats)
        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))

    def test_diel_migration(self):
        num_steps = 168
        num_particles = 4
        start_time = datetime(2013, 4, 1, 0, tzinfo=pytz.utc)

        # Behavior
        behavior_config = open(os.path.normpath(os.path.join(os.path.dirname(__file__), "./resources/files/diel_suncycles.json"))).read()
        lb = LarvaBehavior(json=behavior_config)

        models = [self.transport]
        models.append(lb)

        model = CachingModelController(
            latitude=60.68,
            longitude=-146.42,
            depth=self.start_depth,
            start=start_time,
            step=self.time_step,
            nstep=num_steps,
            npart=num_particles,
            models=models,
            use_bathymetry=True,
            use_shoreline=True,
            time_chunk=24,
            horiz_chunk=2,
            time_method='nearest',
            bathy_path=self.bathy_file)

        model.setup_run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_DAS.nc", cache=self.cache_path)
        model.run(output_path=self.output_path, output_formats=self.output_formats)
        self.assertTrue(os.path.exists(os.path.join(self.output_path, "simple_trackline.geojson")))
        self.draw_trackline(os.path.join(self.output_path, "simple_trackline.geojson"))
