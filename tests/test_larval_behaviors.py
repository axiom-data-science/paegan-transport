import unittest
from paegan.transport.models.transport import Transport
import os
import shutil
import json
import urllib
from paegan.transport.models.behavior import LarvaBehavior
import pytz


class LarvalBehaviourTest(unittest.TestCase):

    def test_behavior_growth_and_settlement(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: test_behavior_growth_and_settlement")

        # 6 days
        num_steps = 144

        num_particles = 2

        # Behavior
        behavior_config = open(os.path.normpath(os.path.join(os.path.dirname(__file__),"./resources/files/behavior_for_run_testing.json"))).read()
        lb = LarvaBehavior(json=behavior_config)

        models = [self.transport]
        models.append(lb)

        model = ModelController(latitude=self.start_lat, longitude=self.start_lon, depth=self.start_depth, start=self.start_time, step=self.time_step, nstep=num_steps, npart=num_particles, models=models, use_bathymetry=True, use_shoreline=True,
            time_chunk=24, horiz_chunk=4, time_method='nearest')

        output_path = os.path.join(self.output_path, "test_behavior_growth_and_settlement")
        shutil.rmtree(output_path, ignore_errors=True)
        os.makedirs(output_path)
        output_formats = ['Shapefile','NetCDF','Trackline','Pickle']

        cache_path = os.path.join(self.cache_path, "test_behavior_growth_and_settlement.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", bathy=self.bathy_file, cache=cache_path, output_path=output_path, output_formats=output_formats)
    """
    """
    def test_quick_settlement(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: test_quick_settlement")

        num_steps = 24

        num_particles = 4

        # Behavior
        behavior_config = open(os.path.normpath(os.path.join(os.path.dirname(__file__),"./resources/files/behavior_quick_settle.json"))).read()
        lb = LarvaBehavior(json=behavior_config)

        models = [self.transport]
        models.append(lb)

        model = ModelController(latitude=self.start_lat, longitude=self.start_lon, depth=self.start_depth, start=self.start_time, step=self.time_step, nstep=num_steps, npart=num_particles, models=models, use_bathymetry=True, use_shoreline=True,
            time_chunk=12, horiz_chunk=2, time_method='nearest')

        output_path = os.path.join(self.output_path, "test_quick_settlement")
        shutil.rmtree(output_path, ignore_errors=True)
        os.makedirs(output_path)
        output_formats = ['Shapefile','NetCDF','Trackline']

        cache_path = os.path.join(self.cache_path, "test_quick_settlement.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_DAS.nc", bathy=self.bathy_file, cache=cache_path, output_path=output_path, output_formats=output_formats)

    def test_kayak_island(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: Kayak Island")

        # 6 days
        num_steps = 1632

        num_particles = 100

        time_step = 3600

        behavior_config = json.loads(urllib.urlopen("http://behaviors.larvamap.asascience.com/library/50ef1bb1cc7b61000700001d.json").read())
        lb = LarvaBehavior(data=behavior_config[u'results'][0])

        #behavior_config = json.loads(open(os.path.normpath(os.path.join(os.path.dirname(__file__),"./resources/files/nick.json"))).read())
        #lb = LarvaBehavior(data=behavior_config[u'results'][0])

        models = [Transport(horizDisp=0.01, vertDisp=0.001)]
        models.append(lb)

        start_time = datetime(2011, 5, 2, 00, tzinfo=pytz.utc)

        start_lat = 59.93517413488866
        start_lon = -144.496213677788
        depth = -1

        shoreline_path = os.path.join(self.shoreline_path, "westcoast", "New_Land_Clean.shp")

        model = ModelController(latitude=start_lat, longitude=start_lon, depth=depth, start=start_time, step=time_step, nstep=num_steps, npart=num_particles, models=models, use_bathymetry=True, use_shoreline=True,
            time_chunk=24, horiz_chunk=5, time_method='interp', shoreline_path=shoreline_path)

        output_path = os.path.join(self.output_path, "kayak_island")
        shutil.rmtree(output_path, ignore_errors=True)
        os.makedirs(output_path)
        output_formats = ['Shapefile','NetCDF','Trackline']

        cache_path = os.path.join(self.cache_path, "kayak_island.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L1_FCST.nc", bathy=self.bathy_file, cache=cache_path, output_path=output_path, output_formats=output_formats)

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

        model = ModelController(latitude=start_lat, longitude=start_lon, depth=depth, start=start_time, step=time_step, nstep=num_steps, npart=num_particles, models=models, use_bathymetry=True, use_shoreline=True,
            time_chunk=24, horiz_chunk=5, time_method='interp', shoreline_path=shoreline_path)

        output_path = os.path.join(self.output_path, "sheep_bay")
        shutil.rmtree(output_path, ignore_errors=True)
        os.makedirs(output_path)
        output_formats = ['Shapefile','NetCDF','Trackline']

        cache_path = os.path.join(self.cache_path, "sheep_bay.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", bathy=self.bathy_file, cache=cache_path, output_path=output_path, output_formats=output_formats)
    """

    """
    def test_diel_migration(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: test_diel_migration")

        num_steps = 168

        num_particles = 4

        start_time = datetime(2013,4,1,0, tzinfo=pytz.utc)

        # Behavior
        behavior_config = open(os.path.normpath(os.path.join(os.path.dirname(__file__),"./resources/files/diel_suncycles.json"))).read()
        lb = LarvaBehavior(json=behavior_config)

        models = [self.transport]
        models.append(lb)

        model = ModelController(latitude=60.68, longitude=-146.42, depth=self.start_depth, start=start_time, step=self.time_step, nstep=num_steps, npart=num_particles, models=models, use_bathymetry=True, use_shoreline=True,
            time_chunk=24, horiz_chunk=2, time_method='nearest')

        output_path = os.path.join(self.output_path, "test_diel_migration")
        shutil.rmtree(output_path, ignore_errors=True)
        os.makedirs(output_path)
        output_formats = ['Shapefile','NetCDF','Trackline']

        cache_path = os.path.join(self.cache_path, "test_diel_migration.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_DAS.nc", bathy=self.bathy_file, cache=cache_path, output_path=output_path, output_formats=output_formats)
