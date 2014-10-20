import time
from datetime import datetime
from paegan.transport.particles.particle import LarvaParticle
from paegan.transport.particles.particle import Particle
from paegan.location4d import Location4D
from paegan.utils.asarandom import AsaRandom
from paegan.transport.utils.asatransport import AsaTransport
from paegan.transport.shoreline import Shoreline
from paegan.cdm.dataset import CommonDataset
from paegan.transport.exceptions import ModelError, BaseDataControllerError
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import cascaded_union
import multiprocessing
from paegan.transport.parallel_manager import CachingDataController, Consumer
from paegan.transport.forcers import CachingForcer, BaseForcer
import paegan.transport.export as ex
import os
import Queue
import pytz
import logging

from paegan.logger import logger

class Runner(object):
    def __call__(self, task):
        return task(None)

class BaseModelController(object):
    """
        Controls the models
    """
    def __init__(self, **kwargs):

        """
            Mandatory named arguments:
            * geometry (Shapely Geometry Object) no default
            * depth (meters) default 0
            * start (DateTime Object) none
            * step (seconds) default 3600
            * npart (number of particles) default 1
            * nstep (number of steps) no default
            * models (list object) no default, so far there is a transport model and a behavior model
            geometry is interchangeable (if it is a point release) with:
            * latitude (DD) no default
            * longitude (DD) no default
            * depth (meters) default 0

            Non-mandatory named arguments:
            * pool (task pool) - defaults to multiprocessing.pool, inject your own for cluster ability
        """

        # Shoreline
        self._use_shoreline         = kwargs.pop('use_shoreline', True)
        self.shoreline_path         = kwargs.get("shoreline_path", None)
        self.shoreline_feature      = kwargs.get("shoreline_feature", None)
        self.shoreline_index_buffer = kwargs.get("shoreline_index_buffer", 0.1)
        self.reverse_distance       = kwargs.get("reverse_distance", 100)

        # Bathy
        self._use_bathymetry = kwargs.pop('use_bathymetry', True)
        self.bathy_path      = kwargs.get("bathy_path", None)

        # SeaSurface
        self._use_seasurface = kwargs.pop('use_seasurface', True)

        self._depth          = kwargs.pop('depth', 0)
        self._npart          = kwargs.pop('npart', 1)
        self._step           = kwargs.pop('step', 3600)
        self.start           = kwargs.get('start', None)
        if self.start is None:
            raise TypeError("must provide a start time to run the model")

        # Always convert to UTC
        if self.start.tzinfo is None:
            self.start = self.start.replace(tzinfo=pytz.utc)
        self.start = self.start.astimezone(pytz.utc)

        self._models = kwargs.pop('models', None)
        self._dirty  = True

        self.particles              = []
        self.time_method            = kwargs.get('time_method', 'interp').lower()
        try:
            assert "interp" == self.time_method or "nearest" == self.time_method
        except:
            raise TypeError("Not a recognized 'time_method' parameter.  Only 'nearest' or 'interp' are allowed.")

        # The model timesteps in datetime objects
        self.datetimes = []

        # Interchangeables
        if "geometry" in kwargs:
            self.geometry = kwargs.pop('geometry')
            if not isinstance(self.geometry, Point) and not isinstance(self.geometry, Polygon) and not isinstance(self.geometry, MultiPolygon):
                raise TypeError("The geometry attribute must be a shapely Point or Polygon")
        elif "latitude" and "longitude" in kwargs:
            self.geometry = Point(kwargs.pop('longitude'), kwargs.pop('latitude'))
        else:
            raise TypeError("must provide a shapely geometry object (point or polygon) or a latitude and a longitude")

        # Errors
        try:
            self._nstep = kwargs.pop('nstep')
        except StandardError:
            logger.exception("Must provide the number of timesteps to the ModelController")
            raise

        self.pool = kwargs.get('pool', None)

    def set_geometry(self, geo):
        # If polygon is passed in, we need to trim it by the coastline
        # so we don't start particles on land
        if isinstance(geo, Polygon) and self._use_shoreline:
            c = geo.centroid
            b = geo.bounds
            spatialbuffer = max(b[2] - b[0], b[3] - b[1])
            shore_geoms = Shoreline(path=self.shoreline_path, feature_name=self.shoreline_feature, point=c, spatialbuffer=spatialbuffer).geoms
            if len(shore_geoms) > 0:
                all_shore = cascaded_union(shore_geoms)
                geo = geo.difference(all_shore)

        self._geometry = geo

    def get_geometry(self):
        return self._geometry
    geometry = property(get_geometry, set_geometry)

    def get_reference_location(self):
        pt = self.geometry.centroid
        return Location4D(latitude=pt.y, longitude=pt.x, depth=self._depth, time=self.start)
    reference_location = property(get_reference_location, None)

    def set_start(self, sta):
        self._start = sta

    def get_start(self):
        return self._start
    start = property(get_start, set_start)

    def set_particles(self, parts):
        self._particles = parts

    def get_particles(self):
        return self._particles
    particles = property(get_particles, set_particles)

    def __str__(self):
        return "*** BaseModelController ***"

    def get_common_variables_from_dataset(self, dataset):

        def getname(names):
            for n in names:
                nm = dataset.get_varname_from_stdname(n)
                if len(nm) > 0:
                    return nm[0]
                else:
                    continue
            return None

        uname = getname(['eastward_sea_water_velocity', 'eastward_current'])
        vname = getname(['northward_sea_water_velocity', 'northward_current'])
        wname = getname('upward_sea_water_velocity')
        temp_name = getname('sea_water_temperature')
        salt_name = getname('sea_water_salinity')

        coords = dataset.get_coord_names(uname)
        xname = coords['xname']
        yname = coords['yname']
        zname = coords['zname']
        tname = coords['tname']
        tname = None  # temporary

        return {
            "u"     :   uname,
            "v"     :   vname,
            "w"     :   wname,
            "temp"  :   temp_name,
            "salt"  :   salt_name,
            "x"     :   xname,
            "y"     :   yname,
            "z"     :   zname,
            "time"  :   tname
        }

    def get_number_of_tasks(self):
        return len(self.particles)

    def listen_for_results(self):
        logger.info("Waiting for %i particle results" % len(self.particles))
        logger.progress((5, "Running model"))

        particles = self.result.get(timeout=300)        # blocks until completion

        # @TODO: smartness with logging partial results
        # waitloop?

        return particles

        """
        logger.warn("Got an unrecognized response from a task.")
        logger.warn("Particle %s has FAILED!!" % tempres.uid)
        logger.info("A zombie process was caught and task was removed from queue")
        logger.info("Particle %d finished" % tempres.uid)
        # We mulitply by 95 here to save 5% for the exporting
        logger.progress((round((retrieved / self.number_of_tasks) * 90., 1), "Particle %d finished" % tempres.uid))
        logger.info("Got a strange result on results queue: %s" % str(tempres))
        logger.info("Retrieved %i/%i results" % (int(retrieved), self.number_of_tasks))
        """

    def cleanup(self):

        # Remove timevar
        del self.timevar

    def start_tasks(self):

        # @TODO: this is more initialization, but need to prevent derived classes from doing this
        if self.pool is None:
            self.pool = multiprocessing.Pool()

        try:
            logger.info('Adding %i particles as tasks' % len(self.particles))
            tasks = []

            for part in self.particles:
                forcer = BaseForcer(self.hydrodataset,
                                    particle=part,
                                    common_variables=self.common_variables,
                                    timevar=self.timevar,
                                    times=self.times,
                                    start_time=self.start,
                                    models=self._models,
                                    release_location_centroid=self.reference_location.point,
                                    usebathy=self._use_bathymetry,
                                    useshore=self._use_shoreline,
                                    usesurface=self._use_seasurface,
                                    reverse_distance=self.reverse_distance,
                                    bathy_path=self.bathy_path,
                                    shoreline_path=self.shoreline_path,
                                    shoreline_feature=self.shoreline_feature,
                                    time_method=self.time_method,
                                    redis_url=self.redis_url,
                                    redis_results_channel=self.redis_results_channel,
                                    shoreline_index_buffer=self.shoreline_index_buffer)
                tasks.append(forcer)

            ar = self.pool.map_async(Runner(), tasks)
            return ar

        except Exception:
            logger.exception("Something didn't start correctly!")
            raise

    def setup_run(self, **kwargs):

        logger.setLevel(logging.PROGRESS)

        self.redis_url             = None
        self.redis_log_channel     = None
        self.redis_results_channel = None
        if "redis" in kwargs.get("output_formats", []):
            from paegan.logger.redis_handler import RedisHandler
            self.redis_url             = kwargs.get("redis_url")
            self.redis_log_channel     = kwargs.get("redis_log_channel")
            self.redis_results_channel = kwargs.get("redis_results_channel")
            rhandler = RedisHandler(self.redis_log_channel, self.redis_url)
            rhandler.setLevel(logging.PROGRESS)
            logger.addHandler(rhandler)

        # Relax.
        time.sleep(0.5)

        # Add ModelController description to logfile
        logger.info(unicode(self))

        # Add the model descriptions to logfile
        for m in self._models:
            logger.info(unicode(m))

        # Calculate the model timesteps
        # We need times = len(self._nstep) + 1 since data is stored one timestep
        # after a particle is forced with the final timestep's data.
        self.times = range(0, (self._step*self._nstep)+1, self._step)
        # Calculate a datetime object for each model timestep
        # This method is duplicated in CachingDataController and CachingForcer
        # using the 'times' variables above.  Will be useful in those other
        # locations for particles released at different times
        # i.e. released over a few days
        self.modelTimestep, self.datetimes = AsaTransport.get_time_objects_from_model_timesteps(self.times, start=self.start)

        logger.progress((1, "Setting up particle start locations"))
        point_locations = []
        if isinstance(self.geometry, Point):
            point_locations = [self.reference_location] * self._npart
        elif isinstance(self.geometry, Polygon) or isinstance(self.geometry, MultiPolygon):
            point_locations = [Location4D(latitude=loc.y, longitude=loc.x, depth=self._depth, time=self.start) for loc in AsaTransport.fill_polygon_with_points(goal=self._npart, polygon=self.geometry)]

        # Initialize the particles
        logger.progress((2, "Initializing particles"))
        for x in xrange(0, self._npart):
            p = LarvaParticle(id=x)
            p.location = point_locations[x]
            # We don't need to fill the location gaps here for environment variables
            # because the first data collected actually relates to this original
            # position.
            # We do need to fill in fields such as settled, halted, etc.
            p.fill_status_gap()
            # Set the inital note
            p.note = p.outputstring()
            p.notes.append(p.note)
            self.particles.append(p)

        # Each particle is a task, plus the CachingDataController
        self.number_of_tasks = self.get_number_of_tasks()

        logger.progress((3, "Initializing and caching hydro model's grid %s" % self.hydrodataset))
        try:
            ds = CommonDataset.open(self.hydrodataset)
        except Exception:
            logger.exception("Failed to access dataset %s" % self.hydrodataset)
            raise BaseDataControllerError("Inaccessible Dataset: %s" % self.hydrodataset)
        # Query the dataset for common variable names
        # and the time variable.
        logger.debug("Retrieving variable information from dataset")
        self.common_variables = self.get_common_variables_from_dataset(ds)

        self.timevar = None
        try:
            assert self.common_variables.get("u") in ds._current_variables
            assert self.common_variables.get("v") in ds._current_variables
            assert self.common_variables.get("x") in ds._current_variables
            assert self.common_variables.get("y") in ds._current_variables

            self.timevar = ds.gettimevar(self.common_variables.get("u"))
        except AssertionError:
            logger.exception("Could not locate variables needed to run model: %s" % unicode(self.common_variables))
            raise BaseDataControllerError("A required data variable was not found in %s" % self.hydrodataset)

        model_start = self.timevar.get_dates()[0]
        model_end   = self.timevar.get_dates()[-1]

        try:
            assert self.start > model_start
            assert self.start < model_end
        except AssertionError:
            raise BaseDataControllerError("Start time for model (%s) is not available in source dataset (%s/%s)" % (self.datetimes[0], model_start, model_end))

        try:
            assert self.datetimes[-1] > model_start
            assert self.datetimes[-1] < model_end
        except AssertionError:
            raise BaseDataControllerError("End time for model (%s) is not available in source dataset (%s/%s)" % (self.datetimes[-1], model_start, model_end))

        ds.closenc()

    def run(self, hydrodataset, **kwargs):

        self.hydrodataset = hydrodataset

        self.setup_run(**kwargs)

        logger.progress((4, "Starting tasks"))
        self.result = self.start_tasks()
        if self.result is None:
            raise BaseDataControllerError("Not all tasks started! Exiting.")

        # This blocks until the tasks are all done.
        self.listen_for_results()

        logger.info('Consumers are all finished!')

        logger.info('Cleaning up')
        self.cleanup()

        if len(self.particles) > 0:
            # If output_formats and path specified,
            # output particle run data to disk when completed
            if "output_formats" in kwargs:

                logger.progress((96, "Exporting results"))

                # Make sure output_path is also included
                if kwargs.get("output_path", None) is not None:
                    formats = kwargs.get("output_formats")
                    output_path = kwargs.get("output_path")
                    if isinstance(formats, list):
                        for format in formats:
                            logger.info("Exporting to: %s" % format)
                            try:
                                self.export(output_path, format=format)
                            except:
                                logger.exception("Failed to export to: %s" % format)
                    else:
                        logger.warn('The output_formats parameter should be a list, not saving any output!')
                else:
                    logger.warn('No output path defined, not saving any output!')
            else:
                logger.warn('No output format defined, not saving any output!')
        else:
            logger.warn("Model didn't actually do anything, check the log.")
            if self.error_code == -2:
                raise BaseDataControllerError("Error in the BaseDataController")
            else:
                raise ModelError("Error in the model")

        logger.progress((97, "Model Run Complete"))
        return self.particles

    def export(self, folder_path, format=None):
        """
            General purpose export method, gets file type
            from filepath extension

            Valid output formats currently are:
                Trackline: trackline or trkl or *.trkl
                Shapefile: shapefile or shape or shp or *.shp
                NetCDF:    netcdf or nc or *.nc
        """

        if format is None:
            raise ValueError("Must export to a specific format, no format specified.")

        format = format.lower()

        if format == "trackline" or format[-4:] == "trkl":
            ex.Trackline.export(folder=folder_path, particles=self.particles, datetimes=self.datetimes)
        elif format == "shape" or format == "shapefile" or format[-3:] == "shp":
            ex.GDALShapefile.export(folder=folder_path, particles=self.particles, datetimes=self.datetimes)
        elif format == "netcdf" or format[-2:] == "nc":
            ex.NetCDF.export(folder=folder_path, particles=self.particles, datetimes=self.datetimes, summary=str(self))
        elif format == "pickle" or format[-3:] == "pkl" or format[-6:] == "pickle":
            ex.Pickle.export(folder=folder_path, particles=self.particles, datetimes=self.datetimes)
        elif format == "redis":
            return


class CachingModelController(BaseModelController):

    def __init__(self, **kwargs):
        super(CachingModelController, self).__init__(**kwargs)
        self.time_chunk  = kwargs.get('time_chunk', 10)
        self.horiz_chunk = kwargs.get('horiz_chunk', 5)

    def get_number_of_tasks(self):
        # Add the CachingDataController to the number of tasks
        return len(self.particles) + 1

    def start_tasks(self):
        try:
            logger.info('Starting CachingDataController')

            # Add data controller to the queue first so that it
            # can get the initial data and is not blocked
            data_controller = CachingDataController(self.hydrodataset, self.common_variables, self.n_run, self.get_data, self.write_lock, self.has_write_lock, self.read_lock, self.read_count,
                                                    self.time_chunk, self.horiz_chunk, self.times, self.start, self.point_get, self.reference_location, cache_path=self.cache_path)
            self.tasks.put(data_controller)
            # Create CachingDataController worker
            self.data_controller_process = Consumer(self.tasks, self.results, self.n_run, self.nproc_lock, self.active, self.get_data, name="CachingDataController")
            self.data_controller_process.start()

            logger.info('Adding %i particles as tasks' % len(self.particles))

            for part in self.particles:
                forcer = CachingForcer(self.cache_path,
                                       particle=part,
                                       common_variables=self.common_variables,
                                       timevar=self.timevar,
                                       times=self.times,
                                       start_time=self.start,
                                       models=self._models,
                                       release_location_centroid=self.reference_location.point,
                                       usebathy=self._use_bathymetry,
                                       useshore=self._use_shoreline,
                                       usesurface=self._use_seasurface,
                                       reverse_distance=self.reverse_distance,
                                       bathy_path=self.bathy_path,
                                       shoreline_path=self.shoreline_path,
                                       shoreline_feature=self.shoreline_feature,
                                       time_method=self.time_method,
                                       redis_url=self.redis_url,
                                       redis_results_channel=self.redis_results_channel,
                                       shoreline_index_buffer=self.shoreline_index_buffer,
                                       get_data=self.get_data,
                                       read_lock=self.read_lock,
                                       has_read_lock=self.has_read_lock,
                                       read_count=self.read_count,
                                       point_get=self.point_get,
                                       data_request_lock=self.data_request_lock,
                                       has_data_request_lock=self.has_data_request_lock
                                      )
                self.tasks.put(forcer)

            # Create workers for the particles.
            self.procs = [Consumer(self.tasks, self.results, self.n_run, self.nproc_lock, self.active, self.get_data, name="CachingForcer-%d" % i)
                          for i in xrange(self.nproc - 1) ]
            for w in self.procs:
                w.start()
                logger.info('Started %s' % w.name)

            return True

        except Exception:
            logger.exception("Something didn't start correctly!")
            return False

    def __str__(self):
        return "*** CachingModelController ***"

    def setup_run(self, **kwargs):

        super(CachingModelController, self).setup_run(**kwargs)

        # Get the number of cores (may take some tuning) and create that
        # many workers then pass particles into the queue for the workers
        self.mgr = multiprocessing.Manager()

        # This tracks if the system is 'alive'.  Most looping whiles will check this
        # and break out if it is False.  This is True until something goes very wrong.
        self.active = self.mgr.Value('bool', True)

        # Either spin up the number of cores, or the number of tasks
        self.nproc = min(multiprocessing.cpu_count() - 1, self.number_of_tasks)

        # Number of tasks that we need to run.  This is decremented everytime something exits.
        self.n_run = self.mgr.Value('int', self.number_of_tasks)

        # The lock that controls access to the 'n_run' variable
        self.nproc_lock = self.mgr.Lock()

        # Create the task queue for all of the particles and the CachingDataController
        self.tasks = multiprocessing.JoinableQueue(self.number_of_tasks)

        # Create the result queue for all of the particles and the CachingDataController
        self.results = self.mgr.Queue(self.number_of_tasks)

        # Should we remove the cache file at the end of the run?
        self.remove_cache        = kwargs.get("remove_cache", False)
        self.cache_path          = kwargs.get("cache_path", None)

        # Create a temp file for the cache if nothing was passed in
        if self.cache_path is None:
            default_cache_dir = os.path.join(os.path.dirname(__file__), "_cache")
            temp_name = AsaRandom.filename(prefix=str(datetime.now().microsecond), suffix=".nc")
            self.cache_path = os.path.join(default_cache_dir, temp_name)

        # Be sure the cache directory exists
        if not os.path.exists(os.path.dirname(self.cache_path)):
            logger.info("Creating cache directory: %s" % self.cache_path)
            os.makedirs(os.path.dirname(self.cache_path))

        # Create the shared state objects

        # Particles use this to tell the Data Controller to "get_data".
        # The CachingDataController sets this to False when it is done writing to the cache file.
        # Particles will wait for this to be False before reading from the cache file.
        # If we are caching, this starts as True so the Particles don't take off.  If we
        # are not caching, this is False so the Particles can start immediatly.
        self.get_data = self.mgr.Value('bool', True)
        # Particles use this to tell the DataContoller which indices to 'get_data' for
        self.point_get = self.mgr.Value('list', [0, 0, 0])

        # This locks access to the 'has_data_request_lock' value
        self.data_request_lock = self.mgr.Lock()
        # This tracks which Particle PID is asking the CachingDataController for data
        self.has_data_request_lock = self.mgr.Value('int', -1)

        # The lock that controls access to modifying 'has_read_lock' and 'read_count'
        self.read_lock = self.mgr.Lock()
        # List of Particle PIDs that are reading from the cache
        self.has_read_lock = self.mgr.list()
        # The number of Particles that are reading from the cache
        self.read_count = self.mgr.Value('int', 0)

        # When something is writing to the cache file
        self.write_lock = self.mgr.Lock()
        # PID of process with lock
        self.has_write_lock = self.mgr.Value('int', -1)

    def listen_for_results(self):
        try:
            # Get results back from queue, test for failed particles
            return_particles = []
            retrieved = 0.
            self.error_code = 0

            logger.info("Waiting for %i particle results" % len(self.particles))
            logger.progress((5, "Running model"))
            while retrieved < self.number_of_tasks:
                try:
                    # Returns a tuple of code, result
                    code, tempres = self.results.get(timeout=240)
                except Queue.Empty:
                    # Poll the active processes to make sure they are all alive and then continue with loop
                    if not self.data_controller_process.is_alive() and self.data_controller_process.exitcode != 0:
                        # Data controller is zombied, kill off other processes.
                        self.get_data.value is False
                        self.results.put((-2, "CachingDataController"))

                    new_procs = []
                    old_procs = []
                    for p in self.procs:
                        if not p.is_alive() and p.exitcode != 0:
                            # Do what the Consumer would do if something finished.
                            # Add something to results queue
                            self.results.put((-3, "ZombieParticle"))
                            # Decrement nproc (CachingDataController exits when this is 0)
                            with self.nproc_lock:
                                self.n_run.value = self.n_run.value - 1

                            # Remove task from queue (so they can be joined later on)
                            self.tasks.task_done()

                            # Start a new Consumer.  It will exit if there are no tasks available.
                            np = Consumer(self.tasks, self.results, self.n_run, self.nproc_lock, self.active, self.get_data, name=p.name)
                            new_procs.append(np)
                            old_procs.append(p)

                            # Release any locks the PID had
                            if p.pid in self.has_read_lock:
                                with self.read_lock:
                                    self.read_count.value -= 1
                                    self.has_read_lock.remove(p.pid)

                            if self.has_data_request_lock.value == p.pid:
                                self.has_data_request_lock.value = -1
                                try:
                                    self.data_request_lock.release()
                                except:
                                    pass

                            if self.has_write_lock.value == p.pid:
                                self.has_write_lock.value = -1
                                try:
                                    self.write_lock.release()
                                except:
                                    pass

                    for p in old_procs:
                        try:
                            self.procs.remove(p)
                        except ValueError:
                            logger.warn("Did not find %s in the list of processes.  Continuing on." % p.name)

                    for p in new_procs:
                        self.procs.append(p)
                        logger.warn("Started a new consumer (%s) to replace a zombie consumer" % p.name)
                        p.start()

                else:
                    # We got one.
                    retrieved += 1
                    if code is None:
                        logger.warn("Got an unrecognized response from a task.")
                    elif code == -1:
                        logger.warn("Particle %s has FAILED!!" % tempres.uid)
                    elif code == -2:
                        self.error_code = code
                        logger.warn("CachingDataController has FAILED!!  Removing cache file so the particles fail.")
                        try:
                            os.remove(self.cache_path)
                        except OSError:
                            logger.debug("Could not remove cache file, it probably never existed")
                            pass
                    elif code == -3:
                        self.error_code = code
                        logger.info("A zombie process was caught and task was removed from queue")
                    elif isinstance(tempres, Particle):
                        logger.info("Particle %d finished" % tempres.uid)
                        return_particles.append(tempres)
                        # We mulitply by 95 here to save 5% for the exporting
                        logger.progress((round((retrieved / self.number_of_tasks) * 90., 1), "Particle %d finished" % tempres.uid))
                    elif tempres == "CachingDataController":
                        logger.info("CachingDataController finished")
                        logger.progress((round((retrieved / self.number_of_tasks) * 90., 1), "CachingDataController finished"))
                    else:
                        logger.info("Got a strange result on results queue")
                        logger.info(str(tempres))

                    logger.info("Retrieved %i/%i results" % (int(retrieved), self.number_of_tasks))

            if len(return_particles) != len(self.particles):
                logger.warn("Some particles failed and are not included in the output")

            # The results queue should be empty at this point
            assert self.results.empty() is True

            # Should be good to join on the tasks now that the queue is empty
            logger.info("Joining the task queue")
            self.tasks.join()

            self.particles = return_particles

        finally:
            # Join all processes
            logger.info("Joining the processes")
            for w in self.procs + [self.data_controller_process]:
                    # Wait 20 seconds
                    w.join(20.)
                    if w.is_alive():
                        # Process is hanging, kill it.
                        logger.info("Terminating %s forcefully.  This should have exited itself." % w.name)
                        w.terminate()

        return self.particles

    def cleanup(self):
        super(CachingModelController, self).cleanup()

        # Remove Manager so it shuts down
        del self.mgr

        # Remove the cache file
        if self.remove_cache is True:
            try:
                os.remove(self.cache_path)
            except OSError:
                logger.debug("Could not remove cache file, it probably never existed")


class DistributedModelController(BaseModelController):

    def __init__(self, **kwargs):
        super(DistributedModelController, self).__init__(**kwargs)

    def setup_run(self, hydrodataset, **kwargs):
        self.hydrodataset = hydrodataset
        kwargs["manager"] = False
        super(DistributedModelController, self).setup_run(**kwargs)

    def run(self, hydrodataset, **kwargs):
        raise NotImplementedError("Distributed models can not be started, only setup.  Use 'setup_run' and start the Forcers manually.")
