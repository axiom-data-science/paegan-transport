import os
import time
import pytz
import logging
import threading
import multiprocessing
import multiprocessing.pool

from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import cascaded_union

from paegan.transport.particles.particle import LarvaParticle
from paegan.location4d import Location4D
from paegan.transport.utils.asatransport import AsaTransport
from paegan.transport.shoreline import Shoreline
from paegan.cdm.dataset import CommonDataset
from paegan.transport.exceptions import ModelError, BaseDataControllerError
from paegan.transport.forcers import BaseForcer
from paegan.logger import logger
import paegan.transport.export as ex


class Runner(object):
    """
    Shim needed for Pool/cluster runnable interface.
    """
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

        # Should we thread the result listener?  Defined in subclasses that need it.
        self.thread_result_listener = False

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

    def total_particle_count(self):
        return len(self.particles)

    def total_task_count(self):
        return self.total_particle_count()

    def listen_for_results(self, output_h5_file, total_particles):
        logger.info("Waiting for %i particle results" % total_particles)
        logger.progress((5, "Running model"))

        particles = []
        retrieved = 0
        timeout = 200

        while retrieved < total_particles:
            try:
                # self.result is an iterator that can timeout on next()
                particle = self.result.next(timeout)
                retrieved += 1
                particles.append(particle)
            except StopIteration:
                assert retrieved >= total_particles
                break
            except:
                logger.exception("Particle has FAILED!!")
                continue

            # We multiply by 90 here to save 10% for the exporting
            logger.progress((round((float(retrieved) / total_particles) * 90., 1), "%s Particle(s) complete" % retrieved))

        logger.info(particles)
        results = ex.ResultsPyTable(output_h5_file)
        for p in particles:
            for x in range(len(p.locations)):
                results.write(p.timestep_index_dump(x))
        results.compute()
        results.close()

        return

    def cleanup(self):

        # Remove timevar
        del self.timevar

    def start_tasks(self, **kwargs):

        if self.pool is None:
            self.pool = multiprocessing.Pool()

        try:
            logger.info('Adding %i particles as tasks' % self.total_particle_count())
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
                                    shoreline_index_buffer=self.shoreline_index_buffer)
                tasks.append(forcer)

            return self.pool.imap(Runner(), tasks)

        except Exception:
            logger.exception("Something didn't start correctly!")
            raise

    def setup_run(self, hydrodataset, **kwargs):

        self.hydrodataset = hydrodataset

        logger.setLevel(logging.PROGRESS)

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

    def run(self, **kwargs):

        logger.progress((4, "Starting tasks"))
        self.result = self.start_tasks(**kwargs)
        if self.result is None:
            raise BaseDataControllerError("Not all tasks started! Exiting.")

        # Store results in hdf5 file for processing later
        output_h5_file = None
        if kwargs.get('output_path') is not None:
            output_h5_file = os.path.join(kwargs.get('output_path'), 'results.h5')

        if self.thread_result_listener is True:
            rl = threading.Thread(name="ResultListener", target=self.listen_for_results, args=(output_h5_file, self.total_particle_count()))
            rl.daemon = True
            rl.start()
            rl.join()  # This blocks until the tasks are all done.
        else:
            self.listen_for_results(output_h5_file, self.total_particle_count())    # This blocks until the tasks are all done.

        logger.info('Tasks are all finished... Cleaning up!!')
        self.cleanup()

        # If output_formats and path specified,
        # output particle run data to disk when completed
        if "output_formats" in kwargs:

            logger.progress((96, "Exporting results"))

            # Make sure output_path is also included
            if kwargs.get("output_path", None) is not None:
                formats = kwargs.get("output_formats")
                output_path = kwargs.get("output_path")
                if isinstance(formats, list):
                    for fmt in formats:
                        logger.info("Exporting to: %s" % fmt)
                        try:
                            # Calls the export function
                            fmt.export(output_path, output_h5_file)
                        except:
                            logger.exception("Failed to export to: %s" % fmt)
                else:
                    logger.warn('The output_formats parameter should be a list, not saving any output!')
            else:
                logger.warn('No output path defined, not saving any output!')
        else:
            logger.warn('No output_formats parameter was defined, not saving any output!')

        logger.progress((97, "Model Run Complete"))

        return
