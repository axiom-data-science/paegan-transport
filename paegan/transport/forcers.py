import os
import time as timer
import math
import time
import json
import logging

import textwrap

import numpy as np

from paegan.location4d import Location4D

from paegan.transport.utils.asatransport import AsaTransport
from paegan.transport.shoreline import Shoreline
from paegan.transport.bathymetry import Bathymetry

from paegan.cdm.dataset import CommonDataset
from paegan.cdm.timevar import date2num

from paegan.logger import logger


class BaseForcer(object):

    def __init__(self, hydrodataset, **kwargs):

        """
        part, common_variables, timevar, times, start_time, models,
        release_location_centroid, usebathy, useshore, usesurface,
        get_data, n_run, read_lock, has_read_lock, read_count,
        point_get, data_request_lock, has_data_request_lock, reverse_distance=None, bathy=None,
        shoreline_path=None, shoreline_feature=None, time_method=None, caching=None, redis_url=None, redis_results_channel=None, shoreline_index_buffer=None):

            This is the task/class/object/job that forces an
            individual particle and communicates with the
            other particles and data controller for local
            cache updates
        """
        assert hydrodataset is not None

        # Common parameters
        self.hydrodataset               = hydrodataset
        self.bathy_path                 = kwargs.get("bathy_path")
        self.release_location_centroid  = kwargs.get("release_location_centroid")
        self.particle                   = kwargs.get("particle")
        self.times                      = kwargs.get("times")
        self.timevar                    = kwargs.get("timevar", None)
        self.start_time                 = kwargs.get("start_time")
        self.models                     = kwargs.get("models", [])
        self.usebathy                   = kwargs.get("usebathy", False)
        self.useshore                   = kwargs.get("useshore", False)
        self.usesurface                 = kwargs.get("usesurface", True)
        self.shoreline_path             = kwargs.get("shoreline_path")
        self.shoreline_feature          = kwargs.get("shoreline_feature", None)
        self.shoreline_index_buffer     = kwargs.get("shoreline_index_buffer", 0.1)
        self.time_method                = kwargs.get("time_method", "nearest")
        self.reverse_distance           = kwargs.get("reverse_distance", 500)

        # Redis for results
        self.redis_url                  = kwargs.get("redis_url", None)
        self.redis_results_channel      = kwargs.get("redis_results_channel", None)

        # Set common variable names
        self.common_variables = kwargs.get("common_variables")
        self.uname      = self.common_variables.get("u", None)
        self.vname      = self.common_variables.get("v", None)
        self.wname      = self.common_variables.get("w", None)
        self.temp_name  = self.common_variables.get("temp", None)
        self.salt_name  = self.common_variables.get("salt", None)
        self.xname      = self.common_variables.get("x", None)
        self.yname      = self.common_variables.get("y", None)
        self.zname      = self.common_variables.get("z", None)
        self.tname      = self.common_variables.get("time", None)

        self.active     = None

    def load_initial_dataset(self):
        """
        Initialize self.dataset, then close it
        A cacher will have to wrap this in locks, while a straight runner will not.
        """
        try:
            self.dataset = CommonDataset.open(self.hydrodataset)
            if self.timevar is None:
                self.timevar = self.dataset.gettimevar(self.common_variables.get("u"))
        except Exception:
            logger.warn("No source dataset: %s.  Particle exiting" % self.hydrodataset)
            raise

    def boundary_interaction(self, **kwargs):
        """
            Returns a list of Location4D objects
        """
        particle = kwargs.pop('particle')
        starting = kwargs.pop('starting')
        ending   = kwargs.pop('ending')

        # shoreline
        if self.useshore:
            intersection_point = self._shoreline.intersect(start_point=starting.point, end_point=ending.point)
            if intersection_point is not None:
                # Set the intersection point.
                hitpoint = Location4D(point=intersection_point['point'], time=starting.time + (ending.time - starting.time))
                particle.location = hitpoint

                # This relies on the shoreline to put the particle in water and not on shore.
                resulting_point = self._shoreline.react(start_point=starting,
                                                        end_point=ending,
                                                        hit_point=hitpoint,
                                                        reverse_distance=self.reverse_distance,
                                                        feature=intersection_point['feature'],
                                                        distance=kwargs.get('distance'),
                                                        angle=kwargs.get('angle'),
                                                        azimuth=kwargs.get('azimuth'),
                                                        reverse_azimuth=kwargs.get('reverse_azimuth'))
                ending.latitude = resulting_point.latitude
                ending.longitude = resulting_point.longitude
                ending.depth = resulting_point.depth
                if logger.isEnabledFor(logging.DEBUG):
                    logger.debug("%s - hit the shoreline at %s.  Setting location to %s." % (particle.logstring(), hitpoint.logstring(),  ending.logstring()))

        # bathymetry
        if self.usebathy:
            if not particle.settled:
                bintersect = self._bathymetry.intersect(start_point=starting, end_point=ending)
                if bintersect:
                    pt = self._bathymetry.react(type='reverse', start_point=starting, end_point=ending)
                    if logger.isEnabledFor(logging.DEBUG):
                        logger.debug("%s - hit the bottom at %s.  Setting location to %s." % (particle.logstring(), ending.logstring(), pt.logstring()))
                    ending.latitude = pt.latitude
                    ending.longitude = pt.longitude
                    ending.depth = pt.depth

        # sea-surface
        if self.usesurface:
            if ending.depth > 0:
                if logger.isEnabledFor(logging.DEBUG):
                    logger.debug("%s - rose out of the water.  Setting depth to 0." % particle.logstring())
                ending.depth = 0

        particle.location = ending

    def get_nearest_data(self, i):
        """ Note: self.dataset.opennc() must be called before calling this function.
            This is because the caching forcer must close it everytime, while a non caching
            forcer can leave the dataset open.
        """
        try:
            # Grab data at time index closest to particle location
            u = np.mean(np.mean(self.dataset.get_values('u', timeinds=[np.asarray([i])], point=self.particle.location )))
            v = np.mean(np.mean(self.dataset.get_values('v', timeinds=[np.asarray([i])], point=self.particle.location )))
            # if there is vertical velocity inthe dataset, get it
            if 'w' in self.dataset.nc.variables:
                w = np.mean(np.mean(self.dataset.get_values('w', timeindsf=[np.asarray([i])], point=self.particle.location )))
            else:
                w = 0.0
            # If there is salt and temp in the dataset, get it
            if self.temp_name is not None and self.salt_name is not None:
                temp = np.mean(np.mean(self.dataset.get_values('temp', timeinds=[np.asarray([i])], point=self.particle.location )))
                salt = np.mean(np.mean(self.dataset.get_values('salt', timeinds=[np.asarray([i])], point=self.particle.location )))

            # Check for nans that occur in the ocean (happens because
            # of model and coastline resolution mismatches)
            if np.isnan(u).any() or np.isnan(v).any() or np.isnan(w).any():
                # Take the mean of the closest 4 points
                # If this includes nan which it will, result is nan
                uarray1 = self.dataset.get_values('u', timeinds=[np.asarray([i])], point=self.particle.location, num=2)
                varray1 = self.dataset.get_values('v', timeinds=[np.asarray([i])], point=self.particle.location, num=2)
                if 'w' in self.dataset.nc.variables:
                    warray1 = self.dataset.get_values('w', timeinds=[np.asarray([i])], point=self.particle.location, num=2)
                    w = warray1.mean()
                else:
                    w = 0.0

                if self.temp_name is not None and self.salt_name is not None:
                    temparray1 = self.dataset.get_values('temp', timeinds=[np.asarray([i])], point=self.particle.location, num=2)
                    saltarray1 = self.dataset.get_values('salt', timeinds=[np.asarray([i])], point=self.particle.location, num=2)
                    temp = temparray1.mean()
                    salt = saltarray1.mean()
                u = uarray1.mean()
                v = varray1.mean()

            if self.temp_name is None:
                temp = np.nan
            if self.salt_name is None:
                salt = np.nan

        except Exception:
            logger.exception("Could not retrieve data.")
            raise

        return u, v, w, temp, salt

    def get_linterp_data(self, i, currenttime):
        """ Note: self.dataset.opennc() must be called before calling this function.
            This is because the caching forcer must close it everytime, while a non caching
            forcer can leave the dataset open.
        """
        try:
            # Grab data at time index closest to particle location
            u = [np.mean(np.mean(self.dataset.get_values('u', timeinds=[np.asarray([i])], point=self.particle.location ))),
                 np.mean(np.mean(self.dataset.get_values('u', timeinds=[np.asarray([i+1])], point=self.particle.location )))]
            v = [np.mean(np.mean(self.dataset.get_values('v', timeinds=[np.asarray([i])], point=self.particle.location ))),
                 np.mean(np.mean(self.dataset.get_values('v', timeinds=[np.asarray([i+1])], point=self.particle.location )))]
            # if there is vertical velocity inthe dataset, get it
            if 'w' in self.dataset.nc.variables:
                w = [np.mean(np.mean(self.dataset.get_values('w', timeinds=[np.asarray([i])], point=self.particle.location ))),
                     np.mean(np.mean(self.dataset.get_values('w', timeinds=[np.asarray([i+1])], point=self.particle.location )))]
            else:
                w = [0.0, 0.0]
            # If there is salt and temp in the dataset, get it
            if self.temp_name is not None and self.salt_name is not None:
                temp = [np.mean(np.mean(self.dataset.get_values('temp', timeinds=[np.asarray([i])], point=self.particle.location ))),
                        np.mean(np.mean(self.dataset.get_values('temp', timeinds=[np.asarray([i+1])], point=self.particle.location )))]
                salt = [np.mean(np.mean(self.dataset.get_values('salt', timeinds=[np.asarray([i])], point=self.particle.location ))),
                        np.mean(np.mean(self.dataset.get_values('salt', timeinds=[np.asarray([i+1])], point=self.particle.location )))]

            # Check for nans that occur in the ocean (happens because
            # of model and coastline resolution mismatches)
            if np.isnan(u).any() or np.isnan(v).any() or np.isnan(w).any():
                # Take the mean of the closest 4 points
                # If this includes nan which it will, result is nan
                uarray1 = self.dataset.get_values('u', timeinds=[np.asarray([i])], point=self.particle.location, num=2)
                varray1 = self.dataset.get_values('v', timeinds=[np.asarray([i])], point=self.particle.location, num=2)
                uarray2 = self.dataset.get_values('u', timeinds=[np.asarray([i+1])], point=self.particle.location, num=2)
                varray2 = self.dataset.get_values('v', timeinds=[np.asarray([i+1])], point=self.particle.location, num=2)
                if 'w' in self.dataset.nc.variables:
                    warray1 = self.dataset.get_values('w', timeinds=[np.asarray([i])], point=self.particle.location, num=2)
                    warray2 = self.dataset.get_values('w', timeinds=[np.asarray([i+1])], point=self.particle.location, num=2)
                    w = [warray1.mean(), warray2.mean()]
                else:
                    w = [0.0, 0.0]

                if self.temp_name is not None and self.salt_name is not None:
                    temparray1 = self.dataset.get_values('temp', timeinds=[np.asarray([i])], point=self.particle.location, num=2)
                    saltarray1 = self.dataset.get_values('salt', timeinds=[np.asarray([i])], point=self.particle.location, num=2)
                    temparray2 = self.dataset.get_values('temp', timeinds=[np.asarray([i+1])], point=self.particle.location, num=2)
                    saltarray2 = self.dataset.get_values('salt', timeinds=[np.asarray([i+1])], point=self.particle.location, num=2)
                    temp = [temparray1.mean(), temparray2.mean()]
                    salt = [saltarray1.mean(), saltarray2.mean()]
                u = [uarray1.mean(), uarray2.mean()]
                v = [varray1.mean(), varray2.mean()]

            # Linear interp of data between timesteps
            currenttime = date2num(currenttime)
            timevar = self.timevar.datenum
            u = self.linterp(timevar[i:i+2], u, currenttime)
            v = self.linterp(timevar[i:i+2], v, currenttime)
            w = self.linterp(timevar[i:i+2], w, currenttime)
            if self.temp_name is not None and self.salt_name is not None:
                temp = self.linterp(timevar[i:i+2], temp, currenttime)
                salt = self.linterp(timevar[i:i+2], salt, currenttime)

            if self.temp_name is None:
                temp = np.nan
            if self.salt_name is None:
                salt = np.nan

        except Exception:
            logger.exception("Could not retrieve data.")
            raise

        return u, v, w, temp, salt

    def linterp(self, setx, sety, x):
        """
            Linear interp of model data values between time steps
        """
        if math.isnan(sety[0]) or math.isnan(setx[0]):
            return np.nan
        return sety[0] + (x - setx[0]) * ( (sety[1]-sety[0]) / (setx[1]-setx[0]) )

    def run(self):

        self.load_initial_dataset()

        redis_connection = None
        if self.redis_url is not None and self.redis_results_channel is not None:
            import redis
            redis_connection = redis.from_url(self.redis_url)

        # Setup shoreline
        self._shoreline = None
        if self.useshore is True:
            self._shoreline = Shoreline(path=self.shoreline_path, feature_name=self.shoreline_feature, point=self.release_location_centroid, spatialbuffer=self.shoreline_index_buffer)
            # Make sure we are not starting on land.  Raises exception if we are.
            self._shoreline.intersect(start_point=self.release_location_centroid, end_point=self.release_location_centroid)

        # Setup Bathymetry
        if self.usebathy is True:
            try:
                self._bathymetry = Bathymetry(file=self.bathy_path)
            except Exception:
                logger.exception("Could not load Bathymetry file: %s, using no Bathymetry for this run!" % self.bathy_path)
                self.usebathy = False

        # Calculate datetime at every timestep
        modelTimestep, newtimes = AsaTransport.get_time_objects_from_model_timesteps(self.times, start=self.start_time)

        if self.time_method == 'interp':
            time_indexs = self.timevar.nearest_index(newtimes, select='before')
        elif self.time_method == 'nearest':
            time_indexs = self.timevar.nearest_index(newtimes)
        else:
            logger.warn("Method for computing u,v,w,temp,salt not supported!")
        try:
            assert len(newtimes) == len(time_indexs)
        except AssertionError:
            logger.exception("Time indexes are messed up. Need to have equal datetime and time indexes")
            raise

        # Keep track of how much time we spend in each area.
        tot_boundary_time = 0.
        tot_model_time    = {}
        tot_read_data     = 0.
        for m in self.models:
            tot_model_time[m.name] = 0.

        # Set the base conditions
        # If using Redis, send the results
        if redis_connection is not None:
            redis_connection.publish(self.redis_results_channel, json.dumps(self.particle.timestep_dump()))

        # loop over timesteps
        # We don't loop over the last time_index because
        # we need to query in the time_index and set the particle's
        # location as the 'newtime' object.
        for loop_i, i in enumerate(time_indexs[0:-1]):

            if self.active and self.active.value is False:
                raise ValueError("Particle exiting due to Failure.")

            newloc = None

            st = time.clock()
            # Get the variable data required by the models
            if self.time_method == 'nearest':
                u, v, w, temp, salt = self.get_nearest_data(i)
            elif self.time_method == 'interp':
                u, v, w, temp, salt = self.get_linterp_data(i, newtimes[loop_i])
            else:
                logger.warn("Method for computing u,v,w,temp,salt is unknown. Only 'nearest' and 'interp' are supported.")
            tot_read_data += (time.clock() - st)

            # Get the bathy value at the particles location
            if self.usebathy is True:
                bathymetry_value = self._bathymetry.get_depth(self.particle.location)
            else:
                bathymetry_value = -999999999999999

            # Age the particle by the modelTimestep (seconds)
            # 'Age' meaning the amount of time it has been forced.
            self.particle.age(seconds=modelTimestep[loop_i])

            # loop over models - sort these in the order you want them to run
            for model in self.models:
                st = time.clock()
                movement = model.move(self.particle, u, v, w, modelTimestep[loop_i], temperature=temp, salinity=salt, bathymetry_value=bathymetry_value)
                newloc = Location4D(latitude=movement['latitude'], longitude=movement['longitude'], depth=movement['depth'], time=newtimes[loop_i+1])
                tot_model_time[m.name] += (time.clock() - st)
                if logger.isEnabledFor(logging.DEBUG):
                    logger.debug("%s - moved %.3f meters (horizontally) and %.3f meters (vertically) by %s with data from %s" % (self.particle.logstring(), movement['distance'], movement['vertical_distance'], model.__class__.__name__, newtimes[loop_i].isoformat()))
                if newloc:
                    st = time.clock()
                    self.boundary_interaction(particle=self.particle, starting=self.particle.location, ending=newloc,
                                              distance=movement['distance'], angle=movement['angle'],
                                              azimuth=movement['azimuth'], reverse_azimuth=movement['reverse_azimuth'],
                                              vertical_distance=movement['vertical_distance'], vertical_angle=movement['vertical_angle'])
                    tot_boundary_time += (time.clock() - st)
                if logger.isEnabledFor(logging.DEBUG):
                    logger.debug("%s - was forced by %s and is now at %s" % (self.particle.logstring(), model.__class__.__name__, self.particle.location.logstring()))

            self.particle.note = self.particle.outputstring()
            # Each timestep, save the particles status and environmental variables.
            # This keep fields such as temp, salt, halted, settled, and dead matched up with the number of timesteps
            self.particle.save()

            # If using Redis, send the results
            if redis_connection is not None:
                redis_connection.publish(self.redis_results_channel, json.dumps(self.particle.timestep_dump()))

        self.dataset.closenc()

        # We won't pull data for the last entry in locations, but we need to populate it with fill data.
        self.particle.fill_gap()

        if self.usebathy is True:
            self._bathymetry.close()

        if self.useshore is True:
            self._shoreline.close()

        logger.info(textwrap.dedent('''Particle %i Stats:
                          Data read: %f seconds
                          Model forcing: %s seconds
                          Boundary intersection: %f seconds''' % (self.particle.uid, tot_read_data, { s: '{:g} seconds'.format(f) for s, f in list(tot_model_time.items()) }, tot_boundary_time)))

        return self.particle

    def __call__(self, active):
        self.active = active
        return self.run()


class CachingForcer(BaseForcer):

    def __init__(self, *args, **kwargs):
        super(CachingForcer, self).__init__(*args, **kwargs)
        self.get_data               = kwargs.get("get_data")
        self.read_lock              = kwargs.get("read_lock")
        self.has_read_lock          = kwargs.get("has_read_lock")
        self.read_count             = kwargs.get("read_count")
        self.point_get              = kwargs.get("point_get")
        self.data_request_lock      = kwargs.get("data_request_lock")
        self.has_data_request_lock  = kwargs.get("has_data_request_lock")

    def load_initial_dataset(self):
        """
        Initialize self.dataset, then close it
        A cacher will have to wrap this in locks, while a straight runner will not.
        """
        try:
            with self.read_lock:
                self.read_count.value += 1
                self.has_read_lock.append(os.getpid())
            self.dataset = CommonDataset.open(self.hydrodataset)
            self.dataset.closenc()
        except Exception:
            logger.warn("No source dataset: %s.  Particle exiting" % self.hydrodataset)
            raise
        finally:
            with self.read_lock:
                self.read_count.value -= 1
                self.has_read_lock.remove(os.getpid())

    def __call__(self, active):

        if active.value is True:
            while self.get_data.value is True:
                logger.info("Waiting for DataController to start...")
                timer.sleep(5)
                pass

        return super(CachingForcer, self).__call__(active)

    def get_nearest_data(self, i):
        self.fill_cache_with_nearest_data(i)

        # Now that the cache is filled, get the actual data and return
        with self.read_lock:
            self.read_count.value += 1
            self.has_read_lock.append(os.getpid())
        try:
            self.dataset.opennc()
            return super(CachingForcer, self).get_nearest_data(i)
        except Exception:
            logger.exception("Could not retrieve data")
        finally:
            self.dataset.closenc()
            with self.read_lock:
                self.read_count.value -= 1
                self.has_read_lock.remove(os.getpid())

    def get_linterp_data(self, i, currenttime):
        self.fill_cache_with_linterp_data(i, currenttime)

        # Now that the cache is filled, get the actual data and return
        with self.read_lock:
            self.read_count.value += 1
            self.has_read_lock.append(os.getpid())
        try:
            self.dataset.opennc()
            return super(CachingForcer, self).get_linterp_data(i, currenttime)
        except Exception:
            logger.exception("Could not retrieve data")
        finally:
            self.dataset.closenc()
            with self.read_lock:
                self.read_count.value -= 1
                self.has_read_lock.remove(os.getpid())

    def need_data(self, i):
        """
            Method to test if cache contains the data that
            the particle needs
        """

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Checking cache for data availability at %s." % self.particle.location.logstring())

        try:
            # Tell the DataController that we are going to be reading from the file
            with self.read_lock:
                self.read_count.value += 1
                self.has_read_lock.append(os.getpid())

            self.dataset.opennc()
            # Test if the cache has the data we need
            # If the point we request contains fill values,
            # we need data
            cached_lookup = self.dataset.get_values('domain', timeinds=[np.asarray([i])], point=self.particle.location)
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug("Type of result: %s" % type(cached_lookup))
                logger.debug("Double mean of result: %s" % np.mean(np.mean(cached_lookup)))
                logger.debug("Type of Double mean of result: %s" % type(np.mean(np.mean(cached_lookup))))
            if type(np.mean(np.mean(cached_lookup))) == np.ma.core.MaskedConstant:
                need = True
                if logger.isEnabledFor(logging.DEBUG):
                    logger.debug("I NEED data.  Got back: %s" % cached_lookup)
            else:
                need = False
                logger.debug("I DO NOT NEED data")
        except Exception:
            # If the time index doesnt even exist, we need
            need = True
            logger.debug("I NEED data (no time index exists in cache)")
        finally:
            self.dataset.closenc()
            with self.read_lock:
                self.read_count.value -= 1
                self.has_read_lock.remove(os.getpid())

        return need  # Returns True if we need data or False if we dont

    def fill_cache_with_nearest_data(self, i):
        """
            Method to streamline request for data from cache,
            Uses nearest time to get u,v,w,temp,salt
        """
        if self.active.value is True:
            while self.get_data.value is True:
                logger.debug("Waiting for DataController to release cache file so I can read from it...")
                timer.sleep(2)
                pass

        if self.need_data(i):
            # Acquire lock for asking for data
            self.data_request_lock.acquire()
            self.has_data_request_lock.value = os.getpid()
            try:
                if self.need_data(i):

                    with self.read_lock:
                        self.read_count.value += 1
                        self.has_read_lock.append(os.getpid())

                    # Open netcdf file on disk from commondataset
                    self.dataset.opennc()
                    # Get the indices for the current particle location
                    indices = self.dataset.get_indices('u', timeinds=[np.asarray([i-1])], point=self.particle.location )
                    self.dataset.closenc()

                    with self.read_lock:
                        self.read_count.value -= 1
                        self.has_read_lock.remove(os.getpid())

                    # Override the time
                    self.point_get.value = [indices[0]+1, indices[-2], indices[-1]]
                    # Request that the data controller update the cache
                    # DATA CONTOLLER STARTS
                    self.get_data.value = True
                    # Wait until the data controller is done
                    if self.active.value is True:
                        while self.get_data.value is True:
                            logger.debug("Waiting for DataController to update cache...")
                            timer.sleep(2)
                            pass
            except Exception:
                raise
            finally:
                self.has_data_request_lock.value = -1
                self.data_request_lock.release()

    def fill_cache_with_linterp_data(self, i, currenttime):
        """
            Method to streamline request for data from cache,
            Uses linear interpolation bewtween timesteps to
            get u,v,w,temp,salt
        """
        if self.active.value is True:
            while self.get_data.value is True:
                logger.debug("Waiting for DataController to release cache file so I can read from it...")
                timer.sleep(2)
                pass

        if self.need_data(i+1):
            # Acquire lock for asking for data
            self.data_request_lock.acquire()
            self.has_data_request_lock.value = os.getpid()
            try:
                # Do I still need data?
                if self.need_data(i+1):

                    # Tell the DataController that we are going to be reading from the file
                    with self.read_lock:
                        self.read_count.value += 1
                        self.has_read_lock.append(os.getpid())

                    # Open netcdf file on disk from commondataset
                    self.dataset.opennc()
                    # Get the indices for the current particle location
                    indices = self.dataset.get_indices('u', timeinds=[np.asarray([i-1])], point=self.particle.location )
                    self.dataset.closenc()

                    with self.read_lock:
                        self.read_count.value -= 1
                        self.has_read_lock.remove(os.getpid())

                    # Override the time
                    # get the current time index data
                    self.point_get.value = [indices[0] + 1, indices[-2], indices[-1]]
                    # Request that the data controller update the cache
                    self.get_data.value = True
                    # Wait until the data controller is done
                    if self.active.value is True:
                        while self.get_data.value is True:
                            logger.debug("Waiting for DataController to update cache with the CURRENT time index")
                            timer.sleep(2)
                            pass

                    # Do we still need to get the next timestep?
                    if self.need_data(i+1):
                        # get the next time index data
                        self.point_get.value = [indices[0] + 2, indices[-2], indices[-1]]
                        # Request that the data controller update the cache
                        self.get_data.value = True
                        # Wait until the data controller is done
                        if self.active.value is True:
                            while self.get_data.value is True:
                                logger.debug("Waiting for DataController to update cache with the NEXT time index")
                                timer.sleep(2)
                                pass
            except Exception:
                logger.warn("Particle failed to request data correctly")
                raise
            finally:
                # Release lock for asking for data
                self.has_data_request_lock.value = -1
                self.data_request_lock.release()
