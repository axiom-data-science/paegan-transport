import os
import sys
import time as timer
import math
import traceback
import Queue
import logging
import multiprocessing

import numpy as np
import netCDF4

from paegan.location4d import Location4D

from paegan.transport.utils.asatransport import AsaTransport
from paegan.transport.shoreline import Shoreline
from paegan.transport.bathymetry import Bathymetry
from paegan.transport.exceptions import CachingDataControllerError

from paegan.cdm.dataset import CommonDataset

from paegan.transport.forcers import BaseForcer

from paegan.logger import logger


class Consumer(multiprocessing.Process):
    def __init__(self, task_queue, result_queue, n_run, nproc_lock, active=True, get_data=None, **kwargs):
        """
            This is the process class that does all the handling of queued tasks
        """
        multiprocessing.Process.__init__(self, **kwargs)
        self.task_queue     = task_queue
        self.result_queue   = result_queue
        self.n_run          = n_run
        self.nproc_lock     = nproc_lock
        self.active         = active
        self.get_data       = get_data

    def run(self):

        while True:

            try:
                next_task = self.task_queue.get(True, 10)
            except Queue.Empty:
                logger.info("No tasks left to complete, closing %s" % self.name)
                break
            else:
                answer = (None, None)
                try:
                    answer = (1, next_task(self.active))
                except Exception:
                    logger.exception("Disabling Error")
                    if isinstance(next_task, CachingDataController):
                        answer = (-2, "CachingDataController")
                        # Tell the particles that the CachingDataController is releasing file
                        self.get_data.value = False
                        # The data controller has died, so don't process any more tasks
                        self.active.value = False
                    elif isinstance(next_task, BaseForcer):
                        answer = (-1, next_task.particle)
                    else:
                        logger.warn("Strange task raised an exception: %s" % str(next_task.__class__))
                        answer = (None, None)
                finally:
                    self.result_queue.put(answer)

                    self.nproc_lock.acquire()
                    self.n_run.value = self.n_run.value - 1
                    self.nproc_lock.release()

                    self.task_queue.task_done()


class CachingDataController(object):
    def __init__(self, hydrodataset, common_variables, n_run, get_data, write_lock, has_write_lock, read_lock, read_count,
                 time_chunk, horiz_chunk, times, start_time, point_get, start, **kwargs):
        """
            The data controller controls the updating of the
            local netcdf data cache
        """
        assert "cache_path" in kwargs
        self.cache_path = kwargs["cache_path"]
        self.caching = kwargs.get("caching", True)
        self.hydrodataset = hydrodataset
        if self.cache_path == self.hydrodataset and self.caching is True:
            raise CachingDataControllerError("Caching is set to True but the cache path and data path are the same.  Refusing to overwrite the data path.")

        self.n_run = n_run
        self.get_data = get_data
        self.write_lock = write_lock
        self.has_write_lock = has_write_lock
        self.read_lock = read_lock
        self.read_count = read_count
        self.inds = None  # np.arange(init_size+1)
        self.time_size = time_chunk
        self.horiz_size = horiz_chunk
        self.point_get = point_get
        self.start_time = start_time
        self.times = times
        self.start = start

        # Set common variable names
        self.uname = common_variables.get("u", None)
        self.vname = common_variables.get("v", None)
        self.wname = common_variables.get("w", None)
        self.temp_name = common_variables.get("temp", None)
        self.salt_name = common_variables.get("salt", None)
        self.xname = common_variables.get("x", None)
        self.yname = common_variables.get("y", None)
        self.zname = common_variables.get("z", None)
        self.tname = common_variables.get("time", None)

    def get_remote_data(self, localvars, remotevars, inds, shape):
        """
            Method that does the updating of local netcdf cache
            with remote data
        """
        # If user specifies 'all' then entire xy domain is
        # grabbed, default is 4, specified in the model controller
        if self.horiz_size == 'all':
            y, y_1 = 0, shape[-2]
            x, x_1 = 0, shape[-1]
        else:
            r = self.horiz_size
            x, x_1 = self.point_get.value[2]-r, self.point_get.value[2]+r+1
            y, y_1 = self.point_get.value[1]-r, self.point_get.value[1]+r+1
            x, x_1 = x[0], x_1[0]
            y, y_1 = y[0], y_1[0]
            if y < 0:
                y = 0
            if x < 0:
                x = 0
            if y_1 > shape[-2]:
                y_1 = shape[-2]
            if x_1 > shape[-1]:
                x_1 = shape[-1]

        # Update domain variable for where we will add data
        domain = self.local.variables['domain']

        if len(shape) == 4:
            domain[inds[0]:inds[-1]+1, 0:shape[1], y:y_1, x:x_1] = np.ones((inds[-1]+1-inds[0], shape[1], y_1-y, x_1-x))
        elif len(shape) == 3:
            domain[inds[0]:inds[-1]+1, y:y_1, x:x_1] = np.ones((inds[-1]+1-inds[0], y_1-y, x_1-x))

        # Update the local variables with remote data
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Filling cache with: Time - %s:%s, Lat - %s:%s, Lon - %s:%s" % (str(inds[0]), str(inds[-1]+1), str(y), str(y_1), str(x), str(x_1)))
        for local, remote in zip(localvars, remotevars):
            if len(shape) == 4:
                local[inds[0]:inds[-1]+1, 0:shape[1], y:y_1, x:x_1] = remote[inds[0]:inds[-1]+1,  0:shape[1], y:y_1, x:x_1]
            else:
                local[inds[0]:inds[-1]+1, y:y_1, x:x_1] = remote[inds[0]:inds[-1]+1, y:y_1, x:x_1]

    def __call__(self, active):
        c = 0

        self.dataset = CommonDataset.open(self.hydrodataset)
        self.remote = self.dataset.nc

        # Calculate the datetimes of the model timesteps like
        # the particle objects do, so we can figure out unique
        # time indices
        modelTimestep, newtimes = AsaTransport.get_time_objects_from_model_timesteps(self.times, start=self.start_time)

        timevar = self.dataset.gettimevar(self.uname)

        # Don't need to grab the last datetime, as it is not needed for forcing, only
        # for setting the time of the final particle forcing
        time_indexs = timevar.nearest_index(newtimes[0:-1], select='before')

        # Have to make sure that we get the plus 1 for the
        # linear interpolation of u,v,w,temp,salt
        self.inds = np.unique(time_indexs)
        self.inds = np.append(self.inds, self.inds.max()+1)

        # While there is at least 1 particle still running,
        # stay alive, if not break
        while self.n_run.value > 1:

            if self.caching is False:
                logger.debug("Caching is False, not doing much.  Just hanging out until all of the particles finish.")
                timer.sleep(10)
                continue

            # If particle asks for data, do the following
            if self.get_data.value is True:
                logger.debug("Particle asked for data!")

                # Wait for particles to get out
                while True:
                    self.read_lock.acquire()

                    logger.debug("Read count: %d" % self.read_count.value)
                    if self.read_count.value > 0:
                        logger.debug("Waiting for write lock on cache file (particles must stop reading)...")
                        self.read_lock.release()
                        timer.sleep(2)
                    else:
                        break

                # Get write lock on the file.  Already have read lock.
                self.write_lock.acquire()
                self.has_write_lock.value = os.getpid()

                if c == 0:
                    logger.debug("Creating cache file")
                    try:
                        # Open local cache for writing, overwrites
                        # existing file with same name
                        self.local = netCDF4.Dataset(self.cache_path, 'w')

                        indices = self.dataset.get_indices(self.uname, timeinds=[np.asarray([0])], point=self.start)
                        self.point_get.value = [self.inds[0], indices[-2], indices[-1]]

                        # Create dimensions for u and v variables
                        self.local.createDimension('time', None)
                        self.local.createDimension('level', None)
                        self.local.createDimension('x', None)
                        self.local.createDimension('y', None)

                        # Create 3d or 4d u and v variables
                        if self.remote.variables[self.uname].ndim == 4:
                            self.ndim = 4
                            dimensions = ('time', 'level', 'y', 'x')
                            coordinates = "time z lon lat"
                        elif self.remote.variables[self.uname].ndim == 3:
                            self.ndim = 3
                            dimensions = ('time', 'y', 'x')
                            coordinates = "time lon lat"
                        shape = self.remote.variables[self.uname].shape

                        # If there is no FillValue defined in the dataset, use np.nan.
                        # Sometimes it will work out correctly and other times we will
                        # have a huge cache file.
                        try:
                            fill = self.remote.variables[self.uname].missing_value
                        except Exception:
                            fill = np.nan

                        # Create domain variable that specifies
                        # where there is data geographically/by time
                        # and where there is not data,
                        #   Used for testing if particle needs to
                        #   ask cache to update
                        domain = self.local.createVariable('domain', 'i', dimensions, zlib=False, fill_value=0)
                        domain.coordinates = coordinates

                        # Create local u and v variables
                        u = self.local.createVariable('u', 'f', dimensions, zlib=False, fill_value=fill)
                        v = self.local.createVariable('v', 'f', dimensions, zlib=False, fill_value=fill)

                        v.coordinates = coordinates
                        u.coordinates = coordinates

                        localvars = [u, v, ]
                        remotevars = [self.remote.variables[self.uname], self.remote.variables[self.vname]]

                        # Create local w variable
                        if self.wname is not None:
                            w = self.local.createVariable('w', 'f', dimensions, zlib=False, fill_value=fill)
                            w.coordinates = coordinates
                            localvars.append(w)
                            remotevars.append(self.remote.variables[self.wname])

                        if self.temp_name is not None and self.salt_name is not None:
                            # Create local temp and salt vars
                            temp = self.local.createVariable('temp', 'f', dimensions, zlib=False, fill_value=fill)
                            salt = self.local.createVariable('salt', 'f', dimensions, zlib=False, fill_value=fill)
                            temp.coordinates = coordinates
                            salt.coordinates = coordinates
                            localvars.append(temp)
                            localvars.append(salt)
                            remotevars.append(self.remote.variables[self.temp_name])
                            remotevars.append(self.remote.variables[self.salt_name])

                        # Create local lat/lon coordinate variables
                        if self.remote.variables[self.xname].ndim == 2:
                            lon = self.local.createVariable('lon', 'f', ("y", "x"), zlib=False)
                            lon[:] = self.remote.variables[self.xname][:, :]
                            lat = self.local.createVariable('lat', 'f', ("y", "x"), zlib=False)
                            lat[:] = self.remote.variables[self.yname][:, :]
                        if self.remote.variables[self.xname].ndim == 1:
                            lon = self.local.createVariable('lon', 'f', ("x"), zlib=False)
                            lon[:] = self.remote.variables[self.xname][:]
                            lat = self.local.createVariable('lat', 'f', ("y"), zlib=False)
                            lat[:] = self.remote.variables[self.yname][:]

                        # Create local z variable
                        if self.zname is not None:
                            if self.remote.variables[self.zname].ndim == 4:
                                z = self.local.createVariable('z', 'f', ("time", "level", "y", "x"), zlib=False)
                                remotez = self.remote.variables[self.zname]
                                localvars.append(z)
                                remotevars.append(remotez)
                            elif self.remote.variables[self.zname].ndim == 3:
                                z = self.local.createVariable('z', 'f', ("level", "y", "x"), zlib=False)
                                z[:] = self.remote.variables[self.zname][:, :, :]
                            elif self.remote.variables[self.zname].ndim == 1:
                                z = self.local.createVariable('z', 'f', ("level",), zlib=False)
                                z[:] = self.remote.variables[self.zname][:]

                        # Create local time variable
                        time = self.local.createVariable('time', 'f8', ("time",), zlib=False)
                        if self.tname is not None:
                            time[:] = self.remote.variables[self.tname][self.inds]

                        if self.point_get.value[0]+self.time_size > np.max(self.inds):
                            current_inds = np.arange(self.point_get.value[0], np.max(self.inds)+1)
                        else:
                            current_inds = np.arange(self.point_get.value[0], self.point_get.value[0] + self.time_size)

                        # Get data from remote dataset and add
                        # to local cache.
                        # Try 20 times on the first attempt
                        current_attempt = 1
                        max_attempts = 20
                        while True:
                            try:
                                assert current_attempt <= max_attempts
                                self.get_remote_data(localvars, remotevars, current_inds, shape)
                            except AssertionError:
                                raise
                            except:
                                logger.warn("CachingDataController failed to get remote data.  Trying again in 20 seconds. %s attemps left." % unicode(max_attempts-current_attempt))
                                logger.exception("Data Access Error")
                                timer.sleep(20)
                                current_attempt += 1
                            else:
                                break

                        c += 1
                    except (StandardError, AssertionError):
                        logger.error("CachingDataController failed to get data (first request)")
                        raise
                    finally:
                        self.local.sync()
                        self.local.close()
                        self.has_write_lock.value = -1
                        self.write_lock.release()
                        self.get_data.value = False
                        self.read_lock.release()
                        logger.debug("Done updating cache file, closing file, and releasing locks")
                else:
                    logger.debug("Updating cache file")
                    try:
                        # Open local cache dataset for appending
                        self.local = netCDF4.Dataset(self.cache_path, 'a')

                        # Create local and remote variable objects
                        # for the variables of interest
                        u = self.local.variables['u']
                        v = self.local.variables['v']
                        time = self.local.variables['time']
                        remoteu = self.remote.variables[self.uname]
                        remotev = self.remote.variables[self.vname]

                        # Create lists of variable objects for
                        # the data updater
                        localvars = [u, v, ]
                        remotevars = [remoteu, remotev, ]
                        if self.salt_name is not None and self.temp_name is not None:
                            salt = self.local.variables['salt']
                            temp = self.local.variables['temp']
                            remotesalt = self.remote.variables[self.salt_name]
                            remotetemp = self.remote.variables[self.temp_name]
                            localvars.append(salt)
                            localvars.append(temp)
                            remotevars.append(remotesalt)
                            remotevars.append(remotetemp)
                        if self.wname is not None:
                            w = self.local.variables['w']
                            remotew = self.remote.variables[self.wname]
                            localvars.append(w)
                            remotevars.append(remotew)
                        if self.zname is not None:
                            remotez = self.remote.variables[self.zname]
                            if remotez.ndim == 4:
                                z = self.local.variables['z']
                                localvars.append(z)
                                remotevars.append(remotez)
                        if self.tname is not None:
                            # remotetime = self.remote.variables[self.tname]
                            time[self.inds] = self.remote.variables[self.inds]

                        if self.point_get.value[0]+self.time_size > np.max(self.inds):
                            current_inds = np.arange(self.point_get.value[0], np.max(self.inds)+1)
                        else:
                            current_inds = np.arange(self.point_get.value[0], self.point_get.value[0] + self.time_size)

                        # Get data from remote dataset and add
                        # to local cache
                        while True:
                            try:
                                self.get_remote_data(localvars, remotevars, current_inds, shape)
                            except:
                                logger.warn("CachingDataController failed to get remote data.  Trying again in 30 seconds")
                                timer.sleep(30)
                            else:
                                break

                        c += 1
                    except StandardError:
                        logger.error("CachingDataController failed to get data (not first request)")
                        raise
                    finally:
                        self.local.sync()
                        self.local.close()
                        self.has_write_lock.value = -1
                        self.write_lock.release()
                        self.get_data.value = False
                        self.read_lock.release()
                        logger.debug("Done updating cache file, closing file, and releasing locks")
            else:
                logger.debug("Particles are still running, waiting for them to request data...")
                timer.sleep(2)

        self.dataset.closenc()

        return "CachingDataController"
