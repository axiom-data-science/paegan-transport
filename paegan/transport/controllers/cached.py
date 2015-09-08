import time
import multiprocessing
import multiprocessing.pool

import os
try:
    import Queue as queue
except ImportError:
    import queue
import multiprocessing
import multiprocessing.pool
from datetime import datetime

from paegan.transport.particles.particle import Particle
from paegan.transport.controllers import BaseModelController
from paegan.transport.forcers import CachingForcer
from paegan.transport.parallel_manager import CachingDataController, Consumer
from paegan.utils.asarandom import AsaRandom
from paegan.logger import logger
import paegan.transport.export as ex


class CachingModelController(BaseModelController):

    def __init__(self, **kwargs):
        super(CachingModelController, self).__init__(**kwargs)
        self.time_chunk  = kwargs.get('time_chunk', 10)
        self.horiz_chunk = kwargs.get('horiz_chunk', 5)

    def total_task_count(self):
        # Add the CachingDataController to the number of tasks
        return self.total_particle_count() + 1

    def start_tasks(self, **kwargs):
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

            logger.info('Adding %i particles as tasks' % self.total_particle_count())

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
                          for i in range(self.nproc - 1) ]
            logger.progress((5, 'Running model'))
            for w in self.procs:
                w.start()
                logger.info('Started %s' % w.name)

            return True

        except Exception:
            logger.exception("Something didn't start correctly!")
            return False

    def __str__(self):
        return "*** CachingModelController ***"

    def setup_run(self, hydrodataset, **kwargs):

        super(CachingModelController, self).setup_run(hydrodataset, **kwargs)

        # Get the number of cores (may take some tuning) and create that
        # many workers then pass particles into the queue for the workers
        self.mgr = multiprocessing.Manager()

        # This tracks if the system is 'alive'.  Most looping whiles will check this
        # and break out if it is False.  This is True until something goes very wrong.
        self.active = self.mgr.Value('bool', True)

        # Either spin up the number of cores, or the number of tasks
        self.nproc = min(multiprocessing.cpu_count() - 1, self.total_task_count())

        # Number of tasks that we need to run.  This is decremented everytime something exits.
        self.n_run = self.mgr.Value('int', self.total_task_count())

        # The lock that controls access to the 'n_run' variable
        self.nproc_lock = self.mgr.Lock()

        # Create the task queue for all of the particles and the CachingDataController
        self.tasks = multiprocessing.JoinableQueue(self.total_task_count())

        # Create the result queue for all of the particles and the CachingDataController
        self.results = self.mgr.Queue(self.total_task_count())

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

    def listen_for_results(self, output_h5_file, total_particles):
        try:
            # Get results back from queue, test for failed particles
            return_particles = []
            retrieved = 0.
            self.error_code = 0

            logger.info("Waiting for %i particle results" % total_particles)
            while retrieved < self.total_task_count():  # One for the CachingDataController

                logger.info("looping in listen_for_results")

                try:
                    # Returns a tuple of code, result
                    code, tempres = self.results.get(timeout=240)
                except queue.Empty:
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
                        logger.progress((round((retrieved / self.total_task_count()) * 90., 1), "Particle %d finished" % tempres.uid))
                    elif tempres == "CachingDataController":
                        logger.info("CachingDataController finished")
                        logger.progress((round((retrieved / self.total_task_count()) * 90., 1), "CachingDataController finished"))
                    else:
                        logger.info("Got a strange result on results queue")
                        logger.info(str(tempres))

                    logger.info("Retrieved %i/%i results" % (int(retrieved), self.total_task_count()))

                # Relax
                time.sleep(1)

            if len(return_particles) != total_particles:
                logger.warn("Some particles failed and are not included in the output")

            # The results queue should be empty at this point
            assert self.results.empty() is True

            # Should be good to join on the tasks now that the queue is empty
            logger.info("Joining the task queue")
            self.tasks.join()
            self.tasks.close()
            self.tasks.join_thread()

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

        if self.error_code == -2:
            raise ValueError("Error in the BaseDataController (error_code was -2)")

        results = ex.ResultsPyTable(output_h5_file)
        for p in return_particles:
            for x in range(len(p.locations)):
                results.write(p.timestep_index_dump(x))
        results.compute()
        results.close()

        return

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
