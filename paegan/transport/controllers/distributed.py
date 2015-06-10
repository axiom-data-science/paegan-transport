import sys
import json
import redis
import logging
import multiprocessing
import traceback
from datetime import datetime
from urlparse import urlparse

from paegan.transport.forcers import BaseForcer
from paegan.transport.controllers import BaseModelController, Runner
from paegan.logger import logger
import paegan.transport.export as ex


def particle_runner(part, model):

    from paegan.logger import logger
    logger.setLevel(logging.PROGRESS)

    from paegan.logger.redis_handler import RedisHandler
    rhandler = RedisHandler(model.redis_log_channel, model.redis_url)
    rhandler.setLevel(logging.PROGRESS)
    logger.addHandler(rhandler)

    try:
        redis_connection = redis.from_url(model.redis_url)
        forcer = BaseForcer(model.hydrodataset,
                            particle=part,
                            common_variables=model.common_variables,
                            times=model.times,
                            start_time=model.start,
                            models=model._models,
                            release_location_centroid=model.reference_location.point,
                            usebathy=model._use_bathymetry,
                            useshore=model._use_shoreline,
                            usesurface=model._use_seasurface,
                            reverse_distance=model.reverse_distance,
                            bathy_path=model.bathy_path,
                            shoreline_path=model.shoreline_path,
                            shoreline_feature=model.shoreline_feature,
                            time_method=model.time_method,
                            redis_url=model.redis_url,
                            redis_results_channel=model.redis_results_channel,
                            shoreline_index_buffer=model.shoreline_index_buffer)
        forcer.run()
    except Exception:
        logger.exception(traceback.format_exc())
        redis_connection.publish(model.redis_results_channel, json.dumps({"status" : "FAILED", "uid" : part.uid }))
    else:
        redis_connection.publish(model.redis_results_channel, json.dumps({"status" : "COMPLETED", "uid" : part.uid }))


class DistributedModelController(BaseModelController):

    def __init__(self, **kwargs):
        super(DistributedModelController, self).__init__(**kwargs)
        self.thread_result_listener = True

    def start_tasks(self, **kwargs):
        logger.progress((5, 'Running model'))
        rc = redis.from_url(self.redis_url)
        if kwargs.get('task_queue_call'):
            for p in self.particles:
                try:
                    kwargs.get('task_queue_call')(func=particle_runner, args=(p, self,))
                except Exception:
                    logger.exception(traceback.format_exc())
                    rc.publish(self.redis_results_channel, json.dumps({"status" : "FAILED", "uid" : p.uid }))
            return True
        else:
            tasks = []
            for p in self.particles:
                f = BaseForcer(self.hydrodataset,
                               particle=p,
                               common_variables=self.common_variables,
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
                tasks.append(f)
            self.pool = multiprocessing.Pool()
            return self.pool.imap(Runner(), tasks)

        return False

    def setup_run(self, hydrodataset, **kwargs):
        from paegan.logger.redis_handler import RedisHandler
        self.redis_url             = kwargs.get("redis_url")
        self.redis_log_channel     = kwargs.get("redis_log_channel")
        self.redis_results_channel = kwargs.get("redis_results_channel")
        rhandler = RedisHandler(self.redis_log_channel, self.redis_url)
        rhandler.setLevel(logging.PROGRESS)
        logger.addHandler(rhandler)

        super(DistributedModelController, self).setup_run(hydrodataset, **kwargs)

    def listen_for_results(self, output_h5_file, total_particles):
            # Create output file (hdf5)
            particles_finished = 0
            results = ex.ResultsPyTable(output_h5_file)

            res = urlparse(self.redis_url)
            redis_pool = redis.ConnectionPool(host=res.hostname, port=res.port, db=res.path[1:])
            r = redis.Redis(connection_pool=redis_pool)

            pubsub = r.pubsub()
            pubsub.subscribe(self.redis_results_channel)
            for msg in pubsub.listen():

                if msg['type'] != "message":
                    continue

                if msg["data"] == "FINISHED":
                    break

                try:
                    json_msg = json.loads(msg["data"])
                    if json_msg.get("status", None):
                        #  "COMPLETED" or "FAILED" when a particle finishes
                        particles_finished += 1
                        percent_complete = 90. * (float(particles_finished) / float(total_particles)) + 5  # Add the 5 progress that was used prior to the particles starting (controller)
                        r.publish(self.redis_log_channel, json.dumps({"time" : datetime.utcnow().isoformat(), "level" : "progress", "value" : percent_complete, "message" : "Particle #%s %s!" % (particles_finished, json_msg.get("status"))}))
                        if particles_finished == total_particles:
                            break
                    else:
                        # Write to HDF file
                        results.write(json_msg)
                except Exception:
                    pass

            pubsub.close()
            results.compute()
            results.close()
            redis_pool.disconnect()
            sys.exit()
