from paegan.transport.controllers import BaseModelController, Runner
from paegan.transport.forcers import BaseForcer
from paegan.logger import logger

import paegan.transport.export as ex


class IPythonClusterModelController(BaseModelController):

    def listen_for_results(self, output_h5_file, total_particles):
        logger.info("Waiting for %i particle results" % total_particles)

        particles = []
        retrieved = 0

        while retrieved < total_particles:
            try:
                # IPython parallel View
                # self.result is an AsyncMapResult
                from IPython.parallel import TimeoutError
                try:
                    new_particles = self.result.get(timeout=1)
                except TimeoutError:
                    pass    # this is fine, get incremental progress below
                else:
                    particles = new_particles

                # progress is absolute, not incremental
                retrieved = self.result.progress
            except StopIteration:
                assert retrieved >= total_particles
                break
            except:
                logger.exception("Particle has FAILED!!")
                continue

            # We multiply by 90 here to save 10% for the exporting
            logger.progress((round((float(retrieved) / total_particles) * 90., 1), "%s Particle(s) complete" % retrieved))

        results = ex.ResultsPyTable(output_h5_file)
        for p in particles:
            for x in range(len(p.locations)):
                results.write(p.timestep_index_dump(x))
        results.compute()
        results.close()

        return

    def start_tasks(self, **kwargs):

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

            logger.progress((5, 'Running model'))
            return self.pool.map_async(Runner(), tasks)

        except Exception:
            logger.exception("Something didn't start correctly!")
            raise
