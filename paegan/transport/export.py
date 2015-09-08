import os
import glob
import math
import zipfile
from shapely.geometry import MultiPoint, Point, mapping
from datetime import datetime

# NetCDF
import netCDF4
import pytz

from fiona import collection

from collections import OrderedDict

from paegan.logger import logger

import geojson

import tables
from tables import *


# Pytables representation of a model run
class ModelResultsTable(IsDescription):
    particle    = UInt8Col()
    time        = Time32Col()
    latitude    = Float32Col()
    longitude   = Float32Col()
    depth       = Float32Col()
    u_vector    = Float32Col()
    v_vector    = Float32Col()
    w_vector    = Float32Col()
    temperature = Float32Col()
    salinity    = Float32Col()
    age         = Float32Col()
    lifestage   = UInt8Col()
    progress    = Float32Col()
    settled     = BoolCol()
    halted      = BoolCol()
    dead        = BoolCol()


class ResultsPyTable(object):
    def __init__(self, output_file):
        self._file  = open_file(output_file, mode="w", title="Model run output")
        self._root  = self._file.create_group("/", "trajectories", "Trajectory Data")
        self._table = self._file.create_table(self._root, "model_results", ModelResultsTable, "Model Results")
        self._table.autoindex = False
        self._table.cols.particle.create_index()
        self._table.cols.time.create_index()
        self._table.cols.latitude.create_index()
        self._table.cols.longitude.create_index()

    def write(self, data):
        record = self._table.row
        for k, v in data.items():
            try:
                record[k] = v
            except Exception:
                # No column named "k", so don't add the data
                pass

        record.append()

    def trackline(self):
        pass

    def metadata(self):
        pass

    def compute(self):
        self.trackline()
        self.metadata()

    def close(self):
        self._table.flush()
        self._table.reindex()
        self._file.close()


class Export(object):
    @classmethod
    def export(cls, **kwargs):
        raise("Please implement the export method of your Export class.")


class H5TracklineWithPoints(Export):
    @classmethod
    def export(cls, folder, h5_file):
        with tables.open_file(h5_file, mode="r") as h5:
            table = h5.root.trajectories.model_results
            timestamps = sorted(list(set([ x["time"] for x in table.iterrows() ])))

            pts = []
            features = []
            for i, ts in enumerate(timestamps):
                points = MultiPoint([ Point(x['longitude'], x['latitude']) for x in table.where("""time == %s""" % ts) if x["latitude"] and x["longitude"] ])
                cp = points.centroid.coords[0]
                geo_pt = geojson.Point(cp)
                pts.append(cp)
                feat = geojson.Feature(id=i, geometry=geo_pt, properties={ "time" : datetime.utcfromtimestamp(ts).replace(tzinfo=pytz.utc).isoformat() })
                features.append(feat)

            geo_ls = geojson.LineString(pts)
            features.append(geojson.Feature(geometry=geo_ls, id='path'))

            fc = geojson.FeatureCollection(features)

            if not os.path.exists(folder):
                os.makedirs(folder)
            filepath = os.path.join(folder, "full_trackline.geojson")
            with open(filepath, "wb") as r:
                r.write(geojson.dumps(fc).encode('utf-8'))


class H5Trackline(Export):
    @classmethod
    def export(cls, folder, h5_file):
        with tables.open_file(h5_file, mode="r") as h5:
            table = h5.root.trajectories.model_results
            timestamps = sorted(list(set([ x["time"] for x in table.iterrows() ])))

            pts = []
            for i, ts in enumerate(timestamps):
                points = MultiPoint([ Point(x['longitude'], x['latitude']) for x in table.where("""time == %s""" % ts) if x["latitude"] and x["longitude"] ])
                pts.append(points.centroid.coords[0])
            geo_ls = geojson.LineString(pts)
            feat = geojson.Feature(geometry=geo_ls, id='path')

            if not os.path.exists(folder):
                os.makedirs(folder)
            filepath = os.path.join(folder, "simple_trackline.geojson")
            with open(filepath, "wb") as r:
                r.write(geojson.dumps(feat).encode('utf-8'))


class H5ParticleTracklines(Export):
    @classmethod
    def export(cls, folder, h5_file):
        with tables.open_file(h5_file, mode="r") as h5:
            table = h5.root.trajectories.model_results
            particles = sorted(list(set([ x["particle"] for x in table.iterrows() ])))

            features = []
            for puid in particles:
                points = [ (x["time"], (x['longitude'], x['latitude'])) for x in table.where("""particle == %s""" % puid) if x["latitude"] and x["longitude"] ]

                geo_ls = geojson.LineString( [ x[1] for x in points ] )
                times  = [ datetime.utcfromtimestamp(x[0]).replace(tzinfo=pytz.utc).isoformat() for x in points ]

                feat = geojson.Feature(geometry=geo_ls, id=puid, properties={ "particle" : puid, "times" : times })
                features.append(feat)

            fc = geojson.FeatureCollection(features)

            if not os.path.exists(folder):
                os.makedirs(folder)
            filepath = os.path.join(folder, "particle_tracklines.geojson")
            with open(filepath, "wb") as r:
                r.write(geojson.dumps(fc).encode('utf-8'))


class H5ParticleMultiPoint(Export):
    @classmethod
    def export(cls, folder, h5_file):
        with tables.open_file(h5_file, mode="r") as h5:
            table = h5.root.trajectories.model_results
            particles = sorted(list(set([ x["particle"] for x in table.iterrows() ])))

            features = []
            for puid in particles:
                points = [ (x["time"], (x['longitude'], x['latitude'])) for x in table.where("""particle == %s""" % puid) if x["latitude"] and x["longitude"] ]

                geo_mp = geojson.MultiPoint( [ x[1] for x in points ] )
                times  = [ x[0] for x in points ]

                feat = geojson.Feature(geometry=geo_mp, id=puid, properties={ "particle" : puid, "time" : times })
                features.append(feat)

            fc = geojson.FeatureCollection(features)

            if not os.path.exists(folder):
                os.makedirs(folder)
            filepath = os.path.join(folder, "particle_multipoint.geojson")
            with open(filepath, "wb") as r:
                r.write(geojson.dumps(fc).encode('utf-8'))


class H5GDALShapefile(Export):
    @classmethod
    def export(cls, folder, h5_file):
        shape_schema = {'geometry': 'Point',
                        'properties': OrderedDict([('particle',     'int'),
                                                   ('date',         'str'),
                                                   ('latitude',     'float'),
                                                   ('longitude',    'float'),
                                                   ('depth',        'float'),
                                                   ('u_vector',     'float'),
                                                   ('v_vector',     'float'),
                                                   ('w_vector',     'float'),
                                                   ('temp',         'float'),
                                                   ('salinity',     'float'),
                                                   ('age',          'float'),
                                                   ('settled',      'str'),
                                                   ('dead',         'str'),
                                                   ('halted',       'str'),
                                                   ('lifestage',    'int'),
                                                   ('progress',     'float')])}
        shape_crs = {'no_defs': True, 'ellps': 'WGS84', 'datum': 'WGS84', 'proj': 'longlat'}

        if not os.path.exists(folder):
            os.makedirs(folder)
        filepath = os.path.join(folder, "particles.shp")

        with tables.open_file(h5_file, mode="r") as h5:
            table = h5.root.trajectories.model_results

            with collection(filepath, "w", driver='ESRI Shapefile', schema=shape_schema, crs=shape_crs) as shape:
                for r in table.iterrows():
                    shape.write({'geometry': mapping(Point(r["longitude"], r["latitude"])),
                                 'properties': OrderedDict([('particle', r["particle"]),
                                                            ('date', str(datetime.utcfromtimestamp(r["time"]).isoformat())),
                                                            ('latitude', float(r["latitude"])),
                                                            ('longitude', float(r["longitude"])),
                                                            ('depth', float(r["depth"])),
                                                            ('temp', float(r["temperature"])),
                                                            ('salinity', float(r["salinity"])),
                                                            ('u_vector', float(r["u_vector"])),
                                                            ('v_vector', float(r["v_vector"])),
                                                            ('w_vector', float(r["w_vector"])),
                                                            ('settled', str(r["settled"])),
                                                            ('dead', str(r["dead"])),
                                                            ('halted', str(r["halted"])),
                                                            ('age', float(r["age"])),
                                                            ('lifestage' , int(r["lifestage"])),
                                                            ('progress' , float(r["progress"]))])})

        # Zip the output
        shpzip = zipfile.ZipFile(os.path.join(folder, "h5shape.shp.zip"), mode='w')
        for f in glob.glob(os.path.join(folder, "particles*")):
            shpzip.write(f, os.path.basename(f))
            os.remove(f)
        shpzip.close()


class NetCDF(Export):
    @classmethod
    def export(cls, folder, particles, datetimes, summary, **kwargs):
        """
            Export particle data to CF trajectory convention
            netcdf file
        """
        time_units = 'seconds since 1990-01-01 00:00:00'

        # Create netcdf file, overwrite existing
        if not os.path.exists(folder):
            os.makedirs(folder)

        filepath = os.path.join(folder, 'trajectories.nc')
        nc = netCDF4.Dataset(filepath, 'w')
        # Create netcdf dimensions
        nc.createDimension('time', None)
        nc.createDimension('particle', None)

        fillvalue = -9999.9

        # Create netcdf variables
        time = nc.createVariable('time', 'i', ('time',))
        part = nc.createVariable('particle', 'i', ('particle',))
        depth = nc.createVariable('depth', 'f', ('time', 'particle'))
        lat = nc.createVariable('lat', 'f', ('time', 'particle'), fill_value=fillvalue)
        lon = nc.createVariable('lon', 'f', ('time', 'particle'), fill_value=fillvalue)
        salt = nc.createVariable('salt', 'f', ('time', 'particle'), fill_value=fillvalue)
        temp = nc.createVariable('temp', 'f', ('time', 'particle'), fill_value=fillvalue)
        u = nc.createVariable('u', 'f', ('time', 'particle'), fill_value=fillvalue)
        v = nc.createVariable('v', 'f', ('time', 'particle'), fill_value=fillvalue)
        w = nc.createVariable('w', 'f', ('time', 'particle'), fill_value=fillvalue)
        settled = nc.createVariable('settled', 'f', ('time', 'particle'), fill_value=fillvalue)
        dead = nc.createVariable('dead', 'f', ('time', 'particle'), fill_value=fillvalue)
        halted = nc.createVariable('halted', 'f', ('time', 'particle'), fill_value=fillvalue)

        # Loop through locations in each particle,
        # add to netcdf file
        for j, particle in enumerate(particles):
            part[j] = particle.uid
            i = 0

            normalized_locations = particle.normalized_locations(datetimes)
            normalized_temps = [x if x is not None and not math.isnan(x) else fillvalue for x in particle.temps]
            normalized_salts = [x if x is not None and not math.isnan(x) else fillvalue for x in particle.salts]
            normalized_u = [x if x is not None and not math.isnan(x) else fillvalue for x in particle.u_vectors]
            normalized_v = [x if x is not None and not math.isnan(x) else fillvalue for x in particle.v_vectors]
            normalized_w = [x if x is not None and not math.isnan(x) else fillvalue for x in particle.w_vectors]
            normalized_settled = [x if x is not None and not math.isnan(x) else fillvalue for x in particle.settles]
            normalized_dead = [x if x is not None and not math.isnan(x) else fillvalue for x in particle.deads]
            normalized_halted = [x if x is not None and not math.isnan(x) else fillvalue for x in particle.halts]

            if len(normalized_locations) != len(normalized_temps):
                logger.info("No temperature being added to netcdf.")
                # Create list of 'fillvalue' equal to the length of locations
                normalized_temps = [fillvalue] * len(normalized_locations)

            if len(normalized_locations) != len(normalized_salts):
                logger.info("No salinity being added to netcdf.")
                # Create list of 'fillvalue' equal to the length of locations
                normalized_salts = [fillvalue] * len(normalized_locations)

            if len(normalized_locations) != len(normalized_u):
                logger.info("No U being added to netcdf.")
                # Create list of 'fillvalue' equal to the length of locations
                normalized_u = [fillvalue] * len(normalized_locations)

            if len(normalized_locations) != len(normalized_v):
                logger.info("No V being added to netcdf.")
                # Create list of 'fillvalue' equal to the length of locations
                normalized_v = [fillvalue] * len(normalized_locations)

            if len(normalized_locations) != len(normalized_w):
                logger.info("No W being added to netcdf.")
                # Create list of 'fillvalue' equal to the length of locations
                normalized_w = [fillvalue] * len(normalized_locations)

            if len(normalized_locations) != len(normalized_settled):
                logger.info("No Settled being added to shapefile.")
                # Create list of 'fillvalue' equal to the length of locations
                normalized_settled = [fillvalue] * len(normalized_locations)

            if len(normalized_locations) != len(normalized_dead):
                logger.info("No Dead being added to shapefile.")
                # Create list of 'fillvalue' equal to the length of locations
                normalized_dead = [fillvalue] * len(normalized_locations)

            if len(normalized_locations) != len(normalized_halted):
                logger.info("No Halted being added to shapefile.")
                # Create list of 'fillvalue' equal to the length of locations
                normalized_halted = [fillvalue] * len(normalized_locations)

            for loc, _temp, _salt, _u, _v, _w, _settled, _dead, _halted in zip(normalized_locations, normalized_temps, normalized_salts, normalized_u, normalized_v, normalized_w, normalized_settled, normalized_dead, normalized_halted):

                if j == 0:
                    time[i] = int(round(netCDF4.date2num(loc.time, time_units)))
                depth[i, j] = loc.depth
                lat[i, j] = loc.latitude
                lon[i, j] = loc.longitude
                salt[i, j] = _salt
                temp[i, j] = _temp
                u[i, j] = _u
                v[i, j] = _v
                w[i, j] = _w
                settled[i, j] = _settled
                dead[i, j] = _dead
                halted[i, j] = _halted
                i += 1

        # Variable attributes
        depth.coordinates = "time particle lat lon"
        depth.standard_name = "depth_below_sea_surface"
        depth.units = "m"
        depth.POSITIVE = "up"
        depth.positive = "up"

        salt.coordinates = "time particle lat lon"
        salt.standard_name = "sea_water_salinity"
        salt.units = "psu"

        temp.coordinates = "time particle lat lon"
        temp.standard_name = "sea_water_temperature"
        temp.units = "degrees_C"

        u.coordinates = "time particle lat lon"
        u.standard_name = "eastward_sea_water_velocity"
        u.units = "m/s"

        v.coordinates = "time particle lat lon"
        v.standard_name = "northward_sea_water_velocity"
        v.units = "m/s"

        w.coordinates = "time particle lat lon"
        w.standard_name = "upward_sea_water_velocity"
        w.units = "m/s"

        settled.coordinates = "time particle lat lon"
        settled.description = "Is the particle settled"
        settled.standard_name = "particle_settled"

        dead.coordinates = "time particle lat lon"
        dead.description = "Is the particle dead"
        dead.standard_name = "particle_dead"

        halted.coordinates = "time particle lat lon"
        halted.description = "Is the particle prevented from being forced by currents"
        halted.standard_name = "particle_halted"

        time.units = time_units
        time.standard_name = "time"

        lat.units = "degrees_north"

        lon.units = "degrees_east"

        part.cf_role = "trajectory_id"

        # Global attributes
        nc.featureType = "trajectory"
        nc.summary = str(summary)
        for key in kwargs:
            nc.__setattr__(key, kwargs.get(key))

        nc.sync()
        nc.close()
