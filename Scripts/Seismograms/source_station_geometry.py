from pathlib import Path

import numpy as np
from obspy import read_events
import pygmt

write = False

# select station to plot
stn = "R33-25"

# define file paths
root_dir = Path.cwd().resolve().parent.parent
data_dir = root_dir / "Scripts" / "Data" / "SPECFEM3D_CARTESIAN" / "DATA"
stations_file = data_dir / "STATIONS"
source_file = data_dir / "CMTSOLUTION"
mesh_par_file = data_dir / "meshfem3D_files" / "Mesh_Par_file"

# read station information

stations_info = np.genfromtxt(stations_file, dtype=str)

station_names = stations_info[:,0]
station_latitudes = np.array([float(i) for i in stations_info[:,2]])
station_longitudes = np.array([float(i) for i in stations_info[:,3]])

mask = np.char.startswith(station_names, 'C')

station_names = station_names[~mask]
station_latitudes = station_latitudes[~mask]/1000
station_longitudes = station_longitudes[~mask]/1000

station_x_coordinate = station_longitudes[station_names == stn]
station_y_coordinate = station_latitudes[station_names == stn]

# read source information
events = read_events(source_file)

event = events[0]
event_origin = event.preferred_origin()
event_latitude = event_origin.latitude / 1000
event_longitude = event_origin.longitude / 1000
event_depth = event_origin.depth / 1000

moment_tensor = event.preferred_focal_mechanism().moment_tensor.tensor
mrr = moment_tensor.m_rr
mtt = moment_tensor.m_tt
mpp = moment_tensor.m_pp
mrt = moment_tensor.m_rt
mrp = moment_tensor.m_rp
mtp = moment_tensor.m_tp

focal_mechanism = np.array([event_longitude, event_latitude, event_depth, mrr, mtt, mpp, mrt, mrp, mtp, 1, 0, 0])

# read mesh information
with open(mesh_par_file, 'r') as fid:
    for line in fid:
        if 'LATITUDE_MIN' in line:
            latitude_min = float(line.split()[2]) / 1000
        if 'LATITUDE_MAX' in line:
            latitude_max = float(line.split()[2]) / 1000
        if 'LONGITUDE_MIN' in line:
            longitude_min = float(line.split()[2]) / 1000
        if 'LONGITUDE_MAX' in line:
            longitude_max = float(line.split()[2]) / 1000

# Plot the source and stations
fig = pygmt.Figure()

fig.basemap(region=[longitude_min, longitude_max, latitude_min, latitude_max],
            projection="X20c/20c", frame=True)

fig.plot(x=station_longitudes, y=station_latitudes, style='t0.125c')
fig.plot(x=station_x_coordinate, y=station_y_coordinate, style='t0.35c', fill='red')
fig.meca(focal_mechanism, convention='mt', scale='1.75c', compressionfill='red')

if write: fig.savefig('source_station_map.png')
fig.show()