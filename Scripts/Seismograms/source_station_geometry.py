import pygmt
import numpy as np

write = False

stations_file = '/scratch/agupta7/specfem/rectangular_grid/Brown_SymGroups/ISO/OUTPUT_FILES/STATIONS_FILTERED'
stations_info = np.genfromtxt(stations_file, dtype=str)

stations = stations_info[:,0]
latitudes = np.array([float(i) for i in stations_info[:,2]])
longitudes = np.array([float(i) for i in stations_info[:,3]])

mask = np.char.startswith(stations, 'C')

stations = stations[~mask]
latitudes = latitudes[~mask]/1000
longitudes = longitudes[~mask]/1000

focal_mechanism = np.array([0, 0, 75000, -4.54, 2.27, 2.27, 0, 0, 0, 16, 0, 0])

fig = pygmt.Figure()

# Set the projection and region
fig.basemap(region=[-156, 156, -156, 156], projection="X10c/10c", frame=True)

fig.plot(x=longitudes, y=latitudes, style='t0.1c')
fig.plot(x=[35], y=[75], style='t0.35c', fill='red')

# Plot the beachball with the focal mechanism
fig.meca(focal_mechanism, convention='mt', scale='6c',compressionfill='red')

if write: fig.savefig('source_station_map.png')

fig.show()