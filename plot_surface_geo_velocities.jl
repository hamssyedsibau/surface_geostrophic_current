using NCDatasets
ds = NCDataset("dataset-duacs-rep-global-merged-allsat-phy-l4_1633328415183.nc")
# Latitude is from 45.875 -to- 25.875
lat = ds["latitude"]
# Longitude is form -4.125 -to- -84.875
lon = ds["longitude"]
# Time is from 2018-07-01 -to- 2020-06-03 with the stepsize of 01 days
time = ds["time"]
# Absolute geostrophic velocity: zonal component
u = ds["ugos"]
# Absolute geostrophic velocity: meridian component
v = ds["vgos"]
# Absolute dynamic topography is the sea surface height above geoid
adt = ds["adt"]
# Sea level anomaly is the sea surface height above mean sea surface
sla = ds["sla"]

nx = length(lon) # no. of grid: longitude
ny = length(lat) # no. of grid: latitude
nt = length(time)# no. of time steps: per day

dx = (lon[2]-lon[1]) / 2 # grid thickness along longitude
dy = (lat[2]-lat[1]) / 2 # grid thickness along latitude

# grid points: generate the 2d matrix containing grid point
XG = repeat(lon,1, length(lat))
YG = repeat(lat', length(lon),1)

# grid centers: generate the 2d matrix containing grid centers
XC = XG[1:end-1,1:end-1] .+ dx
YC = YG[1:end-1,1:end-1] .+ dy


# Plot subset of the velocity fields
i = 1 # at time (In day)
nx_sub = 1:4:nx
ny_sub = 1:2:ny

using Plots
gr()
U = u[nx_sub,ny_sub,i]
V = v[nx_sub,ny_sub,i]
GR.setlinewidth(0)
p1 = contour(lon[nx_sub], lat[ny_sub], U', fill = true, linewidth=0,levels=10)
p2 = contour(lon[nx_sub], lat[ny_sub], V', fill = true, linewidth=0,levels=10)
plot(p1,p2)

x = XG[nx_sub,ny_sub]|>vec
y = YG[nx_sub,ny_sub]|>vec
vx = U |> vec
vy = V |> vec
quiver(x, y, quiver=(vx, vy), arrow=:closed, arrowsize=0.1)


# average flow
#U = u[nx_sub,ny_sub,:]
#V = v[nx_sub,ny_sub,:]
#u_mean = sum(U,dims=3)/nt
#v_mean = sum(V,dims=3)/nt

#using Plots
#gr()
#U = u_mean[:,:,1]
#V = v_mean[:,:,1]
#GR.setlinewidth(0)
#p1 = contour(lon[nx_sub], lat[ny_sub], U', fill = true, linewidth=0,levels=10)
#p2 = contour(lon[nx_sub], lat[ny_sub], V', fill = true, linewidth=0,levels=10)
#plot(p1,p2)
