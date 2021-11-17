using NCDatasets
using Dates
ds = NCDataset("dataset-duacs-rep-global-merged-allsat-phy-l4_1633328415183.nc");
# Latitude is from 45.875 -to- 25.875
lat = ds["latitude"][:];
# Longitude is form -4.125 -to- -84.875
lon = ds["longitude"][:];
# Time is from 2018-07-01 -to- 2020-06-03 with the stepsize of 01 days
time = ds["time"][:];
# Absolute geostrophic velocity: zonal component
u = ds["ugos"][:,:,:];
# Absolute geostrophic velocity: meridian component
v = ds["vgos"][:,:,:];
# Absolute dynamic topography is the sea surface height above geoid
adt = ds["adt"][:,:,:];
# Sea level anomaly is the sea surface height above mean sea surface
sla = ds["sla"][:,:,:];

# get the length and thickness in each direction
nx = length(lon) # no. of grid points in longitude direction
ny = length(lat) # no. of grid points in latitude direction
nt = length(time)# no. of time steps: per day

dx = (lon[2]-lon[1]) # grid thickness along longitude
dy = (lat[2]-lat[1]) # grid thickness along latitude

# reshape the velocity component from shape (nx,ny,nt) to (ny,nx,nt)
U = zeros(ny,nx,nt)
V = zeros(ny,nx,nt)
for i = 1:nt
    U[:,:,i] = u[:,:,i]'
    V[:,:,i] = v[:,:,i]'
end

# grid points: generate the 2d matrix containing grid point
XG = repeat(lon',length(lat),1)
YG = repeat(lat,1, length(lon))

# compute grid thickness
DX = dx*ones(ny,nx)
DY = dy*ones(ny,nx)

# grid centers: generate the 2d matrix containing grid centers
XC = XG .+ DX./2
YC = YG .+ DY./2

# compute the divervence of the velocity field
r = 6371.0   # radius of the earth
DUDX = zeros(ny,nx-1,nt)
DVDY = zeros(ny-1,nx,nt)
DIV = zeros(ny-1,nx-1,nt)
for i in 1:nt
    DUDX[:,:,i] = diff(U[1:end,1:end,i],dims=2) ./ DX[1:end,1:end-1]
    DVDY[:,:,i] = diff(V[1:end,1:end,i].*(cos.(DY)),dims=1) ./ DY[1:end-1,1:end]
    DIV[:,:,i] = (1 ./ (r*cos.(DY[1:end-1,1:end-1]))).*(DUDX[1:end-1,1:end,i] + DVDY[1:end,1:end-1,i])
end

# compute surface area
A = r*r*cos.(YG[1:end-1,1:end-1]).*DY[1:end-1,1:end-1].*DX[1:end-1,1:end-1]
sumA = sum(sum(A))

# Plot the divergence of velocity
using Plots
gr()
# subsets of the grid
nx_sub = 1:1:nx
ny_sub = 1:1:ny
xc = lon[1:end-1] .+ dx/2 # longitude grid centers
yc = lat[1:end-1] .+ dy/2 # latitude grid centers
begin
    anim = @animate for i=1:nt # at time (In day)
    #GR.setlinewidth(0)
    contourf(xc,yc,DIV[:,:,i],linewidth=0,levels=50)
    title!(string(Date(time[i])))
    end
    gif(anim, "/tmp/div.gif", fps = 1)
end

# Plot subset of the velocity field
x = lon[nx_sub] # longitude grid points
y = lat[ny_sub] # latitude grid points
begin
    anim = @animate for i=1:nt # at time (In day)
    U1 = U[ny_sub,nx_sub,i]
    U2 = V[ny_sub,nx_sub,i]
    W = sqrt.(U1.*U1 + U2.*U2)
    #GR.setlinewidth(0)
    contourf(x,y,W,linewidth=0,levels=50)
    title!(string(Date(time[i])))
    end
    gif(anim, "/tmp/vel.gif", fps = 1)
end
