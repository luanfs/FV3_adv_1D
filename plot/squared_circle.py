import numpy as np
import matplotlib.pyplot as plt

# Python script to plot the squared circle

# directory
graphdir='../graphs/' # must exist
fig_format = 'png'

# some constants
nbfaces = 4 # number of faces in a square
rad2deg = 180.0/np.pi
deg2rad = 1.0/rad2deg
pi      = np.pi
pio2    = pi*0.5
pio4    = pio2*0.5
R       = 1.0 # circle radius
a       = R/np.sqrt(2.0)

gridtype = 2
if gridtype==2:
  aref = pio4
  Rref = 1.0
elif gridtype==0:
  aref = np.arcsin(1.0/np.sqrt(3.0))
  Rref = np.sqrt(2.0)
elif gridtype==1:
  aref = 1.0
  Rref = 1.0

N = 9
ng = 3
dx = 2*aref/N

ksd = 0
ks  = ng
ke  = ks+N
ked = ke+ng

xc = np.linspace(-aref-ng*dx,aref+ng*dx,N+1+2*ng)
if gridtype != 1:
   xc = Rref*np.tan(xc)
yc = np.ones(np.shape(xc))

xa = (xc[:ked]+xc[1:])*0.5
ya = np.ones(np.shape(xa))

xa = a*xa
ya = a*ya

xc = a*xc
yc = a*yc

rc = np.hypot(xc, a)
Xc, Yc = R*xc/rc, R*a/rc

ra = np.hypot(xa, a)
Xa, Ya = R*xa/ra, R*a/ra

# a grid
pa_X, pa_Y, pa_x, pa_y = \
np.zeros((np.shape(Xa)[0],nbfaces)), np.zeros((np.shape(Ya)[0],nbfaces)), np.zeros((np.shape(Xa)[0],nbfaces)), np.zeros((np.shape(Ya)[0],nbfaces))
pa_X[:,0], pa_Y[:,0], pa_x[:,0], pa_y[:,0]  =  Xa,  Ya,  xa,  ya
pa_X[:,1], pa_Y[:,1], pa_x[:,1], pa_y[:,1]  =  Ya, -Xa,  ya, -xa
pa_X[:,2], pa_Y[:,2], pa_x[:,2], pa_y[:,2]  = -Xa, -Ya, -xa, -ya
pa_X[:,3], pa_Y[:,3], pa_x[:,3], pa_y[:,3]  = -Ya,  Xa, -ya,  xa

# c grid
pc_X, pc_Y, pc_x, pc_y = \
np.zeros((np.shape(Xc)[0],nbfaces)), np.zeros((np.shape(Yc)[0],nbfaces)), np.zeros((np.shape(Xc)[0],nbfaces)), np.zeros((np.shape(Yc)[0],nbfaces))
pc_X[:,0], pc_Y[:,0], pc_x[:,0], pc_y[:,0]  =  Xc,  Yc,  xc,  yc
pc_X[:,1], pc_Y[:,1], pc_x[:,1], pc_y[:,1]  =  Yc, -Xc,  yc, -xc
pc_X[:,2], pc_Y[:,2], pc_x[:,2], pc_y[:,2]  = -Xc, -Yc, -xc, -yc
pc_X[:,3], pc_Y[:,3], pc_x[:,3], pc_y[:,3]  = -Yc,  Xc, -yc,  xc

#------------------------------------------------------------------
# Figure quality
dpi = 100
plt.figure(figsize=(1000/dpi, 1000/dpi), dpi=dpi)

#------------------------------------------------------------------
Colors = ('blue', 'green', 'red', 'orange' )
for p in range(0,nbfaces):
   Color = Colors[p]
   plt.plot(pc_X[ks:ke+1,p], pc_Y[ks:ke+1,p], marker='s',  markersize=10, color = Color)
   plt.plot(pa_X[ks:ke,p], pa_Y[ks:ke,p], marker='o', markersize=5, color = Color)
   plt.plot(pc_x[ks:ke+1,p], pc_y[ks:ke+1,p], marker='s', markersize=5, color = Color)

# Drawing lines connecting corresponding points
p = 0
Color = Colors[p]
for i in range(len(xc[ks:ke+1])):
    plt.plot([pc_X[ks+i,p], pc_x[ks+i,p]], [pc_Y[ks+i,p], pc_y[ks+i,p]], color=Color, linestyle='--')
    plt.plot([0, pc_x[ks+i,p]], [0, pc_y[ks+i,p]], color=Color, linestyle='--')

# plot ghost cells for p=0
for k in range(ksd,ks):
   plt.plot(pa_X[k,p], pa_Y[k,p], marker='*',  markersize=7, color = Color)
   plt.plot(pc_X[k,p], pc_Y[k,p], marker='*',  markersize=10, color = Color)

for k in range(ke,ked):
   plt.plot(pa_X[k,p], pa_Y[k,p], marker='*',  markersize=7, color = Color)
   plt.plot(pc_X[k+1,p], pc_Y[k+1,p], marker='*',  markersize=10, color = Color)

#------------------------------------------------------------------
plt.xlim(-1.10*R, 1.10*R) 
plt.xlim(-1.10*R, 1.10*R)

plt.axhline(0, color='black')
plt.axvline(0, color='black')

# Title
title = 'g'+str(gridtype)+' - sc'+str(N)
plt.title(title)

# Save the figure
plt.savefig(graphdir+'g'+str(gridtype)+'_sc'+str(N)+'.'+fig_format, format=fig_format)

#plt.show()
