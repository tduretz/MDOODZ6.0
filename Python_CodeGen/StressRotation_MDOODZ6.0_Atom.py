from sympy import *
import numpy as np
init_printing()

# %% ANALYTICAL STRESS ROTATION

# Definitions
angle, txx, txz, Nx, Nz = symbols( 'angle, txx, txz, Nx, Nz' )
R   = Matrix([[cos(angle), -sin(angle)], [sin(angle), cos(angle)]]);
T   = Matrix([ [txx, txz], [txz, -txx] ])
N   = Matrix([ [Nx], [Nz] ])

N1  = R * N

# Rotate
T1  = R * T * R.transpose()

print ('nx = ' + ccode(N1[0] ) + ';')
print ('nz = ' + ccode(N1[1] ) + ';')
for i in range(1):
    print ('particles->sxxd[k] = ' + ccode(T1[0,0] ) + ';')
    print ('particles->sxz[k]  = ' + ccode(T1[0,1] ) + ';')
# %% UPPER CONVECTED STRESS DERIVATIVE
dudx,dvdx,dudz,dvdz = symbols( 'dudx,dvdx,dudz,dvdz' )
L   = Matrix([[dudx,dudz], [dvdx,dvdz]]); # (grad u) = du_i / dx_j
duc = - L*T - T*L.transpose()
print(duc[0,0])
print(duc[1,0])

# %% UPPER CONVECTED STRESS DERIVATIVE
d_dx,d_dz, u, v = symbols('d_dx,d_dz, u, v')
grad = Matrix([d_dx, d_dz])
V   = Matrix([u, v])
grad*V.transpose()
# %%
# Update formula for the director vector from Mülhaus
w12,e11,e12,nx,nz = symbols('w12,e11,e12,nx,nz')
W = Matrix([[  0, w12], [-w12,   0]]);   # vorticity tensor
D = Matrix([[e11, e12], [e12, -e11]]);   # deviatoric strain rate tensor
Director      = Matrix([[nx], [nz]]);
sca1          = (Director.transpose()*Director)
Director_dot  = W*Director - np.multiply(sca1,D*Director)
sca2        = (Director.transpose()*( D*Director));
Director_dot  = Director_dot + np.multiply(sca2,Director)
print( 'ndotx = ' + ccode(Director_dot[0]) + ';')
print( 'ndotz = ' + ccode(Director_dot[1]) + ';')

# Director_dot = W*Director - D*Director*(Director'*Director) + (Director'*(D*Director))*Director;


# %%
print(Director.transpose()*Director)
# %%
