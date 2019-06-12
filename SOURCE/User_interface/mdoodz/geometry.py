#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 11:40:39 2019

@author: abauville
"""



import numpy as np
from mdoodz.scaling import Scaling
from mdoodz.utils import use_numba, maybe_numba

   
# =============================================================================
#
#                           Class Rotation
#
# =============================================================================
class Rotation():
    def __init__(self,xc=0.0,zc=0.0,angle=0.0):
        self.xc = xc
        self.zc = zc
        self.angle = angle
        
    def apply_to_arrays(self,x,z,reverse=False):
        xtemp = x-self.xc
        ztemp = z-self.zc
        
        if reverse:
            angle = -self.angle
        else:
            angle = self.angle                
        
        x[:] = np.cos(angle)*xtemp - np.sin(angle)*ztemp + self.xc
        z[:] = np.sin(angle)*xtemp + np.cos(angle)*ztemp + self.zc
        
    def apply_to_polygon(self,polygon,reverse=False):
        if not isinstance(polygon,Polygon):
            raise ValueError("'Polygon' must be an instance of 'Geometry.Polygon'")
            
        self.apply_to_arrays(polygon.x,polygon.z,reverse)

    def copy(self):
        return Rotation(xc=self.xc,zc=self.zc,angle=self.angle)



# =============================================================================
#
#                           Class Polygon
#
# =============================================================================
class Polygon():
    def __init__(self,x,z,phase=0,scaled=False):
        if isinstance(x,list) or isinstance(x,tuple):
            self.x = np.array(x)
        elif isinstance(x,np.ndarray):
            self.x = x
        else:
            raise TypeError("x must be a list, tuple or numpy.ndarray")
            
        if isinstance(z,list) or isinstance(z,tuple):
            self.z = np.array(z)
        elif isinstance(x,np.ndarray):
            self.z = z
        else:
            raise TypeError("z must be a list, tuple or numpy.ndarray")


        self.phase = phase
        self.scaled = False

        self.rotation = Rotation()
        
    def rotate(self,angle,xc=0.0,zc=0.0,reverse=False):
        """
            This function rotates the polygon by a angle "angle" (in radians) around
            the point of coordinates (xc,zc)
            
        """
        self.rotation = Rotation(xc=xc,zc=zc,angle=angle)
        self.rotation.apply_to_polygon(self,reverse)
        
    def derotate(self):
        """
            This function rotates the polygon by a angle -"angle" (in radians) around
            the point of coordinates (xc,zc)
            
        """
        self.rotation.apply_to_polygon(self,reverse=True)
        self.rotation = Rotation()
        

    def scale(self,scaling):
        if not isinstance(scaling,Scaling):
            raise TypeError("'scale' must be an instance of Scaling")
            
        if self.scaled == True:
            raise ValueError("Trying to scale an already scaled Polygon")
        self.x /= scaling.L
        self.z /= scaling.L
        self.rotation.xc /= scaling.L
        self.rotation.zc /= scaling.L
        self.scaled = True
        
    def unscale(self,scaling):
        if not isinstance(scaling,Scaling):
            raise TypeError("'scale' must be an instance of Scaling")
            
        if self.scaled == False:
            raise ValueError("Trying to unscale a non-scaled Polygon")
        self.x *= scaling.L
        self.z *= scaling.L
        self.rotation.xc *= scaling.L
        self.rotation.zc *= scaling.L
        self.scaled = False
        
    def assign_phase_to_particles(self,particles,atol=1e-6):
        """ ---------------------------------------------------------------------------
        Fast detection points inside a polygonal region.
        
        Originally written as a MATLAB mexFunction by:
        A. David Redish      email: adr@nsma.arizona.edu
        Guillaume Jacquenot  email: guillaume.jacquenot@gmail.com
        
        Output is simplified to skip distinguishing points IN and ON the polygon.
        Coordinates storage format is changed from strided to interleaved.
        Bounding box computation is separated from the main test function.
        #---------------------------------------------------------------------------
        """
        
        xPoints = particles.x
        zPoints = particles.z
        xVertices = self.x
        zVertices = self.z
        
        # get bounding box
        xmin = np.min(xVertices)
        xmax = np.max(xVertices)
        zmin = np.min(zVertices)
        zmax = np.max(zVertices)    
    
        nPoints = xPoints.size        
        nVertices = xVertices.size
        
        @maybe_numba(use_numba)
        def main_loop(particles_phase, polygon_phase):

        	# test whether each point is in polygon
            for ip in range(nPoints):
    
                # get point coordinates            
                zp = zPoints[ip]    
        		# check bounding box
                if(zp < zmin): continue
                if(zp > zmax): continue
                xp = xPoints[ip]
                if(xp < xmin): continue
                if(xp > xmax): continue
                
        
        		# count the number of intersections
                nIntersect = 0.0;
        
                for iv in range(nVertices):
        			# does the line PQ intersect the line AB?
                    if iv == nVertices-1:
                        ax = xVertices[(nVertices-1) ];
                        ay = zVertices[(nVertices-1)];
                        bx = xVertices[0         ];
                        by = zVertices[0         ];
                    else:
                        ax = xVertices[iv];
                        ay = zVertices[iv];
                        bx = xVertices[iv+1];
                        by = zVertices[iv+1];
                    # end if iV == nVertices-1:
        
                    if ax == bx:			
                        # vertical points
                        if xp == ax:
                            # ensure order correct
                            if ay > by:
                                tmp = ay; ay = by; by = tmp;
                                
                            if zp >= ay and zp <= by:
                                # point_on   = 1;
                                nIntersect = 0.0;
                                break;
                        # end if xp == ax:
                    else:
        				# non-vertical points
                        if xp < min(ax, bx) or max(ax, bx) < xp: continue;
                        
                        intersecty = ay + (xp - ax)/(bx - ax)*(by - ay);
        
                        if np.abs(intersecty - zp) < atol:
                            # point_on   = 1;
                            nIntersect = 0.0;
                            break;
                        elif intersecty < zp and (ax == xp or bx == xp):				
                            if ax == xp:
                                if iv == 0:
                                    ind = nVertices-1;
                                else:
                                    ind = iv-1;
        
                                xvind = xVertices[ind];
        
                                if min(bx, xvind) < xp and xp < max(bx, xvind):
                                    nIntersect += 1.0;
                            # end if    				
                        elif (intersecty < zp):
                            nIntersect += 1.0;
                        # end if
                    # end if ax == bx:
                # end for iV
                   
                # check if the contour polygon is closed
                point_in = np.int((nIntersect - 2.0*np.floor(nIntersect/2.0)));
                if point_in:
                    particles_phase[ip] = polygon_phase
        	# end for iP
        # end function main_loop
        main_loop(particles.phase, self.phase)
#     end assign_phase_to_particles

   
    

# =============================================================================
#
#                               Class Box
#
# =============================================================================     
class Box(Polygon):
    """
    Box(box=[x0,x1,z0,z1],phase=0)
    """
    def __init__(self,box,phase=0,scaled=False):
        if not isinstance(box,(list, tuple, np.ndarray)):
            raise TypeError("x must be a list, tuple or ndarray")
        
        
        self.x = np.array([box[0], box[1], box[1], box[0]])
        self.z = np.array([box[2], box[2], box[3], box[3]])
        self.phase = phase
        self.scaled = False

        self.rotation = Rotation()
    def assign_phase_to_particles(self,particles,atol=1e-6):
        """ ---------------------------------------------------------------------------
        Fast detection points inside a polygonal region.
        
        Originally written as a MATLAB mexFunction by:
        A. David Redish      email: adr@nsma.arizona.edu
        Guillaume Jacquenot  email: guillaume.jacquenot@gmail.com
        
        Output is simplified to skip distinguishing points IN and ON the polygon.
        Coordinates storage format is changed from strided to interleaved.
        Bounding box computation is separated from the main test function.
        #---------------------------------------------------------------------------
        """
        xPoints = particles.x.copy()
        zPoints = particles.z.copy()
        nPoints = xPoints.size 
        
        if self.rotation.angle==0.0:
            xmin = np.min(self.x)
            xmax = np.max(self.x)
            zmin = np.min(self.z)
            zmax = np.max(self.z)            
        else:
            # get bounding box
            rotation = self.rotation.copy()
            self.derotate()
            xmin = np.min(self.x)
            xmax = np.max(self.x)
            zmin = np.min(self.z)
            zmax = np.max(self.z)
            self.rotate(xc=rotation.xc,zc=rotation.zc,angle=rotation.angle)
            rotation.apply_to_arrays(xPoints,zPoints,reverse=True)
        # end if rotation==0.0
        
        @maybe_numba(use_numba)
        def main_loop(particles_phase, polygon_phase):
            for ip in range(nPoints):
                # get point coordinates            
                zp = zPoints[ip]
        		# check bounding box
                if(zp < zmin): continue
                if(zp > zmax): continue
                xp = xPoints[ip]
                if(xp < xmin): continue
                if(xp > xmax): continue
            
                particles_phase[ip] = polygon_phase    
#            return 
            
        main_loop(particles.phase, self.phase)
        
        
        

# =============================================================================
#
#                           Class BinaryNoise
#
# ============================================================================= 
class BinaryNoise():
    def __init__(self,proportion,N=500, phase_ids=(0,1),rL=1.0,h=1.0,clx=0.1,cly=0.1,angle=0.0):
        from scipy.optimize import minimize     
        self.proportion = proportion
        self.N = N
        self.rL = rL
        self.h = h
        self.clx = clx
        self.cly = cly
        self.angle = angle
        self.phase_ids = phase_ids
        
        def error_function(threshold,f,N,targetPhaseProportion):
            phaseProportion = np.sum(f<threshold)/N**2
            return (phaseProportion - targetPhaseProportion)**2
    

        f = self._create_random_surface(N=N,rL=rL,h=h,clx=clx,cly=cly,angle=angle)
           
        iniThreshold = 0.0
        idealThreshold = minimize(error_function, x0=iniThreshold, args=(f,self.N,self.proportion),method='Nelder-Mead').x[0]
        
        #print(errorFunction(idealThreshold.x[0],f,N,targetPhaseProportion))
        
        self.phase = phase_ids[0]*np.ones(f.shape,dtype=int)
        self.phase[f<idealThreshold] = phase_ids[1]
        phaseProportion = np.sum(f<idealThreshold)/N**2
        print("phaseProportion = %.2f %%" % (phaseProportion*100.0))
        
        
        
        


    def _create_random_surface(self,N,rL=1.0,h=1.0,clx=0.1,cly=0.1,angle=0.0):
        from numpy import exp, sqrt, pi
        from numpy.fft import fft2, ifft2
        # 2D Gaussian random rough surface with Gaussian autocovariance function
        x = np.linspace(-1.0,1.0,N)
        y = np.linspace(-1.0,1.0,N)
        
        X, Y =np.meshgrid(x,y)
        X = X.T
        Y = Y.T
    
        Z = np.random.randn(N,N)
        
        # if angle == 0.0:
        F = exp(-(X**2/(0.5*clx**2)+Y**2/(0.5*cly**2))); # Gaussian filter
        f = 2.0/sqrt(pi)*rL/N/sqrt(clx)/sqrt(cly)*ifft2(fft2(Z)*fft2(F)); # correlated surface generation including convolution (faltning) and inverse Fourier transform
        #   else:
        #       #tilted Gaussian
        #       a = ((cosd(angle)^2) / (2*clx^2)) + ((sind(angle)^2) / (2*cly^2));
        #       b = -((sind(2*angle)) / (4*clx^2)) + ((sind(2*angle)) / (4*cly^2));
        #       c = ((sind(angle)^2) / (2*clx^2)) + ((cosd(angle)^2) / (2*cly^2));
        #        
        #       F = exp(-(a*(X).^2 + 2*b*(X).*(Y) + c*(Y).^2)); % tilted Gaussian filter
        #       f = 2/sqrt(pi)*rL/N/sqrt(clx)/sqrt(cly)*ifft2(fft2(Z).*fft2(F)); % correlated surface generation including convolution (faltning) and inverse Fourier transform    
        import warnings
#            warnings.filterwarnings('ignore')
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Casting complex values to real discards the imaginary part")
            f = f.astype("float64")
            
        return f
        
    def assign_phase_to_particles(self,particles,atol=1e-6):
        
        def main_loop(particles_phase, noise_phase, particles_x, particles_z, xmin, xmax, zmin, zmax, N):

            for ip in range(particles_phase.size):
                x = (particles.x[ip] - xmin)/(xmax-xmin)
                z = (particles.z[ip] - zmin)/(zmax-zmin)
                
                # find the index of the closest node
                # /!\ maybe that's good or maybe cell center are better (important for periodicity)
                ix = int(np.floor(x*N))
                iz = int(np.floor(z*N))
                
                particles_phase[ip] = noise_phase[ix,iz]
                

        main_loop(particles.phase, self.phase, particles.x, particles.z, particles._xmin, particles._xmax, particles._zmin, particles._zmax, self.N)
                
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        