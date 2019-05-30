#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 11:40:39 2019

@author: abauville
"""

""" ---------------------------------------------------------------------------
Fast detection points inside a polygonal region.

Originally written as a MATLAB mexFunction by:
A. David Redish      email: adr@nsma.arizona.edu
Guillaume Jacquenot  email: guillaume.jacquenot@gmail.com

Modified to be callable directly from C by:
Anton A. Popov       email: popov@uni-mainz.de

Output is simplified to skip distinguishing points IN and ON the polygon.
Coordinates storage format is changed from strided to interleaved.
Bounding box computation is separated from the main test function.

This function has been copied from the open source code of LaMEM - Lithosphere 
and Mantle Evolution Model, commit number: 1dfe60b, Dec. 12, 2018
ref: Kaus, Boris JP, et al. "Forward and inverse modelling of lithospheric 
deformation on geological timescales." NIC Symposium 2016-Proceedings, NIC Series, Jülich, DE. 2016.
https:#bitbucket.org/bkaus/lamem/
#---------------------------------------------------------------------------
"""

import numpy as np
import InputClassDef as Input

try:
    from numba import jit #, prange
    useNumba = True
except:
    print("Warning: the package numba was not found. just-in-time compilation is off")
    useNumba = False
    
def maybe_numba(useNumba):
    return jit(nopython=True) if useNumba else lambda x:x



#import math
class Rotation():
    def __init__(self,xc=0.0,zc=0.0,angle=0.0):
        self.xc = xc
        self.zc = zc
        self.angle = angle
        
    def applyToArrays(self,x,z,reverse=False):
        xtemp = x-self.xc
        ztemp = z-self.zc
        
        if reverse:
            angle = -self.angle
        else:
            angle = self.angle                
        
        x[:] = np.cos(angle)*xtemp - np.sin(angle)*ztemp + self.xc
        z[:] = np.sin(angle)*xtemp + np.cos(angle)*ztemp + self.zc
        
    def applyToPolygon(self,polygon,reverse=False):
        if not isinstance(polygon,Polygon):
            raise ValueError("'Polygon' must be an instance of 'Geometry.Polygon'")
            
        self.applyToArrays(polygon.x,polygon.z,reverse)

    def copy(self):
        return Rotation(xc=self.xc,zc=self.zc,angle=self.angle)




class Polygon():
    def __init__(self,x,z,phase=0,scaled=False):
        if type(x)==list or type(x)==tuple:
            self.x = np.array(x)
        elif type(x)==np.ndarray:
            self.x = x
        else:
            raise TypeError("x must be a list, tuple or ndarray")
            
        if type(z)==list or type(z)==tuple:
            self.z = np.array(z)
        elif type(z)==np.ndarray:
            self.z = z
        else:
            raise TypeError("x must be a list, tuple or ndarray")


        self.phase = phase
        self.scaled = False

        self.rotation = Rotation()
        
    def rotate(self,angle,xc=0.0,zc=0.0,reverse=False):
        """
            This function rotates the polygon by a angle "angle" (in radians) around
            the point of coordinates (xc,zc)
            
        """
        self.rotation = Rotation(xc=xc,zc=zc,angle=angle)
        self.rotation.applyToPolygon(self,reverse)
        
    def derotate(self):
        """
            This function rotates the polygon by a angle -"angle" (in radians) around
            the point of coordinates (xc,zc)
            
        """
        self.rotation.applyToPolygon(self,reverse=True)
        self.rotation = Rotation()
        

    def scale(self,scaling):
        if type(scaling)!=Input.Scaling:
            raise TypeError("'scale' must be an instance of Input.Scaling")
            
        if self.scaled == True:
            raise ValueError("Trying to scale an already scaled Polygon")
        self.x /= scaling.L
        self.z /= scaling.L
        self.rotation.xc /= scaling.L
        self.rotation.zc /= scaling.L
        self.scaled = True
        
    def unscale(self,scaling):
        if type(scaling)!=Input.Scaling:
            raise TypeError("'scale' must be an instance of Input.Scaling")
            
        if self.scaled == False:
            raise ValueError("Trying to unscale a non-scaled Polygon")
        self.x *= scaling.L
        self.z *= scaling.L
        self.rotation.xc *= scaling.L
        self.rotation.zc *= scaling.L
        self.scaled = False
        
    def assignPhaseToParticles(self,particles,atol=1e-6):
        
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
        
        @maybe_numba(useNumba)
        def mainLoop(particles_phase, polygon_phase):

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
        # end function mainLoop
        mainLoop(particles.phase, self.phase)
#     end assignPhaseToParticles

        
class Box(Polygon):
    """
    Box(box=[x0,x1,z0,z1],phase=0)
    """
    def __init__(self,box,phase=0,scaled=False):
        if type(box) not in (list, tuple, np.ndarray):
            raise TypeError("x must be a list, tuple or ndarray")
        
        
        self.x = np.array([box[0], box[1], box[1], box[0]])
        self.z = np.array([box[2], box[2], box[3], box[3]])
        self.phase = phase
        self.scaled = False

        self.rotation = Rotation()
    def assignPhaseToParticles(self,particles,atol=1e-6):
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
            rotation.applyToArrays(xPoints,zPoints,reverse=True)
        # end if rotation==0.0
        
        @maybe_numba(useNumba)
        def mainLoop(particles_phase, polygon_phase):
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
            
        mainLoop(particles.phase, self.phase)