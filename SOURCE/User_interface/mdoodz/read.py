#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 15:15:28 2019

@author: abauville
"""
import h5py
import numpy as np

def printFileTree(filename):
    with h5py.File(filename, 'r') as f:
        for key in f.keys():
            print("\n" + key ) #Names of the groups in HDF5 file.
            for key2 in f[key].keys():
                print("  " + key2 + " - " + str(f[key][key2].shape))    
                
class Grid():
    def __init__(self,filename):
        with h5py.File(filename, 'r') as f:
            params = f['Model/Params'][()]
            self.width  = params[1]
            self.height = params[2]
            self.nx     = int(params[3])
            self.nz     = int(params[4])
            self.ncx    = int(self.nx-1)
            self.ncz    = int(self.nz-1)
            self.dx     = params[5]
            self.dz     = params[6]
            
            energy = f['Model/Energy'][()]
            self.shortening = energy[1]
            self.Ut         = energy[2]
            self.Ue         = energy[3]
            self.W          = energy[4]
            
            
            self.xc_coord  = f['Model/xc_coord'][()]
            self.zc_coord  = f['Model/zc_coord'][()]
            self.xg_coord  = f['Model/xg_coord'][()]
            self.zg_coord  = f['Model/zg_coord'][()]
            self.xvz_coord = f['Model/xvz_coord'][()]
            self.zvx_coord = f['Model/zvx_coord'][()]
            
            self.xmin = np.min(self.xg_coord)
            self.xmax = np.max(self.xg_coord)
            self.zmin = np.min(self.zg_coord)
            self.zmax = np.max(self.zg_coord)
            
class Topo():
    def __init__(self,filename):
        with h5py.File(filename, 'r') as f:            
            self.height = f['Topo/height'][()]
            self.phase  = f['Topo/phase'][()]
            self.vx     = f['Topo/vx'][()]
            self.vz     = f['Topo/vz'][()]
            self.x      = f['Topo/x'][()]
            self.z      = f['Topo/z'][()]

class Time():
    def __init__(self,filename):
        with h5py.File(filename, 'r') as f:
            params = f['Model/Params'][()]
                    
            time = params[0]
            self.dt    = params[7]
            
            if time < 60:
                timeLabel = ' t = %.3f s ' % time
            elif time < 3600:
                timeLabel = ' t = %.3f mn ' % (time/60)
            elif time < 24*3600:
                timeLabel = ' t = %.3f h ' % (time/3600) + ' h '
            elif time < 365.25*24*3600:
                timeLabel = ' t = %.3f j ' % (time/24/3600) + ' j '
            elif time < 1e3*365.25*24*3600:
                timeLabel = ' t = %.3f a ' % (time/365/24/3600) + ' a '
            elif time < 1e6*365.25*24*3600:
                timeLabel = ' t = %.3f ka ' % (time/1e3/365.25/24/3600) + ' ka '
            else:
                timeLabel = ' t = %.3f Ma ' % (time/1e6/365.25/24/3600) + ' Ma '
            
            self.time  = time
            self.timeLabel = timeLabel
           
    

        
def get_center_data(filename,key):
        with h5py.File(filename, 'r') as f:
            params = f['Model/Params'][()]
            ncx   = int(params[3])-1
            ncz   = int(params[4])-1
            
            return np.reshape(f['Centers'][key][()],(ncz,ncx))
        
def get_vertex_data(filename,key):
        with h5py.File(filename, 'r') as f:
            params = f['Model/Params'][()]
            nx   = int(params[3])
            nz   = int(params[4])
            
            return np.reshape(f['Vertices'][key][()],(nz,nx))
            
def get_viz_grid_data(filename,key):
        with h5py.File(filename, 'r') as f:
            params = f['Model/Params'][()]
            ncx   = int(params[3])-1
            ncz   = int(params[4])-1
            if key[-3:] == '_hr':
                ncx *= 2
                ncz *= 2
            return np.reshape(f['VizGrid'][key][()],(ncz,ncx))
        
def get_vx_nodes_data(filename,key):
        with h5py.File(filename, 'r') as f:
            params = f['Model/Params'][()]
            nx   = int(params[3])
            nz   = int(params[4])+1
            return np.reshape(f['VxNodes'][key][()],(nz,nx))

def get_vz_nodes_data(filename,key):
        with h5py.File(filename, 'r') as f:
            params = f['Model/Params'][()]
            nx   = int(params[3])+1
            nz   = int(params[4])
            return np.reshape(f['VzNodes'][key][()],(nz,nx))
        
        
        
def get_data(filename,key,average="none"):
    # average = "none" (default), "center" or "vertex"
    # note: unfinished, only "vertex" is implemented for vx and vy
    found = False
    with h5py.File(filename, 'r') as f:
        for superKey in f.keys():
            print("\n" + superKey ) #Names of the groups in HDF5 file.
            for thisKey in f[superKey].keys():
                print("  " + thisKey + " - " + str(f[superKey][thisKey].shape))    
                if thisKey==key:
                    if superKey == 'Centers':
                        data = get_center_data(filename,key)
                    elif superKey == 'Vertices':
                        data = get_vertex_data(filename,key)
                    elif superKey == 'VizGrid':
                        data = get_viz_grid_data(filename,key)
                    elif superKey == 'VxNodes':
                        data = get_vx_nodes_data(filename,key)
                        if average.lower()=="vertex":
                            data = (data[:-1:,:] + data[1::,:])/2.0                            
                    elif superKey == 'VzNodes':
                        data = get_vz_nodes_data(filename,key)
                        if average.lower()=="vertex":
                            data = (data[:,:-1:] + data[:,1::])/2.0
                    else:
                        raise ValueError("Unknown superKey %s" % superKey)
                    found = True
                    break
                # end if thisKey=key
            # end for thisKey
            if found:
                break
        # end for superKey
    return data
            