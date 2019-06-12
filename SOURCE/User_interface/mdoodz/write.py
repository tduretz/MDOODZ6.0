#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 14:34:37 2019

@author: abauville
"""
import numpy as np
import json
from mdoodz.model import Model

class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return json.JSONEncoder.default(self, obj)


# File writing functions 
# =============================================================================
def JsonFile(topo_chain,
                  filename="input.json"):
    
    myJsonFile = dict(TopoChain = vars(topo_chain))   
    with open(filename, "w") as write_file:
        json.dump(myJsonFile, write_file , indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False, cls=NumpyEncoder)





# Write files
# ==========================================
def text_file(model, scaling, particles, materials_list, filename="input.txt"):
    
    particle_dict = {"Nx_part":particles.Nx_part,
                     "Nz_part":particles.Nz_part,
                     "min_part_cell":particles.min_part_cell}
    model_dict = vars(model)
    scaling_dict = vars(scaling)
    
    mergeDict = {}
    for key in particle_dict:
        if key[0] != '_':
            mergeDict[key] = particle_dict[key]
    for key in model_dict:
        if key[0] != '_':
            mergeDict[key] = model_dict[key]
    for key in scaling_dict:
        if key[0] != '_':
            mergeDict[key] = scaling_dict[key]
        
        
    
    
#    with open(filename, "w") as write_file:
#        json.dump(myJsonFile, write_file , indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False, cls=NumpyEncoder)
        
    myFile = json.dumps(mergeDict, indent=4, separators=('', ' = '), ensure_ascii=False, cls=NumpyEncoder)
    
    myFile += "\n\n/**** MATERIALS ****/"
    
    if len(materials_list) != model.Nb_phases:
        raise ValueError("model.Nb_phases=%i but len(materials_list)=%i. Be sure to pass all defined materials to write the file." % (model.Nb_phases, len(materials_list)))
#    else:
#        myFile += "\nNb_phases = %i\n" % model.Nb_phases
        
    for iMat in range(len(materials_list)):
        myFile += "\n\n"
        myFile += json.dumps(vars(materials_list[iMat]), indent=4, separators=('', ' = '), ensure_ascii=False, cls=NumpyEncoder)
        
    
    myFile = myFile.replace("{\n","")
    myFile = myFile.replace("\n}","")
    myFile = myFile.replace("    ","")
    myFile = myFile.replace('"','')
#    
    
#    print(myFile)
    
    with open(filename, "w") as write_file:
        write_file.write(myFile)

    
      
def ini_particles_file(model, particles, topo_chain=None, filename = "IniParticles.dat"):

    if model.isScaled:
        raise ValueError("the model should not be scaled to write the input file. (This error raising should disappear when development is mature)")
        
    s1 = 4 # sizeof(int)
    s2 = 8 # sizeof(double)
    s3 = 8 # sizeof(DoodzFP)
    s4 = 1 # sizeof(char)
    
    # Global variables
    s1Type = np.int32
    s2Type = np.float64
    s3Type = np.float64
    s4Type = np.str
    
#    print("Writing file...")
    
    with open(filename, "w") as file:
        myArray = np.array([s1,s2,s3,s4],dtype=np.int8)
        myArray.tofile(file)
        #fwrite( &s1, 1, 1, file);
        #fwrite( &s2, 1, 1, file);
        #fwrite( &s3, 1, 1, file);
        #fwrite( &s4, 1, 1, file);    
        
        if model.free_surf == 1:
            myArray = np.array([topo_chain.Nb_part],dtype=s1Type)   
            myArray.tofile(file)
            #    fwrite( &topo_chain->Nb_part,  s1,                   1, file );
                
            topo_chain.x.astype(s3Type).tofile(file)
            topo_chain.z.astype(s3Type).tofile(file)
#            topo_chain.Vx.astype(s3Type).tofile(file)
#            topo_chain.Vz.astype(s3Type).tofile(file)
            #    fwrite( topo_chain->x,         s3, topo_chain->Nb_part, file );
            #    fwrite( topo_chain->z,         s3, topo_chain->Nb_part, file );
            #    fwrite( topo_chain->Vx,        s3, topo_chain->Nb_part, file );
            #    fwrite( topo_chain->Vz,        s3, topo_chain->Nb_part, file );
    
    
    
    
        myArray = np.array([particles.Nb_part],dtype=s1Type)   
        myArray.tofile(file)
        #fwrite( &particles->Nb_part,  s1, 1, file);
        particles.x.astype(s3Type).tofile(file)
        particles.z.astype(s3Type).tofile(file)
#        particles.P.astype(s3Type).tofile(file)
#        particles.Vx.astype(s3Type).tofile(file)
#        particles.Vz.astype(s3Type).tofile(file)
#        particles.phi.astype(s3Type).tofile(file)
#        particles.X.astype(s3Type).tofile(file)
        particles.phase.astype(s1Type).tofile(file)
        #fwrite( particles->x,     s3, particles->Nb_part, file);
        #fwrite( particles->z,     s3, particles->Nb_part, file);
        #fwrite( particles->P,     s3, particles->Nb_part, file);
        #fwrite( particles->Vx,    s3, particles->Nb_part, file);
        #fwrite( particles->Vz,    s3, particles->Nb_part, file);
        #fwrite( particles->phi,   s3, particles->Nb_part, file);
        #fwrite( particles->X  ,   s3, particles->Nb_part, file);
        #fwrite( particles->phase, s1, particles->Nb_part, file);
        #
    
    
#    print("file writing took %.2f s\n" % (time.time()-tic))
        
