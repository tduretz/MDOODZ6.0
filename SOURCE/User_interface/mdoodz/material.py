#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 11:53:41 2019

@author: abauville
"""

from mdoodz.model import Model
from copy import deepcopy

class Material():
    def __init__(self,
                 model,
         # Read general parameters
                 ID=0,
                 rho=2700.0,
                 mu=1.0e10,
                 Cv=1.0e3,
                 k=1.0e-6,
                 Qr=1.0e-30,
                 C=1.0e7,
                 phi=30.0,
                 Slim=1.0e10,
                 alp=0.0,
                 bet=1.0e-40,
                 drho=0.0,
                 
         # Read flow law settings
                 cstv=1.0,
                 pwlv=0.0,
                 linv=0.0,
                 gbsv=0.0,
                 expv=0.0,
                 gsel=0.0,
                 eta0=1.0e20,
                 npwl=1.0,
                 Qpwl=0.0,
                 pref_pwl=1.0,
                 gs=0.0,
                 gsref=2.0e-3,
                 
         # Strain softening
                 Ce=1.0e-7,
                 phie=30.0,
                 plss=1.0e6,
                 plse=1.0e6,
                 
         # Density models
                 density_model=1,
                 phase_diagram=-1
                 ):
        
        
        if not isinstance(model, Model):
            raise TypeError("'model' must be an instance of Model")
        
        model.Nb_phases += 1 
        model._materials_list.append(self)
# Read general parameters
        self.ID=ID
        self.rho=rho
        self.mu=mu
        self.Cv=Cv
        self.k=k
        self.Qr=Qr
        self.C=C
        self.phi=phi
        self.Slim=Slim
        self.alp=alp
        self.bet=bet
        self.drho=drho
         
 # Read flow law settings
        self.cstv=cstv
        self.pwlv=pwlv
        self.linv=linv
        self.gbsv=gbsv
        self.expv=expv
        self.gsel=gsel
        self.eta0=eta0
        self.npwl=npwl
        self.Qpwl=Qpwl
        self.pref_pwl=pref_pwl
        self.gs=gs
        self.gsref=gsref
         
 # Strain softening
        self.Ce=Ce
        self.phie=phie
        self.plss=plss
        self.plse=plse
        
 # Density models
        self.density_model=density_model
        self.phase_diagram=phase_diagram
        
    def copy_from(self,original_material):
        # Copy of everything except ID
        if not isinstance(original_material, self.__class__):
            raise ValueError("original_material should be an instance of mdoodz.material.Material")

        
        self.rho=original_material.rho
        self.mu=original_material.mu
        self.Cv=original_material.Cv
        self.k=original_material.k
        self.Qr=original_material.Qr
        self.C=original_material.C
        self.phi=original_material.phi
        self.Slim=original_material.Slim
        self.alp=original_material.alp
        self.bet=original_material.bet
        self.drho=original_material.drho
         
 # Read flow law settings
        self.cstv=original_material.cstv
        self.pwlv=original_material.pwlv
        self.linv=original_material.linv
        self.gbsv=original_material.gbsv
        self.expv=original_material.expv
        self.gsel=original_material.gsel
        self.eta0=original_material.eta0
        self.npwl=original_material.npwl
        self.Qpwl=original_material.Qpwl
        self.pref_pwl=original_material.pref_pwl
        self.gs=original_material.gs
        self.gsref=original_material.gsref
         
 # Strain softening
        self.Ce=original_material.Ce
        self.phie=original_material.phie
        self.plss=original_material.plss
        self.plse=original_material.plse
        
 # Density models
        self.density_model=original_material.density_model
        self.phase_diagram=original_material.phase_diagram
        
        