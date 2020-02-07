#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
\file parameters.py

\brief This file contains all information about fovea structure

\author
Petrucio Ricardo Tavares de Medeiros \n
Universidade Federal do Rio Grande do Norte \n 
Departamento de Computacao e Automacao Industrial \n
petrucior at ufrn (dot) edu (dot) br

\version 0.1
\date January 2020

\note
This code demands the instalation of PyYAML library:
pip install pyyaml

This file is part of projectFovea software.
This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 2 of the License, or (at your option) any later
version. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import yaml

class Parameters:
    '''
    \class Parameters
    
    \brief This class implements the Parameters TAD of fovea
    '''
    
    def __init__( self, file ):
        '''
        \fn Parameters()

        \brief Constructor default

        \param file - File that contains all information about the fovea
        '''
        with open( file, 'r' ) as stream :
            data = yaml.full_load( stream )

            # Loading parameters to attributes ( For more details visit the projectFovea page )
            self.m = data['numberOfLevels'] # numberOfLevels = m
            self.w = data['smallestLevel'] # smallestLevel = [ wx, wy ]
            self.f = data['foveas'] # foveas = [ fx1, fy1, fx2, fy2, ..., fxk, fyk ]
            self.bvector = data['bvector'] # bvector = [ b1, b2, ..., bk ]
            self.etavector = data['etavector'] # etavector = [ e1, e2, ..., ek ]
            self.levelvector = data['levelvector'] # levelvector = [ l1, l2, ..., lk ]
            self.nOctaveLayers = data['nOctaveLayers'] # nOctaveLayers = x
            self.hessianThreshold = data['hessianThreshold'] # hessianThreshold = x
            self.growthfactor = data['growthfactor'] # growthfactor = x
            self.typeShape = data['typeShape'] # typeShape = Blocks (0) | Polygons (1)
            self.typeFovea = data['typeFovea'] # typeFovea = MRMF (0) | MMF (1)
            self.typeMultifovea = data['typeMultifovea'] # typeMultifovea = REEXECUTION (0) | PIXELBYPIXEL (1) | BITMAP (2) | SENDINGBLOCKS (3)

    def updateParameterSizeImage( self, u ):
        '''
        \fn updateParameterSizeImage( u )

        \brief Updating the U = ( ux, uy ) that corresponds the image size

        \param u - Size of image
        '''
        self.u = u

    def updateParameters( self, m, w, u, f, bvector, etavector, levelvector, nOctaveLayers, hessianThreshold, growthfactor, typeShape, typeFovea, typeMultifovea ):
        '''
        \fn updateParameters( m, w, u, f, bvector, etavector, levelvector, nOctaveLayers, hessianThreshold, growthfactor, typeShape, typeFovea, typeMultifovea )

        \brief Updating all parameters to fovea structure
        
        \param m - Number levels of fovea
        \param w - Size of levels
        \param u - Size of image
        \param f - Position (x, y) to build the fovea
        \param bvector: [b1, b2, ..., bn] - a vector where bi is 0 if the feature extraction step number i should be discarded or 1, otherwise
        \param etavector: [e1, e2, ..., en] - a vector where ei is the octave (> 0) for which the feature extraction step number i should be performed
        \param levelvector: [l1, l2, ..., ln] - a vector where li is the foveated model level (>= 0 and < numberOfLevels) for which the feature extraction step number should be performed
        \param nOctaveLayers - Number of octaves associated the SURF feature
        \param hessianThreshold - Threshold associated the SURF feature
        \param growthfactor - Growth factor associated the SURF feature 
        \param typeShape - Feature specification configured, where is selected blocks or polygons
        \param typeFovea - This code indicates which method to foveation will be used. If code is zero, then MRMF is chosen, otherwise MMF
        \param typeMultifovea - identify the approach: reexecution, pixel-by-pixel, bitmap or block-based
        '''
        self.m = m
        self.w = w
        self.u = u
        self.f = f
        self.bvector = bvector
        self.etavector = etavector
        self.levelvector = levelvector
        self.nOctaveLayers = nOctaveLayers
        self.hessianThreshold = hessianThreshold
        self.growthfactor = growthfactor
        self.typeShape = typeShape
        self.typeFovea = typeFovea
        self.typeMultifovea = typeMultifovea


#How to instantiate and use this class
#params = Parameters('params.yaml')
#print( params.w )
