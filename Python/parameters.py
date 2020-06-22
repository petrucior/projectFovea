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
            self.colors = data['colors'] # colors = [ B1, G1, R1, B2, G2, R2, ..., Bk, Gk, Rk ] - Blue, Green and Red channels
            self.typeShape = data['typeShape'] # typeShape = Blocks (0) | Polygons (1)
            self.typeFovea = data['typeFovea'] # typeFovea = MRMF (0) | MMF (1)
            self.typeMultifovea = data['typeMultifovea'] # typeMultifovea = REEXECUTION (0) | PIXELBYPIXEL (1) | BITMAP (2) | SENDINGBLOCKS (3)

            # MMF in foveatedHessian
            self.bvector = data['bvector'] # bvector = [ b1, b2, ..., bk ]
            self.etavector = data['etavector'] # etavector = [ e1, e2, ..., ek ]
            self.levelvector = data['levelvector'] # levelvector = [ l1, l2, ..., lk ]
            self.growthfactor = data['growthfactor'] # growthfactor = x
            
            # Features
            self.features = data['features'] # features = [ ORB, KAZE, SURF, AKAZE, BRISK, SIFT ] ( binary vector )
            self.thresholdFeatures = data['thresholdFeatures']
            
            # ORB configuration
            self.orb_nfeatures = data['orb_nfeatures']
            self.orb_scaleFactor = data['orb_scaleFactor']
            self.orb_nlevels = data['orb_nlevels']
            self.orb_edgeThreshold = data['orb_edgeThreshold']
            self.orb_firstLevel = data['orb_firstLevel']
            self.orb_WTA_K = data['orb_WTA_K']
            self.orb_scoreType = data['orb_scoreType']
            self.orb_patchSize = data['orb_patchSize']
            self.orb_fastThreshold = data['orb_fastThreshold']

            # KAZE Configuration
            self.kaze_extended = data['kaze_extended']
            self.kaze_upright = data['kaze_upright']
            self.kaze_threshold = data['kaze_threshold']
            self.kaze_nOctaves = data['kaze_nOctaves']
            self.kaze_nOctaveLayers = data['kaze_nOctaves']
            self.kaze_diffusivity = data['kaze_diffusivity']

            # SURF Configuration
            self.surf_hessianThreshold = data['surf_hessianThreshold']
            self.surf_nOctaves = data['surf_nOctaves']
            self.surf_nOctaveLayers = data['surf_nOctaveLayers']
            self.surf_extended = data['surf_extended']
            self.surf_upright = data['surf_upright']
  
            # AKAZE Configuration
            self.akaze_descriptor_type = data['akaze_descriptor_type']
            self.akaze_descriptor_size = data['akaze_descriptor_size']
            self.akaze_descriptor_channels = data['akaze_descriptor_channels']
            self.akaze_threshold = data['akaze_threshold']
            self.akaze_nOctaves = data['akaze_nOctaves']
            self.akaze_nOctaveLayers = data['akaze_nOctaveLayers']
            self.akaze_diffusivity = data['akaze_diffusivity']

            # BRISK Configuration
            self.brisk_thresh = data['brisk_thresh']
            self.brisk_octaves = data['brisk_octaves']
            self.brisk_patternScale = data['brisk_patternScale']

            # SIFT Configuration
            self.sift_nfeatures = data['sift_nfeatures']
            self.sift_nOctaveLayers = data['sift_nOctaveLayers']
            self.sift_contrastThreshold = data['sift_contrastThreshold']
            self.sift_edgeThreshold = data['sift_edgeThreshold']
            self.sift_sigma = data['sift_sigma']
            
            
    def updateParameterSizeImage( self, u ):
        '''
        \fn updateParameterSizeImage( u )

        \brief Updating the U = ( ux, uy ) that corresponds the image size

        \param u - Size of image
        '''
        self.u = u
        # checking parameters
        self.checkParameters( self.m, self.w, self.u, self.f, self.bvector, self.etavector, self.levelvector, \
                              self.growthfactor, self.typeShape, self.typeFovea, self.typeMultifovea )

    def updateParameterFoveas( self, f, colors ):
        '''
        \fn updateParameterFoveas( f, colors )

        \brief Updating the foveas

        \param f - Image points where we will add foveas
        '''
        self.f = f
        self.colors = colors
        # checking parameters
        self.checkParameters( self.m, self.w, self.u, self.f, self.bvector, self.etavector, self.levelvector, \
                              self.growthfactor, self.typeShape, self.typeFovea, self.typeMultifovea )


    def fixFovea( self ):
        '''
        \fn fixFovea()

        \brief Fix the fovea position: if fovea is outsite image domain, snap it to the closest valid position independently for each coordinate
        '''
        self.f[0] = min((self.u[0] - self.w[0])/2, self.f[0])
        self.f[0] = max((self.w[0] - self.u[0])/2, self.f[0])
        self.f[1] = min((self.u[1] - self.w[1])/2, self.f[1])
        self.f[1] = max((self.w[1] - self.u[1])/2, self.f[1])
        
    
    def setFovea( self, px ):
        '''
        \fn setFovea( px )

        \brief Convert image points to axis fovea.

        \param px - Image points ( x, y )
        '''
        self.f[0] = px[0] - (self.u[0]/2)
        self.f[1] = px[1] - (self.u[1]/2)
        self.fixFovea()
        
        # checking parameters
        self.checkParameters( self.m, self.w, self.u, self.f, self.bvector, self.etavector, self.levelvector, \
                              self.growthfactor, self.typeShape, self.typeFovea, self.typeMultifovea )


    def updateParameters( self, m, w, u, f, bvector, etavector, levelvector, growthfactor, typeShape, typeFovea, typeMultifovea ):
        '''
        \fn updateParameters( m, w, u, f, bvector, etavector, levelvector, growthfactor, typeShape, typeFovea, typeMultifovea )

        \brief Updating all parameters to fovea structure
        
        \param m - Number levels of fovea
        \param w - Size of levels
        \param u - Size of image
        \param f - Position (x, y) to build the fovea
        \param bvector: [b1, b2, ..., bn] - a vector where bi is 0 if the feature extraction step number i should be discarded or 1, otherwise
        \param etavector: [e1, e2, ..., en] - a vector where ei is the octave (> 0) for which the feature extraction step number i should be performed
        \param levelvector: [l1, l2, ..., ln] - a vector where li is the foveated model level (>= 0 and < numberOfLevels) for which the feature extraction step number should be performed
        \param growthfactor - Growth factor associated the SURF feature 
        \param typeShape - Feature specification configured, where is selected blocks or polygons
        \param typeFovea - This code indicates which method to foveation will be used. If code is zero, then MRMF is chosen, otherwise MMF
        \param typeMultifovea - identify the approach: reexecution, pixel-by-pixel, bitmap or block-based
        '''
        # Checking parameters
        self.checkParameters( m, w, u, f, bvector, etavector, levelvector, growthfactor, typeShape, typeFovea, typeMultifovea )
        
        self.m = m
        self.w = w
        self.u = u
        self.f = f
        self.bvector = bvector
        self.etavector = etavector
        self.levelvector = levelvector
        self.growthfactor = growthfactor
        self.typeShape = typeShape
        self.typeFovea = typeFovea
        self.typeMultifovea = typeMultifovea


    def checkParameters( self, m, w, u, f, bvector, etavector, levelvector, growthfactor, typeShape, typeFovea, typeMultifovea ):
        '''
        \fn checkParameters( m, w, u, f, bvector, etavector, levelvector, growthfactor, typeShape, typeFovea, typeMultifovea )

        \brief This method check the parameters to build a fovea
        
        \param m - Number levels of fovea
        \param w - Size of levels
        \param u - Size of image
        \param f - Position (x, y) to build the fovea
        \param bvector: [b1, b2, ..., bn] - a vector where bi is 0 if the feature extraction step number i should be discarded or 1, otherwise
        \param etavector: [e1, e2, ..., en] - a vector where ei is the octave (> 0) for which the feature extraction step number i should be performed
        \param levelvector: [l1, l2, ..., ln] - a vector where li is the foveated model level (>= 0 and < numberOfLevels) for which the feature extraction step number should be performed
        \param growthfactor - Growth factor associated the SURF feature 
        \param typeShape - Feature specification configured, where is selected blocks or polygons
        \param typeFovea - This code indicates which method to foveation will be used. If code is zero, then MRMF is chosen, otherwise MMF
        \param typeMultifovea - identify the approach: reexecution, pixel-by-pixel, bitmap or block-based
        '''
        assert( m >= 1 )
        # Verify if w is bigger than zero and lower than image size
        assert( ( w[0] > 0 ) and ( w[0] < u[0] ) )
        assert( ( w[1] > 0 ) and ( w[1] < u[1] ) )
        # Verify if u is bigger than zero
        assert( ( u[0] > 0 ) and ( u[1] > 0 ) )
        # Verify if bvector, etavector and levelvector have size equal to m
        assert( ( len( bvector ) == m ) and ( len( etavector ) == m ) and ( len( levelvector ) == m ) )


#How to instantiate and use this class
#params = Parameters('params.yaml')
#print( params.w )

