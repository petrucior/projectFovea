#!/usr/bin/pythonOA
# -*- coding: utf-8 -*-

'''
\file features.py

\brief This file contains the prototype of features.

\author
Petrucio Ricardo Tavares de Medeiros \n
Universidade Federal do Rio Grande do Norte \n 
Departamento de Computacao e Automacao Industrial \n
petrucior at ufrn (dot) edu (dot) br

\version 0.1
\date February 2020
 
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
import cv2
from parameters import Parameters
from level import Level
from matplotlib import pyplot
import matplotlib.pyplot as plt
import numpy as np

class Features :
    '''
    \class Features
    
    \brief This class implements the Features TAD

    \note https://docs.opencv.org/master/d0/d13/classcv_1_1Feature2D.html
          https://python.docow.com/8279/plotar-imagens-lado-a-lado-usando-o-matplotlib.html
    '''
    
    def __init__( self, img, levels, parameters ):
        '''
        \fn Features( levels, parameters )
        
        \brief Constructor default
        
        \param img - image to be foveated
        \param levels - set of levels to build the fovea
        \param parameters - Parameters of fovea structure
        '''
        assert( len(levels) > 0 ) # Feature extraction only from a fovea structure
        p = parameters
        # ORB feature configured
        self.orb = cv2.ORB_create( p.orb_nfeatures, p.orb_scaleFactor, p.orb_nlevels, p.orb_edgeThreshold, p.orb_firstLevel,
                              p.orb_WTA_K, cv2.ORB_FAST_SCORE, p.orb_patchSize, p.orb_fastThreshold )
                                      
        # KAZE feature configured
        self.kaze = cv2.KAZE_create( p.kaze_extended, p.kaze_upright, p.kaze_threshold, p.kaze_nOctaves, p.kaze_nOctaveLayers,
                                cv2.KAZE_DIFF_PM_G2 )

        # SURF feature configured
        self.surf = cv2.xfeatures2d.SURF_create( p.surf_hessianThreshold, p.surf_nOctaves, p.surf_nOctaveLayers, p.surf_extended, p.surf_upright )

        # AKAZE feature configured
        self.akaze = cv2.AKAZE_create( cv2.AKAZE_DESCRIPTOR_MLDB, p.akaze_descriptor_size, p.akaze_descriptor_channels,
                                  p.akaze_threshold, p.akaze_nOctaves, p.akaze_nOctaveLayers, cv2.KAZE_DIFF_PM_G2 )

        # BRISK feature configured
        self.brisk = cv2.BRISK_create( p.brisk_thresh, p.brisk_octaves, p.brisk_patternScale )

        # SIFT feature configured
        self.sift = cv2.xfeatures2d.SIFT_create( p.sift_nfeatures, p.sift_nOctaveLayers, p.sift_contrastThreshold, p.sift_edgeThreshold, p.sift_sigma )

        self.kpORB = []; self.kpKAZE = []; self.kpSURF = []; self.kpAKAZE = []; self.kpBRISK = []; self.kpSIFT = []; # cleaning keypoints
        self.desORB = []; self.desKAZE = []; self.desSURF = []; self.desAKAZE = []; self.desBRISK = []; self.desSIFT = []; # cleaning descriptors
        
        # features = [ ORB, KAZE, SURF, AKAZE, BRISK, SIFT ]
        levelsKp = []
        for k in range( 0, p.m + 1 ):
            imgLevel = levels[k].getLevel( img, p )
            levelsKp.append( imgLevel )
            if ( p.features[0] == 1 ): # ORB feature activated
                #print('ORB feature activated')
                kpAux, desAux = self.orb.detectAndCompute( imgLevel, None )
                self.kpORB.append( kpAux )
                self.desORB.append( desAux )
                imgLevelORB = cv2.drawKeypoints( imgLevel, self.kpORB[k], None, (255,0,0), 4 )
                levelsKp.append( imgLevelORB )
            if ( p.features[1] == 1 ): # KAZE feature activated
                #print('KAZE feature activated')
                kpAux, desAux = self.kaze.detectAndCompute( imgLevel, None )
                self.kpKAZE.append( kpAux )
                self.desKAZE.append( desAux )
                imgLevelKAZE = cv2.drawKeypoints( imgLevel, self.kpKAZE[k], None, (0,255,0), 4 )
                levelsKp.append( imgLevelKAZE )                
            if ( p.features[2] == 1 ): # SURF feature activated
                #print('SURF feature activated')
                kpAux, desAux = self.surf.detectAndCompute( imgLevel, None )
                self.kpSURF.append( kpAux )
                self.desSURF.append( desAux )
                imgLevelSURF = cv2.drawKeypoints( imgLevel, self.kpSURF[k], None, (0,0,255), 4 )
                levelsKp.append( imgLevelSURF )
            if ( p.features[3] == 1 ): # AKAZE feature activated
                #print('AKAZE feature activated')
                kpAux, desAux = self.akaze.detectAndCompute( imgLevel, None )
                self.kpAKAZE.append( kpAux )
                self.desAKAZE.append( desAux )
                imgLevelAKAZE = cv2.drawKeypoints( imgLevel, self.kpAKAZE[k], None, (255,255,0), 4 )
                levelsKp.append( imgLevelAKAZE )
            if ( p.features[4] == 1 ): # BRISK feature activated
                #print('BRISK feature activated')
                kpAux, desAux = self.brisk.detectAndCompute( imgLevel, None )
                self.kpBRISK.append( kpAux )
                self.desBRISK.append( desAux )
                imgLevelBRISK = cv2.drawKeypoints( imgLevel, self.kpBRISK[k], None, (255,0,255), 4 )
                levelsKp.append( imgLevelBRISK )
            if ( p.features[5] == 1 ): # SIFT feature activated
                #print('SIFT feature activated')
                kpAux, desAux = self.sift.detectAndCompute( imgLevel, None )
                self.kpSIFT.append( kpAux )
                self.desSIFT.append( desAux )
                imgLevelSIFT = cv2.drawKeypoints( imgLevel, self.kpSIFT[k], None, (0,255,255), 4 )
                levelsKp.append( imgLevelSIFT )

        # Figure shows distribution of the features ORB, KAZE, SURF, AKAZE, BRISK and SIFT by levels
        '''
        fig, axs = plt.subplots(p.m+1, int(len(levelsKp)/((p.m)+1)) )#, figsize=(3, 3))
        fig.suptitle('Distribuição das features ORB, KAZE, SURF, AKAZE, BRISK e SIFT')
        axs = axs.flatten()
        for img, ax in zip(levelsKp, axs):
            ax.imshow(img)
            ax.set(xlabel='', ylabel='')
            ax.label_outer()
        plt.show()
        '''

    def matchingFeatures( self, model, parameters ):
        '''
        \fn matchingFeatures( model, parameters )
        
        \brief Maching between multiresolution domain and model
        
        \param model - Template that represent the visual stimulli
        \param parameters - Parameters of fovea structure
        '''
        # Matching happens only if features were extracted
        assert( len(self.kpORB) > 0 or len(self.kpKAZE) > 0 or len(self.kpSURF) > 0 or len(self.kpAKAZE) > 0 or len(self.kpBRISK) > 0 or len(self.kpSIFT) > 0)
        p = parameters
        # Model keypoints and descriptors 
        #self.kpModelORB = []; self.kpModelKAZE = []; self.kpModelSURF = []; self.kpModelAKAZE = []; self.kpModelBRISK = []; self.kpModelSIFT = []; # cleaning keypoints
        #self.desModelORB = []; self.desModelKAZE = []; self.desModelSURF = []; self.desModelAKAZE = []; self.desModelBRISK = []; self.desModelSIFT = []; # cleaning descriptors

        # Matches vectors
        self.matchORB = []; self.matchKAZE = []; self.matchSURF= []; self.matchAKAZE = []; self.matchBRISK = []; self.matchSIFT = []; # cleaning matches

        # Inliers rate
        self.inlierRateORB = []; self.inlierRateKAZE = []; self.inlierRateSURF = []; self.inlierRateAKAZE = []; self.inlierRateBRISK = []; self.inlierRateSIFT = [];

        # features = [ ORB, KAZE, SURF, AKAZE, BRISK, SIFT ]
        if ( p.features[0] == 1 ): # ORB feature activated
            #print('ORB feature activated')
            self.kpModelORB, self.desModelORB = self.orb.detectAndCompute( model, None )
            
        if ( p.features[1] == 1 ): # KAZE feature activated
            #print('KAZE feature activated')
            self.kpModelKAZE, self.desModelKAZE = self.kaze.detectAndCompute( model, None )
                        
        if ( p.features[2] == 1 ): # SURF feature activated
            #print('SURF feature activated')
            self.kpModelSURF, self.desModelSURF = self.surf.detectAndCompute( model, None )
                        
        if ( p.features[3] == 1 ): # AKAZE feature activated
            #print('AKAZE feature activated')
            self.kpModelAKAZE, self.desModelAKAZE = self.akaze.detectAndCompute( model, None )
            
        if ( p.features[4] == 1 ): # BRISK feature activated
            #print('BRISK feature activated')
            self.kpModelBRISK, self.desModelBRISK = self.brisk.detectAndCompute( model, None )
            
        if ( p.features[5] == 1 ): # SIFT feature activated
            #print('SIFT feature activated')
            self.kpModelSIFT, self.desModelSIFT = self.sift.detectAndCompute( model, None )
        
        for k in range( 0, p.m + 1 ):
            if ( p.features[0] == 1 ): # ORB feature activated
                if ( len(self.kpModelORB) > 0 and len(self.kpORB[k]) > 0 ):
                    # create BFMatcher object
                    bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=True)
                    # Match descriptors
                    matches = bf.match(self.desModelORB, self.desORB[k])
                    # Sort them in the order of their distance
                    matches = sorted(matches, key = lambda x:x.distance)
                    matches = matches[:50]
                    # Saving matches
                    self.matchORB.append( matches )

                    # Localize the object
                    src_pts = np.float32([self.kpModelORB[m.queryIdx].pt for m in matches]).reshape(-1, 1, 2)
                    dst_pts = np.float32([self.kpORB[k][m.trainIdx].pt for m in matches]).reshape(-1, 1, 2)
                    # Find Homography and do perspective transform
                    H, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC, 5.0)
                    inliers = 0
                    for i in range(0, len(mask)):
                        if ( mask[i] == 1 ):
                            inliers += 1
                    #if ( inliers > parameters.thresholdFeatures ):
                    if ( len(mask) > parameters.thresholdFeatures ):
                        self.inlierRateORB.append( float(inliers/len(mask)) )
                    else:
                        self.inlierRateORB.append( 0.0 )
                else:
                    self.matchORB.append( 0 )
                    self.inlierRateORB.append( 0 )
            
            if ( p.features[1] == 1 ): # KAZE feature activated
                if ( len(self.kpModelKAZE) > 0 and len(self.kpKAZE[k]) > 0 ):                
                    # create BFMatcher object
                    bf = cv2.BFMatcher(cv2.NORM_L2, crossCheck=True)
                    # Match descriptors
                    matches = bf.match(self.desModelKAZE, self.desKAZE[k])
                    # Sort them in the order of their distance
                    matches = sorted(matches, key = lambda x:x.distance)
                    matches = matches[:50]
                    # Saving matches
                    self.matchKAZE.append( matches )
                
                    # Localize the object 
                    src_pts = np.float32([self.kpModelKAZE[m.queryIdx].pt for m in matches]).reshape(-1, 1, 2)
                    dst_pts = np.float32([self.kpKAZE[k][m.trainIdx].pt for m in matches]).reshape(-1, 1, 2)
                    # Find Homography and do perspective transform
                    H, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC, 5.0)
                    inliers = 0
                    for i in range(0, len(mask)):
                        if ( mask[i] == 1 ):
                            inliers += 1
                    #if ( inliers > parameters.thresholdFeatures ):
                    if ( len(mask) > parameters.thresholdFeatures ):
                        self.inlierRateKAZE.append( float(inliers/len(mask)) )
                    else:
                        self.inlierRateKAZE.append( 0.0 )
                else:
                    self.matchKAZE.append( 0 )
                    self.inlierRateKAZE.append( 0 )
            
            if ( p.features[2] == 1 ): # SURF feature activated
                if ( len(self.kpModelSURF) > 0 and len(self.kpSURF[k]) > 0 ):
                    # create BFMatcher object
                    bf = cv2.BFMatcher(cv2.NORM_L2, crossCheck=True)
                    # Match descriptors
                    matches = bf.match(self.desModelSURF, self.desSURF[k])
                    # Sort them in the order of their distance
                    matches = sorted(matches, key = lambda x:x.distance)
                    #matches = matches[:50]
                    # Saving matches
                    self.matchSURF.append( matches )
                
                    # Localize the object 
                    src_pts = np.float32([self.kpModelSURF[m.queryIdx].pt for m in matches]).reshape(-1, 1, 2)
                    dst_pts = np.float32([self.kpSURF[k][m.trainIdx].pt for m in matches]).reshape(-1, 1, 2)
                    # Find Homography and do perspective transform
                    H, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC, 5.0)
                    inliers = 0
                    for i in range(0, len(mask)):
                        if ( mask[i] == 1 ):
                            inliers += 1
                    #print( "matches: "+ str(len(matches))+", inliers: "+ str(inliers)+", ratio: "+str(float(inliers/len(mask))) )
                    #if ( inliers > parameters.thresholdFeatures ):
                    if ( len(mask) > parameters.thresholdFeatures ):
                        self.inlierRateSURF.append( float(inliers/len(mask)) )
                    else:
                        self.inlierRateSURF.append( 0.0 )
                else:
                    self.matchSURF.append( 0 )
                    self.inlierRateSURF.append( 0 )                    
            
            if ( p.features[3] == 1 ): # AKAZE feature activated
                if ( len(self.kpModelAKAZE) > 0 and len(self.kpAKAZE[k]) > 0 ):                
                    # create BFMatcher object
                    bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=True)
                    # Match descriptors
                    matches = bf.match(self.desModelAKAZE, self.desAKAZE[k])
                    # Sort them in the order of their distance
                    matches = sorted(matches, key = lambda x:x.distance)
                    matches = matches[:50]
                    # Saving matches
                    self.matchAKAZE.append( matches )
                    
                    # Localize the object 
                    src_pts = np.float32([self.kpModelAKAZE[m.queryIdx].pt for m in matches]).reshape(-1, 1, 2)
                    dst_pts = np.float32([self.kpAKAZE[k][m.trainIdx].pt for m in matches]).reshape(-1, 1, 2)
                    # Find Homography and do perspective transform
                    H, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC, 5.0)
                    inliers = 0
                    for i in range(0, len(mask)):
                        if ( mask[i] == 1 ):
                            inliers += 1
                    #if ( inliers > parameters.thresholdFeatures ):
                    if ( len(mask) > parameters.thresholdFeatures ):
                        self.inlierRateAKAZE.append( float(inliers/len(mask)) )
                    else:
                        self.inlierRateAKAZE.append( 0.0 )
                else:
                    self.matchAKAZE.append( 0 )
                    self.inlierRateAKAZE.append( 0 )
                    
            if ( p.features[4] == 1 ): # BRISK feature activated
                if ( len(self.kpModelBRISK) > 0 and len(self.kpBRISK[k]) > 0 ):
                    # create BFMatcher object
                    bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=True)
                    # Match descriptors
                    matches = bf.match(self.desModelBRISK, self.desBRISK[k])
                    # Sort them in the order of their distance
                    matches = sorted(matches, key = lambda x:x.distance)
                    matches = matches[:50]
                    # Saving matches
                    self.matchBRISK.append( matches )
                    
                    # Localize the object 
                    src_pts = np.float32([self.kpModelBRISK[m.queryIdx].pt for m in matches]).reshape(-1, 1, 2)
                    dst_pts = np.float32([self.kpBRISK[k][m.trainIdx].pt for m in matches]).reshape(-1, 1, 2)
                    # Find Homography and do perspective transform
                    H, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC, 5.0)
                    
                    inliers = 0;
                    for i in range(0, len(mask)):
                        if ( mask[i] == 1 ):
                            inliers += 1
                    #if ( inliers > parameters.thresholdFeatures ):
                    if ( len(mask) > parameters.thresholdFeatures ):
                        self.inlierRateBRISK.append( float(inliers/len(mask)) )
                    else:
                        self.inlierRateBRISK.append( 0.0 )
                else:
                    self.matchBRISK.append( 0 )
                    self.inlierRateBRISK.append( 0 )
                    
            if ( p.features[5] == 1 ): # SIFT feature activated
                if ( len(self.kpModelSIFT) > 0 and len(self.kpSIFT[k]) > 0 ):                
                    # create BFMatcher object
                    bf = cv2.BFMatcher(cv2.NORM_L2, crossCheck=True)
                    # Match descriptors
                    matches = bf.match(self.desModelSIFT, self.desSIFT[k])
                    # Sort them in the order of their distance
                    matches = sorted(matches, key = lambda x:x.distance)
                    matches = matches[:50]
                    # Saving matches
                    self.matchSIFT.append( matches )
                    
                    # Localize the object 
                    src_pts = np.float32([self.kpModelSIFT[m.queryIdx].pt for m in matches]).reshape(-1, 1, 2)
                    dst_pts = np.float32([self.kpSIFT[k][m.trainIdx].pt for m in matches]).reshape(-1, 1, 2)
                    # Find Homography and do perspective transform
                    H, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC, 5.0)
                    inliers = 0
                    for i in range(0, len(mask)):
                        if ( mask[i] == 1 ):
                            inliers += 1
                    #if ( inliers > parameters.thresholdFeatures ):
                    if ( len(mask) > parameters.thresholdFeatures ):
                        self.inlierRateSIFT.append( float(inliers/len(mask)) )
                    else:
                        self.inlierRateSIFT.append( 0.0 )
                else:
                    self.matchSIFT.append( 0 )
                    self.inlierRateSIFT.append( 0 )
                    
                

#How to instantiate and use this class
#if __name__ == '__main__':
#    params = Parameters('params.yaml')
#    img = cv2.imread('../../box_in_scene.png')
#    rows,cols = img.shape[:2]
#    u = [ cols, rows ]
#    fovea = Fovea( u, params )
#    color = (255, 255, 255) # B G R
#    output = fovea.foveatedImage( img, params, color )
#    cv2.imshow( "foveated image", output )
#    cv2.waitKey( 0 )    
