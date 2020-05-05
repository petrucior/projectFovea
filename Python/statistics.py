#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
\file statistics.py

\brief This file contains the prototype of statistics from fovea.

\author
Petrucio Ricardo Tavares de Medeiros \n
Universidade Federal do Rio Grande do Norte \n 
Departamento de Computacao e Automacao Industrial \n
petrucior at ufrn (dot) edu (dot) br

\version 0.1
\date January 2020
 
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
import numpy
from parameters import Parameters
from level import Level
from fovea import Fovea
from multifovea import Multifovea
#from feature import Features

from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np


class Statistics :
    '''
    \class Statistics
    
    \brief This class implements the Statistics TAD
    '''
    
    def __init__( self ):
        '''
        \fn Statistics( )
        
        \brief Constructor default
        '''
        print("Welcome to class of data manipulation")
        

    def plotProportion( self, model, scene, fovea, parameters ):
        '''
        \fn plotProportion( model, scene, fovea, parameters )
        
        \brief Plot proportion

        \param model - Template that represent the visual stimulli
        \param scene - Scene to be foveated
        \param fovea - Structure that contain the features ( keypoints, descriptors and matches )
        \param parameters - Parameters of fovea structure
        '''
        file = open("data/graph3d.dat", "w")
        stepx = 60; stepy = 60;
        complete = False
        savedValue = 0.0
        # loop image
        for i in range( 0, parameters.u[1], stepx ): # rows
            for j in range( 0, parameters.u[0], stepy ): # cols
                file.write(str(i) + '\t' + str(j) + '\t')
                # level
                for k in range( 0, parameters.m + 1 ):
                    pixel = [i, j]
                    parameters.setFovea( pixel )
                    fovea.updateFovea( parameters )
                    #output = fovea.foveatedImage( scene, params )
                    #cv2.imshow( "foveated image", output )
                    #cv2.waitKey( 0 )
                    fovea.foveatedFeatures( scene, params )
                    fovea.foveatedMatching( model, params )
                    file.write(str(fovea.features.inlierRateSURF[k]) + '\t')
                file.write('\n')


    def plotFunction( self, model, scene, fovea, parameters ):
        '''
        \fn plotFunction( model, scene, fovea, parameters )
        
        \brief Plot function

        \param model - Template that represent the visual stimulli
        \param scene - Scene to be foveated
        \param fovea - Structure that contain the features ( keypoints, descriptors and matches )
        \param parameters - Parameters of fovea structure
        '''
        file = open("data/function.dat", "w")
        stepx = 20; stepy = 20;
        complete = False
        savedValue = 0.0
        # loop image
        for i in range( 0, parameters.u[1], stepx ): # rows
            for j in range( 0, parameters.u[0], stepy ): # cols
                file.write(str(i) + '\t' + str(j) + '\t')
                pixel = [i, j]
                parameters.setFovea( pixel )
                fovea.updateFovea( parameters )
                #output = fovea.foveatedImage( scene, params )
                #cv2.imshow( "foveated image", output )
                #cv2.waitKey( 0 )
                fovea.foveatedFeatures( scene, params )
                fovea.foveatedMatching( model, params )
                #file.write(str(fovea.features.inlierRateSURF[k]) + '\t')
                file.write( str(self.weightedFunctionFovea( fovea, parameters )) + '\t' )
                file.write('\n')

    
    def weightedFunctionFovea( self, fovea, parameters ):
        '''
        \fn weightedFunctionFovea( self, fovea, parameters )
        
        \brief Function to calculate the function associated to fovea
        
        \param fovea - Structure that contain the features ( keypoints, descriptors and matches )
        \param parameters - Parameters of fovea structure
        '''
        p = parameters
        self.weightedFunction = []
        # Adding all regions
        self.regionSum = [0, 0]
        for r in range( 0, p.m + 1 ):
            self.regionSum = numpy.add( self.regionSum, fovea.levels[r].getSize( parameters ) )
        # Weighting levels
        for i in range( 0, p.m + 1 ):
            den = fovea.levels[p.m - i].getSize( parameters )
            self.weightedFunction.append( numpy.divide( den , self.regionSum ) )
        #print( self.weightedFunction )
        # Weighting fovea
        self.fpdf = 0.0
        for k in range( 0, p.m + 1 ):
            inlierRate = fovea.features.inlierRateSURF[k]
            self.fpdf += self.weightedFunction[k][0] * inlierRate + self.weightedFunction[k][1] * inlierRate
        #print( self.fpdf )
        return self.fpdf
        

    def weightedFunctionMultifovea( self, multifovea, parameters ):
        
        '''
        \fn weightedFunctionMultifovea( self, multifovea, parameters )
        
        \brief Function to calculate the function associated to multifovea
        
        \param multifovea - Structure that contain the features ( keypoints, descriptors and matches )
        \param parameters - Parameters of fovea structure
        '''
        print( parameters.m + 1 )
        
        
#How to instantiate and use this class
if __name__ == '__main__':
    params = Parameters('params.yaml')
    model = cv2.imread('../../box.png')
    scene = cv2.imread('../../box_in_scene.png')
    rows,cols = scene.shape[:2]
    u = [ cols, rows ]
    fovea = Fovea( u, params )
    #output = fovea.foveatedImage( scene, params )
    #cv2.imshow( "foveated image", output )
    #cv2.waitKey( 0 )
    fovea.foveatedFeatures( scene, params )
    fovea.features.matchingFeatures( model, params )
    statistics = Statistics()
    #statistics.plotProportion( model, scene, fovea, params )
    #statistics.plotFunction( model, scene, fovea, params )
    #statistics.weightedFunctionFovea( fovea, params )
    cv2.destroyAllWindows()
    
'''
    params = Parameters('params.yaml')
    model = cv2.imread('../../box.png')
    scene = cv2.imread('../../box_in_scene.png')
    rows,cols = scene.shape[:2]
    u = [ cols, rows ]
    #print( u )
    fovea = Fovea( u, params )
    #output = fovea.foveatedImage( scene, params )
    #cv2.imshow( "foveated image", output )
    #cv2.waitKey( 0 )
    fovea.foveatedFeatures( scene, params )
    fovea.features.matchingFeatures( model, params )
    statistics = Statistics()
    #statistics.plotProportion( model, scene, fovea, params )
    statistics.weightedFunctionFovea( fovea, params )
    cv2.destroyAllWindows()
'''
