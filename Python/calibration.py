#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
\file calibration.py

\brief This file changes the parameters fovea

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
import sys, os
import numpy as np
from parameters import Parameters
from fovea import Fovea
from multifovea import Multifovea
from statistics import Statistics
from functools import partial

class Structure:
    def __init__( self, parameters, multifovea ):
        self.parameters = parameters
        self.multifovea = multifovea
        self.indexFovea = 0
    
    def updateFovea( self ):
        if ( self.indexFovea + 1 <= int( len(self.parameters.f)/2) ):
            self.indexFovea += 1
        else:
            self.indexFovea = 0
    
    def update( self, position ):
        #print( str(position[0]) + ", " + str(position[1]) )
        self.parameters.f[2*self.indexFovea] = position[0]
        self.parameters.f[2*self.indexFovea+1] = position[1]
        self.parameters.fixFovea()
        self.multifovea.updateFovea( self.indexFovea, self.parameters )
    

def callback( event, x, y, flags, structure ):
    '''
    \fn callback( event, x, y, flags, params )
    
    \brief This function is responsible for controlling all actions with the user's mouse
    '''
    #print( str(x) + ", " + str(y) )
    position = [ x - (structure.parameters.u[0]/2), y - (structure.parameters.u[1]/2) ]
    structure.update( position )
    
def filter( value, parameters, index, multifovea ):
    parameters.thresholdFeatures = value
    multifovea.updateFovea( index, parameters )
    
def main():    
    arguments = sys.argv[1:]
    if ( len( arguments ) != 3 ):
        print ('python main.py ''< model >'' ''< dataset >'' ''<yaml>'' ')
        print ('Example: python main.py ~/<model> ~/<scene> params.yaml')
    else:
        model = cv2.imread( arguments[0] )
        scene = cv2.imread( arguments[1] )
        params = Parameters( arguments[2] )
        
        # Starting Statistics
        statistics = Statistics()
        
        winname = 'image foveated'
        cv2.namedWindow( winname )
        
        # Updating size of scene
        rows, cols = scene.shape[:2]
        u = [ cols, rows ]
        params.updateParameterSizeImage( u )
        
        index = 0
        multifovea = Multifovea( u, params )
        structure = Structure( params, multifovea )
        
        trackbar_name = 'Filter'
        cv2.createTrackbar(trackbar_name, winname , 0, 100, partial(filter, parameters=params, index=index, multifovea=multifovea))
        
        while( True ):
            # bind the callback function to window
            cv2.setMouseCallback( winname, callback, structure )
            output = multifovea.multifoveatedImage( scene, params )
            cv2.imshow( winname, output )
            key = cv2.waitKey( 2 )
            if ( key == ord('q') ):
                break
            # Controls from fovea
            if ( key == ord('d') ):
                params.w[0] = min( params.u[0], params.w[0] + 10 )
            if ( key == ord('a') ):
                params.w[0] = max( 1, params.w[0] - 10 )
            if ( key == ord('c') ):
                params.w[1] = min( params.u[1] - 1, params.w[1] + 10 )
            if ( key == ord('z') ):
                params.w[1] = max( 1, params.w[1] - 10 )
            if ( key == ord('m') ): # changes the index of the fovea
                index += 1
                index = index % int( len(params.f)/2 )
            if ( key == ord('g') ):
                statistics.plotProportion( model, scene, fovea, parameters )
                statistics.plotFunction( model, scene, fovea, params )
            
        # Configuration with 4 foveas ( config 0 )
        #params = statistics.poseFoveas( params, 0 )
        #parameters = statistics.findMaximum( model, scene, params, [ [0, 1, 2], [1, 2, 3], [0, 2, 3], [0, 1, 3] ] )
        #pose = statistics.trilaterationEstimator( model, scene, parameters ) 
        #poseMLE0 = statistics.maximumLikelihoodEstimator( model, scene, params, 0.1, 0 )
        #statistics.displayConfigurations( model, scene, params, pose )
        
        # Configuration with 5 foveas ( config 1 )
        #params = statistics.poseFoveas( params, 1 )
        #parameters = statistics.findMaximum( model, scene, params, [ [0, 1, 2], [0, 2, 3], [1, 2, 4], [2, 3, 4] ] )
        #pose = statistics.trilaterationEstimator( model, scene, parameters ) 
        #poseMLE1 = statistics.maximumLikelihoodEstimator( model, scene, params, 0.1, 0 )
        #statistics.displayConfigurations( model, scene, params, pose )   
        
        # Configuration with 6 foveas ( config 2 )
        #params = statistics.poseFoveas( params, 2 )
        #parameters = statistics.findMaximum( model, scene, params, [ [0, 1, 2], [0, 2, 3], [1, 2, 4], [2, 3, 4] ] )
        #pose = statistics.trilaterationEstimator( model, scene, parameters ) 
        #poseMLE2 = statistics.maximumLikelihoodEstimator( model, scene, params, 0.1, 0 )
        #statistics.displayConfigurations( model, scene, params, [-1, -1] )   
        
        # Configuration with 9 foveas ( config 3 )
        #params = statistics.poseFoveas( params, 3 )
        #poseMLE3 = statistics.maximumLikelihoodEstimator( model, scene, params, 0.1, 0 )
        
        #print("MLE: "+str(poseMLE0)+" "+str(poseMLE1)+" "+str(poseMLE2)+" "+str(poseMLE3))
        
        #multifovea = Multifovea( u, params )
        #output = multifovea.multifoveatedImage( scene, params )
        #cv2.imshow( winname, output )
        
    cv2.destroyAllWindows()

if __name__ == '__main__':
    main()

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

