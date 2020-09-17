#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
\file main2.py

\brief This file executable to multifovea structure in a dataset

\author
Petrucio Ricardo Tavares de Medeiros \n
Universidade Federal do Rio Grande do Norte \n 
Departamento de Computacao e Automacao Industrial \n
petrucior at ufrn (dot) edu (dot) br

\version 0.1
\date May 2020
 
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
import math
from parameters import Parameters
from fovea import Fovea
from multifovea import Multifovea
from statistics import Statistics

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


def calcError( point, groundtruth, index):
    '''
    \fn calcError( point, file )

    \brief This function calculates the error between the center of marked region and point

    \param point - Point (x, y) estimated
    \param groundtruth - Points defining the region of detection
    \param index - The line index of file
    '''
    x = groundtruth[index][0]; sx = groundtruth[index][0]+groundtruth[index][2];
    y = groundtruth[index][1]; sy = groundtruth[index][1]+groundtruth[index][3];
    center = ( (x + sx)/2, (y + sy)/2 )
    error = math.sqrt( pow( center[0] - point[0], 2.0 ) + pow( center[1] - point[1], 2.0 ) )
    #print( point, center, error )
    return center, error

def main():
    arguments = sys.argv[1:]
    if ( len( arguments ) != 5 ):
        print ('python main.py ''< model >'' ''< dataset >'' ''<yaml>'' ''<file_save_dates>'' ''<groundtruth>'' ')
        print ('Example: python main.py ~/<model>.png ~/<dataset>/ fovea.yaml data/graph.dat ~/groundtruth.txt')
    else:
        model = cv2.imread( arguments[0] )
        params = Parameters( arguments[2] )
        dir_path = arguments[1]
        fileName = arguments[3]
        groundtruth = arguments[4]
        
        lines = []
        file = open(groundtruth, 'r')
        for line in file:
            #lines.append( list( map(int, line.split(",")) ) )
            lines.append( list( map(int, line.split()) ) )
        
        # Starting Statistics
        statistics = Statistics()

        # Open File
        file = open(fileName, "w")

        # Creating Window
        #winname = 'image foveated'
        #cv2.namedWindow( winname )
        count = 0
        for _, _, arquivo in os.walk( dir_path ):
            arquivo.sort()
            for i in range(0, len(arquivo)):
                print("processing: "+str(int((i*100)/len(arquivo)))+"%")
                scene = cv2.imread( dir_path +'/'+arquivo[i] )
                
                '''
                # Ground truth
                #imgGroundTruth = cv2.rectangle( scene, ( lines[i][1], lines[i][2]), (lines[i][1]+lines[i][3], lines[i][2]+lines[i][4]), (255,0,0), 2 )
                imgGroundTruth = cv2.rectangle( scene, ( lines[i][0], lines[i][1]), (lines[i][0]+lines[i][2], lines[i][1]+lines[i][3]), (255,0,0), 2 )
                cv2.namedWindow( "ground truth", cv2.WINDOW_NORMAL)
                cv2.imshow("ground truth", imgGroundTruth)
                cv2.waitKey( 0 )
                '''

                # Updating size of scene
                rows, cols = scene.shape[:2]
                u = [ cols, rows ]
                #thresholdFeatures = 16
                #w = [90, 90]
                thresholdFeatures = 15
                w = [80, 80]
                params.updateParameterSizes( u, w, thresholdFeatures ) 

                
                '''
                for c in range(0, 4):
                    # Configuration with 4, 5, 6 and 7 foveas ( config 0, 1, 2, 3 )
                    params = statistics.poseFoveas( params, c )
                    #statistics.displayConfigurations( model, scene, params, [-1, -1] )
                    
                    ## Display multifovea
                    #multifovea = Multifovea( u, params )
                    #output = multifovea.multifoveatedImage( scene, params )
                    #cv2.imshow( winname, output ) 
                    
                    ############################################
                    ## MLE, Multi, Tril and Coord Barycentric ##
                    ############################################
                    # MLE
                    poseMLE = statistics.maximumLikelihoodEstimator( model, scene, params, 0.06, 0 )
                    #statistics.displayConfigurations( model, scene, params, poseMLE )
                    pose, errorMLE = calcError( poseMLE, lines, i )

                    # Multilateration
                    poseMult = statistics.multilaterationEstimator( model, scene, params )
                    #statistics.displayConfigurations( model, scene, params, poseMult )
                    pose, errorMult = calcError( poseMult, lines, i )

                    # Parameters configuration 
                    parameters = statistics.findMaximum( model, scene, params, [ [0, 1, 2], [1, 2, 3], [0, 2, 3], [0, 1, 3] ] )
                    if ( c == 1 ):
                        parameters = statistics.findMaximum( model, scene, params, [ [0, 1, 2], [1, 4, 2], [4, 3, 2], [3, 0, 2] ] )
                    if ( c == 2 ):
                        parameters = statistics.findMaximum( model, scene, params, [ [0, 1, 3], [1, 3, 4], [0, 3, 4], [0, 1, 4], [1, 2, 4], [2, 4, 5], [1, 4, 5], [1, 2, 5] ] )
                    if ( c == 3 ):
                        parameters = statistics.findMaximum( model, scene, params, [ [0, 1, 3], [1, 3, 4], [0, 3, 4], [0, 1, 4], [1, 2, 4], [2, 4, 5], [1, 4, 5], [1, 2, 5], [3, 4, 6], [4, 6, 7], [3, 6, 7], [3, 4, 7], [4, 5, 7], [5, 7, 8], [4, 7, 8], [4, 5, 8] ])
                    
                    # Trilateration
                    poseTril = statistics.trilaterationEstimator( model, scene, parameters )
                    #statistics.displayConfigurations( model, scene, params, poseTril )
                    pose, errorTril = calcError( poseTril, lines, i )
                
                    # Barycentric Coordinates 
                    poseBC = statistics.barycentricCoordinates( model, scene, parameters )
                    #statistics.displayConfigurations( model, scene, params, poseBC )
                    pose, errorBC = calcError( poseBC, lines, i )
                    
                    # Writing on the file
                    if ( c == 0 ):
                        file.write( str( count ) + '\t' +
                                    str( "{0:.2f}".format(round(errorMLE,2)) ) + '\t' + \
                                    str( "{0:.2f}".format(round(errorMult,2)) ) + '\t' + \
                                    str( "{0:.2f}".format(round(errorTril,2)) ) + '\t' + \
                                    str( "{0:.2f}".format(round(errorBC,2)) ) + '\t' )
                    else:
                        file.write( str( "{0:.2f}".format(round(errorMLE,2)) ) + '\t' + \
                                    str( "{0:.2f}".format(round(errorMult,2)) ) + '\t' + \
                                    str( "{0:.2f}".format(round(errorTril,2)) ) + '\t' + \
                                    str( "{0:.2f}".format(round(errorBC,2)) ) + '\t' )
                    
                '''


                for c in range(0, 3):
                
                    ######################################################
                    ## Potential Field - Reduce Region and Intersection ##
                    ######################################################

                    # Configurations
                    iterations = 10
                    jump = 10
                    delta = [ int(params.w[0]/2) + jump, int(params.w[1]/2) + jump ]
                    size = [ params.u[0] - int(params.w[0]/2) - jump, params.u[1] - int(params.w[1]/2) - jump ]
                    region = [ delta, size ]

                    # Reduce Region by Local Gradient
                    poseRRLG = statistics.reduceRegionByLocalGradient( model, scene, params, region, iterations, jump, c )
                    #statistics.displayConfigurations( model, scene, params, poseRRLG )
                    pose, errorRRLG = calcError( poseRRLG, lines, i )

                    # Intersection Local Gradient
                    poseILG = ( 0.0, 0.0 )
                    if ( c != 0 ):
                        poseILG = statistics.intersectionLocalGradient( model, scene, params, jump, c )
                    #statistics.displayConfigurations( model, scene, params, poseILG )
                    pose, errorILG = calcError( poseILG, lines, i )
                    

                    # Writing on the file
                    if ( c == 0 ):
                        file.write( str( count ) + '\t' +
                                    str( "{0:.2f}".format(round(errorRRLG,2)) ) + '\t' + \
                                    str( "{0:.2f}".format(round(errorILG,2)) ) + '\t' )
                    else:
                        file.write( str( "{0:.2f}".format(round(errorRRLG,2)) ) + '\t' + \
                                    str( "{0:.2f}".format(round(errorILG,2)) ) + '\t' )
               
                    
                    key = cv2.waitKey( 1 )
                    if ( key == ord('q') ):
                        break

                # Line break
                file.write('\n')
                count += 1
                
                
    cv2.destroyAllWindows()


if __name__ == '__main__':
    main()
