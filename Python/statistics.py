#!/usr/bin/pythoAn
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
import math
import numpy
import random
import copy
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
        multiFPDF = []
        pArray = []
        for f in range( 0, len(parameters.f), 2 ):
            p = copy.copy(parameters)
            p.updateParameterFoveas( [ parameters.f[f], parameters.f[f+1] ], p.colors )
            pArray.append( p )
        for f in range( 0, int(len(parameters.f) / 2) ):
            multiFPDF.append( self.weightedFunctionFovea( multifovea.foveas[f], pArray[f] )  )
        return multiFPDF


    def localGradient( self, referencePoint, model, scene, parameters, jump, config ):
        '''
        \fn localGradient( referencePoint, model, scene, parameters, jump, config )
        
        \brief Function that computes the direction of local gradient
        
        \param referencePoint - Reference fovea location (x, y)
        \param model - Template that represent the visual stimulli
        \param scene - Scene to be foveated
        \param parameters - Parameters of fovea structure
        \param jump - Distance between neighborhood foveas
        \param config - Indicates how the potentials were extracted, both in the clockwise direction
        - If config == 0, then potentials extracted from up (North), right (EAST), down (SOUTH) and left (WEST)
        - If config == 1, then potentials extracted from northeast, southeast, southwest and northwest
        - if config == 2, then potentials extracted from 8 regions
        
        #        vt
        #  vtl   ^   vtr
        #     `  | ´
        # vl <--   --> vr
        #     ,  | . 
        # vld    v   vrd
        #        vd
        
        \return Index related to increased local potential
        '''
        # Fixing the reference
        referencePoint = [ referencePoint[0] - int(parameters.u[0]/2), referencePoint[1] - int(parameters.u[1]/2) ]
        # Creating foveas
        foveas = []
        foveas.append( referencePoint[0] ); foveas.append( referencePoint[1] );
        if ( config == 0 ):
            '''
            # Debug
            # reference ( red )
            cv2.circle(scene, (referencePoint[0] + int(parameters.u[0]/2), referencePoint[1] + int(parameters.u[1]/2)), 1, (0, 0, 255), -1)
            # top ( blue )
            cv2.circle(scene, (referencePoint[0] + int(parameters.u[0]/2), referencePoint[1] - jump + int(parameters.u[1]/2)), 1, (255, 0, 0), -1)
            # right ( green )
            cv2.circle(scene, (referencePoint[0] + jump + int(parameters.u[0]/2), referencePoint[1] + int(parameters.u[1]/2)), 1, (0, 255, 0), -1)
            # down ( yellow )
            cv2.circle(scene, (referencePoint[0] + int(parameters.u[0]/2), referencePoint[1] + jump + int(parameters.u[1]/2)), 1, (0, 255, 255), -1)
            # left ( pink )
            cv2.circle(scene, (referencePoint[0] - jump + int(parameters.u[0]/2), referencePoint[1] + int(parameters.u[1]/2)), 1, (255, 0, 255), -1)
            cv2.imshow('posicao', scene)
            cv2.waitKey( 0 )
            '''
            foveas.append( referencePoint[0] ); foveas.append( referencePoint[1] - jump ); # Top
            foveas.append( referencePoint[0] + jump ); foveas.append( referencePoint[1] ); # Right
            foveas.append( referencePoint[0] ); foveas.append( referencePoint[1] + jump ); # Down
            foveas.append( referencePoint[0] - jump ); foveas.append( referencePoint[1] ); # Left
        if ( config == 1 ):
            '''
            # Debug
            # reference ( red )
            cv2.circle(scene, (referencePoint[0] + int(parameters.u[0]/2), referencePoint[1] + int(parameters.u[1]/2)), 1, (0, 0, 255), -1)
            # northeast ( blue )
            cv2.circle(scene, (referencePoint[0] + jump + int(parameters.u[0]/2), referencePoint[1] - jump + int(parameters.u[1]/2)), 1, (255, 0, 0), -1)
            # southeast ( green )
            cv2.circle(scene, (referencePoint[0] + jump + int(parameters.u[0]/2), referencePoint[1] + jump + int(parameters.u[1]/2)), 1, (0, 255, 0), -1)
            # southwest ( yellow )
            cv2.circle(scene, (referencePoint[0] - jump + int(parameters.u[0]/2), referencePoint[1] + jump + int(parameters.u[1]/2)), 1, (0, 255, 255), -1)
            # northwest ( pink )
            cv2.circle(scene, (referencePoint[0] - jump + int(parameters.u[0]/2), referencePoint[1] - jump + int(parameters.u[1]/2)), 1, (255, 0, 255), -1)
            cv2.imshow('posicao', scene)
            cv2.waitKey( 0 )
            '''
            foveas.append( referencePoint[0] + jump ); foveas.append( referencePoint[1] - jump ); # Northeast
            foveas.append( referencePoint[0] + jump ); foveas.append( referencePoint[1] + jump ); # Southeast
            foveas.append( referencePoint[0] - jump ); foveas.append( referencePoint[1] + jump ); # Southwest
            foveas.append( referencePoint[0] - jump ); foveas.append( referencePoint[1] - jump ); # Northwest
        if ( config == 2 ):
            foveas.append( referencePoint[0] ); foveas.append( referencePoint[1] - jump ); # Top
            foveas.append( referencePoint[0] + jump ); foveas.append( referencePoint[1] - jump ); # Northeast
            foveas.append( referencePoint[0] + jump ); foveas.append( referencePoint[1] ); # Right
            foveas.append( referencePoint[0] + jump ); foveas.append( referencePoint[1] + jump ); # Southeast
            foveas.append( referencePoint[0] ); foveas.append( referencePoint[1] + jump ); # Down
            foveas.append( referencePoint[0] - jump ); foveas.append( referencePoint[1] + jump ); # Southwest
            foveas.append( referencePoint[0] - jump ); foveas.append( referencePoint[1] ); # Left
            foveas.append( referencePoint[0] - jump ); foveas.append( referencePoint[1] - jump ); # Northwest
        
        params = copy.copy( parameters )
        colors = []
        for i in range( 0, int(len( foveas )/2) ):
            colors.append( random.randint(0, 255) ) # blue
            colors.append( random.randint(0, 255) ) # green
            colors.append( random.randint(0, 255) ) # red           
        params.updateParameterFoveas( foveas, colors )
        multifovea = Multifovea( params.u, params )
        multifovea.multifoveatedFeatures( scene, params )
        multifovea.multifoveatedMatching( model, params )
        # Display fovea
        #output = multifovea.multifoveatedImage( scene, params )
        #cv2.imshow( "multifoveated image", output )
        #cv2.waitKey( 0 )
        #cv2.destroyAllWindows()        

        # FPDF
        potentialFromMultifovea = self.weightedFunctionMultifovea( multifovea, params )
        #print( potentialFromMultifovea )
        referencePotential = potentialFromMultifovea[0]
        potentials = potentialFromMultifovea[1:]

        # Returning the index related to increased local potential
        maxPotential = 0.0; index = 0;
        for p in range( 0, len(potentials) ):
            if ( maxPotential < ( potentials[p] - referencePotential ) ):
                maxPotential = ( potentials[p] - referencePotential )
                index = p
        return index
        

    def reduceRegionByLocalGradient( self, model, scene, parameters, regionUnderAnalysis, iterations, jump, config ):
        '''
        \fn reduceRegionByLocalGradient( model, scene, parameters, regionUnderAnalysis, iterations, jump, config )
        
        \brief Function to calculate the potential local approach
        
        \param model - Template that represent the visual stimulli
        \param scene - Scene to be foveated
        \param parameters - Parameters of fovea structure
        \param regionUnderAnalysis - Region that needs to be investigated [ [xi, yi], [xf, yf] ]
        \param iterations - Number of search iterations
        \param jump - Distance between neighborhood foveas
        \param config - Indicates how the potentials were extracted, both in the clockwise direction
        - If config == 0, then potentials extracted from up (North), right (EAST), down (SOUTH) and left (WEST)
        - If config == 1, then potentials extracted from northeast, southeast, southwest and northwest
        - if config == 2, then potentials extracted from 8 regions
        '''
        # This method will work until iterations to be null
        if ( iterations == 0 ):
            print( regionUnderAnalysis )
            return regionUnderAnalysis;
        else :
            print( regionUnderAnalysis )
            xmin = regionUnderAnalysis[0][0]; ymin = regionUnderAnalysis[0][1];
            xmax = regionUnderAnalysis[1][0]; ymax = regionUnderAnalysis[1][1];
            # It will work until regionUnderAnalysis is bigger than 25 x 25
            if ( xmax - xmin > int(parameters.w[0]/2) and ymax - ymin > int(parameters.w[1]/2) ):
                # Chosen point
                x = random.randint(xmin, xmax)
                y = random.randint(ymin, ymax)
                print( x, y )
                cv2.rectangle(scene,(xmin, ymin),(xmax, ymax),(0,255,0),3)
                cv2.rectangle(scene,(0, 30), (20, 60 ), (100, 100, 0), 2 )
                cv2.circle(scene, (x, y), 3, (0, 0, 255), -1)
                cv2.imshow('region', scene)
                cv2.waitKey( 0 )
                # Detecting bigger potential
                index = self.localGradient( [x, y], model, scene, parameters, jump, config )
                if ( config == 0 ): # index related to up (North), right (EAST), down (SOUTH) and left (WEST)
                    if ( index == 0 ): # top
                        ymax = y
                    if ( index == 1 ): # right
                        xmin = x
                    if ( index == 2 ): # down
                        ymin = y;
                    if ( index == 3 ): # left
                        xmax = x
                        
                if ( config == 1 ): # index related to northeast, southeast, southwest and northwest
                    if ( index == 0 ): # notheast
                        xmin = x; ymax = y;
                    if ( index == 1 ): # southeast
                        xmin = x; ymin = y;
                    if ( index == 2 ): # southwest
                        xmax = x; ymin = y;
                    if ( index == 3 ): # northwest
                        xmax = x; ymax = y;

                if ( config == 2 ):
                    if ( index == 0 ): # top
                        ymax = y
                    if ( index == 1 ): # notheast
                        xmin = x; ymax = y;
                    if ( index == 2 ): # right
                        xmin = x
                    if ( index == 3 ): # southeast
                        xmin = x; ymin = y;
                    if ( index == 4 ): # down
                        ymin = y;
                    if ( index == 5 ): # southwest
                        xmax = x; ymin = y;
                    if ( index == 6 ): # left
                        xmax = x
                    if ( index == 7 ): # northwest
                        xmax = x; ymax = y;

                cv2.rectangle(scene,(xmin, ymin),(xmax, ymax),(244,255,0),3)
                cv2.imshow('region', scene)
                cv2.waitKey( 0 )
                regionUnderAnalysis = [ [xmin, ymin], [xmax, ymax] ]
                self.reduceRegionByLocalGradient( model, scene, parameters, regionUnderAnalysis, iterations - 1, jump, config )

            else :
                return regionUnderAnalysis
            
    
    def intersectionLocalGradient( self, model, scene, parameters, jump, config):
        '''
        \fn intersectionLocalGradient( model, scene, parameters, jump, config )
        
        \brief Function that computes the direction of local gradient
        
        \param model - Template that represent the visual stimulli
        \param scene - Scene to be foveated
        \param parameters - Parameters of fovea structure
        \param jump - Distance between neighborhood foveas
        \param config - Indicates how the potentials were extracted, both in the clockwise direction
        - If config == 0, then potentials extracted from up (North), right (EAST), down (SOUTH) and left (WEST)
        - If config == 1, then potentials extracted from northeast, southeast, southwest and northwest
        - if config == 2, then potentials extracted from 8 regions
        '''
        # Chosen points a and c
        xa = random.randint(0 + parameters.w[0]/2 + jump, parameters.u[0] - parameters.w[0]/2 - jump)
        ya = random.randint(0 + parameters.w[1]/2 + jump, parameters.u[1]  - parameters.w[1]/2 - jump)
        xc = random.randint(0 + parameters.w[0]/2 + jump, parameters.u[0]  - parameters.w[0]/2 - jump)
        yc = random.randint(0 + parameters.w[1]/2 + jump, parameters.u[1]  - parameters.w[1]/2 - jump)
        # Finding the potential local
        indexA = self.localGradient( [xa, ya], model, scene, parameters, jump, config )
        indexC = self.localGradient( [xc, yc], model, scene, parameters, jump, config )
        #print( indexA, indexC )
        # Defining points b and d (Neighbors of the chosen points)
        xb = 0; yb = 0;
        xd = 0; yd = 0;
        if ( config == 0 ):
            # Point B
            if ( indexA == 0 ): # Top
                xb = xa; yb = ya - jump;
            if ( indexA == 1 ): # Right
                xb = xa + jump; yb = ya;
            if ( indexA == 2 ): # Down
                xb = xa; yb = ya + jump;
            if ( indexA == 3 ): # Left
                xb = xa - jump; yb = ya;
            # Point D
            if ( indexC == 0 ): # Top
                xd = xc; yd = yc - jump;
            if ( indexC == 1 ): # Right
                xd = xc + jump; yd = yc;
            if ( indexC == 2 ): # Down
                xd = xc; yd = yc + jump;
            if ( indexC == 3 ): # Left
                xd = xc - jump; yd = yc;
        if ( config == 1 ):
            # Point B
            if ( indexA == 0 ): # Northeast
                xb = xa + jump; yb = ya - jump;
            if ( indexA == 1 ): # Southeast
                xb = xa + jump; yb = ya + jump;
            if ( indexA == 2 ): # Southwest 
                xb = xa - jump; yb = ya + jump;
            if ( indexA == 3 ): # Northwest
                xb = xa - jump; yb = ya - jump;
            # Point D
            if ( indexC == 0 ): # Northeast
                xd = xc + jump; yd = yc - jump;
            if ( indexC == 1 ): # Southeast
                xd = xc + jump; yd = yc + jump;
            if ( indexC == 2 ): # Southwest
                xd = xc - jump; yd = yc + jump;
            if ( indexC == 3 ): # Northwest
                xd = xc - jump; yd = yc - jump;
        if ( config == 2 ):
            # Point B
            if ( indexA == 0 ): # Top
                xb = xa; yb = ya - jump;
            if ( indexA == 1 ): # Northeast
                xb = xa + jump; yb = ya - jump;
            if ( indexA == 2 ): # Right
                xb = xa + jump; yb = ya;
            if ( indexA == 3 ): # Southeast
                xb = xa + jump; yb = ya + jump;
            if ( indexA == 4 ): # Down
                xb = xa; yb = ya + jump;
            if ( indexA == 5 ): # Southwest 
                xb = xa - jump; yb = ya + jump;
            if ( indexA == 6 ): # Left
                xb = xa - jump; yb = ya;
            if ( indexA == 7 ): # Northwest
                xb = xa - jump; yb = ya - jump;
            # Point D
            if ( indexC == 0 ): # Top
                xd = xc; yd = yc - jump;
            if ( indexC == 1 ): # Northeast
                xd = xc + jump; yd = yc - jump;
            if ( indexC == 2 ): # Right
                xd = xc + jump; yd = yc;
            if ( indexC == 3 ): # Southeast
                xd = xc + jump; yd = yc + jump;
            if ( indexC == 4 ): # Down
                xd = xc; yd = yc + jump;
            if ( indexC == 5 ): # Southwest
                xd = xc - jump; yd = yc + jump;
            if ( indexC == 6 ): # Left
                xd = xc - jump; yd = yc;
            if ( indexC == 7 ): # Northwest
                xd = xc - jump; yd = yc - jump;

        '''
        # Display in the image the intersection between lines
        #xa = 333; ya = 303; xb = 338; yb = 298;
        #xc = 219; yc = 283; xd = 214; yd = 278;
        print( xa, ya ); print( xb, yb );
        print( xc, yc ); print( xd, yd );
        b = random.randint(0, 255); g = random.randint(0, 255); r = random.randint(0, 255);
        cv2.line(scene,(xa, ya),(xb, yb),(b,g,r),2)
        cv2.line(scene,(xc, yc),(xd, yd),(b,g,r),2)
        cv2.imshow('intersection', scene)
        cv2.waitKey( 0 )
        '''
        
        '''
              x    y    1
        A =   xa   ya   1
              xb   yb   1
        det( A ) = 0
        ( x * ya * 1 ) + ( y * 1 * xb ) + ( xa * yb * 1 ) - ( xb * ya * 1 ) - ( x * yb * 1 ) - ( y * xa * 1 ) = 0
        xya + yxb + xayb - xbya - xyb - yxa = 0
        eq1:    x ( ya - yb ) + y ( xb - xa ) = xbya - xayb

              x    y   1
        C =   xc   yc  1
              xd   yd  1
        det( C ) = 0
        ( x * yc * 1 ) + ( y * 1 * xd ) + ( xc * yd * 1 ) - ( xd * yc * 1 ) - ( x * yd * 1 ) - ( y * xc * 1 ) = 0
        xyc + yxd + xcyd - xdyc - xyd - yxc = 0
        eq2:    x ( yc - yd ) + y ( xd - xc ) = xdyc - xcyd

        if ( ya == yb )
           y = ( xbya - xayb ) / ( xb - xa )
           x = ( ( xdyc - xcyd ) - y ( xd - xc ) ) / ( yc - yd )

        if ( xb == xa )
           x = ( xbya - xayb ) / ( ya - yb )
           y = ( xdyc - xcyd ) - x ( yc - yd ) / ( xd - xc )

        if ( yc == yd )
           y = ( xdyc - xcyd ) / ( xd - xc )
           x = ( xbya - xayb ) - y ( xb - xa ) / ( ya - yb )

        if ( xd == xc )
           x = ( xdyc - xcyd ) / ( yc - yd )
           y = ( xbya - xayb ) - x ( ya - yb ) / ( xb - xa )
        
        '''
        
        x = int(parameters.u[0]/2); y = int(parameters.u[1]/2);
        if ( ( xb - xa != 0.0 ) and ( xd - xc != 0.0 ) ):
            # angular coefficient
            mA = (yb - ya)/(xb - xa)
            mB = (yd - yc)/(xd - xc)
            #print( mA, mB )
            # linear coefficient
            # y = mx - mx_{1} + y_{1} => y = mx + c
            cA = mA * xa + ya
            cB = mB * xc + yc
            if ( mA == mB ):
                '''
                print("Angular coefficients are equals, then parallel ")
                if ( cA != cB ):
                    print("and distintes lines")
                if ( cA == cB ):
                    print("and coincident lines")
                '''
                # recalculating
                return self.intersectionLocalGradient( model, scene, parameters, jump, config );
            else:
                if ( mA * mB == -1 ):
                    #print("Perpendicular lines")
                    x = ( yc - ya + (mA * xa) - (mB * xc) ) / (mA - mB)
                    y = mA * ( x - xa ) + ya
                    return x, y
                
                if ( ( mA == 0 ) or ( mB == 0 ) ):
                    #print("Vertical parallel line")
                    x = ( yc - ya + (mA * xa) - (mB * xc) ) / (mA - mB)
                    y = mA * ( x - xa ) + ya
                    return x, y
                    # recalculating
                    #self.intersectionLocalGradient(  model, scene, parameters, jump, config )
                    
        else:
            #print("Horizontal parallel lines")
            if ( ya == yb ):
                y = ( xb * ya - xa * yb ) / ( xb - xa )
                x = ( ( xd * yc - xc * yd ) - y * ( xd - xc ) ) / ( yc - yd )
                
            if ( xb == xa ):
                x = ( xb * ya - xa * yb ) / ( ya - yb )
                y = ( ( xd * yc - xc * yd ) - x * ( yc - yd ) ) / ( xd - xc )
            
            if ( yc == yd ):
                y = ( xd * yc - xc * yd ) / ( xd - xc )
                x = ( ( xb * ya - xa * yb ) - y * ( xb - xa ) ) / ( ya - yb )
                
            if ( xd == xc ):
                x = ( xd * yc - xc * yd ) / ( yc - yd )
                y = ( ( xb * ya - xa * yb ) - x * ( ya - yb ) ) / ( xb - xa )
        
        return x, y
        
    
    def maximumLikelihoodEstimator( self, model, scene, parameters, threshold,  method ):
        '''
        \fn maximumLikelihoodEstimator( model, scene, parameters, threshold, method )
        
        \brief Calculate the Maximum Likelihood Estimator (MLE)
        
        \param model - Template that represent the visual stimulli
        \param scene - Scene to be foveated
        \param parameters - Parameters of fovea structure
        \param threshold - Threshold is a potential limiter
        \param method - Arithmetic mean (0) and weighted average (1)
        '''
        print("MLE")
        multifovea = Multifovea( parameters.u, parameters )
        multifovea.multifoveatedFeatures( scene, params )
        multifovea.multifoveatedMatching( model, params )
        potentials = self.weightedFunctionMultifovea( multifovea, params )
        # Transforming points to cartesian domain
        pose = []
        for f in range( 0, len(parameters.f), 2 ):
            pose.append( [ parameters.f[f] + int( parameters.u[0]/2 ), parameters.f[f+1] + int( parameters.u[1]/2 ) ] )
        x = 0.0; y = 0.0;
        if ( method == 0 ): # Arithmetic mean
            counter = 0
            for i in range( 0, len( potentials ) ):
                if ( potentials[i] >= threshold ):
                    x += pose[i][0]
                    y += pose[i][1]
                    counter += 1
            if ( counter == 0 ):
                x = 0.0; y = 0.0;
            else:
                x /= counter
                y /= counter
        
        if ( method == 1 ): # Weighted average
            for i in range( 0, len( potentials ) ):
                x += pose[i][0] * potentials[i]
                y += pose[i][1] * potentials[i]
            x /= len( potentials )
            y /= len( potentials )
        
        return x, y
    
    
    def trilaterationEstimator( self, model, scene, parameters ):
        '''
        \fn trilaterationEstimator( model, scene, parameters )
        
        \brief Calculate the trilateration Estimator
        
        \param model - Template that represent the visual stimulli
        \param scene - Scene to be foveated
        \param parameters - Parameters of fovea structure
        '''
        print("Trilateration")
        multifovea = Multifovea( parameters.u, parameters )
        multifovea.multifoveatedFeatures( scene, params )
        multifovea.multifoveatedMatching( model, params )
        potentials = self.weightedFunctionMultifovea( multifovea, params )
        # Descending ordering of potential
        potentialsAux = sorted( potentials, reverse=True )
        potential = potentialsAux[:3]
        first = 0; second = 0; third = 0;
        for p in range( 0, len(potentials)):
            if ( potentials[p] == potential[0] ):
                first = p * 2
            if ( potentials[p] == potential[1] ):
                second = p * 2
            if ( potentials[p] == potential[2] ):
                third = p * 2
        
        # Points
        x1 = parameters.f[first]; y1 = parameters.f[first+1];
        x2 = parameters.f[second]; y2 = parameters.f[second+1];
        x3 = parameters.f[third]; y3 = parameters.f[third+1];

        # Transforming points to cartesian domain
        x1 += int( parameters.u[0]/2 ); y1 += int( parameters.u[1]/2 );
        x2 += int( parameters.u[0]/2 ); y2 += int( parameters.u[1]/2 );
        x3 += int( parameters.u[0]/2 ); y3 += int( parameters.u[1]/2 );
        #print( x1, y1, x2, y2, x3, y3 )
        
        # Inverse of detection rate
        r1 = ((1.0 - potential[0]) * 0.001 )/ 0.0001
        r2 = ((1.0 - potential[1]) * 0.001 )/ 0.0001
        r3 = ((1.0 - potential[2]) * 0.001 )/ 0.0001
        
        a = (-2 * x1) + (2 * x2)
        b = (-2 * y1) + (2 * y2)
        c = (r1*r1) - (r2*r2) - (x1*x1) + (x2*x2) - (y1*y1) + (y2*y2)
        d = (-2 * x2) + (2 * x3)
        e = (-2 * y2) + (2 * y3)
        f = (r2*r2) - (r3*r3) - (x2*x2) + (x3*x3) - (y2*y2) + (y3*y3)

        den = (a * e) - (b * d)
        x = ((c * e) - (b * f)) / den
        y = ((a * f) - (c * d)) / den

        return x, y

    
    def multilaterationEstimator( self, model, scene, parameters ):
        '''
        \fn multilaterationEstimator( model, scene, parameters )
        
        \brief Calculate the multilateration Estimator
        
        \param model - Template that represent the visual stimulli
        \param scene - Scene to be foveated
        \param parameters - Parameters of fovea structure
        '''
        print("Multilateration")
        multifovea = Multifovea( parameters.u, parameters )
        multifovea.multifoveatedFeatures( scene, params )
        multifovea.multifoveatedMatching( model, params )
        potentials = self.weightedFunctionMultifovea( multifovea, params )
        # Inverse of detection rate
        invpotential = []
        for i in range( 0, len( potentials ) ):
            invpotential.append( ((1.0 - potentials[i]) * 0.001 )/ 0.0001 )
        #print( invpotential )
        x = 0.0; y = 0.0;
        if ( len(parameters.f) == 0 ): return x, y
        m = int(len(parameters.f)/2)
        b = numpy.zeros( (m - 1, 1), dtype=numpy.float64 )
        A = numpy.zeros( (m - 1, 2), dtype=numpy.float64 )
        #print( parameters.f )
        # Transforming points to cartesian domain
        pose = []
        for f in range( 0, len(parameters.f), 2 ):
            pose.append( [ parameters.f[f] + int( parameters.u[0]/2 ), parameters.f[f+1] + int( parameters.u[1]/2 ) ] )
        #print( pose )
        for i in range( 0, m - 1 ):
            A[i][0]  = ( -2 * pose[i][0] ) + ( 2 * pose[i+1][0] )            
            A[i][1]  = ( -2 * pose[i][1] ) + ( 2 * pose[i+1][1] )
            b[i][0] = pow( invpotential[i], 2.0 ) - pow( invpotential[i+1], 2.0 ) - pow( pose[i][0], 2.0 ) + pow( pose[i+1][0], 2.0 ) - pow( pose[i][1], 2.0 ) + pow( pose[i+1][1], 2.0 )
        #print( A, b )

        # Transposed A
        at = A.transpose()
        # Multiplication between transposed A and A
        mult = numpy.dot( at, A )
        # Inverse of the multiplication
        inv = numpy.linalg.inv( mult )
        # Multiplication between inv and transposed A
        multInv = numpy.dot( inv, at )
        # Multiplication between multInv and b
        estimated = numpy.dot( multInv, b )

        # Debug
        '''
        print("Matriz b")
        print( b )
        print("Matrix A")
        print( A )
        print("Transposed A")
        print( at )
        print("Multiplication between transposed A and A")
        print( mult )
        print("Inverse of the multiplication")
        print( inv )
        print("Multiplication between inv and transposed A")
        print( multInv )
        print("Multiplication between multInv and b")
        print( estimated )
        '''
        x = (estimated[0])[0]; y = (estimated[1])[0]
        return x, y


    def baricentricCoordinates( self, model, scene, parameters ):
        '''
        \fn baricentricCoordinates( model, scene, parameters )
        
        \brief Calculate the baricentric coordinates estimator
        
        \param model - Template that represent the visual stimulli
        \param scene - Scene to be foveated
        \param parameters - Parameters of fovea structure
        '''
        print("Weighted Barycentric Coordinates")
        multifovea = Multifovea( parameters.u, parameters )
        multifovea.multifoveatedFeatures( scene, params )
        multifovea.multifoveatedMatching( model, params )
        potentials = self.weightedFunctionMultifovea( multifovea, params )
        detectionRateTotal = 0.0
        for d in range( 0, len(potentials) ):
            detectionRateTotal += potentials[d]
        #print( detectionRateTotal )
        # Transforming points to cartesian domain
        pose = []
        for f in range( 0, len(parameters.f), 2 ):
            pose.append( [ parameters.f[f] + int( parameters.u[0]/2 ), parameters.f[f+1] + int( parameters.u[1]/2 ) ] )
        x = 0.0; y = 0.0;
        for f in range( 0, len(potentials) ):
            x += ( potentials[f] * pose[f][0] ) / detectionRateTotal
            y += ( potentials[f] * pose[f][1] ) / detectionRateTotal
        #print( x, y )
        return x, y

    
        
#How to instantiate and use this class
if __name__ == '__main__':
    '''
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
    
    params = Parameters('params.yaml')
    model = cv2.imread('../../box.png')
    scene = cv2.imread('../../box_in_scene.png')
    rows,cols = scene.shape[:2]
    u = [ cols, rows ]
    params.updateParameterSizeImage( u )
    #multifovea = Multifovea( u, params )
    #multifovea.multifoveatedFeatures( scene, params )
    #multifovea.multifoveatedMatching( model, params )
    statistics = Statistics()
    #print(statistics.weightedFunctionMultifovea( multifovea, params ))
    iterations = 10
    jump = 10
    config = 2
    delta = [ int(params.w[0]/2) + jump, int(params.w[1]/2) + jump ]
    size = [ params.u[0] - int(params.w[0]/2) - jump, params.u[1] - int(params.w[1]/2) - jump ]
    region = [ delta, size ]
    #print(statistics.localGradient( [ params.u[0] - int(params.w[0]/2) - jump, params.u[1] - int(params.w[1]/2) - jump ], model, scene, params, jump, config))
    #print(statistics.localGradient( [ 256, 192 ], model, scene, params, jump, config))
    #print( statistics.reduceRegionByLocalGradient( model, scene, params, region, iterations, jump, config ) )
    #print( statistics.intersectionLocalGradient( model, scene, params, jump, config ) )
    #print( statistics.maximumLikelihoodEstimator( model, scene, params, 0.3, 0 ) )
    #print( statistics.trilaterationEstimator( model, scene, params ) )
    #print( statistics.multilaterationEstimator( model, scene, params ) )
    print( statistics.baricentricCoordinates( model, scene, params ) )
