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
                    fovea.foveatedFeatures( scene, parameters )
                    fovea.foveatedMatching( model, parameters )
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
                fovea.foveatedFeatures( scene, parameters )
                fovea.foveatedMatching( model, parameters )
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
        # Adding all weight
        sumWeight = 0.0
        for k in range( 0, p.m + 1 ):
            sumWeight += self.weightedFunction[k][0] + self.weightedFunction[k][1]
        # Weighting fovea
        self.fpdf = 0.0
        for k in range( 0, p.m + 1 ):
            inlierRate = ( fovea.features.inlierRateSURF[k] )
#            print("taxa de inlier: " +str(inlierRate)+", pesox: "+str(self.weightedFunction[k][0])+", pesoy: "+str(self.weightedFunction[k][1]))
            #self.fpdf += 1.0 * inlierRate
            self.fpdf += (self.weightedFunction[k][0]/sumWeight) * inlierRate + (self.weightedFunction[k][1]/sumWeight) * inlierRate
#        print("fpdf"+str(self.fpdf))
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
            x = regionUnderAnalysis[0][0]; y = regionUnderAnalysis[0][1];
            sx = regionUnderAnalysis[1][0]; sy = regionUnderAnalysis[1][1];
            center = ( (x + sx)/2, (y + sy)/2 )
            return center
        else:
            xmin = regionUnderAnalysis[0][0]; ymin = regionUnderAnalysis[0][1];
            xmax = regionUnderAnalysis[1][0]; ymax = regionUnderAnalysis[1][1];
            # It will work until regionUnderAnalysis is bigger than 25 x 25
            if ( (xmax - xmin) > int(parameters.w[0]/2) and (ymax - ymin) > int(parameters.w[1]/2) ):
                # Chosen point
                x = random.randint(xmin, xmax)
                y = random.randint(ymin, ymax)
                
                '''
                cv2.rectangle(scene,(xmin, ymin),(xmax, ymax),(0,255,0),3)
                cv2.rectangle(scene,(0, 30), (20, 60 ), (100, 100, 0), 2 )
                cv2.circle(scene, (x, y), 3, (0, 0, 255), -1)
                cv2.imshow('region', scene)
                cv2.waitKey( 0 )
                '''
                
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

                    '''
                    cv2.rectangle(scene,(xmin, ymin),(xmax, ymax),(244,255,0),3)
                    cv2.imshow('region', scene)
                    cv2.waitKey( 0 )
                    '''
            regionUnderAnalysis = [ [xmin, ymin], [xmax, ymax] ]
            return self.reduceRegionByLocalGradient( model, scene, parameters, regionUnderAnalysis, iterations - 1, jump, config )


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
            return self.intersectionLocalGradient( model, scene, parameters, jump, config )
        '''
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
        '''
        
    
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
        #print("MLE")
        multifovea = Multifovea( parameters.u, parameters )
        multifovea.multifoveatedFeatures( scene, parameters )
        multifovea.multifoveatedMatching( model, parameters )
        potentials = self.weightedFunctionMultifovea( multifovea, parameters )
        # Transforming points to cartesian domain
        pose = []
        for f in range( 0, len(parameters.f), 2 ):
            pose.append( [ parameters.f[f] + int( parameters.u[0]/2 ), parameters.f[f+1] + int( parameters.u[1]/2 ) ] )
        x = 0.0; y = 0.0;
        if ( method == 0 ): # Arithmetic mean
            counter = 0
            for i in range( 0, len( potentials ) ):
                if ( potentials[i] > threshold ):
                    x += pose[i][0]
                    y += pose[i][1]
                    counter += 1
            if ( counter == 0 ):
                x = 0.0; y = 0.0;
            else:
                x /= counter
                y /= counter
        
        if ( method == 1 ): # Weighted average
            sumPotentials = 0.0
            for i in range( 0, len( potentials ) ):
                x += pose[i][0] * potentials[i]
                y += pose[i][1] * potentials[i]
                sumPotentials += potentials[i]
            x /= sumPotentials
            y /= sumPotentials
        
        return x, y


    def distanceEstimator( self, points, potential, parameters ):
        '''
        \fn distanceEstimator( points, potential )

        \param points - Foveas position
        \param potential - Detection rate

        \return A vector of estimated distances
        '''
        # Intersection between secant circumferences
        components = numpy.zeros(int(len(points)/2))
        radius = [];
        for i in range( len(potential) ):
            radius.append(1 - potential[i])
        found = False
        while ( found == False ):
            found = True
            for i in range( len(potential) ):
                for j in range( i, len(potential) ):
                    if ( i != j ):

                        d = math.sqrt( pow(points[2*i] - points[2*j], 2.0) + pow(points[2*i+1] - points[2*j+1], 2.0) )
                        if ( ( d > radius[i] - radius[j] ) and ( d < radius[i] + radius[j] ) ):
                            found = found and True
                        else:
                            found = False
                            radius[i] += ( 1 - potential[i] )
                            radius[j] += ( 1 - potential[j] )
        #print( radius )
        components = radius
        return components
    
    
    def trilaterationEstimator( self, model, scene, parameters ):
        '''
        \fn trilaterationEstimator( model, scene, parameters )
        
        \brief Calculate the trilateration Estimator
        
        \param model - Template that represent the visual stimulli
        \param scene - Scene to be foveated
        \param parameters - Parameters of fovea structure
        '''
        #print("Trilateration")
        multifovea = Multifovea( parameters.u, parameters )
        multifovea.multifoveatedFeatures( scene, parameters )
        multifovea.multifoveatedMatching( model, parameters )
        potential = self.weightedFunctionMultifovea( multifovea, parameters )
        
        # Points
        points = []
        for p in range( 0, len(parameters.f) ):
            points.append( parameters.f[p] )

        # Transforming points to cartesian domain
        for p in range( 0, len(points) ):
            if ( p % 2 == 0 ):
                points[p] += int( parameters.u[0]/2 )
            else:
                points[p] += int( parameters.u[1]/2 )        
        #print( points )

        # Distance Estimator
        r = self.distanceEstimator( points, potential, parameters )
        
        '''
        output = multifovea.multifoveatedImage( scene, parameters )
        for i in range( 0, len(potential) ):
            cv2.circle( output, (int(points[2*i]), int(points[2*i+1])), int(r[i]), (parameters.colors[3*i], parameters.colors[3*i+1], parameters.colors[3*i+2]), 1 )
        cv2.imshow( "trilateration", output )
        cv2.waitKey( 0 )
        '''

        x1 = points[0]; y1 = points[1];
        x2 = points[2]; y2 = points[3];
        x3 = points[4]; y3 = points[5];

        r1 = r[0]; r2 = r[1]; r3 = r[2];

        a = (-2 * x1) + (2 * x2)
        b = (-2 * y1) + (2 * y2)
        c = (r1*r1) - (r2*r2) - (x1*x1) + (x2*x2) - (y1*y1) + (y2*y2)
        d = (-2 * x2) + (2 * x3)
        e = (-2 * y2) + (2 * y3)
        f = (r2*r2) - (r3*r3) - (x2*x2) + (x3*x3) - (y2*y2) + (y3*y3)

        den = (a * e) - (b * d)
        if ( den != 0 ):
            x = ((c * e) - (b * f)) / den
            y = ((a * f) - (c * d)) / den
        else:
            x = 0
            y = 0
        
        return x, y    
    
    def calcArea( self, points ):
        '''
        \fn calcArea( points )
        
        \brief Calculates the triangle area
        
        \param points - The triangle vertex points [ p0x, p0y, p1x, p1y, p2x, p2y ]
        '''
        a = 0; b = 1; c = 2;
        return ( ( (points[2*b] - points[2*a])*(points[(2*c)+1] - points[(2*a)+1]) ) - ( (points[2*c] - points[2*a])*(points[(2*b)+1] - points[(2*a)+1]) ) ) 

    
    def multilaterationEstimator( self, model, scene, parameters ):
        '''
        \fn multilaterationEstimator( model, scene, parameters )
        
        \brief Calculate the multilateration Estimator
        
        \param model - Template that represent the visual stimulli
        \param scene - Scene to be foveated
        \param parameters - Parameters of fovea structure
        '''
        #print("Multilateration")
        multifovea = Multifovea( parameters.u, parameters )
        multifovea.multifoveatedFeatures( scene, parameters )
        multifovea.multifoveatedMatching( model, parameters )
        potentials = self.weightedFunctionMultifovea( multifovea, parameters )
        
        # Transforming points to cartesian domain
        pose = []
        for f in range( 0, len(parameters.f), 2 ):
            pose.append( [ parameters.f[f] + int( parameters.u[0]/2 ), parameters.f[f+1] + int( parameters.u[1]/2 ) ] )
        #print( pose )

        # Points
        points = []
        for p in range( 0, len(parameters.f) ):
            points.append( parameters.f[p] )
        
        # Transforming points to cartesian domain
        for p in range( 0, len(points) ):
            if ( p % 2 == 0 ):
                points[p] += int( parameters.u[0]/2 )
            else:
                points[p] += int( parameters.u[1]/2 )
        
        # Distance Estimator
        r = self.distanceEstimator( points, potentials, parameters )

        '''
        output = multifovea.multifoveatedImage( scene, parameters )
        for p in range( 0, len(potentials) ):
            cv2.circle( output, (int(pose[p][0]), int(pose[p][1])), int(r[p]), (parameters.colors[3*p], parameters.colors[3*p+1], parameters.colors[3*p+2]), 1 )
        cv2.imshow( "testando", output )
        cv2.waitKey( 0 )
        '''
        
        x = 0.0; y = 0.0;
        if ( len(parameters.f) == 0 ): return x, y
        m = int(len(parameters.f)/2)
        b = numpy.zeros( (m - 1, 1), dtype=numpy.float64 )
        A = numpy.zeros( (m - 1, 2), dtype=numpy.float64 )
        #print( parameters.f )
        for i in range( 0, m - 1 ):
            A[i][0]  = ( -2 * pose[i][0] ) + ( 2 * pose[i+1][0] )            
            A[i][1]  = ( -2 * pose[i][1] ) + ( 2 * pose[i+1][1] )
            b[i][0] = pow( r[i], 2.0 ) - pow( r[i+1], 2.0 ) - pow( pose[i][0], 2.0 ) + pow( pose[i+1][0], 2.0 ) - pow( pose[i][1], 2.0 ) + pow( pose[i+1][1], 2.0 )
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

    
    def barycentricCoordinates( self, model, scene, parameters ):
        '''
        \fn barycentricCoordinates( model, scene, parameters )
        
        \brief Calculate the barycentric coordinates estimator
        
        \param model - Template that represent the visual stimulli
        \param scene - Scene to be foveated
        \param parameters - Parameters of fovea structure
        '''
        #print("Weighted Barycentric Coordinates")
        multifovea = Multifovea( parameters.u, parameters )
        multifovea.multifoveatedFeatures( scene, parameters )
        multifovea.multifoveatedMatching( model, parameters )
        potentials = self.weightedFunctionMultifovea( multifovea, parameters )
        # reoorganizing the potentials
        for d in range( 0, len(potentials) ):
            potentials[d] = 1.0 - (0.57 - potentials[d])
        #print( potentials )
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
    

    def displayConfigurations( self, model, scene, parameters, poseEstimated ):
        '''
        \fn displayConfigurations( model, scene, parameters, poseEstimated )
        
        \brief Displaying the configurations
        
        \param model - Template that represent the visual stimulli
        \param scene - Scene to be foveated
        \param parameters - Parameters of fovea structure
        \param poseEstimated - Pose of localization estimated
        '''
        multifovea = Multifovea( parameters.u, parameters )
        multifovea.multifoveatedFeatures( scene, parameters )
        multifovea.multifoveatedMatching( model, parameters )
        fpdf = self.weightedFunctionMultifovea( multifovea, parameters )
        #print( fpdf )
        output = multifovea.multifoveatedImage( scene, parameters )
        for i in range( 0, int(len( parameters.f )/2) ):
            cv2.putText( output, "f"+str(i)+": "+str(round(fpdf[i], 2)), ( int(parameters.f[i*2] + parameters.u[0]/2 - 30), int(parameters.f[i*2 + 1] + parameters.u[1]/2) ), cv2.FONT_HERSHEY_PLAIN, 1.0, [parameters.colors[i*3], parameters.colors[i*3+1], parameters.colors[i*3 + 2]], 2 )
        if (( poseEstimated[0] != -1 ) and ( poseEstimated[1] != -1 ) ):
            # Radius of circle 
            radius = 5
            # Blue color in BGR 
            color = (0, 0, 255)
            # Line thickness of 1 px 
            thickness = 2
            # put the point estimated
            cv2.circle( output, (int(poseEstimated[0]), int(poseEstimated[1])), radius, color, thickness )
        winname = 'Estimation multifoveated'
        cv2.namedWindow( winname, cv2.WINDOW_NORMAL )
        cv2.imshow( winname, output )
        cv2.waitKey( 1 )
        #cv2.destroyAllWindows()

    
    def findMaximum(  self, model, scene, parameters, array ):
        '''
        \fn findMaximum( array )
        
        \brief Finding the higher fpdf value of foveas
        
        \param model - Template that represent the visual stimulli
        \param scene - Scene to be foveated
        \param parameters - Parameters of fovea structure
        \param array - Contains in each position 3 foveas that are distributed in a triangle
        
        \return Position of maximum configuration
        '''
        params = copy.copy(parameters)
        multifovea = Multifovea( parameters.u, params )
        multifovea.multifoveatedFeatures( scene, params )
        multifovea.multifoveatedMatching( model, params )
        fpdf = self.weightedFunctionMultifovea( multifovea, params )
        #print( fpdf )
        indexSave = 0; biggerSum = 0.0;
        for index in range( 0, len( array ) ):
            sumAux = 0.0;
            for indexAux in range( 0, len( array[index] ) ):
                sumAux += fpdf[ array[index][indexAux] ]
            if ( biggerSum < sumAux ):
                biggerSum = sumAux
                indexSave = index
        #print( array[indexSave] )
        # Using only 3 foveas
        foveas = []; colors = [];
        for i in range( 0, len( array[indexSave] ) ):
            foveas.append( params.f[ 2 * array[indexSave][i] ] )
            foveas.append( params.f[ 2 * array[indexSave][i] + 1 ] )
            colors.append( params.colors[ 3 * array[indexSave][i] ] )
            colors.append( params.colors[ 3 * array[indexSave][i] + 1 ] )
            colors.append( params.colors[ 3 * array[indexSave][i] + 2 ] )
        #print( foveas )
        params.updateParameterFoveas( foveas, colors )
        return params
                
    def poseFoveas(  self, parameters, config ):
        '''
        \fn poseFoveas( parameters, config )
        
        \brief Defining the position of foveas based on defined configuration
        
        \param parameters - Parameters of fovea structure
        \param config - Configuration (0) with 4 foveas, (1) 5 foveas, (2) 6 foveas and (3) 9 foveas
        '''
        scale = [ [4, 4],
                  [4, 4],
                  [6, 4],
                  [6, 6] ]
        jump = [ [2, 2],
                 [2, 1],
                 [2, 2],
                 [2, 2] ]
        distpointsx = int(parameters.u[0]/scale[config][0]); distpointsy = int(parameters.u[1]/scale[config][1]);
        foveas = []
        for j in range( distpointsy, parameters.u[1], jump[config][1] * distpointsy ):
            for i in range( distpointsx, parameters.u[0], jump[config][0] * distpointsx ):
                if ( ( config == 1 ) and ( j == parameters.u[1]/2 ) ):
                    i = int(parameters.u[0]/2)
                foveas.append( (i - math.floor( parameters.u[0]/2 )) )
                foveas.append( (j - math.floor( parameters.u[1]/2 )) )
            if ( ( config == 1 ) and ( j == parameters.u[1]/2 ) ):
                foveas.pop(); foveas.pop(); # removing x and y
               
        #print( foveas )
        colors = []
        for i in range( 0, int(len( foveas )/2) ):
            colors.append( random.randint(0, 255) ) # blue
            colors.append( random.randint(0, 255) ) # green
            colors.append( random.randint(0, 255) ) # red
        parameters.updateParameterFoveas( foveas, colors )
        
        return parameters
        
        
    
    
#How to instantiate and use this class
#if __name__ == '__main__':
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
    print( statistics.barycentricCoordinates( model, scene, params ) )
    '''

    '''
    # First code
    # 0. Display Model and Scene
    # 1. Display Multirresolution Domain
    # 2. Features Detection
    # 3. Detecting Problems In Step 2
    # 4. Minimum Fovea Size To Solve Step 2
    
    # Getting parameters, model and scene
    params = Parameters('params.yaml')
    model = cv2.imread('../../box.png')
    scene = cv2.imread('../../box_in_scene.png')
    # Updateing size of scene
    rows,cols = scene.shape[:2]
    u = [ cols, rows ]
    params.updateParameterSizeImage( u )
    # Creating fovea
    fovea = Fovea( u, params )
    
    # Display Model and Scene
    cv2.imshow( "Model", model )
    cv2.imshow( "Scene", scene )
    cv2.waitKey( 0 )

    # Display Multirresolution Domain
    output = fovea.foveatedImage( scene, params )
    cv2.imshow( "fovea", output ) 
    cv2.waitKey( 0 )
    fovea.saveLevels( "nivel", scene, params )
    cv2.waitKey( 0 )

    # Features Detection
    fovea.foveatedFeatures( scene, params )
    
    
    # Second Code
    # 1. Proportion of Inliers by Level
    # 2. Detection Probability Density Function
    
    # Features Matches
    fovea.features.matchingFeatures( model, params )

    # Starting Statistics
    statistics = Statistics()

    # Proportion of Inliers by Level
    print("Ploting proportion")
    statistics.plotProportion( model, scene, fovea, params )

    # Detection Probability Density Function
    print("Building FPDF")
    statistics.plotFunction( model, scene, fovea, params )

    
    # Destroy all windows
    cv2.destroyAllWindows()
    '''
    
    '''
    # Third Code
    # 1. Calculating and Displaying the Results of reduceRegion and intersectionPotential
    
    # Getting parameters, model and scene
    params = Parameters('params.yaml')
    model = cv2.imread('../../box.png')
    scene = cv2.imread('../../box_in_scene.png')
    # Updateing size of scene
    rows,cols = scene.shape[:2]
    u = [ cols, rows ]
    params.updateParameterSizeImage( u )

    # Starting Statistics
    statistics = Statistics()
    
    # Configurations
    iterations = 10
    jump = 15
    config = 2
    delta = [ int(params.w[0]/2) + jump, int(params.w[1]/2) + jump ]
    size = [ params.u[0] - int(params.w[0]/2) - jump, params.u[1] - int(params.w[1]/2) - jump ]
    region = [ delta, size ]

    #print( statistics.reduceRegionByLocalGradient( model, scene, params, region, iterations, jump, config ) )
    #print( statistics.intersectionLocalGradient( model, scene, params, jump, config ) )
    '''

    '''
    # Fourth Code
    # 1. Displaying the 4 configurations
    # 2. Calculating and Displaying the MLE
    # 3. Calculating and Displaying the trilateration and multilateration
    # 4. Calculating and Displaying the Barycentric coordinates

    # Getting parameters, model and scene
    params = Parameters('params.yaml')
    model = cv2.imread('../../box.png')
    scene = cv2.imread('../../box_in_scene.png')
    #cv2.imshow( "Model", model )
    #cv2.imshow( "Scene", scene )
    #cv2.waitKey( 0 )

    # Updating size of scene
    rows,cols = scene.shape[:2]
    u = [ cols, rows ]
    params.updateParameterSizeImage( u )

    # Starting Statistics
    statistics = Statistics()
    
    # ---------------------------
    # Configuracao com 4 foveas
    # ---------------------------
    params = statistics.poseFoveas( params, 0 )
    #foveas = [  -128.0, -96.0, 128.0, -96.0, 128.0, 96.0, -128.0, 96.0 ]
    #colors = []
    #for i in range( 0, int(len( foveas )/2) ):
    #    colors.append( random.randint(0, 255) ) # blue
    #    colors.append( random.randint(0, 255) ) # green
    #    colors.append( random.randint(0, 255) ) # red
    #params.updateParameterFoveas( foveas, colors )
    
    #pose = statistics.maximumLikelihoodEstimator( model, scene, params, 0.1, 0 )
    #pose = statistics.maximumLikelihoodEstimator( model, scene, params, 0.12, 1 )
    parameters = statistics.findMaximum( model, scene, params, [ [0, 1, 3], [1, 2, 3], [0, 2, 3], [0, 1, 2] ] )
    pose = statistics.trilaterationEstimator( model, scene, parameters )
    #pose = statistics.multilaterationEstimator( model, scene, params )
    #parameters = statistics.findMaximum( model, scene, params, [ [0, 1, 3], [1, 2, 3], [0, 2, 3], [0, 1, 2] ] )
    #pose = statistics.barycentricCoordinates( model, scene, parameters )
    statistics.displayConfigurations( model, scene, params, pose )

    # ---------------------------
    # Configuracao com 5 foveas
    # ---------------------------
    params = statistics.poseFoveas( params, 1 )
    #foveas = [  -128.0, -96.0, 128.0, -96.0, 128.0, 96.0, -128.0, 96.0, 0.0, 0.0 ]
    #colors = []
    #for i in range( 0, int(len( foveas )/2) ):
    #    colors.append( random.randint(0, 255) ) # blue
    #    colors.append( random.randint(0, 255) ) # green
    #    colors.append( random.randint(0, 255) ) # red
    #params.updateParameterFoveas( foveas, colors )
    
    #pose = statistics.maximumLikelihoodEstimator( model, scene, params, 0.1, 0 )
    #pose = statistics.maximumLikelihoodEstimator( model, scene, params, 0.12, 1 )
    parameters = statistics.findMaximum( model, scene, params, [ [0, 1, 4], [1, 2, 4], [2, 3, 4], [0, 3, 4] ] )
    pose = statistics.trilaterationEstimator( model, scene, parameters )
    #pose = statistics.multilaterationEstimator( model, scene, params )
    #parameters = statistics.findMaximum( model, scene, params, [ [0, 1, 4], [1, 2, 4], [2, 3, 4], [0, 3, 4] ] )
    #pose = statistics.barycentricCoordinates( model, scene, parameters )
    statistics.displayConfigurations( model, scene, params, pose )
    
    # ---------------------------
    # Configuracao com 6 foveas
    # ---------------------------
    params = statistics.poseFoveas( params, 2 )
    #foveas = [ -170.0, -96.0, 0.0, -96.0, 170.0, -96.0, 170.0, 96.0, 0.0, 96.0, -170.0, 96.0 ]
    #colors = []
    #for i in range( 0, int(len( foveas )/2) ):
    #    colors.append( random.randint(0, 255) ) # blue
    #    colors.append( random.randint(0, 255) ) # green
    #    colors.append( random.randint(0, 255) ) # red
    #params.updateParameterFoveas( foveas, colors )
    
    #pose = statistics.maximumLikelihoodEstimator( model, scene, params, 0.1, 0 )
    #pose = statistics.maximumLikelihoodEstimator( model, scene, params, 0.12, 1 )
    parameters = statistics.findMaximum( model, scene, params, [ [0, 1, 4], [0, 4, 5], [1, 4, 5], [0, 1, 5], [1, 2, 3], [1, 3, 4], [1, 2, 4], [2, 3, 4] ] )
    pose = statistics.trilaterationEstimator( model, scene, parameters )
    #pose = statistics.multilaterationEstimator( model, scene, params )
    #parameters = statistics.findMaximum( model, scene, params, [ [0, 1, 4], [0, 4, 5], [1, 4, 5], [0, 1, 5], [1, 2, 3], [1, 3, 4], [1, 2, 4], [2, 3, 4] ] )
    #pose = statistics.barycentricCoordinates( model, scene, parameters )
    statistics.displayConfigurations( model, scene, params, pose )
    
    # ---------------------------
    # Configuracao com 9 foveas
    # ---------------------------
    params = statistics.poseFoveas( params, 3 )
    #foveas = [ -170.0, -128.0, 0.0, -128.0, 170.0, -128.0, 170.0, 0.0, 0.0, 0.0, -170.0, 0.0, -170.0, 128.0, 0.0, 128.0, 170.0, 128.0 ]
    #colors = []
    #for i in range( 0, int(len( foveas )/2) ):
    #    colors.append( random.randint(0, 255) ) # blue
    #    colors.append( random.randint(0, 255) ) # green
    #    colors.append( random.randint(0, 255) ) # red
    #params.updateParameterFoveas( foveas, colors )
    
    #pose = statistics.maximumLikelihoodEstimator( model, scene, params, 0.1, 0 )
    #pose = statistics.maximumLikelihoodEstimator( model, scene, params, 0.12, 1 )
    parameters = statistics.findMaximum( model, scene, params, [ [0, 1, 4], [0, 4, 5], [1, 4, 5], [0, 1, 5], [1, 2, 3], [1, 3, 4], [1, 2, 4], [2, 3, 4], [5, 4, 7], [5, 6, 7], [4, 6, 7], [4, 5, 6], [3, 4, 8], [4, 7, 8], [3, 4, 7], [3, 7, 8] ] )
    pose = statistics.trilaterationEstimator( model, scene, parameters )
    #pose = statistics.multilaterationEstimator( model, scene, params )
    #parameters = statistics.findMaximum( model, scene, params, [ [0, 1, 4], [0, 4, 5], [1, 4, 5], [0, 1, 5], [1, 2, 3], [1, 3, 4], [1, 2, 4], [2, 3, 4], [5, 4, 7], [5, 6, 7], [4, 6, 7], [4, 5, 6], [3, 4, 8], [4, 7, 8], [3, 4, 7], [3, 7, 8] ] )
    #pose = statistics.barycentricCoordinates( model, scene, parameters )
    statistics.displayConfigurations( model, scene, params, pose )
    
#    foveas = [ -60, +40 ]
#    #foveas = [ -128, +96 ]
#    #foveas = [0, 0] # 22
#    colors = [ 0, 0, 250 ]
#    params.updateParameterFoveas( foveas, colors )
#    pose = [-1, -1]
#    statistics.displayConfigurations( model, scene, params, pose )

    cv2.destroyAllWindows()
    '''
    
