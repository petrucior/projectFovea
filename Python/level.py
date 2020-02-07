#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
\file level.py

\brief This file contains the prototype of levels.

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

from shape import Shape
from blocks import Blocks
from parameters import Parameters

class Level :
    '''
    \class Level
    
    \brief This class implements the Level TAD to represent levels
    '''
    
    def __init__( self, k, parameters ):
        '''
        \fn Level(u, parameters)
        
        \brief Constructor default
    
        \param k - Level of fovea
        \param parameters - Parameters of fovea structure
        '''
        # Updating the image size parameter ( U ) was done in configuration of fovea         
        
        if ( parameters.typeShape == 0 ): # Blocks
            self.boundingBox = [] # Clear boundingBox
            self.shapeLevel = [] # Clear shapeLevel
            self.boundingBox.append( self.getDelta( k, parameters ) )
            self.boundingBox.append( self.getSize( k, parameters ) )
            block = Blocks( self.boundingBox )
            block.printVertices()
            self.shapeLevel.append( block )
        else: # Polygons
            print( "Shape wasn't configured" )


    #def __init__( self, k, u, parameters ):
    #    '''
    #    \fn Level(u, parameters)
    #    
    #    \brief Constructor default
    #
    #    \param k - Level of fovea
    #    \param u - Size of image
    #    \param parameters - Parameters of fovea structure
    #    '''
    #    # Updating the image size parameter ( U )
    #    parameters.updateParameterSizeImage( u )
    #    
    #    if ( parameters.typeShape == 0 ): # Blocks
    #        self.boundingBox = [] # Clear boundingBox
    #        self.shapeLevel = [] # Clear shapeLevel
    #        self.boundingBox.append( self.getDelta( k, parameters ) )
    #        self.boundingBox.append( self.getSize( k, parameters ) )
    #        block = Blocks( self.boundingBox )
    #        block.printVertices()
    #        self.shapeLevel.append( block )
    #    else: # Polygons
    #        print( "Shape wasn't configured" )

        

    def getDelta( self, k, parameters ):
        '''
        \fn getDelta( k, parameters )

        \brief Calculates the initial pixel to build the boundingBox

        \param k - Level of fovea
        \param parameters - Parameters of fovea structure
        '''
        p = parameters
        dx = int( k * ( p.u[0] - p.w[0] + ( 2 * p.f[0] ) ) )/ ( 2 * p.m );
        dy = int( k * ( p.u[1] - p.w[1] + ( 2 * p.f[1] ) ) )/ ( 2 * p.m );
        return [ dx, dy ]

    def getSize( self, k, parameters ):
        '''
        \fn getDelta( k, parameters )

        \brief Calculates the initial pixel to build the boundingBox

        \param k - Level of fovea
        \param parameters - Parameters of fovea structure
        '''
        p = parameters
        sx = int( ((p.m * p.u[0]) + (p.w[0] * k) - (k * p.u[0])) / p.m );
        sy = int( ((p.m * p.u[1]) + (p.w[1] * k) - (k * p.u[1])) / p.m );
        return [ sx, sy ]

#How to instantiate and use this class
if __name__ == '__main__':
    params = Parameters('params.yaml')
    level = Level( 0, [ 250, 250 ], params )
