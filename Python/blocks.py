#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
\file blocks.py

\brief This file contains the prototype of shape blocks.

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

from shape import Shape # Abstract class

class Blocks( Shape ) :
    '''
    \class Blocks
    
    \brief This class implements the Blocks TAD to represent the block shape
    '''

    def __init__( self, boundingBox ):
        '''
        \fn Blocks()

        \brief Constructor default

        \param params - Parameters of fovea structure
        '''
        self.boundingBoxShape = boundingBox

    def printVertices( self ):
        '''
        \fn printVertices( self )
        
        \brief Pure virtual method caracterize this class like abstract.
        This method print all vertices of the shape.
        '''
        for vertex in self.boundingBoxShape:
            print( "(" + str( vertex[0] ) + ", " + str( vertex[1] ) + ")" ) 
            

#How to instantiate and use this class
#if __name__ == '__main__':
#    boundingBox = [ [0, 0], [3, 2] ]
#    block = Blocks( boundingBox )
#    block.printVertices()
