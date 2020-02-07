#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
\file shape.py

\brief This file contains the prototype of shapes to build the 
fovea structure

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

from abc import ABC, abstractmethod # Abstract Base Classes (ABC) library

class Shape( ABC ):
    '''
    \class Shape
    
    \brief This class implements the Shape TAD to represent abstract structure of fovea
    '''
    vertices = []
    boundingBoxShape = []
    
    @abstractmethod
    def printVertices( self ):
        '''
        \fn printVertices( self )
        
        \brief Pure virtual method caracterize this class like abstract.
        This method print all vertices of the shape.
        '''
        pass
    
        
#Note that you can't to instantiate this class
#shape = Shape()
