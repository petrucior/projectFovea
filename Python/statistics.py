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
from parameters import Parameters
from level import Level
from fovea import Fovea
from feature import Features

class Statistics :
    '''
    \class Statistics
    
    \brief This class implements the Statistics TAD
    '''
    
    def __init__( self, ):
        '''
        \fn Statistics( )
        
        \brief Constructor default
    
        \param u - Size of image
        \param parameters - Parameters of fovea structure
        '''
        # Updating the image size parameter ( U )
        parameters.updateParameterSizeImage( u )

        self.levels = [] # cleaning levels
        for k in range( 0, parameters.m + 1 ):
            level = Level( k, parameters )
            self.levels.append( level )
        
