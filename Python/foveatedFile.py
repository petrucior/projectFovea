#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
\file foveatedFile.py

\brief This file executable creates foveated files.

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
import sys
import os
import numpy as np
from parameters import Parameters
from fovea import Fovea

def main():
    arguments = sys.argv[1:]
    print ('python foveatedFile.py ''dir'' ''file <yaml>'' ')
    print ('Example: python foveatedFile.py ~/Download/flower_photos fovea.yaml')

    path = arguments[0]
    params = Parameters( arguments[1] )
    cont = 1
    for _, _, files in os.walk( path ):
        for file in files:
            print( 'processing '+ str( cont ) +' de '+ str(len(files)) )
            img = cv2.imread( path + file )
            rows, cols = img.shape[:2]
            u = [ cols, rows ]
            fovea = Fovea( u, params )
            fovea.saveLevels( file, img, params )
            cont += 1



        
if __name__ == '__main__':
    main()
