#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
\file mean_sd.py

\brief This file calculates mean and standard desviation of the .dat file

\author
Petrucio Ricardo Tavares de Medeiros \n 
Universidade Federal do Rio Grande do Norte \n
Departamento de Computacao e Automacao Industrial \n
petrucior at ufrn (dot) edu (dot) br

\version 0.1
\date Aug 2020

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

import sys
from io import StringIO
import numpy

def main():
    arguments = sys.argv[1:]
    if ( len( arguments ) != 1 ):
        print ('python mean_sd.py ''< data_generated (.dat) >'' ')
        print ('Example: python main.py dataset.dat')
    else:
        data = numpy.genfromtxt( arguments[0],
                                 names=True,
                                 dtype=None,
                                 delimiter='\t')
        #print( data )
        matrix = []
        for l in range( len( data ) ):
            line = list( data[l] )
            matrix.append( line )

        print( "Mean and Standard Desviation For First Approach" )
        print( numpy.mean( matrix, axis = 0 ) )
        print( numpy.std( matrix, axis = 0 ) )
        
        #print("\nmean of matrix, axis = 0 : ", numpy.mean(matrix, axis = 0)) 
        

if __name__ == '__main__':
    main()
