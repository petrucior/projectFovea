/**
 * \file graph.hpp
 *
 * \brief This file contains the prototype and implementation of graph class
 * that is used to plot graphics using gnuplot.
 *
 * \author
 * Petrucio Ricardo Tavares de Medeiros \n
 * Universidade Federal do Rio Grande do Norte \n 
 * Departamento de Computacao e Automacao Industrial \n
 * petrucior at ufrn (dot) edu (dot) br
 *
 * \version 0.1
 * \date November 2019
 *
 * This file is part of projectFovea software.
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 2 of the License, or (at your option) any later
 * version. This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details. You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <iostream> //std::cout, std::endl
#include <string>

/**
 * \defgroup ProjectFovea Project Fovea
 * @{
 */

/**
 * \class Graph
 *
 * \brief This class implements the Graph TAD
 */
class Graph{
public:
  //
  // Methods
  //
  /**
   * \fn Graph()
   *
   * \brief Constructor default of graph class
   */
  Graph();
  
  /**
   * \fn ~Graph()
   *
   * \brief Destructor default of graph class.
   */
  ~Graph();
  
  /**
   * \fn void operator () ( const string & command ) 
   *
   * \brief Defining operator () that will be used to input the command
   */
  void operator () ( const std::string & command );
  
  
protected:
  //
  // Attributes
  //
  FILE *gnuFile; ///< File pointer
  
};

#endif

/**
 * \fn Graph()
 *
 * \brief Constructor default of graph class
 */
Graph::Graph(){
  gnuFile = popen( "gnuplot -persist", "w" ); // Open gnuplot
  if ( !gnuFile )
    std::cerr << ( "Gnuplot not found!" );
}

/**
 * \fn ~Graph()
 *
 * \brief Destructor default of graph class.
 */
Graph::~Graph(){
  fprintf( gnuFile, "exit\n" );
  pclose( gnuFile );
}

/**
 * \fn void operator () ( const string & command ) 
 *
 * \brief Defining operator () that will be used to input the command
 */
void 
Graph::operator () ( const std::string & command ){
  fprintf( gnuFile, "%s\n", command.c_str() );
  fflush( gnuFile );
}

/** @} */ //end of group class.
