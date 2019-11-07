#include <iostream>
#include <fstream>
#include <stdio.h>
#include "../graph.hpp"

void menu(){
  std::cout << "Information about this code: " << std::endl;
  std::cout << "compile: g++ testegraph.cpp -o <exec>" << std::endl;
  std::cout << "run: ./<exec> /dir/<name_file> <option_plot" << std::endl;
  std::cout << "<name_file> : file .dat without extension" << std::endl;
  std::cout << "<option_plot> : 2 ( 2D plot ) or 3 ( 3D plot )" << std::endl;
}

int main( int argc, char** argv ){
  menu();
  Graph g;
  g("set term postscript eps");
  const std::string file = argv[1]; // File should be passed like: /dir/<name_file>
  g("set output \""+ file +".eps\" ");
  int dimension = atoi(argv[2]);
  if ( dimension == 2 )
    g("plot \'"+ file +".dat\' u 1:2 w l ");
  else{
    // http://lowrank.net/gnuplot/plot3d2-e.html
    g("set dgrid3d 30,30");
    g("set hidden3d");
    g("splot \'"+ file +".dat\' u 1:2:3 with lines ");
  }

  return 0;
}
