#
# Executing file:
# $ gnuplot -c graphicError.gnu graphError.dat test.eps
#
set term postscript landscape color solid  
set xlabel 'Imagens'
set ylabel 'Error'
set key top left
set style function linespoints
set style line 1 lw 4 lc rgb '#990042' ps 0.7 pt 1# pi 5
set style line 2 lw 3 lc rgb '#31f120' ps 0.5 pt 2# pi 3
set style line 3 lw 3 lc rgb '#0044a5' ps 0.5 pt 6# pi 5
set style line 4 lw 4 lc rgb '#888888' ps 0.3 pt 4# pi 4

dataset = ARG1
file = ARG2

# Configuration 0
set output file
#plot dataset using ($1):($2) with lines title "MLE" ls 1, \
#     dataset using ($1):($3) with lines title "Multilateracao" ls 2, \
#     dataset using ($1):($4) with lines title "Trilateracao" ls 3, \
#     dataset using ($1):($5) with lines title "Coordenadas Baricentricas" ls 4

plot dataset using ($1):($2) with points title "MLE" ls 1, \
     dataset using ($1):($3) with points title "Multilateracao" ls 2, \
     dataset using ($1):($4) with points title "Trilateracao" ls 3, \
     dataset using ($1):($5) with points title "Coordenadas Baricentricas" ls 4




