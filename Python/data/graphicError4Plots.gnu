#
# Executing file:
# $ gnuplot -c graphicError4Plots.gnu graphError.dat test1.(gif/eps) test2.(gif/eps)
#
# GIF
#set terminal gif size 1200,600 enhanced
# EPS
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
file1 = ARG2
file2 = ARG3

###############
# First Plot  #
###############
set output file1
set multiplot layout 2,2

# Configuration 0
set title "Configuracao 0"
plot dataset using ($1):($2) with points title "MLE" ls 1, \
     dataset using ($1):($3) with points title "Multilateracao" ls 2, \
     dataset using ($1):($4) with points title "Trilateracao" ls 3, \
     dataset using ($1):($5) with points title "Coordenadas Baricentricas" ls 4

# Configuration 1
set title "Configuracao 1"
plot dataset using ($1):($6) with points title "MLE" ls 1, \
     dataset using ($1):($7) with points title "Multilateracao" ls 2, \
     dataset using ($1):($8) with points title "Trilateracao" ls 3, \
     dataset using ($1):($9) with points title "Coordenadas Baricentricas" ls 4


# Configuration 2
set title "Configuracao 2"
plot dataset using ($1):($10) with points title "MLE" ls 1, \
     dataset using ($1):($11) with points title "Multilateracao" ls 2, \
     dataset using ($1):($12) with points title "Trilateracao" ls 3, \
     dataset using ($1):($13) with points title "Coordenadas Baricentricas" ls 4


# Configuration 3
set title "Configuracao 3"
plot dataset using ($1):($14) with points title "MLE" ls 1, \
     dataset using ($1):($15) with points title "Multilateracao" ls 2, \
     dataset using ($1):($16) with points title "Trilateracao" ls 3, \
     dataset using ($1):($17) with points title "Coordenadas Baricentricas" ls 4

# Showing mean and standard deviation
#f(x) = mean_y
#fit f(x) dataset u ($1):($2) via mean_y
#stddev_y = sqrt(FIT_WSSR / (FIT_NDF + 1 ))

unset multiplot

###############
# Second Plot #
###############
set output file2
set multiplot layout 2,2

# Configurations of MLE
set title "MLE"
plot dataset using ($1):($2) with lines title "Config 0" ls 1, \
     dataset using ($1):($6) with lines title "Config 1" ls 2, \
     dataset using ($1):($10) with lines title "Config 2" ls 3, \
     dataset using ($1):($14) with lines title "Config 3" ls 4

# Configurations of Mult
set title "Multilateracao"
plot dataset using ($1):($3) with lines title "Config 0" ls 1, \
     dataset using ($1):($7) with lines title "Config 1" ls 2, \
     dataset using ($1):($11) with lines title "Config 2" ls 3, \
     dataset using ($1):($15) with lines title "Config 3" ls 4


# Configurations of Trilateracao
set title "Trilateracao"
plot dataset using ($1):($4) with lines title "Config 0" ls 1, \
     dataset using ($1):($8) with lines title "Config 1" ls 2, \
     dataset using ($1):($12) with lines title "Config 2" ls 3, \
     dataset using ($1):($16) with lines title "Config 3" ls 4


# Configurations of Coordenadas Baricentricas
set title "Coordenadas Baricentricas"
plot dataset using ($1):($5) with lines title "Config 0" ls 1, \
     dataset using ($1):($9) with lines title "Config 1" ls 2, \
     dataset using ($1):($13) with lines title "Config 2" ls 3, \
     dataset using ($1):($17) with lines title "Config 3" ls 4

# Showing mean and standard deviation
#f(x) = mean_y
#fit f(x) dataset u ($1):($2) via mean_y
#stddev_y = sqrt(FIT_WSSR / (FIT_NDF + 1 ))

unset multiplot
quit
