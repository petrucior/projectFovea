#
# Executing file:
# $ gnuplot -c graphicError4Plots.gnu graphError1.dat graphError2.dat test1A.(gif/eps) test2A.(gif/eps) test1B.(gif/eps) test2B.(gif/eps)
# graphError1.dat -> RRLG and ILG
# graphError2.dat -> MLE, Mult, Tril and CP
#
# GIF
#set terminal gif size 1200,600 enhanced
# EPS
set term postscript landscape color solid  
set xlabel 'Imagens'
set ylabel 'Erro (pixels)'
set key top left
set style function linespoints
set style line 1 lw 4 lc rgb '#990042' ps 0.7 pt 1# pi 5
set style line 2 lw 3 lc rgb '#31f120' ps 0.5 pt 2# pi 3
set style line 3 lw 3 lc rgb '#0044a5' ps 0.5 pt 6# pi 5
set style line 4 lw 4 lc rgb '#888888' ps 0.3 pt 4# pi 4

# Setting without legends
#set nokey

datasetA = ARG1
datasetB = ARG2
file1A = ARG3
file2A = ARG4
file1B = ARG5
file2B = ARG6

################
# First Plot A #
################
set output file1A
set multiplot layout 2,2

# Abordagem 0
set title "Abordagem 0"
plot datasetA using ($1):($2) with points title "Reducao de Regiao" ls 1 #, \
     #datasetA using ($1):($3) with points title "Interseccao de Potenciais" ls 2

# Abordagem 1
set title "Abordagem 1"
plot datasetA using ($1):($4) with points title "Reducao de Regiao" ls 1, \
     datasetA using ($1):($5) with points title "Intersecao de Potenciais" ls 2


# Abordagem 2
set title "Abordagem 2"
plot datasetA using ($1):($6) with points title "Reducao de Regiao" ls 1, \
     datasetA using ($1):($7) with points title "Intersecao de Potenciais" ls 2


# Showing mean and standard deviation
#f(x) = mean_y
#fit f(x) dataset u ($1):($2) via mean_y
#stddev_y = sqrt(FIT_WSSR / (FIT_NDF + 1 ))

unset multiplot

#################
# Second Plot A #
#################
set output file2A
set multiplot layout 2,2

# Abordagens of Reduction of Region
set title "Reducao de Regiao"
plot datasetA using ($1):($2) with lines title "Abordagem 0" ls 1, \
     datasetA using ($1):($4) with lines title "Abordagem 1" ls 2, \
     datasetA using ($1):($6) with lines title "Abordagem 2" ls 3

# Configurations of Intersection of Potentials
set title "Intersecao de Potenciais"
plot datasetA using ($1):($5) with lines title "Abordagem 1" ls 2, \
     datasetA using ($1):($7) with lines title "Abordagem 2" ls 3


# Showing mean and standard deviation
#f(x) = mean_y
#fit f(x) dataset u ($1):($2) via mean_y
#stddev_y = sqrt(FIT_WSSR / (FIT_NDF + 1 ))

unset multiplot


################
# First Plot B #
################
set output file1B
set multiplot layout 2,2

# Configuration 0
set title "Configuracao 0"
plot datasetB using ($1):($2) with points title "MLE" ls 1, \
     datasetB using ($1):($3) with points title "Multilateracao" ls 2, \
     datasetB using ($1):($4) with points title "Trilateracao" ls 3, \
     datasetB using ($1):($5) with points title "Coordenadas Baricentricas" ls 4

# Configuration 1
set title "Configuracao 1"
plot datasetB using ($1):($6) with points title "MLE" ls 1, \
     datasetB using ($1):($7) with points title "Multilateracao" ls 2, \
     datasetB using ($1):($8) with points title "Trilateracao" ls 3, \
     datasetB using ($1):($9) with points title "Coordenadas Baricentricas" ls 4


# Configuration 2
set title "Configuracao 2"
plot datasetB using ($1):($10) with points title "MLE" ls 1, \
     datasetB using ($1):($11) with points title "Multilateracao" ls 2, \
     datasetB using ($1):($12) with points title "Trilateracao" ls 3, \
     datasetB using ($1):($13) with points title "Coordenadas Baricentricas" ls 4


# Configuration 3
set title "Configuracao 3"
plot datasetB using ($1):($14) with points title "MLE" ls 1, \
     datasetB using ($1):($15) with points title "Multilateracao" ls 2, \
     datasetB using ($1):($16) with points title "Trilateracao" ls 3, \
     datasetB using ($1):($17) with points title "Coordenadas Baricentricas" ls 4

# Showing mean and standard deviation
#f(x) = mean_y
#fit f(x) dataset u ($1):($2) via mean_y
#stddev_y = sqrt(FIT_WSSR / (FIT_NDF + 1 ))

unset multiplot

#################
# Second Plot B #
#################
set output file2B
set multiplot layout 2,2

# Configurations of MLE
set title "MLE"
plot datasetB using ($1):($2) with lines title "Config 0" ls 1, \
     datasetB using ($1):($6) with lines title "Config 1" ls 2, \
     datasetB using ($1):($10) with lines title "Config 2" ls 3, \
     datasetB using ($1):($14) with lines title "Config 3" ls 4

# Configurations of Mult
set title "Multilateracao"
plot datasetB using ($1):($3) with lines title "Config 0" ls 1, \
     datasetB using ($1):($7) with lines title "Config 1" ls 2, \
     datasetB using ($1):($11) with lines title "Config 2" ls 3, \
     datasetB using ($1):($15) with lines title "Config 3" ls 4


# Configurations of Trilateracao
set title "Trilateracao"
plot datasetB using ($1):($4) with lines title "Config 0" ls 1, \
     datasetB using ($1):($8) with lines title "Config 1" ls 2, \
     datasetB using ($1):($12) with lines title "Config 2" ls 3, \
     datasetB using ($1):($16) with lines title "Config 3" ls 4


# Configurations of Coordenadas Baricentricas
set title "Coordenadas Baricentricas"
plot datasetB using ($1):($5) with lines title "Config 0" ls 1, \
     datasetB using ($1):($9) with lines title "Config 1" ls 2, \
     datasetB using ($1):($13) with lines title "Config 2" ls 3, \
     datasetB using ($1):($17) with lines title "Config 3" ls 4

# Showing mean and standard deviation
#f(x) = mean_y
#fit f(x) dataset u ($1):($2) via mean_y
#stddev_y = sqrt(FIT_WSSR / (FIT_NDF + 1 ))

unset multiplot
quit
