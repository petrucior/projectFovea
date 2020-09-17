set term postscript eps
#set view 75, 30, 1, 1
set view 75, 110, 1, 1
#set view 90, 0, 1, 1
set xlabel 'x'
set ylabel 'y'
set dgrid3d 30, 30
set zrange [0 : 1.2]
set hidden3d
# Level 0
set output "graph3d_level0.eps"
splot 'graph3d.dat' u 1:2:3 with lines title "Nivel 0", 'pose.dat' u 1:2:3 with lines lc rgb "blue" title "posicao do estimulo"
# Level 1
set output "graph3d_level1.eps"
splot 'graph3d.dat' u 1:2:4 with lines title "Nivel 1", 'pose.dat' u 1:2:3 with lines lc rgb "blue" title "posicao do estimulo"
# Level 2
set output "graph3d_level2.eps"
splot 'graph3d.dat' u 1:2:5 with lines title "Nivel 2", 'pose.dat' u 1:2:3 with lines lc rgb "blue" title "posicao do estimulo"
# Level 3
set output "graph3d_level3.eps"
splot 'graph3d.dat' u 1:2:6 with lines title "Nivel 3", 'pose.dat' u 1:2:3 with lines lc rgb "blue" title "posicao do estimulo"
# Level 4
#set output "graph3d_level4.eps"
#splot 'graph3d.dat' u 1:2:7 with lines title "Nivel 4", 'pose.dat' u 1:2:3 with lines lc rgb "blue" title "posicao do estimulo"

# Plot function
set output "plotFunction.eps"
splot 'function.dat' u 1:2:3 with lines title "Funcao da fovea"
