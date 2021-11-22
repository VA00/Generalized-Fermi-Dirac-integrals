reset
set terminal png enhanced truecolor font "/usr/X11R6/lib/fonts/TTF/Vera.ttf" 14 size 1200, 1200
#set terminal postscript
set format y "%.16e"


set arrow from 256,graph(0,0) to 256, graph(1,1) nohead lc rgb "blue"
set label "Sommerfeld\nexpansion" at graph 0.85,0.5
set label "Double exp.\nquadrature" at graph 0.05,0.5
set label "k=1, theta = 256 \261 128 ULP's, theta=4" at graph 0.5,0.05
set x2tics

set output "continuity0.png"
set title "F"
plot 'continuity.dat' u 12:1 every 2 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 12:1 every 3 with points lc rgb "red" title "Mathematica", 'continuity_fxt.dat' u 12:1 every 5 with points lc rgb "green" title "dfermi", '' u 11:(1/0) ax x2y1 notitle

set output "continuity1.png"
set title "dF/dtheta"
plot 'continuity.dat' u 12:2 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 12:2 with points lc rgb "red" title "Mathematica", 'continuity_fxt.dat' u 12:2 with points lc rgb "green" title "dfermi", '' u 11:(1/0) ax x2y1 notitle

set output "continuity2.png"
set title "d^2F/dtheta^2"
plot 'continuity.dat' u 12:3 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 12:3 with points lc rgb "red" title "Mathematica", 'continuity_fxt.dat' u 12:3 with points lc rgb "green" title "dfermi", '' u 11:(1/0) ax x2y1 notitle


set output "continuity3.png"
set title "d^3F/dtheta^3"
plot 'continuity.dat' u 12:4 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 12:4 with points lc rgb "red" title "Mathematica", '' u 11:(1/0) ax x2y1 notitle

set output "continuity4.png"
set title "dF/deta"
plot 'continuity.dat' u 12:5 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 12:5 with points lc rgb "red"   title "Mathematica", 'continuity_fxt.dat' u 12:5 with points lc rgb "green" title "dfermi", '' u 11:(1/0) ax x2y1 notitle

set output "continuity5.png"
set title "d^2F/deta/dtheta"                                                                                                        
plot 'continuity.dat' u 12:6 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 12:6 with points lc rgb "red"   title "Mathematica", 'continuity_fxt.dat' u 12:6 with points lc rgb "green" title "dfermi", '' u 11:(1/0) ax x2y1 notitle

set output "continuity6.png"                                                                                                        
set title "d^3F/deta/dtheta^2"
plot 'continuity.dat' u 12:7 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 12:7 with points lc rgb "red"   title "Mathematica", '' u 11:(1/0) ax x2y1 notitle

set output "continuity7.png"                                                                                                        
set title "d^2F/deta^2"
plot 'continuity.dat' u 12:8 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 12:8 with points lc rgb "red"   title "Mma", 'continuity_fxt.dat' u 12:8 with points lc rgb "green" title "dfermi", '' u 11:(1/0) ax x2y1 notitle

set output "continuity8.png"
set title "d^3F/deta^2/dtheta"
plot 'continuity.dat' u 12:9 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 12:9 with points lc rgb "red" title "Mma", '' u 11:(1/0) ax x2y1 notitle

set output "continuity9.png"
set title "d^3F/deta^3"
plot 'continuity.dat' u 12:10 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 12:10 with points lc rgb "red" title "Mma", '' u 11:(1/0) ax x2y1 notitle
