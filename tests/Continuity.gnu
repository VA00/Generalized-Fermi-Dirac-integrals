reset
set terminal png enhanced truecolor font "/usr/X11R6/lib/fonts/TTF/Vera.ttf" 14 size 1200, 1200
#set terminal postscript
set format y "%.16e"

set output "continuity2.png"
set title "F"
plot 'continuity.dat' u 1:2 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 1:2 with points lc rgb "red" title "Mathematica"

set output "continuity3.png"
set title "dF/dtheta"
plot 'continuity.dat' u 1:3 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 1:3 with points lc rgb "red" title "Mathematica"

set output "continuity4.png"
set title "d^2F/dtheta^2"
plot 'continuity.dat' u 1:4 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 1:4 with points lc rgb "red" title "Mathematica"


set output "continuity5.png"
set title "d^3F/dtheta^3"
plot 'continuity.dat' u 1:5 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 1:5 with points lc rgb "red" title "Mathematica"

set output "continuity6.png"
set title "dF/deta"
plot 'continuity.dat' u 1:6 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 1:6 with points lc rgb "red"   title "Mathematica"

set output "continuity7.png"
set title "d^2F/deta/dtheta"                                                                                                        
plot 'continuity.dat' u 1:7 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 1:7 with points lc rgb "red"   title "Mathematica"

set output "continuity8.png"                                                                                                        
set title "d^3F/deta/dtheta^2"
plot 'continuity.dat' u 1:8 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 1:8 with points lc rgb "red"   title "Mathematica"

set output "continuity9.png"                                                                                                        
set title "d^2F/deta^2"
plot 'continuity.dat' u 1:9 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 1:9 with points lc rgb "red"   title "Mma"

set output "continuity10.png"
set title "d^3F/deta^2/dtheta"
plot 'continuity.dat' u 1:10 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 1:10 with points lc rgb "red" title "Mma"

set output "continuity11.png"
set title "d^3F/deta^3"
plot 'continuity.dat' u 1:11 with points lc rgb "black" title "libfermidirac", 'continuity_Mma.dat' u 1:11 with points lc rgb "red" title "Mma"
