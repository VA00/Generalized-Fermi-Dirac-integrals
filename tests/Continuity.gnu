reset
set terminal png enhanced truecolor font "/usr/X11R6/lib/fonts/TTF/Vera.ttf" 14 size 1200, 1200
#set terminal postscript
set format y "%.16e"
plot 'continuity1.dat' u 3 with points lc rgb "black" title "icc", 'continuity1_gcc.dat' u 3 with points lc rgb "red" title "gcc"