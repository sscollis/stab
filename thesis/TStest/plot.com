set xrange [0:15]
set xlabel "y"
set ylabel "u_r, u_i"
plot "space.1" u 1:4 w l title "u_r", "space.1" u 1:5 w l title "u_i"
