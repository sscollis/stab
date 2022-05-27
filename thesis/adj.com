set xrange [0:15]
set xlabel "y"
set ylabel "Real, Imag"
plot "efun.out" u 1:4 w l title "u_r", "efun.out" u 1:5 w l title "u_i", "adj.out" u 1:4 w l title "v_r", "adj.out" u 1:5 w l title "v_i"
