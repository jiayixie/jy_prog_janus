
set order = 4
foreach factor ( 2 5 8 )
./get_BS $factor $order
gnuplot << EOF
set term postscript color enhanced
set output "BS_f${factor}_o${order}.ps"
set title "factor${factor} Order${order}"
plot "Bs.0.dat" w l, "Bs.1.dat" w l, "Bs.2.dat" w l, "Bs.3.dat" w l, "Bs.4.dat" w l
q
EOF
end

set order = 3
foreach factor ( 1 2 3 )
./get_BS $factor $order
gnuplot << EOF
set term postscript color enhanced
set output "BS_f${factor}_o${order}.ps"
set title "factor${factor} Order${order}"
plot "Bs.0.dat" w l, "Bs.1.dat" w l, "Bs.2.dat" w l, "Bs.3.dat" w l, "Bs.4.dat" w l
q
EOF


end
