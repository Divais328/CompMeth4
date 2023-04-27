set terminal jpeg size 700,600 font "Arials,16"

set xlabel "w"
set ylabel "Number of iterations"
set grid

set output "SORM_W_32.jpeg"
set title "Successive OverRelaxation method: w optimization.\nN = M = 32" 
plot [1:2][400:800] "SORM_W_32.txt"  w l lw 3



set xlabel "N"
set ylabel "Error"
set logscale
set output "Errors.jpeg"
set title "Gauss-Seidel method Error" 
plot  [4:256][1.e-7:1]"Errors_GSM.txt"  w l lw 3 t "5-points method", "Errors_GS9PM.txt"  w l lw 3 t "9-points method", 1./x**2 lw 3 t "1/N^2", 1./x**4 lw 3 t "1/N^4"



set xlabel "X"
set ylabel "Y"
set zlabel "U"
set pm3d
set hidden3d 
unset logscale 

set output "Exact_Solution_128.jpeg"
set title "Exact Solution: N = M = 128" 
splot "Exact_Solution_128.txt"  w l 

set output "JM_32.jpeg"
set title "Jacobi method: N = M = 32" 
splot "JM_32.txt"  w l 

set output "GSM_32.jpeg"
set title "Gauss-Seidel method: N = M = 32" 
splot "GSM_32.txt"  w l 

set output "SORM_32_1.69.jpeg"
set title "Successive OverRelaxation method: N = M = 32, w = 1.69" 
splot "SORM_32_1.69.txt"  w l 


