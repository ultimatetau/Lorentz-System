# Lorentz-System

Very simple Lorentz System utilizing Runge-Kutta 4 in C++. 

Compile with: "g++ lorenz-System-RK4.cpp -o lorenz"
run with: "./lorenz"

open gnuplot with "gnuplot". Output will be saved as "LorenzHistory.txt"

To plot, type: 'splot "LorenzHistory.txt" using 2:3:4 lt 1 lw 2 with lines'

