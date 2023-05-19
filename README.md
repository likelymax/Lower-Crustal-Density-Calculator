# Lower Crustal Density Calculator

This repository contains code for calculating lower crustal density by comparing theoretical and actual P-to-Ps transmission coefficients. The code is written in C and comprises different components for various calculations. The theoretical coefficient is based on Aki & Richards (2002) Eqns 5.39 and 5.40.

__Code Structure__

RSVRTmatrix: Computes the theoretical P-to-Ps transmission coefficient.  
Trancoeff: Compares the actual and theoretical P-to-Ps transmission coefficients.  
reverb_Trans: Reads all the input files, and return the outputs.  

__Getting Started__

First, compile the code using the provided Makefile to generate the executable file reverb_Trans:

make reverb_Trans

__Usage__

Execute the reverb_Trans program with the following command:

./reverb_Trans -V$vel -C$ratio -D$diff -R$limit/$mdensity/$bound/$inv -S$file -Ooutput.txt

__Parameters__

-V : The input velocity model.  
-C : The actual P-to-Ps transmission coefficient.  
-D : The threshold that can be used to decide the output.  
-R : Lower band of lower crustal density / Upper mantle density / The upper bound of lower crustal density / Interval.  
-S : Input receiver function.  
-O : Output result.  

__Output__

The output file output.txt contains the estimated lower crustal density.

__Dependencies__

This project requires a Unix-like system with a C compiler.

__Contributing__

If you'd like to contribute, please fork the repository and use a feature branch. Pull requests are warmly welcome.

__Contact__

Please reach out if you have any questions or encounter issues using the tool. You can reach me via email at yitanwang@ufl.edu.

__License__

The code in this project is licensed under the MIT license.
