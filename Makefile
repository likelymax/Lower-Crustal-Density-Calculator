SUBS = fft.o Complex.o Intp2D.o sacio.o



Trancoeff.o: Trancoeff.c
	gcc -c Trancoeff.c

RSVRTmatrix.o: RSVRTmatrix.c
	gcc -c RSVRTmatrix.c
