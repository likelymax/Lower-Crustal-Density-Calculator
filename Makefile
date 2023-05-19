SUBS = fft.o Complex.o Intp2D.o sacio.o

reverb_Trans: reverb_Trans.o Trancoeff.o RSVRTmatrix.o $(SUBS)
	$(LINK.c) -o $@ $@.o Trancoeff.o  RSVRTmatrix.o $(SUBS)
	gcc $@.o -o $@ RSVRTmatrix.o Trancoeff.o fft.o Complex.o Intp2D.o sacio.o GMT_grd_init.o  

reverb_Trans.o: reverb_Trans.c
	gcc -c reverb_Trans.c

Trancoeff.o: Trancoeff.c
	gcc -c Trancoeff.c

RSVRTmatrix.o: RSVRTmatrix.c
	gcc -c RSVRTmatrix.c
	
clean:
	rm -f *.o
