CC=mpicc
CFLAGS=-O2 -Wno-unused-result
INCLUDEDIR   = -I./ 
LDFLAGS      = -lgsl -lgslcblas

diagonalization.o: diagonalization.c  
	${CC} ${CFLAGS} ${INCLUDEDIR}  $*.c -c

davidson.o: diagonalization.o davidson.c  
	${CC} ${CFLAGS} ${INCLUDEDIR}  $*.c -c

davidson_drv: davidson_drv.c davidson.o diagonalization.o
	${CC} ${CFLAGS} ${INCLUDEDIR}   diagonalization.o davidson.o davidson_drv.c -o $@ $(CPPFLAGS) $(LDFLAGS)  $(LAPACKLIB) $(BLASLIB) $(F2CLIB) -lm
