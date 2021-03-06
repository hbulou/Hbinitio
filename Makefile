PLPLOTINC=-I/home/bulou/ownCloud/src/PLplot/plplot-5.10.0/examples/c -I/usr/include/plplot
PLPLOTLIB=-lplplotd 
CC=mpicc
CFLAGS=-O2 -Wno-unused-result
INCLUDEDIR   = -I./ 
LDFLAGS      = -lgsl -lgslcblas

diagonalization.o: diagonalization.c  
	${CC} ${CFLAGS} ${INCLUDEDIR}  $*.c -c

GramSchmidt.o: GramSchmidt.c dot.o  
	${CC} ${CFLAGS} ${INCLUDEDIR}  $*.c -c


davidson.o: diagonalization.o GramSchmidt.o dot.o davidson.c  
	${CC} ${CFLAGS} ${INCLUDEDIR}  $*.c -c

davidson2D.o: diagonalization.o GramSchmidt.o dot.o davidson2D.c  
	${CC} ${CFLAGS} ${INCLUDEDIR}  $*.c -c

steepest_descent.o: steepest_descent.c
	${CC} ${CFLAGS} ${INCLUDEDIR}  $*.c -c

steepest_descent_drv: steepest_descent_drv.c steepest_descent.o
	${CC} ${CFLAGS} ${INCLUDEDIR}   steepest_descent_drv.c steepest_descent.o -o $@ $(CPPFLAGS) $(LDFLAGS)  $(LAPACKLIB) $(BLASLIB) $(F2CLIB) -lm

conjugate_gradient.o: conjugate_gradient.c
	${CC} ${CFLAGS} ${INCLUDEDIR}  $*.c -c

conjugate_gradient_drv: conjugate_gradient_drv.c conjugate_gradient.o
	${CC} ${CFLAGS} ${INCLUDEDIR}   conjugate_gradient_drv.c conjugate_gradient.o -o $@ $(CPPFLAGS) $(LDFLAGS)  $(LAPACKLIB) $(BLASLIB) $(F2CLIB) -lm

dot.o: dot.c  
	${CC} ${CFLAGS} ${INCLUDEDIR}  $*.c -c

davidson_drv: davidson_drv.c davidson.o GramSchmidt.o diagonalization.o dot.o
	${CC} ${CFLAGS} ${INCLUDEDIR}   diagonalization.o davidson.o dot.o GramSchmidt.o davidson_drv.c -o $@ $(CPPFLAGS) $(LDFLAGS)  $(LAPACKLIB) $(BLASLIB) $(F2CLIB) -lm

plot2D.o : plot2D.c
	${CC} ${CFLAGS} ${INCLUDEDIR} $(PLPLOTINC) $*.c -c

davidson2D_drv: davidson2D_drv.c davidson2D.o GramSchmidt.o diagonalization.o dot.o 
	${CC} ${CFLAGS} ${INCLUDEDIR}     diagonalization.o davidson2D.o dot.o GramSchmidt.o davidson2D_drv.c -o $@ $(CPPFLAGS) $(LDFLAGS)  $(LAPACKLIB) $(BLASLIB) $(F2CLIB) $(PLPLOTLIB) -lm
