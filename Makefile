# -*- coding: utf-8; mode: makefile -*-

# mpirun -np 2 -mca btl tcp,self SEM_p

par = 1
ifeq ($(par),1)
	CC = mpic++
    MPISRC= mpidata.cpp
	PROJECT= SEM_p
	CFLAGS1= -DMY_MPI -DMY_OVERLAP
else
 # Serial
	CC = g++
	PROJECT= SEM_s
endif

SRCS= SEM_parallel.cpp fedata.cpp vtk.cpp aux.cpp transient.cpp \
	${MPISRC}
SRCS_C= gauss_legendre.c
objects= $(patsubst %.cpp, %.o,$(SRCS)) $(patsubst %.c, %.o,$(SRCS_C))
LINKFLAGS=
CFLAGS= ${CFLAGS1} -Wno-deprecated -g -O2

########
# BLAS #
########
LIBBLAS = -lblas
INCBLAS = -I/usr/include/

#######
# VTK #
#######
LIBVTK = -L/usr/lib/ -lvtkCommon -lvtkIO -lvtkFiltering
INCVTK = -I/usr/include/vtk-5.8


LIBS= ${LIBVTK} ${LIBBLAS} -Llib/ -lelement
INCLUDE = ${INCVTK} ${INCBLAS} -Ilib/


.PHONY: all
all:${PROJECT}
	$(shell cd src; etags $(SRCS) $(SRCS_C))

#$(shell ctags -Re)


# regel for exe-filen, dvs link filer
${PROJECT}: $(objects)
	$(CC) $(objects) $(LINKFLAGS) $(LIBS) -o ${PROJECT}


%.o: %.cpp
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@
%.o: src/%.cpp
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@
%.o: src/%.c
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@


# special rule for the file containing exodus_ref
.PHONY: clean
clean:
	rm -f *.o
	rm -f $(PROJECT)
	rm -f TAGS
