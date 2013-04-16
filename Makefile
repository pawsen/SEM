# -*- coding: utf-8; mode: makefile -*-

# mpirun -np 2 -mca btl tcp,self SEM_p

par = 1
ifeq ($(par),1)
	CC = mpic++
	SRCS = SEM_parallel.cpp fedata.cpp vtk.cpp aux.cpp transient.cpp \
	gauss_legendre.c mpidata.cpp
	PROJECT= SEM_p
	CFLAGS1= -DMY_MPI
else
 # Serial
	CC = g++
	SRCS = SEM_parallel.cpp fedata.cpp vtk.cpp aux.cpp transient.cpp \
	gauss_legendre.c
	PROJECT= SEM_s
endif

objects = $(patsubst %.cpp, %.o,$(SRCS))
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
	$(shell etags $(SRCS))

#$(shell ctags -Re)


# regel for exe-filen, dvs link filer
${PROJECT}: $(objects)
	$(CC) $(objects) $(LINKFLAGS) $(LIBS) -o ${PROJECT}


%.o: %.cpp
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@


# special rule for the file containing exodus_ref
.PHONY: clean
clean:
	rm -f *.o
	rm -f $(PROJECT)
	rm -f TAGS
