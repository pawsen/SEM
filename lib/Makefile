# -*- coding: utf-8; mode: makefile -*-

# source files.
SRC =  libelement.cpp

OBJ = $(SRC:.cpp=.o)

OUT = libelement.a

# include directories
INCLUDES = -I. -I../

# C++ compiler flags (-g -O2 -Wall)
CCFLAGS = -O2

# compiler
CCC = g++

.SUFFIXES: .cpp

default: dep $(OUT)

.cpp.o:
	$(CCC) $(INCLUDES) $(CCFLAGS) -c $< -o $@

$(OUT): $(OBJ)
	ar rcs $(OUT) $(OBJ)

depend: dep

dep:
	makedepend -- $(CFLAGS) -- $(INCLUDES) $(SRC)

clean:
	rm -f $(OBJ) $(OUT)
# DO NOT DELETE

libelement.o: libelement.h ../fedata.h /usr/include/math.h
libelement.o: /usr/include/features.h /usr/include/sys/cdefs.h
libelement.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
libelement.o: /usr/include/gnu/stubs-64.h /usr/include/bits/huge_val.h
libelement.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
libelement.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
libelement.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
libelement.o: ../matrix/CD.c ../matrix/B3.c ../matrix/B4.c ../matrix/B5.c
libelement.o: ../matrix/B6.c ../matrix/B7.c ../matrix/Me3.c ../matrix/Me4.c
libelement.o: ../matrix/Me5.c ../matrix/Me6.c ../matrix/Me7.c
