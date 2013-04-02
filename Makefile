# -*- coding: utf-8 -*-

par = 0
ifeq ($(par),1)
	CC = mpic++
	SRCS = Game_Of_Life_Parallel.cc
	PROJECT= GOL_P
else
 # Serial
	CC = g++
	SRCS = SEM_serial.cpp fedata.cpp
	PROJECT= SEM_s
endif

objects = $(patsubst %.cpp, %.o,$(SRCS))

LINKFLAGS=
LIBS=
CFLAGS= -Wno-write-strings -g
INCLUDE =


.PHONY: all
all:${PROJECT}
	$(shell etags $(SRCS))

#$(shell ctags -Re)


# regel for exe-filen, dvs link filer
${PROJECT}: $(objects)
	$(CC) $(objects) $(LINKFLAGS) $(LIBS) -o ${PROJECT}


%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@


# special rule for the file containing exodus_ref
.PHONY: clean
clean:
	rm -f *.o
	rm -f $(PROJECT)
	rm -f TAGS
