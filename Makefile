# -*- coding: utf-8 -*-

par = 0
ifeq ($(par),1)
	CC = mpic++
	SRCS = Game_Of_Life_Parallel.cc
	PROJECT= GOL_P
else
 # Serial
	CC = g++
	SRCS = SEM_serial.cpp
	PROJECT= SEM_s
# EXODUS SPECIFIC STUFF
	INC_EXODUS =
	LIBEXODUS =
endif

objects = $(patsubst %.cc, %.o,$(SRCS))

LINKFLAGS=
LIBS= $(LIBEXODUS)
CFLAGS= -Wno-write-strings -g
INCLUDE = $(INC_EXODUS)


.PHONY: all
all:${PROJECT}
	$(shell etags $(SRCS))

#$(shell ctags -Re)


# regel for exe-filen, dvs link filer
${PROJECT}: $(objects)
	$(CC) $(objects) $(LINKFLAGS) $(LIBS) -o ${PROJECT}


%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

# special rule for the file containing exodus_ref
exodus.o : exodus.cc
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@ 


# special rule for the file containing exodus_ref
.PHONY: clean
clean:
	rm -f *.o
	rm -f $(PROJECT)
	rm -f TAGS
