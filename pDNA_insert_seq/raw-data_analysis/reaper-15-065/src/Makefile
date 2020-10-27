



CC = gcc
COPTIM = -Wall -O3 -std=gnu99 -DHAVE_ZLIB=1 -DTRACING_ON=0
	#
	# -march=native -mtune=native
	# icc options for optimisation.  COPTIM = -xHOST -O3 -no-prec-div -static -DHAVE_ZLIB=1

CDEBUG = -g -Wall
CARRRR = -Wall -O3 -fPIC -std=gnu99
CBIG   = -Wall -O2 -g -DMAXFIELDSIZE=255
CSTATIC= -Wall -O3 -std=gnu99 -DHAVE_ZLIB=1 -static
CSMALL = -Wall -Os -std=gnu99 -DHAVE_ZLIB=1 -fpack-struct

ifndef R_INCLUDE_DIR
   # R_INCLUDE_DIR = /nfs/research2/enright/local/lib64/R/include
	# T = optim
endif

C=$(CDEBUG)
L = -lz -lm
RI = -I${R_INCLUDE_DIR} -I${R_INCLUDE_DIR}/${R_ARCH}
RI = -I${R_INCLUDE_DIR}

ifndef T
	T = optim
endif

ifdef R_INCLUDE_DIR
	T = reaper.so
	.DEFAULT_GOAL := reaper.so
else
	#T = optim
	#.DEFAULT_GOAL := .DUMMY
endif



ifeq ("$(T)", "optim")
	C = $(COPTIM)
else ifeq ("$(T)", "big")
	C = $(CBIG)
else ifeq ("$(T)", "lean")
	C = $(CSMALL)
else ifeq ("$(T)", "debug")
	C = $(CDEBUG)
else ifeq ("$(T)", "static")
	C = $(CSTATIC)
else ifeq ("$(T)", "reaper.so")
	C = $(CARRRR)
else
	C = $(COPTIM)
endif


all : reaper tally minion swan
.PHONY : all


reaper: reaper.o sw.o slib.o table.o trint.o version.o dna.o
	$(CC) -fstrict-aliasing $(C) -o reaper trint.o version.o dna.o reaper.o sw.o table.o slib.o $(L)

swan: swan.o sw.o slib.o dna.o
	$(CC) -fstrict-aliasing $(C) -o swan slib.o swan.o sw.o dna.o $(L)

tally: tally.o slib.o trint.o version.o dna.o sw.o
	$(CC) -fstrict-aliasing $(C) -o tally trint.o dna.o tally.o version.o slib.o $(L)

minion: minion.o slib.o version.o trint.o sw.o
	$(CC) -fstrict-aliasing $(C) -o minion version.o trint.o slib.o minion.o sw.o $(L)

reaper.so: .clean bindings.o tally-R.o reaper-R.o sw.o slib.o table.o trint.o version.o
	R CMD SHLIB -o reaper.so bindings.o reaper-R.o tally-R.o dna.o sw.o slib.o table.o trint.o version.o
	# $(CC) -std=gnu99 $(C) -shared -o reaper.so bindings.o sw.o slib.o table.o trint.o


tmptst:
	$(CC) -DRREAPERR -fstrict-aliasing -Wall -O2 -fPIC -c -o bindings.o reaper.c -I/nfs/research2/enright/local/lib64/R/include -I/nfs/research2/enright/local/lib64/R/include/

mytest:
	echo $(.TARGET)
	echo $(C)
	echo $(T) x=$(x) y=$(y)


bindings.o: bindings.c reaper.h tally.h
	$(CC) -fstrict-aliasing $(C) -c -o bindings.o bindings.c $(RI)

chair: chair.o table.o slib.o
	$(CC) -g -fstrict-aliasing $(C) -o chair chair.o table.o slib.o $(L)

trint.o: trint.c trint.h
	$(CC) -fstrict-aliasing $(C) -c -o trint.o trint.c

sw.o: sw.c sw.h
	$(CC) -fstrict-aliasing $(C) -c -o sw.o sw.c

version.o: version.c version.h
	$(CC) -fstrict-aliasing $(C) -c -o version.o version.c

slib.o: slib.c slib.h
	$(CC) -fstrict-aliasing $(C) -c -o slib.o slib.c

minion.o: minion.c
	$(CC) -fstrict-aliasing $(C) -c -o minion.o minion.c

swan.o: swan.c
	$(CC) -fstrict-aliasing $(C) -c -o swan.o swan.c

dna.o: dna.c dna.h
	$(CC) -g -fstrict-aliasing $(C) -c -o dna.o dna.c

table.o: table.c table.h slib.h
	$(CC) -g -fstrict-aliasing $(C) -c -o table.o table.c

#reaper.o: reaper.c slib.h sw.c slib.c trint.c table.c
#	$(CC) -fstrict-aliasing $(C) -c -o reaper.o reaper.c

reaper.o: reaper.c sw.h slib.h trint.h table.h
	$(CC) -fstrict-aliasing $(C) -c -o reaper.o reaper.c

reaper-R.o: reaper.c sw.h slib.h trint.h table.h
	$(CC) -DBUILD_R_BINDINGS -fstrict-aliasing $(C) -c -o reaper-R.o reaper.c


tally.o: tally.c
	$(CC) -fstrict-aliasing $(C) -c -o tally.o tally.c

matecare4.o: matecare4.c
	$(CC) -fstrict-aliasing $(C) -c -o matecare4.o matecare4.c

tally-R.o: tally.c
	$(CC) -DBUILD_R_BINDINGS -fstrict-aliasing $(C) -c -o tally-R.o tally.c


wc:
	wc trint.c reaper.c sw.c table.c slib.c tally.c swan.c

clean .clean:
	rm -f {bindings,trint,reaper,reaper-R,slib,sw,table,tally,tally-R,version,minion,dna,swan}.o



matecare4: matecare4.o slib.o trint.o version.o dna.o sw.o
	$(CC) -fstrict-aliasing $(C) -o matecare4 sw.o trint.o dna.o matecare4.o version.o slib.o $(L)
