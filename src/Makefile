#compiler
CC=g++

#flags
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64

BAMTOOLS_ROOT=../bamtools

LIBS = -L./ -L$(BAMTOOLS_ROOT)/lib -lz -lm
INCLUDE = -I$(BAMTOOLS_ROOT)/src 

all: ../bin/stream_sort_pairs

$(BAMTOOLS_ROOT)/lib/libbamtools.a:
	cd $(BAMTOOLS_ROOT) && mkdir -p build && cd build && cmake .. && $(MAKE)
	
stream_sort_pairs ../bin/stream_sort_pairs: stream_sort_pairs.o 
	$(CC) $(CFLAGS) $(INCLUDE) -v stream_sort_pairs.o -o ../bin/stream_sort_pairs $(LIBS)

stream_sort_pairs.o: stream_sort_pairs.cc stream_sort_pairs.h $(BAMTOOLS_ROOT)/lib/libbamtools.a
	$(CC) $(CFLAGS) $(INCLUDE) -c $^

clean:
	rm ./*.o ./*.gch
