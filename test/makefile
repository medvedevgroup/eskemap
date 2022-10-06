CC = g++
CFLAGS = -std=c++2a -O3
LINKS = -pthread -lgtest -lgtest_main -lz -lm -lpthread
SRCFILES = gtest_main.cpp Sketch_unittest.cpp IO_unittest.cpp Measures_unittest.cpp Thomology_unittest.cpp Index_unittest.cpp ../../software/minimap2/libminimap2.a

all:
	$(CC) $(CFLAGS) -o MyTest $(SRCFILES) $(LINKS) -I../src
