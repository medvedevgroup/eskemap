CC = g++
CFLAGS = -std=c++2a -O3
LINKS = -pthread -lgtest -lgtest_main
SRCFILES = gtest_main.cpp Sketch_unittest.cpp IO_unittest.cpp Measures_unittest.cpp Thomology_unittest.cpp

all:
	$(CC) $(CFLAGS) -o MyTest $(SRCFILES) $(LINKS) -I../src
