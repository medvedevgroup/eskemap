CC = g++

#Parasail
# CFLAGS = -std=c++2a -O3 -Wall -lz -lpthread -lparasail

# all:
# 	$(CC) FindSimSeqs.cpp $(CFLAGS) -o FindSimSeqs

CFLAGS = -std=c++2a -O3 -Wall -L edlib -lz -lpthread -ledlib -lparasail

all:
	$(CC) CalcGlobEditDistance.cpp $(CFLAGS) -o CalcGlobEditDistance