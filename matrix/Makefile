#! Makefile

CC=g++
CCFLAGS=-Wall -Werror -std=c++11
LDFLAGS=-lm

MAIN=test.cpp

DEPS= vector.h 

SRC= vector.cpp $(MAIN)

OBJS=$(SRC:%.cpp=%.o)

BIN=test

all: $(BIN)

$(BIN): $(OBJS)
	$(CC) $(LDFLAGS) $^  -o $@

$(OBJS):$(SRC) $(DEPS)
	$(CC)  $(CCFLAS) -c $(SRC)

clean:
	rm -f *.o
