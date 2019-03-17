#! Makefile

CC=g++
CCFLAGS=-Wall -Werror -std=c++11
LDFLAGS=-lm

MAIN=task2.cpp

DEPS=linalg.h linalg.cpp

SRC= $(MAIN)

OBJS=$(SRC:%.cpp=%.o)

BIN=task2

all: $(BIN)

$(BIN): $(OBJS)
	$(CC) $^  $(LDFLAGS)  -o $@

$(OBJS):$(SRC) $(DEPS)
	$(CC)  $(SRC) -c $(CCFLAS) 

clean:
	rm -f *.o
