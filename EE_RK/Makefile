#! Makefile

CC=g++
CCFLAGS=-Wall -Werror -std=c++11
LDFLAGS=-lm

MAIN=611_Soshin.cpp

DEPS=linalg.h linalg.cpp

SRC= $(MAIN)

OBJS=$(SRC:%.cpp=%.o)

BIN=task

all: $(BIN)

$(BIN): $(OBJS)
	$(CC) $^  $(LDFLAGS)  -o $@

$(OBJS):$(SRC) $(DEPS)
	$(CC)  $(SRC) -c $(CCFLAS) 

clean:
	rm -f *.o
