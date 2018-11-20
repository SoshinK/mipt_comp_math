#! Makefile

CC=gcc
CCFLAGS=-Wall -Werror
LDFLAGS=-lm

MAIN=611_soshin.c

DEPS=mtrx.h

SRC=mtrx.c $(MAIN)

OBJS=$(SRC:%.c=%.o)

BIN=SOR

all: $(BIN)

$(BIN): $(OBJS)
	$(CC) $^  $(LDFLAGS)  -o $@

$(OBJS):$(SRC) $(DEPS)
	$(CC)  $(SRC) -c $(CCFLAS) 

clean:
	rm -f *.o
