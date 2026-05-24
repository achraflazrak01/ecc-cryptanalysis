CC = gcc
CFLAGS = -std=c99 -Wall -I./include
SRC = $(wildcard src/*.c)
TEST_SRC = $(wildcard tests/test_*.c)
# Choose only one test file (change the name as needed)
MAIN_TEST = tests/test_bn.c

all: ecc-test

ecc-test: $(SRC) $(MAIN_TEST)
	$(CC) $(CFLAGS) -o $@ $^

clean:
	del *.exe *.o 2>nul || rm -f *.exe *.o

.PHONY: all clean