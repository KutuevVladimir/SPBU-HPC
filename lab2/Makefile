CC = gcc
CFLAGS = -fopenmp -Wall -Wextra -pedantic -O3 -ffast-math
LD = gcc 
LDFLAGS = -fopenmp -lm
DBGCFLAGS = -ggdb -fopenmp -lm

PROGRAM = lu

.PHONY: all debug clean

all: $(PROGRAM)

$(PROGRAM) : $(PROGRAM).o
	$(LD) $(LDFLAGS) $< -o $(PROGRAM)

debug : $(PROGRAM).c
	$(CC) $(DBGCFLAGS) $< -o $(PROGRAM)_dbg

$(PROGRAM).o : $(PROGRAM).c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(PROGRAM) $(PROGRAM).o $(PROGRAM)_dbg
