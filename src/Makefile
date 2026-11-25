# Build the planar binary square-well DFT solver.
#
#   make                 # uses $(CC) with OpenMP if available
#   make CC=gcc-15       # e.g. Homebrew GCC on macOS
#   make OPENMP=0        # serial build
#
# OpenMP is enabled automatically if the compiler supports -fopenmp
# (Apple clang does not; use Homebrew/MacPorts GCC for a parallel build).
CC      ?= cc
CFLAGS  ?= -O3 -march=native
OPENMP  ?= $(shell echo 'int main(void){return 0;}' | $(CC) -fopenmp -x c - -o /dev/null 2>/dev/null && echo 1 || echo 0)
WARN    = -Wall -Wextra -std=c11 -D_POSIX_C_SOURCE=200809L

ifeq ($(OPENMP),1)
  OMPFLAGS = -fopenmp
endif

BIN  = sw_interface
OBJS = sw_interface.o white_bear.o

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) $(OMPFLAGS) $(OBJS) -lm -o $@

%.o: %.c white_bear.h
	$(CC) $(CFLAGS) $(OMPFLAGS) $(WARN) -c $< -o $@

clean:
	rm -f $(OBJS) $(BIN)

.PHONY: clean
