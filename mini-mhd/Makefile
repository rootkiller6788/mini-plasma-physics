# Makefile for mini-mhd
# Build: make        -- compile library and tests
#        make test   -- run tests
#        make examples -- build examples
#        make clean  -- remove binaries

CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99 -I./include
LDFLAGS = -lm

SRCDIR = src
SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(SRCS:.c=.o)

TEST_SRC = tests/test_mhd.c
TEST_BIN = test_mhd

EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(EXAMPLE_SRCS:.c=.exe)

.PHONY: all test examples clean

all: $(OBJS)
	@echo "Library compiled (objects in src/)"

test: $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): $(TEST_SRC) $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(TEST_SRC) $(filter-out src/mhd.lean,$(OBJS)) $(LDFLAGS)

examples: $(EXAMPLE_BINS)
	@echo "Examples built."
	@for e in $(EXAMPLE_BINS); do echo "Running $$e:"; ./$$e; echo; done

examples/%.exe: examples/%.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ $< $(filter-out src/mhd.lean,$(OBJS)) $(LDFLAGS)

$(SRCDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(TEST_BIN) $(EXAMPLE_BINS)
