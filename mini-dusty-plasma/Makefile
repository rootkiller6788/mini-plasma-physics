# Makefile -- mini-dusty-plasma
# Builds dusty plasma C library and test suite.
# Usage: make | make test | make examples | make clean

CC       = gcc
CFLAGS   = -std=c99 -Wall -Wextra -Wno-unused-parameter -O2 -g -Iinclude
LDFLAGS  = -lm
BUILDDIR = build

SRCS = src/dusty_plasma_core.c src/dusty_charging.c src/dusty_waves.c
SRCS += src/dusty_crystal.c src/dusty_forces.c src/dusty_transport.c
SRCS += src/dusty_dynamics.c src/dusty_potential.c
OBJS = $(SRCS:src/%.c=$(BUILDDIR)/%.o)
TESTBIN = $(BUILDDIR)/test_dusty

.PHONY: all lib test examples clean

all: lib test

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BUILDDIR)/%.o: src/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(TESTBIN): $(OBJS) tests/test_dusty.c | $(BUILDDIR)
	$(CC) $(CFLAGS) tests/test_dusty.c $(OBJS) -o $(TESTBIN) $(LDFLAGS)

lib: $(OBJS)

test: $(TESTBIN)
	./$(TESTBIN)

$(BUILDDIR)/%: examples/%.c $(OBJS) | $(BUILDDIR)
	$(CC) $(CFLAGS) $< $(OBJS) -o $@ $(LDFLAGS)

examples: $(BUILDDIR)/example_daw $(BUILDDIR)/example_charging
examples: $(BUILDDIR)/example_crystal $(BUILDDIR)/example_sheath

clean:
	rm -rf $(BUILDDIR)
