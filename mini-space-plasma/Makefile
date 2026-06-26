# mini-space-plasma Makefile
# Reference: Kivelson & Russell, Goldston & Rutherford, MIT 22.611
CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
AR = ar
ARFLAGS = rcs

SRCDIR = src
INCDIR = include
TESTDIR = tests
EXAMPLEDIR = examples
BUILDDIR = build
LIBDIR = lib

# Source files
SRCS = $(SRCDIR)/plasma_parameters.c \
       $(SRCDIR)/mhd_core.c \
       $(SRCDIR)/plasma_waves.c \
       $(SRCDIR)/solar_wind.c \
       $(SRCDIR)/magnetosphere.c

OBJS = $(SRCS:$(SRCDIR)/%.c=$(BUILDDIR)/%.o)
LIB  = $(LIBDIR)/libspaceplasma.a
TEST = $(BUILDDIR)/test_plasma

# Examples
EXAMPLES = $(BUILDDIR)/example_solar_wind \
           $(BUILDDIR)/example_magnetosphere \
           $(BUILDDIR)/example_plasma_waves

.PHONY: all clean test lib examples

all: lib test examples

# ---- Directories ----
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(LIBDIR):
	mkdir -p $(LIBDIR)

# ---- Library ----
$(BUILDDIR)/%.o: $(SRCDIR)/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(LIB): $(OBJS) | $(LIBDIR)
	$(AR) $(ARFLAGS) $@ $^

lib: $(LIB)

# ---- Test ----
$(TEST): $(TESTDIR)/test_plasma.c $(LIB) | $(BUILDDIR)
	$(CC) $(CFLAGS) $< -L$(LIBDIR) -lspaceplasma -o $@ $(LDFLAGS)

test: $(TEST)
	./$(TEST)

# ---- Examples ----
$(BUILDDIR)/example_solar_wind: $(EXAMPLEDIR)/example_solar_wind.c $(LIB) | $(BUILDDIR)
	$(CC) $(CFLAGS) $< -L$(LIBDIR) -lspaceplasma -o $@ $(LDFLAGS)

$(BUILDDIR)/example_magnetosphere: $(EXAMPLEDIR)/example_magnetosphere.c $(LIB) | $(BUILDDIR)
	$(CC) $(CFLAGS) $< -L$(LIBDIR) -lspaceplasma -o $@ $(LDFLAGS)

$(BUILDDIR)/example_plasma_waves: $(EXAMPLEDIR)/example_plasma_waves.c $(LIB) | $(BUILDDIR)
	$(CC) $(CFLAGS) $< -L$(LIBDIR) -lspaceplasma -o $@ $(LDFLAGS)

examples: $(EXAMPLES)

# ---- Clean ----
clean:
	rm -rf $(BUILDDIR) $(LIBDIR)

# ---- Quick compile test (no linking) ----
compile-test: $(OBJS)

# ---- Help ----
help:
	@echo "Targets:"
	@echo "  all      - Build library, tests, and examples"
	@echo "  lib      - Build static library libspaceplasma.a"
	@echo "  test     - Build and run tests"
	@echo "  examples - Build example programs"
	@echo "  clean    - Remove build artifacts"
	@echo "  help     - Show this help"
