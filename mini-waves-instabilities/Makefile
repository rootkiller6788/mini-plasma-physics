# Makefile -- mini-waves-instabilities
# Plasma waves and instabilities library, tests, and examples.
#
# References:
#   Stix (1992), Swanson (2003), Mikhailovskii (1974), Goldston & Rutherford (1995)
#   MIT 22.611, Princeton PHY 521, Berkeley PHYS 242

CC       := gcc
CFLAGS   := -Wall -Wextra -std=c11 -O2 -g -pedantic
INC_DIR  := include
SRC_DIR  := src
TST_DIR  := tests
EXA_DIR  := examples
DM_DIR   := demos
BNC_DIR  := benches

INCLUDES := -I$(INC_DIR)
LIBS     := -lm

SRCS     := $(wildcard $(SRC_DIR)/*.c)
OBJS     := $(SRCS:.c=.o)
LIB      := libwaves.a

TEST_SRCS := $(wildcard $(TST_DIR)/*.c)
TESTS     := $(TEST_SRCS:.c=.exe)

EXAMPLE_SRCS := $(wildcard $(EXA_DIR)/*.c)
EXAMPLES     := $(EXAMPLE_SRCS:.c=.exe)

DEMO_SRCS := $(wildcard $(DM_DIR)/*.c)
DEMOS     := $(DEMO_SRCS:.c=.exe)

BENCH_SRCS := $(wildcard $(BNC_DIR)/*.c)
BENCHES    := $(BENCH_SRCS:.c=.exe)

.PHONY: all clean test examples demos benches lib

all: lib test examples demos benches

lib: $(LIB)

$(LIB): $(OBJS)
	ar rcs $@ $^

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(INC_DIR)/*.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(TST_DIR)/%.exe: $(TST_DIR)/%.c $(LIB)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L. -lwaves $(LIBS) -o $@

$(EXA_DIR)/%.exe: $(EXA_DIR)/%.c $(LIB)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L. -lwaves $(LIBS) -o $@

$(DM_DIR)/%.exe: $(DM_DIR)/%.c $(LIB)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L. -lwaves $(LIBS) -o $@

$(BNC_DIR)/%.exe: $(BNC_DIR)/%.c $(LIB)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L. -lwaves $(LIBS) -o $@

test: $(TESTS)
	@echo "=== Running Tests ==="
	@for t in $(TESTS); do echo "--- $$t ---"; ./$$t || exit 1; done
	@echo "=== All Tests Passed ==="

examples: $(EXAMPLES)
	@echo "=== Running Examples ==="
	@for e in $(EXAMPLES); do echo "--- $$e ---"; ./$$e; done

demos: $(DEMOS)
	@for d in $(DEMOS); do echo "--- $$d ---"; ./$$d; done

benches: $(BENCHES)
	@for b in $(BENCHES); do echo "--- $$b ---"; ./$$b; done

clean:
	rm -f $(SRC_DIR)/*.o $(LIB)
	rm -f $(TST_DIR)/*.exe $(EXA_DIR)/*.exe
	rm -f $(DM_DIR)/*.exe $(BNC_DIR)/*.exe
