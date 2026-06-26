# Makefile -- mini-laser-plasma
# MIT 22.611 / Kruer / Gibbon / Goldston & Rutherford

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -Iinclude -D_USE_MATH_DEFINES
LDFLAGS  = -lm

SRC_DIR  = src
INC_DIR  = include
TST_DIR  = tests
EXM_DIR  = examples
BLD_DIR  = build

SRCS = $(SRC_DIR)/plasma_params.c        $(SRC_DIR)/laser_propagation.c        $(SRC_DIR)/wakefield.c        $(SRC_DIR)/instabilities.c        $(SRC_DIR)/absorption.c        $(SRC_DIR)/ionization.c        $(SRC_DIR)/particle_motion.c

OBJS = $(patsubst $(SRC_DIR)/%.c,$(BLD_DIR)/%.o,$(SRCS))

TEST_SRC = $(TST_DIR)/test_laser_plasma.c
TEST_BIN = $(BLD_DIR)/test_laser_plasma

EX1_SRC = $(EXM_DIR)/example_critical_density.c
EX1_BIN = $(BLD_DIR)/example_critical_density
EX2_SRC = $(EXM_DIR)/example_wakefield.c
EX2_BIN = $(BLD_DIR)/example_wakefield
EX3_SRC = $(EXM_DIR)/example_instabilities.c
EX3_BIN = $(BLD_DIR)/example_instabilities

.PHONY: all test examples clean

all: test examples

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.c | $(BLD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_BIN): $(OBJS) $(TEST_SRC) | $(BLD_DIR)
	$(CC) $(CFLAGS) $(TEST_SRC) $(OBJS) $(LDFLAGS) -o $@

test: $(TEST_BIN)
	$(TEST_BIN)

$(EX1_BIN): $(SRC_DIR)/plasma_params.c $(EX1_SRC) | $(BLD_DIR)
	$(CC) $(CFLAGS) $(EX1_SRC) $(SRC_DIR)/plasma_params.c $(LDFLAGS) -o $@

$(EX2_BIN): $(SRC_DIR)/plasma_params.c $(SRC_DIR)/wakefield.c $(EX2_SRC) | $(BLD_DIR)
	$(CC) $(CFLAGS) $(EX2_SRC) $(SRC_DIR)/plasma_params.c $(SRC_DIR)/wakefield.c $(LDFLAGS) -o $@

$(EX3_BIN): $(SRC_DIR)/plasma_params.c $(SRC_DIR)/instabilities.c $(SRC_DIR)/laser_propagation.c $(EX3_SRC) | $(BLD_DIR)
	$(CC) $(CFLAGS) $(EX3_SRC) $(SRC_DIR)/plasma_params.c $(SRC_DIR)/instabilities.c $(SRC_DIR)/laser_propagation.c $(LDFLAGS) -o $@

examples: $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)
	$(EX1_BIN)
	$(EX2_BIN)
	$(EX3_BIN)

clean:
	rm -rf $(BLD_DIR)
