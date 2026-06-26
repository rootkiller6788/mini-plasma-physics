CC       = gcc
CFLAGS   = -std=c99 -Wall -Wextra -O2 -I. -D_DEFAULT_SOURCE
LDFLAGS  = -lm

SRC_DIR  = src
INC_DIR  = include
TST_DIR  = tests
EXA_DIR  = examples
BEN_DIR  = benches

LIBSRCS  = $(SRC_DIR)/fusion_core.c $(SRC_DIR)/fusion_confinement.c \
           $(SRC_DIR)/fusion_mhd.c $(SRC_DIR)/fusion_heating.c \
           $(SRC_DIR)/fusion_equilibrium.c $(SRC_DIR)/fusion_reactor.c

LIBOBJS  = $(LIBSRCS:.c=.o)
LIBNAME  = libfusionplasma.a

TESTSRCS = $(wildcard $(TST_DIR)/*.c)
EXASRCS  = $(wildcard $(EXA_DIR)/*.c)

.PHONY: all clean test examples lib count audit

all: lib

lib: $(LIBNAME)

$(LIBNAME): $(LIBOBJS)
	ar rcs $@ $^

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(wildcard $(INC_DIR)/*.h)
	$(CC) $(CFLAGS) -c $< -o $@

test: lib
	@for f in $(TESTSRCS); do \
		name=$$(basename $$f .c); \
		echo "=== Testing $$name ==="; \
		$(CC) $(CFLAGS) $$f -L. -lfusionplasma $(LDFLAGS) -o $(TST_DIR)/$$name.exe; \
		$(TST_DIR)/$$name.exe; \
	done

examples: lib
	@for f in $(EXASRCS); do \
		name=$$(basename $$f .c); \
		echo "=== Running $$name ==="; \
		$(CC) $(CFLAGS) $$f -L. -lfusionplasma $(LDFLAGS) -o $(EXA_DIR)/$$name.exe; \
		$(EXA_DIR)/$$name.exe; \
	done

clean:
	rm -f $(SRC_DIR)/*.o
	rm -f $(LIBNAME)
	rm -f $(TST_DIR)/*.exe $(EXA_DIR)/*.exe $(BEN_DIR)/*.exe

count:
	@echo Header files: && wc -l $(INC_DIR)/*.h
	@echo Source files: && wc -l $(SRC_DIR)/*.c
	@echo ---
	@find . \( -path "./include/*.h" -o -path "./src/*.c" \) -exec cat {} + | wc -l | xargs echo Total:

audit:
	@echo === Filler Detection ===
	@grep -rn _fn[0-9] $(SRC_DIR)/*.c $(INC_DIR)/*.h || echo PASS
	@grep -rn _aux[0-9] $(SRC_DIR)/*.c $(INC_DIR)/*.h || echo PASS
	@grep -rn "algorithm variant" $(SRC_DIR)/*.c || echo PASS
	@echo === TODO/FIXME Check ===
	@grep -rn "TODO\|FIXME\|stub\|placeholder" $(SRC_DIR)/*.c $(INC_DIR)/*.h || echo PASS
