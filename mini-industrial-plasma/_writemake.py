
import os
base = r"F:\nano-everything\mini-pure-physics\8. mini-plasma-physics\mini-industrial-plasma"
mf = """CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99 -Iinclude
LDFLAGS = -lm
SRC = src/sheath.c src/discharge.c src/eedf.c src/surface.c src/plasma_params.c src/plasma_chemistry.c
OBJ = build/sheath.o build/discharge.o build/eedf.o build/surface.o build/plasma_params.o build/plasma_chemistry.o

.PHONY: all clean test examples dirs

all: dirs libindustrialplasma.a

dirs:
\tmkdir -p build

libindustrialplasma.a: $(OBJ)
\tar rcs build/libindustrialplasma.a $(OBJ)

build/%.o: src/%.c
\t$(CC) $(CFLAGS) -c $< -o $@

test: dirs $(TEST_BIN)
\t@echo === Running tests ===
\t./build/test_sheath
\t./build/test_discharge
\t./build/test_eedf
\t./build/test_surface
\t./build/test_plasma
\t./build/test_chemistry
\t@echo === All tests passed ===

examples: dirs build/ex1_ccp_etch build/ex2_pecvd build/ex3_paschen

build/test_sheath: tests/test_sheath.c libindustrialplasma.a
\t$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/test_discharge: tests/test_discharge.c libindustrialplasma.a
\t$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/test_eedf: tests/test_eedf.c libindustrialplasma.a
\t$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/test_surface: tests/test_surface.c libindustrialplasma.a
\t$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/test_plasma: tests/test_plasma.c libindustrialplasma.a
\t$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/test_chemistry: tests/test_chemistry.c libindustrialplasma.a
\t$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/ex1_ccp_etch: examples/ex1_ccp_etch.c libindustrialplasma.a
\t$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/ex2_pecvd: examples/ex2_pecvd_deposition.c libindustrialplasma.a
\t$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/ex3_paschen: examples/ex3_paschen_curve.c libindustrialplasma.a
\t$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

clean:
\trm -rf build

bench: dirs
\t$(CC) $(CFLAGS) -O3 benches/bench_eedf.c src/eedf.c src/plasma_params.c src/sheath.c -lm -o build/bench_eedf
\t./build/bench_eedf
"""
with open(os.path.join(base, "Makefile"), "w") as f:
    f.write(mf)
print("Makefile written")
