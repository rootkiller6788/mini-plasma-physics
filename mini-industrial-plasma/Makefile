CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99 -Iinclude

.PHONY: all clean test examples dirs

all: dirs libindustrialplasma.a

dirs:
	mkdir -p build

OBJ = build/sheath.o build/discharge.o build/eedf.o build/surface.o build/plasma_params.o build/plasma_chemistry.o

libindustrialplasma.a: $(OBJ)
	ar rcs build/libindustrialplasma.a $(OBJ)

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: dirs build/test_sheath build/test_discharge build/test_eedf build/test_surface build/test_plasma build/test_chemistry
	@echo === Running tests ===
	./build/test_sheath && ./build/test_discharge && ./build/test_eedf && ./build/test_surface && ./build/test_plasma && ./build/test_chemistry && echo ALL TESTS PASSED

examples: dirs build/ex1_ccp_etch build/ex2_pecvd build/ex3_paschen

build/test_sheath: tests/test_sheath.c libindustrialplasma.a
	$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/test_discharge: tests/test_discharge.c libindustrialplasma.a
	$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/test_eedf: tests/test_eedf.c libindustrialplasma.a
	$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/test_surface: tests/test_surface.c libindustrialplasma.a
	$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/test_plasma: tests/test_plasma.c libindustrialplasma.a
	$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/test_chemistry: tests/test_chemistry.c libindustrialplasma.a
	$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/ex1_ccp_etch: examples/ex1_ccp_etch.c libindustrialplasma.a
	$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/ex2_pecvd: examples/ex2_pecvd_deposition.c libindustrialplasma.a
	$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

build/ex3_paschen: examples/ex3_paschen_curve.c libindustrialplasma.a
	$(CC) $(CFLAGS) $< -Lbuild -lindustrialplasma -lm -o $@

clean:
	rm -rf build

bench: dirs
	$(CC) $(CFLAGS) -O3 benches/bench_eedf.c src/eedf.c src/plasma_params.c src/sheath.c -lm -o build/bench_eedf
	./build/bench_eedf
