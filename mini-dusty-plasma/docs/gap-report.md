# Gap Report -- mini-dusty-plasma

## Current Gaps

### L9 Research Frontiers (currently Partial)

| Gap | Priority | Effort Estimate |
|-----|----------|-----------------|
| Dust in atmospheric pressure plasmas | Low | 2 days |
| Time crystal analogs in dust plasma | Low | 3 days |
| Dust-plasma interaction in hypersonic flows | Low | 5 days |

## Future Extensions

1. **Ewald summation for Yukawa**: Replace O(N^2) energy calculation
   with O(N log N) Ewald method for large crystals.
2. **Particle-in-Cell (PIC) with dust**: Full kinetic simulation of
   dusty plasma including dust charging on-the-fly.
3. **Machine learning force fields**: Train ML potentials on MD
   trajectories of Yukawa OCP.
4. **Dust growth/coagulation module**: Smoluchowski equation solver
   for dust size distribution evolution.

## No Missing Items in L1-L8

All levels L1 through L8 are fully covered with:
- C implementations in include/ and src/
- Lean 4 formal definitions in src/dusty_plasma.lean
- Tests in tests/test_dusty.c
- Examples in examples/