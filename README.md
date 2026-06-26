# Mini Plasma Physics

A collection of **from-scratch, zero-dependency C implementations** of university-level plasma physics and controlled fusion. Each module maps to MIT, Princeton, Stanford, Cambridge, and other top-tier university courses, covering the full spectrum from fundamental plasma theory (Debye shielding, kinetic theory, MHD) to specialized applications in fusion energy, space physics, laser-plasma interactions, industrial processing, dusty (complex) plasmas, and wave-instability theory.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|--------|--------|-------------|
| [mini-plasma-fundamental](mini-plasma-fundamental/) | Debye shielding, plasma frequency/parameters, kinetic theory (Vlasov–Fokker–Planck), single-fluid/two-fluid MHD, particle simulation (Boris pusher, PIC), plasma wave dispersion | MIT 22.611, Princeton PHY 521, Berkeley PHYS 242 |
| [mini-mhd](mini-mhd/) | Ideal/resistive MHD equations, MHD equilibria (Grad–Shafranov, Z-pinch, theta-pinch), MHD instabilities (kink, sausage, interchange, ballooning), MHD wave modes (Alfvén, fast/slow magnetosonic), MHD numerics (Riemann solvers, constrained transport, ∇·B cleaning) | MIT 22.611/22.615, Princeton AST 554, Cambridge Part III MHD |
| [mini-waves-instabilities](mini-waves-instabilities/) | Dispersion relation solvers (Brent, complex Newton, Nyquist), kinetic dispersion (Landau damping, cyclotron damping, Bernstein modes), nonlinear waves (three-wave coupling, Zakharov equations, solitons), instability catalog (two-stream, Buneman, ion-acoustic, drift-wave) | MIT 22.611, Princeton PHY 521, Berkeley PHYS 242 |
| [mini-dusty-plasma](mini-dusty-plasma/) | Dust grain charging (OML theory), forces on dust (ion drag, neutral drag, thermophoretic, Lorentz), dust crystal formation (Coulomb coupling, phase transitions), dusty plasma wave modes (DAW, DIAW, DLW), dusty transport coefficients | MIT 22.611/22.612, Princeton PHY 535 |
| [mini-fusion-plasma](mini-fusion-plasma/) | Magnetic confinement (tokamak geometry, flux surfaces, Shafranov shift), MHD equilibrium/stability (Grad–Shafranov, tearing modes, disruptions), heating & current drive (NBI, ICRF, ECRH, LH), fusion transport (neoclassical, turbulent), Lawson criterion and power balance, ITER/DEMO/SPARC operating points | MIT 22.611/22.615, Princeton PHY 563, Cambridge Part III Fusion |
| [mini-space-plasma](mini-space-plasma/) | Planetary magnetospheres (Chapman–Ferraro, Dungey reconnection), solar wind (Parker model, heliospheric structure), MHD in space contexts, plasma waves in space (whistler, chorus, ULF), space weather parameters | MIT 22.611/22.616, Cambridge Part III Plasma Physics |
| [mini-laser-plasma](mini-laser-plasma/) | Laser absorption mechanisms (inverse bremsstrahlung, resonance absorption, Brunel heating), parametric instabilities (SRS, SBS, TPD, filamentation), field/collisional ionization (BSI, ADK, avalanche), relativistic self-focusing, wakefield acceleration (LWFA), single-particle motion in laser fields (Boris integrator) | MIT 22.611, Princeton PHY 525, Stanford PHYSICS 370 |
| [mini-industrial-plasma](mini-industrial-plasma/) | RF discharge models (CCP/ICP), electron energy distribution functions (Maxwellian, Druyvesteyn), plasma sheath (Bohm criterion, Child–Langmuir), plasma–surface interactions (etching, PECVD, sputtering), collision cross-sections and rate coefficients | MIT 22.611, Stanford EE 414, Berkeley EECS 245 |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `Makefile`, `include/`, `src/`, `examples/`, `demos/`, `tests/`
- **Theory-to-code mapping** — every module includes `docs/` with course-alignment notes and reference derivations
- **Practical demos** — Grad–Shafranov equilibrium solvers, CMA diagrams, particle pushers, sheath models, wave dispersion roots, instability growth-rate calculators, and more

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-plasma-fundamental
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-plasma-physics/
├── mini-plasma-fundamental/    # Fundamental Plasma Physics (parameters, kinetic, MHD, particles, waves)
├── mini-mhd/                   # Magnetohydrodynamics (equations, equilibria, instabilities, numerics, waves)
├── mini-waves-instabilities/   # Plasma Waves & Instabilities (dispersion solvers, kinetic, nonlinear, catalog)
├── mini-dusty-plasma/          # Dusty (Complex) Plasma Physics (charging, forces, crystals, waves, transport)
├── mini-fusion-plasma/         # Fusion Plasma Physics (confinement, equilibrium, heating, MHD, transport)
├── mini-space-plasma/          # Space Plasma Physics (magnetosphere, solar wind, MHD, waves)
├── mini-laser-plasma/          # Laser Plasma Interactions (absorption, instabilities, ionization, wakefield)
└── mini-industrial-plasma/     # Industrial Plasma Processing (discharges, EEDF, sheath, surface)

```

## License

MIT
