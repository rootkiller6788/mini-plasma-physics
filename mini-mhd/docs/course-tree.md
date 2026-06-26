# Course Dependency Tree — mini-mhd

## Prerequisites (What You Need Before MHD)

```
mini-classical-mechanics (0)
  ├── Lagrangian/Hamiltonian mechanics
  ├── Conservation laws (Noether)
  └── Phase space

mini-electromagnetism (1)
  ├── Maxwell's equations
  ├── Poynting vector
  ├── Maxwell stress tensor
  └── Ohm's law

mini-fluid-dynamics (6)
  ├── Navier-Stokes equations
  ├── Compressible flow
  ├── Sound waves
  └── Shocks and discontinuities

mini-thermodynamics-statphys (2)
  ├── Ideal gas law
  ├── Adiabatic processes
  └── Temperature and pressure
```

## MHD Knowledge Map

```
MHD Fundamentals
├── MHD Equations (L4)
│   ├── Continuity: d_t(rho) + div(rho v) = 0
│   ├── Momentum:  d_t(rho v) + div(...) = J x B - grad p
│   ├── Induction:  d_t(B) = curl(v x B) + eta Lap(B)
│   ├── Energy:     d_t(E) + div(F_E) = eta J^2
│   └── div(B) = 0 constraint
├── MHD Approximations (L2)
│   ├── Quasineutrality (lambda_D << L)
│   ├── Non-relativistic (v << c, no displacement current)
│   ├── Small gyroradius (rho_i << L)
│   └── Collisional (MHD closure)
├── MHD Equilibria (L6)
│   ├── Force balance: J x B = grad p
│   ├── 1D pinches (Z, theta, screw)
│   ├── 2D axisymmetric (Grad-Shafranov)
│   └── Force-free fields (curl B = alpha B)
├── MHD Waves (L2, L6)
│   ├── Alfven wave (incompressible, transverse)
│   ├── Fast magnetosonic (compressible, fast)
│   └── Slow magnetosonic (compressible, slow)
├── MHD Stability (L4, L6, L8)
│   ├── Energy principle (delta_W)
│   ├── Ideal instabilities (kink, sausage, interchange)
│   ├── Resistive instabilities (tearing)
│   └── Ballooning modes
├── MHD Invariants (L4)
│   ├── Magnetic helicity H_M = int(A.B dV)
│   ├── Cross helicity H_C = int(v.B dV)
│   └── Total energy
└── Computational MHD (L5)
    ├── Finite volume method
    ├── Riemann solvers (HLL, HLLD, Roe)
    ├── Constrained transport (div(B)=0)
    ├── Time integration (RK2/3/4)
    └── Boundary conditions
```

## Postrequisites (What Builds On MHD)

```
mini-plasma-physics (8)
  ├── mini-mhd (this module) ─────────┐
  ├── kinetic theory                   │
  └── waves and instabilities          │
                                       ▼
Extended MHD ─────────────────────────────
├── Hall MHD (d_i effects)
├── Two-fluid MHD
├── Braginskii equations
└── Gyrokinetics

Fusion Plasma Physics ────────────────────
├── Tokamak physics
├── Stellarator optimization
└── ITER/DEMO

Astrophysical MHD ────────────────────────
├── Solar corona and wind
├── Accretion disks
├── Jets and outflows
└── ISM turbulence

Research Frontiers (L9) ─────────────────
├── Magnetic reconnection (fast)
├── Dynamo theory
├── Relativistic MHD
└── Laboratory astrophysics
```
