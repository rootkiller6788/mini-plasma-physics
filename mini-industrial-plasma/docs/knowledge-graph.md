# Knowledge Graph - mini-industrial-plasma
## L1: Definitions (Complete)
- PlasmaState, PlasmaSpecies, PlasmaReaction, CrossSectionModel
- DischargeGeometry, DischargeParams
- TransportCoeffs, RateMatrix, EEDFState
- SheathSolution, PaschenCurve
- EtchModel, DepositionModel, ElectronegativeParams
- Maxwellian EEDF, Druyvesteyn EEDF
- Sputtering yield definitions (Sigmund, Bohdansky, Yamamura)

## L2: Core Concepts (Complete)
- Bohm criterion for sheath formation
- Debye shielding and plasma quasineutrality
- Townsend discharge and avalanche breakdown
- Paschen curve and minimum breakdown
- DC discharge regimes (Townsend, glow, arc)
- Capacitively Coupled Plasma (CCP)
- Inductively Coupled Plasma (ICP)
- Ion-enhanced etching (Coburn-Winters synergism)
- Bi-Maxwellian EEDF (bulk + tail populations)
- Plasma beta and magnetization

## L3: Mathematical Structures (Complete)
- Rate matrix for chemical kinetics
- Cross section parametrization models
- Poisson equation for plasma sheaths
- Child-Langmuir space-charge law
- Paschen breakdown condition
- EEDF normalization integrals (Gamma functions)
- Generalized EEDF with shape parameter

## L4: Fundamental Laws (Complete)
- Child-Langmuir law: J = (4/9)*eps0*sqrt(2e/m)*V^{3/2}/d^2
- Bohm criterion: v_i >= c_s = sqrt(e*T_e/m_i)
- Floating potential: V_f = -(T_e/2)*ln(m_i/(2*pi*m_e))
- Paschen law: V_br = B*pd/ln(A*pd/ln(1+1/gamma))
- Townsend criterion: gamma*(exp(alpha*d)-1) = 1
- Ion impact energy with collisions
- Sheath capacitance nonlinearity

## L5: Computational Methods (Complete)
- RK4 ODE solver for planar sheath
- Simpson integration for rate coefficients
- Two-term Boltzmann equation iterative solver
- Euler integration for 0-D chemical kinetics
- Global (volume-averaged) discharge model
- Adaptive step size control for stiff ODEs
- Matrix (uniform ion density) sheath model

## L6: Canonical Systems (Complete)
- DC glow discharge characteristic
- CCP asymmetric discharge with self-bias
- ICP with transformer model and E-H transition
- Argon plasma chemistry (7 reactions, 5 species)
- SF6 etching chemistry (5 reactions, 10 species)
- SiH4/H2 PECVD chemistry (5 reactions, 7 species)
- O2 electronegative plasma chemistry (6 reactions, 7 species)
- Paschen curve computation and analysis

## L7: Applications (Complete)
- RIE process model for semiconductor etching
- PECVD a-Si:H deposition with film quality prediction
- Wafer-scale uniformity analysis
- Industrial process parameter database (RIE, PECVD, ICP, sputter)
- RF discharge density estimation from engineering parameters
- Scaling laws for discharge similarity

## L8: Advanced Topics (Partial)
- Electronegative plasma chemistry (O2)
- Atomic Layer Etching (ALE) model
- Non-local EEDF effects in low-pressure ICP
- Generalized EEDF formalism with shape parameter

## L9: Research Frontiers (Partial)
- ALE synergy parameter analysis
- Plasma medicine concepts (documented)
- CO2 conversion in plasmas (documented)
