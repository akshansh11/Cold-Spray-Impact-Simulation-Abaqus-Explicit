# Cold Spray Impact Simulation — Abaqus/Explicit

A fully automated Abaqus Python script that models a single aluminium particle impacting an aluminium substrate at cold spray velocities. The simulation uses a coupled thermomechanical explicit formulation with Johnson-Cook plasticity and a Mie-Gruneisen (Us-Up) equation of state.


<img width="1702" height="793" alt="cold" src="https://github.com/user-attachments/assets/2f09b67d-24ed-45df-881f-f785f372e73b" />

---

## What It Simulates

```
        Particle  (R = 40 um, Al)
        |  v = -700 m/s (downward)
        v
 +--------------------------+
 |   Substrate              |  (diameter 500 um x 250 um deep, Al)
 |   (fixed base)           |
 +--------------------------+
```

- High-velocity particle deformation and bonding zone
- Adiabatic shear instability at the particle-substrate interface
- Temperature rise driven by plastic dissipation
- Contact mechanics with Coulomb friction (mu = 0.3)

---

## Repository Structure

```
.
├── cold_spray_simulation.py   # Main Abaqus script
└── README.md
```

---

## Model Parameters

| Parameter             | Value                  |
|-----------------------|------------------------|
| Particle radius       | 40 um                  |
| Substrate diameter    | 500 um                 |
| Substrate depth       | 250 um                 |
| Impact velocity       | -700 m/s (Y-direction) |
| Initial temperature   | 298 K (both parts)     |
| Step duration         | 60 ns                  |
| Output frames         | 30                     |
| CPUs                  | 4                      |

---

## Material — Aluminium (Both Parts)

### Mechanical and Thermal Properties

| Property             | Value         |
|----------------------|---------------|
| Density              | 2700 kg/m3    |
| Shear modulus        | 27 GPa        |
| Thermal conductivity | 237.2 W/m·K   |
| Specific heat        | 898.2 J/kg·K  |

### Johnson-Cook Plasticity

```
sigma = (A + B * eps^n) * (1 + C * ln(eps_dot*)) * (1 - T*^m)
```

| A (MPa) | B (MPa) | n     | m     | T_melt (K) | T_ref (K) | C     | eps_dot_ref |
|---------|---------|-------|-------|------------|-----------|-------|-------------|
| 148.4   | 345.5   | 0.183 | 0.895 | 916        | 298       | 0.001 | 1.0         |

### Mie-Gruneisen Equation of State (Us-Up)

| c0 (m/s) | s     | Gamma_0 |
|----------|-------|---------|
| 5386     | 1.339 | 1.97    |

---

## Mesh Strategy

### Particle

- Partitioned into 8 octants using XY, YZ, and XZ datum planes
- Uniform seed size: 2.5 um
- Element type: C3D8RT (thermomechanical hex, Explicit)

### Substrate

- Circular refinement zone on top face aligned with particle footprint
- Global seed size: 35 um; impact-zone edge seed: 3.5 um
- Element type: C3D8RT (thermomechanical hex, Explicit)

---

## Contact Definition

| Setting              | Value                          |
|----------------------|--------------------------------|
| Formulation          | Surface-to-surface, kinematic  |
| Sliding              | Finite                         |
| Normal behaviour     | Hard contact, no separation    |
| Tangential behaviour | Penalty friction, mu = 0.3     |

---

## Boundary Conditions and Predefined Fields

| Name                | Region            | Type                  | Value    |
|---------------------|-------------------|-----------------------|----------|
| BC-1                | Substrate base/OD | Encastre              | —        |
| Predefined Field-1  | Particle          | Initial velocity (Y)  | -700 m/s |
| Predefined Field-2  | Particle          | Uniform temperature   | 298 K    |
| Predefined Field-3  | Substrate         | Uniform temperature   | 298 K    |

---

## Field Output Variables

```
S  SVAVG  PE  PEVAVG  PEEQ  PEEQVAVG  LE  U  V  A  RF
CSTRESS  NT  TEMP  HFL  RFL  EVF  STATUS
```

---

## Requirements

- Abaqus 2020 or later (tested on 2025)
- Abaqus/Explicit solver
- 4 CPU cores recommended

---

## Usage

Run without the GUI from the command line:

```bash
abaqus cae noGUI=cold_spray_simulation.py
```

Or from within Abaqus CAE:

```
File > Run Script > cold_spray_simulation.py
```

The script builds the model, submits the job, and waits for completion. Results are written to `cs1.odb` in the working directory.

---

## Viewing Results

Open the results database in Abaqus/Viewer:

```bash
abaqus viewer database=cs1.odb
```

Recommended contour plots:

| Variable  | Description                              |
|-----------|------------------------------------------|
| PEEQ      | Equivalent plastic strain (bonding zone) |
| TEMP      | Temperature (adiabatic shear bands)      |
| S, Mises  | Von Mises stress                         |
| STATUS    | Element deletion / failure               |

---

## Implementation Notes

- ALE adaptive meshing is applied to the particle and inner substrate cells to handle large deformations without excessive element distortion.
- `elemDeletion=ON` allows failed elements to be removed if distortion limits are exceeded.
- `waitForCompletion()` ensures the script does not exit before the solver finishes.
- Both parts share the same Aluminum material definition; section assignment is handled in a loop.

---

## License

[![CC BY-NC 4.0](https://licensebuttons.net/l/by-nc/4.0/88x31.png)](https://creativecommons.org/licenses/by-nc/4.0/)

This project is licensed under the
[Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0)](https://creativecommons.org/licenses/by-nc/4.0/).

You are free to:
- Share — copy and redistribute the material in any medium or format
- Adapt — remix, transform, and build upon the material

Under the following terms:
- Attribution — you must give appropriate credit to the original author
- NonCommercial — you may not use the material for commercial purposes

Copyright (c) 2026 Akshansh Mishra ([@akshansh11](https://github.com/akshansh11))
