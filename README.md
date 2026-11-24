<div align="center">

# Hi there, I'm Emre Sakarya 👋

[![Typing SVG](https://readme-typing-svg.demolab.com?font=Fira+Code&size=18&pause=1000&color=C0392B&center=true&vCenter=true&width=650&lines=Nuclear+Engineering+%40+Hacettepe+University;Reactor+Physics+%7C+Thermal-Hydraulics+%7C+CFD;Python+%7C+Fortran+%7C+ANSYS+Fluent+%7C+Star-CCM%2B;Computational+Engineering+%26+Simulation)](https://git.io/typing-svg)

I'm a senior **Nuclear Engineering** student at Hacettepe University working on computational engineering, simulation and numerical analysis.
I write the solvers myself and check them against benchmarks, from reactor kinetics and core physics to full-bundle CFD.

[![Portfolio](https://img.shields.io/badge/Portfolio-emresakarya.github.io-C0392B?style=for-the-badge&logo=githubpages&logoColor=white)](https://emresakarya.github.io/)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-0077B5?style=for-the-badge&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/emre-sakarya-ab8967300/)
[![Email](https://img.shields.io/badge/Email-D14836?style=for-the-badge&logo=gmail&logoColor=white)](mailto:emresakarya@hacettepe.edu.tr)

</div>

---

## ⭐ Flagship: CoreForge

**Interactive 2-D / 3-D multigroup reactor core simulator in the browser** · [▶ Live demo](https://coreforge.streamlit.app)

A core-physics code I built and validated. A Fortran engine solves the multigroup neutron diffusion eigenvalue problem for a whole core (2-D or full x-y-z) and returns k<sub>eff</sub>, flux and power maps, peaking factors, rod worths and critical boron; a Streamlit web app sits on top.

| | |
|:--|:--|
| **Engine** | Fortran · mesh-centred finite differences · power iteration + red-black SOR · OpenMP |
| **Cross sections** | Published few-group sets, a built-in fuel designer, or DRAGON5 lattice tables (172 groups, JEF-2.2) at 11 enrichments, 1.6–4.95 % |
| **Validation** | Seven published cores (IAEA-2D/3D, BIBLIS-2D, KOEBERG, LRA, TWIGL, C5G7), 3 analytic cases, 39 automated checks; IAEA-3D within **−6.6 pcm** |
| **Measured reactor** | BEAVRS Cycle 1 HZP: **+592 pcm** against the measured critical boron, nothing tuned |
| **Transients** | Point kinetics with Doppler and moderator feedback: rod ejection, ramp, scram |
| **Interface** | Bilingual (TR/EN) web app: core wizard, live rod control, standalone HTML reports |

Since v12 the scope is neutronics only: depletion and thermal-hydraulics were removed because no benchmark backed them. The companion [C5G7 S<sub>N</sub> transport solver](https://github.com/EmreSakarya/c5g7-2d-transport-benchmark) covers the heterogeneous transport side.

---

## 🛠️ Technical Stack

<table>
<tr>
<td valign="top" width="33%">

**Languages**

![Python](https://img.shields.io/badge/Python-3776AB?style=flat&logo=python&logoColor=white) Advanced\
![Fortran](https://img.shields.io/badge/Fortran-734F96?style=flat&logo=fortran&logoColor=white) Intel ifx / gfortran\
![C++](https://img.shields.io/badge/C%2B%2B-00599C?style=flat&logo=c%2B%2B&logoColor=white) Intermediate\
![C](https://img.shields.io/badge/C-A8B9CC?style=flat&logo=c&logoColor=black) Intermediate

</td>
<td valign="top" width="33%">

**Simulation Tools**

![ANSYS](https://img.shields.io/badge/ANSYS_Fluent-FFB71B?style=flat&logo=ansys&logoColor=black)
![Star-CCM+](https://img.shields.io/badge/Star--CCM%2B-CC0000?style=flat&logoColor=white)\
![DRAGON5](https://img.shields.io/badge/DRAGON5-734F96?style=flat&logoColor=white)
![KORIGEN](https://img.shields.io/badge/KORIGEN-734F96?style=flat&logoColor=white)\
![Streamlit](https://img.shields.io/badge/Streamlit-FF4B4B?style=flat&logo=streamlit&logoColor=white)
![ImageJ](https://img.shields.io/badge/ImageJ%20%2F%20Fiji-57606a?style=flat&logoColor=white)\
![Linux](https://img.shields.io/badge/Linux%2FBash-FCC624?style=flat&logo=linux&logoColor=black)

</td>
<td valign="top" width="33%">

**Numerical Methods & Libraries**

🧮 NumPy · SciPy · Pandas\
📊 Matplotlib · Plotly\
💧 PyXSteam (IAPWS-97)\
⚛️ Monte Carlo · FDM · FVM · S<sub>N</sub> · RANS

</td>
</tr>
</table>

**Domain Expertise:** Neutron Diffusion & Transport · Core Benchmarks · Point Reactor Kinetics · Source Term & Atmospheric Dispersion · Two-Phase Flow · CFD Turbulence Modeling (k-ε, k-ω SST) · Thermodynamic Cycle Analysis

---

## ⚛️ Reactor Physics & Neutronics · Fortran / Python

> Diffusion and transport solvers written from scratch and checked against published benchmarks or exact solutions.

<table>
<tr>
<td width="50%" valign="top">

#### [IAEA-3D PWR Benchmark · Neutron Diffusion Solver](https://github.com/EmreSakarya/iaea-3d-pwr-benchmark)
`Fortran` `OpenMP` `Finite Difference` `Eigenvalue` `Red-Black SOR`

From-scratch two-group neutron diffusion solver for the 3-D IAEA PWR core benchmark (ANL-7416). Mesh-centered FDM with harmonic-mean interface coefficients, Robin vacuum BCs, and power iteration + red-black SOR. Reproduces the reference **k_eff = 1.02903** to within **−7.9 pcm** and matches the 1977 VENTURE mesh-convergence sequence.

</td>
<td width="50%" valign="top">

#### [C5G7 · 2-D 7-Group Neutron Transport Benchmark (S<sub>N</sub>)](https://github.com/EmreSakarya/c5g7-2d-transport-benchmark)
`Fortran` `OpenMP` `Discrete Ordinates (S_N)` `7-Group Transport` `Eigenvalue`

From-scratch discrete-ordinates (S<sub>N</sub>) neutron transport solver for the OECD/NEA **C5G7 MOX** benchmark: heterogeneous pin-by-pin UO₂ + MOX geometry with no spatial homogenisation. Product angular quadrature, diamond-difference sweep, volume-preserving "digital-disk" pin discretisation, reflective/vacuum BCs, power iteration with OpenMP. Reproduces the MCNP reference **k_eff = 1.18655** to within **−182 pcm**; max/min pin powers match to ~6%.

</td>
</tr>
<tr>
<td width="50%" valign="top">

#### 🆕 [2-D One-Speed Neutron Diffusion by Finite Differences](https://github.com/EmreSakarya/2d-neutron-diffusion-fdm)
`Python` `SciPy Sparse` `Finite Difference` `Ghost Nodes` `Analytic Verification`

Five-point FDM for steady one-speed diffusion in a 30 × 24 cm non-multiplying region with vacuum and reflective (ghost-node) boundaries. Within **0.11%** of a double Fourier series solution, clean second-order convergence, and a factor of **3.75** rise in peak flux from two reflective sides. [Report](https://github.com/EmreSakarya/2d-neutron-diffusion-fdm/blob/main/docs/NEM493_Project1_Report.pdf)

</td>
<td width="50%" valign="top">

#### [Point Reactor Kinetics Analysis](https://github.com/EmreSakarya/point-reactor-kinetics-analysis)
`Python` `RK4` `Heun's Method` `Stiff Systems`

Point reactor kinetics with delayed neutrons under reactivity insertions. Implemented **RK4** and **Heun's Predictor-Corrector** integrators; stiffness-ratio analysis against analytical solutions.

</td>
</tr>
<tr>
<td width="50%" valign="top">

#### [U-238 Decay Chain Modeling](https://github.com/EmreSakarya/uranium-decay-chain-analysis)
`Fortran 90` `Implicit Euler` `Stiff ODEs`

14-stage decay chain (U-238 → Pb-206) with half-lives spanning 10⁹ years to seconds, a severely **stiff system**. Generic Implicit Euler solver in Fortran, mirroring legacy HPC nuclear-code logic.

</td>
<td width="50%" valign="top"></td>
</tr>
</table>

---

## 🎲 Monte Carlo & Particle Transport

<table>
<tr>
<td width="50%" valign="top">

#### 🆕 [One-Speed Transport in a Three-Region Slab (Analytical + Monte Carlo)](https://github.com/EmreSakarya/three-region-slab-monte-carlo)
`Python` `Monte Carlo` `Collision Estimator` `Exponential Integrals`

Analog Monte Carlo with 10⁶ histories, verified against an E<sub>n</sub> solution: probabilities within 1σ, flux within **1.32%** in all 20 intervals. Isotropic scattering halves the total leakage, from 0.08692 to **0.04333**. [Report](https://github.com/EmreSakarya/three-region-slab-monte-carlo/blob/main/docs/NEM493_Project2_Report.pdf)

</td>
<td width="50%" valign="top">

#### [Stochastic Neutron Transport Simulation](https://github.com/EmreSakarya/neutron-transport-monte-carlo)
`Monte Carlo` `Statistical Analysis` `Slab Geometry`

Stochastic neutron-trajectory algorithm from scratch in slab geometry. Validated against analytical **Exponential Integral** functions; demonstrated 1/√N error scaling.

</td>
</tr>
<tr>
<td width="50%" valign="top">

#### [Charged Particle Energy Loss Simulation](https://github.com/EmreSakarya/bethe-bloch-analysis)
`Bethe-Bloch Formula` `Adaptive Integration`

Energy deposition of protons and alpha particles in water / hydrogen media. Simulated **phase-transition effects** on stopping power using adaptive numerical integration.

</td>
<td width="50%" valign="top"></td>
</tr>
</table>

---

## ☢️ Source Term & Environmental Dispersion · CAREM-25 SMR

> From the end-of-cycle core inventory to off-site air concentration: the KORIGEN result is the source term of the plume study.

<table>
<tr>
<td width="50%" valign="top">

#### [CAREM-25 End-of-Cycle Activity (KORIGEN)](https://github.com/EmreSakarya/carem-25-korigen-activity)
`KORIGEN` `Depletion & Decay` `Source Term` `Accident-Tolerant Fuel`

End-of-cycle radionuclide inventory of the CAREM-25 SMR from the VVER-1000 KORIGEN deck adapted to CAREM geometry (6588 rods, 26.04 MW/tHM, 24 MWd/kg). Core total **2.12 × 10¹⁹ Bq**, plus five variations: low burnup, 20 w/o enrichment, Th-232, U₃Si₂ fuel and FeCrAl cladding.

</td>
<td width="50%" valign="top">

#### 🆕 [Gaussian Plume Dispersion of a CAREM-25 Release](https://github.com/EmreSakarya/carem-25-gaussian-plume)
`Python` `Gaussian Plume` `Pasquill–Gifford` `Pandas`

1% of the CAREM-25 EOC inventory released over one day, followed through 24 hourly weather conditions on a 41 × 41 grid with and without ground reflection. 24-hour time-integrated peak of **4.45 × 10¹¹ Bq·s/m³**, 1.12 km SSE of the reactor. [Report](https://github.com/EmreSakarya/carem-25-gaussian-plume/blob/main/docs/NEM467_Application2_Report.pdf)

</td>
</tr>
</table>

---

## 🌡️ Heat Transfer & Thermal-Hydraulics · Python

<table>
<tr>
<td width="50%" valign="top">

#### [2D Transient Thermal Analysis of a Rectangular Plate](https://github.com/EmreSakarya/2d-transient-thermal-analysis)
`Python` `Control Volume FDM` `Explicit Euler` `Natural Convection`

Two-stage heat treatment on a 3×5 nodal grid using **CVFDM** with temperature-dependent properties. Heating under adiabatic boundaries → **1000 K**; cooling via Ra-based natural convection in a liquid-metal bath → **600 K** (~1920 s).

</td>
<td width="50%" valign="top">

#### [Thermal-Hydraulic Analysis (HEM)](https://github.com/EmreSakarya/vertical-pipe-flow-analysis)
`Python` `Finite Volume Method` `Two-Phase Flow` `IAPWS-97`

1D solver for downward water flow under constant heat flux using the **Homogeneous Equilibrium Model**. Predicts boiling inception and void-fraction profiles with IAPWS-97 steam tables.

</td>
</tr>
</table>

---

## 🌊 CFD Simulations · ANSYS Fluent & Star-CCM+

> CFD projects (NEM 358) cross-validated between two commercial solvers and analytical correlations.

<table>
<tr>
<td width="50%" valign="top">

#### [CFD Simulation of a 3×3 PWR Rod Bundle](https://github.com/EmreSakarya/pwr-rod-bundle-cfd)
`ANSYS Fluent` `Star-CCM+` `Conjugate Heat Transfer` `k-ω SST`

Full conjugate heat transfer simulation of a 3×3 PWR fuel rod bundle. Applied a **non-uniform radial power distribution** (up to 2.51×10⁸ W/m³) and cross-validated ANSYS Fluent (1/4-symmetry) against Star-CCM+ (full geometry) and analytical correlations (Petukhov, Dittus-Boelter).

</td>
<td width="50%" valign="top">

#### [Flow Analysis in a PWR Fuel Rod Subchannel](https://github.com/EmreSakarya/pwr-subchannel-flow-cfd)
`ANSYS Fluent` `Star-CCM+` `RANS` `k-ε` `k-ω SST`

Turbulent coolant flow in a 3.66 m bare-rod subchannel (D_h = 11.78 mm). Two independent CFD codes converged to within **~1%** of each other and of the Darcy-Weisbach / Filonenko analytical reference.

</td>
</tr>
<tr>
<td width="50%" valign="top">

#### [Steady-State Heat Conduction in a TRISO Particle](https://github.com/EmreSakarya/triso-heat-conduction-cfd)
`Star-CCM+` `Spherical Conduction` `Mesh Independence`

3D Star-CCM+ model of a 5-layer TRISO fuel particle (2 GW/m³ source) validated against a 1D analytical solution (T_max = **1106.35 K**). Mesh error < **0.0084%** across all refinement levels.

</td>
<td width="50%" valign="top">

#### [Single-Phase CFD of a VVER-1200 Rod Bundle](https://github.com/EmreSakarya/vver-1200-rod-bundle-cfd)
`Star-CCM+` `Conjugate Heat Transfer` `SST k-ω γ-Reθ` `Hexagonal Bundle`

Conjugate heat transfer simulation of a 19-rod (18 fuel + guide tube) **VVER-1200** hexagonal mini-bundle over the bottom 1/10 of the core with a sinusoidal axial power source (~1.8M cells). CFD matched analytical references to **0.18%** on peak fuel temperature, **1.9%** on coolant ΔT, and **8.1%** on pressure drop.

</td>
</tr>
</table>

---

## ⚙️ Thermodynamic Cycles & Reactor Systems

<table>
<tr>
<td width="50%" valign="top">

#### [Combined Cycle Power Plant Optimization](https://github.com/EmreSakarya/combined-cycle-analysis)
`Python` `Thermodynamics` `SciPy` `Optimization`

Gas-Steam Combined Cycle (Brayton + Rankine) with variable specific heat. Optimized compression ratios via iterative solvers, reaching a peak efficiency of **50.1%**.

</td>
<td width="50%" valign="top">

#### [iPWR Startup Transient Analysis](https://github.com/EmreSakarya/iPWR-Reactor-Simulation-Analysis)
`IAEA Simulator` `Reactivity Management` `Transient Analysis`

Controlled power ascension from 0% to 100% using the IAEA iPWR Simulator. Studied **SUR** fluctuations and verified safety margins against SCRAM setpoints.

</td>
</tr>
<tr>
<td width="50%" valign="top">

#### [Regenerative Rankine Cycle Analysis](https://github.com/EmreSakarya/regenerative-rankine-cycle)
`Thermodynamics` `Reheat` `OFWH` `IAPWS-97`

Trade-off analysis between turbine inlet temperature and steam quality using IAPWS-97 standards. Integrated reheat and open feedwater heating logic.

</td>
<td width="50%" valign="top"></td>
</tr>
</table>

---

## 🔭 Measurement & Image Analysis · ImageJ

<table>
<tr>
<td width="50%" valign="top">

#### 🆕 [Dating an Undated Moon Image with ImageJ](https://github.com/EmreSakarya/moon-image-dating-imagej)
`ImageJ / Fiji` `Macro Language` `Image Analysis` `Curve Fitting`

ImageJ plugin that dates an undated Moon image from six dated ones using size-corrected luminance, a Lambert-sphere phase model and a sunlit-pixel count. Estimate: **23 Oct 2026, 21:43 ± 2.9 h UTC**, within an hour of the NASA Dial-a-Moon check. [Report](https://github.com/EmreSakarya/moon-image-dating-imagej/blob/main/docs/NEM469_Project1_Report.pdf)

</td>
<td width="50%" valign="top"></td>
</tr>
</table>

---

<div align="center">

*Always open to collaborating on computational physics and engineering simulation projects.*

</div>
