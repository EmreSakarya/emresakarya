<div align="center">

# Hi there, I'm Emre Sakarya 👋

[![Typing SVG](https://readme-typing-svg.demolab.com?font=Fira+Code&size=18&pause=1000&color=C0392B&center=true&vCenter=true&width=650&lines=Nuclear+Engineering+%40+Hacettepe+University;Reactor+Physics+%7C+Thermal-Hydraulics+%7C+CFD;Python+%7C+Fortran+%7C+ANSYS+Fluent+%7C+Star-CCM%2B;Computational+Engineering+%26+Simulation)](https://git.io/typing-svg)

I'm a senior **Nuclear Engineering** student passionate about **Computational Engineering**, **Simulation**, and **Numerical Analysis**.
I bridge the gap between theoretical physics and software to solve complex engineering problems — from reactor kinetics to full-bundle CFD simulations.

[![LinkedIn](https://img.shields.io/badge/LinkedIn-0077B5?style=for-the-badge&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/emre-sakarya-ab8967300/)
[![Email](https://img.shields.io/badge/Email-D14836?style=for-the-badge&logo=gmail&logoColor=white)](mailto:emresakarya@hacettepe.edu.tr)

</div>

---

## 🛠️ Technical Stack

<table>
<tr>
<td valign="top" width="33%">

**Languages**

![Python](https://img.shields.io/badge/Python-3776AB?style=flat&logo=python&logoColor=white) Advanced\
![Fortran](https://img.shields.io/badge/Fortran-734F96?style=flat&logo=fortran&logoColor=white) Intel ifx / gfortran\
![C++](https://img.shields.io/badge/C%2B%2B-00599C?style=flat&logo=c%2B%2B&logoColor=white) Learning

</td>
<td valign="top" width="33%">

**CFD & Simulation Tools**

![ANSYS](https://img.shields.io/badge/ANSYS_Fluent-FFB71B?style=flat&logo=ansys&logoColor=black)\
![Star-CCM+](https://img.shields.io/badge/Star--CCM%2B-CC0000?style=flat&logoColor=white)\
![Linux](https://img.shields.io/badge/Linux%2FBash-FCC624?style=flat&logo=linux&logoColor=black)

</td>
<td valign="top" width="33%">

**Numerical Methods & Libraries**

🧮 NumPy · SciPy · Pandas\
📊 Matplotlib\
💧 PyXSteam (IAPWS-97)\
⚛️ Monte Carlo · FDM · FVM · RANS

</td>
</tr>
</table>

**Domain Expertise:** Point Reactor Kinetics · Transient Heat Transfer · Two-Phase Flow · CFD Turbulence Modeling (k-ε, k-ω SST) · Neutron Transport · Thermodynamic Cycle Analysis

---

## 📊 GitHub Stats

<div align="center">

[![GitHub Stats](https://github-readme-stats.vercel.app/api?username=EmreSakarya&show_icons=true&theme=github_dark&hide_border=true&title_color=C0392B&icon_color=C0392B)](https://github.com/EmreSakarya)
&nbsp;&nbsp;
[![Top Languages](https://github-readme-stats.vercel.app/api/top-langs/?username=EmreSakarya&layout=compact&theme=github_dark&hide_border=true&title_color=C0392B)](https://github.com/EmreSakarya)

</div>

---

## 🌊 CFD Simulations · ANSYS Fluent & Star-CCM+

> Industrial-grade computational fluid dynamics projects (NEM 358), cross-validated between two commercial solvers and analytical correlations.

#### ⚛️ [CFD Simulation of a 3×3 PWR Rod Bundle](https://github.com/EmreSakarya/pwr-rod-bundle-cfd)
> **Keywords:** *ANSYS Fluent, Star-CCM+, Conjugate Heat Transfer, k-ω SST, Rod Bundle, PWR*
* Performed a full **conjugate heat transfer (CHT)** simulation of a 3×3 PWR fuel rod bundle with a central guide tube (D_NFR = 12.243 mm), using parameters from a typical 4-loop PWR at 15.51 MPa.
* Applied a **non-uniform radial power distribution** (corner rods 0.85×, edge rods 1.00×) with volumetric heat sources up to **2.51×10⁸ W/m³**.
* ANSYS Fluent used 1/4-symmetry; Star-CCM+ solved the full geometry — cross-validated against analytical correlations (Petukhov friction factor, Dittus-Boelter Nusselt).
* Simultaneously resolved 3D velocity, temperature, and pressure fields across all subchannels.

#### 💧 [Flow Analysis in a PWR Fuel Rod Subchannel](https://github.com/EmreSakarya/pwr-subchannel-flow-cfd)
> **Keywords:** *ANSYS Fluent, Star-CCM+, Turbulent Flow, Darcy-Weisbach, RANS, k-ε, k-ω SST*
* Simulated turbulent coolant flow in a **3.66 m bare-rod PWR subchannel** (D = 9.5 mm, pitch = 12.6 mm, D_h = 11.78 mm) using two independent CFD codes.
* Applied **realizable k-ε** (Star-CCM+) and **k-ω SST** (ANSYS Fluent) closures; benchmarked the pressure drop against the **Darcy-Weisbach / Filonenko** analytical solution.
* Cross-validated ANSYS Fluent (~873 k cells) and Star-CCM+ (~920 k cells) — both converged to within **~1%** of each other and of the analytical reference.

#### 🔥 [Steady-State Heat Conduction in a TRISO Particle](https://github.com/EmreSakarya/triso-heat-conduction-cfd)
> **Keywords:** *Star-CCM+, CFD, Spherical Conduction, Mesh Independence, UO₂*
* Modeled steady-state conduction in a 5-layer **TRISO fuel particle** (kernel → buffer → IPyC → SiC → OPyC) with a **2 GW/m³** source in the UO₂ kernel.
* Validated a 3D Star-CCM+ model against a 1D analytical solution (T_max = **1106.35 K**) through a rigorous mesh-independence study (max error **< 0.0084%**).
* Extended the model with **temperature-dependent UO₂ conductivity** k(T) and quantified the resulting peak-temperature shift.

---

## 🔬 Numerical Methods & Heat Transfer · Python / Fortran

> Solvers built from scratch for non-linear, stiff, and multi-physics heat-transfer problems.

#### 🌡️ [2D Transient Thermal Analysis of a Rectangular Plate](https://github.com/EmreSakarya/2d-transient-thermal-analysis)
> **Keywords:** *Python, Control Volume FDM, Explicit Euler, Natural Convection*
* Solved a two-stage heat-treatment problem on a 3×5 nodal grid using the **Control Volume Finite Difference Method (CVFDM)** with temperature-dependent properties updated every step.
* **Heating:** variable corner heat sources (28–37 kW/m³) under adiabatic boundaries → all nodes exceed **1000 K**.
* **Cooling:** liquid-metal bath with Ra-based **natural convection** → all nodes drop below **600 K** (~1920 s).

#### 📈 [Point Reactor Kinetics Analysis](https://github.com/EmreSakarya/point-reactor-kinetics-analysis)
> **Keywords:** *Python, Runge-Kutta (RK4), Heun's Method, Stiff Systems*
* Solved the point reactor kinetics equations with **delayed neutrons** to analyze transient behavior under reactivity insertions.
* Implemented **RK4** and **Heun's Predictor-Corrector** integrators; conducted stiffness-ratio analysis against analytical solutions.

#### ☢️ [U-238 Decay Chain Modeling](https://github.com/EmreSakarya/uranium-decay-chain-analysis)
> **Keywords:** *Fortran 90, Implicit Euler, Stiff Differential Equations*
* Solved a 14-stage decay chain (U-238 → Pb-206) with half-lives spanning $10^9$ years to seconds — a severely **stiff system**.
* Developed a generic **Implicit Euler solver in Fortran** for stability where explicit methods fail, mirroring legacy HPC nuclear-code logic.

#### 🌊 [Thermal-Hydraulic Analysis (HEM)](https://github.com/EmreSakarya/vertical-pipe-flow-analysis)
> **Keywords:** *Python, Finite Volume Method, Two-Phase Flow*
* Built a **1D numerical solver** for downward water flow under constant heat flux.
* Used the **Homogeneous Equilibrium Model (HEM)** and IAPWS-97 standards to predict boiling inception and void-fraction profiles.

---

## ⚛️ Monte Carlo & Particle Transport

#### 🎲 [Stochastic Neutron Transport Simulation](https://github.com/EmreSakarya/neutron-transport-monte-carlo)
> **Keywords:** *Monte Carlo Method, Statistical Analysis, Slab Geometry*
* Developed a stochastic neutron-trajectory algorithm from scratch in slab geometry.
* Validated against analytical **Exponential Integral** functions and demonstrated $1/\sqrt{N}$ error scaling.

#### ⚡ [Charged Particle Energy Loss Simulation](https://github.com/EmreSakarya/bethe-bloch-analysis)
> **Keywords:** *Bethe-Bloch Formula, Adaptive Step-Size Integration*
* Modeled energy deposition of protons and alpha particles in water / hydrogen media.
* Simulated **phase-transition effects** on stopping power using adaptive numerical integration.

---

## ⚙️ Thermodynamic Cycles & Reactor Systems

#### 🏭 [Combined Cycle Power Plant Optimization](https://github.com/EmreSakarya/combined-cycle-analysis)
> **Keywords:** *Python, Thermodynamics, SciPy, Optimization*
* Modeled a Gas-Steam Combined Cycle (Brayton + Rankine) with **variable specific heat** for realistic gas behavior.
* Optimized compression ratios via iterative solvers (`fsolve`), reaching a peak efficiency of **50.1%**; analyzed heat-exchanger pinch points.

#### ♻️ [Regenerative Rankine Cycle Analysis](https://github.com/EmreSakarya/regenerative-rankine-cycle)
> **Keywords:** *Thermodynamics, Reheat, OFWH, IAPWS-97*
* Analyzed trade-offs between turbine inlet temperature and steam quality using IAPWS-97 standards.
* Integrated reheat and open feedwater heating logic for thermodynamic efficiency gains.

#### 🛰️ [iPWR Startup Transient Analysis](https://github.com/EmreSakarya/iPWR-Reactor-Simulation-Analysis)
> **Keywords:** *IAEA Simulator, Reactivity Management, Transient Analysis*
* Analyzed a controlled power ascension from 0% to 100% using the **IAEA iPWR Simulator**.
* Studied **Start-Up Rate (SUR)** fluctuations and verified safety margins against SCRAM setpoints.

---

<div align="center">

*Always open to collaborating on computational physics and engineering simulation projects.*

</div>
