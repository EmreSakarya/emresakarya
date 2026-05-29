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

## 🚀 Top Featured Projects

#### 📈 [Point Reactor Kinetics Analysis](https://github.com/EmreSakarya/point-reactor-kinetics-analysis)
> **Keywords:** *Python, Runge-Kutta (RK4), Heun's Method, Stiff Systems*
* Solved point reactor kinetics equations with **delayed neutrons** to analyze transient behavior under reactivity insertions.
* Implemented **Runge-Kutta (RK4)** and **Heun's Predictor-Corrector** methods for high-precision numerical integration.
* Conducted stiffness ratio analysis to evaluate numerical stability against analytical solutions.

#### ☢️ [U-238 Decay Chain Modeling](https://github.com/EmreSakarya/uranium-decay-chain-analysis)
> **Keywords:** *Fortran 90, Implicit Euler, Stiff Differential Equations*
* Solved a 14-stage decay chain (U-238 to Pb-206) containing **stiff systems** with half-lives ranging from $10^9$ years to seconds.
* Developed a generic **Implicit Euler Solver in Fortran** to ensure numerical stability where explicit methods (Python) fail.
* Demonstrated proficiency in **High-Performance Computing (HPC)** logic used in legacy nuclear codes.

#### ⚙️ [Combined Cycle Power Plant Optimization](https://github.com/EmreSakarya/combined-cycle-analysis)
> **Keywords:** *Python, Thermodynamics, SciPy, Optimization*
* Modeled a Gas-Steam Combined Cycle (Brayton + Rankine) using **Variable Specific Heat** assumptions for realistic gas behavior.
* Implemented iterative solvers (`fsolve`) to optimize compression ratios, achieving a peak efficiency of **50.1%**.
* Analyzed pinch-point constraints in heat exchangers.

#### 🌊 [Thermal-Hydraulic Analysis (HEM)](https://github.com/EmreSakarya/vertical-pipe-flow-analysis)
> **Keywords:** *Python, Finite Volume Method, Two-Phase Flow*
* Developed a **1D Numerical Solver** for downward water flow subjected to constant heat flux.
* Utilized the **Homogeneous Equilibrium Model (HEM)** and IAPWS-97 standards to predict boiling inception and void fraction profiles.

#### 🌡️ [2D Transient Thermal Analysis of a Rectangular Plate](https://github.com/EmreSakarya/2d-transient-thermal-analysis)
> **Keywords:** *Python, Control Volume FDM, Explicit Euler, Transient Heat Transfer, Natural Convection*
* Solved a two-stage heat treatment problem for a 6×16 cm rectangular plate on a 3×5 nodal grid using the **Control Volume Finite Difference Method (CVFDM)**.
* **Heating Phase:** Applied variable volumetric heat sources (28–37 kW/m³) at four corner nodes under adiabatic boundaries; all nodes exceeded **1000 K** in ~28.6 seconds (simulated ~8 hours).
* **Cooling Phase:** Immersed the plate in a liquid metal bath with temperature-dependent **natural convection** (Ra-based correlations); all nodes cooled below **600 K** in ~1920 seconds.

---

## 🌊 CFD Simulations — NEM 358 (ANSYS Fluent & Star-CCM+)

#### 🔥 [Steady-State Heat Conduction in a TRISO Particle](https://github.com/EmreSakarya/triso-heat-conduction-cfd)
> **Keywords:** *Star-CCM+, CFD, Spherical Conduction, Mesh Independence, UO₂*
* Modeled steady-state heat conduction in a 5-layer **TRISO fuel particle** (kernel → buffer → IPyC → SiC → OPyC) with a volumetric heat source of **2 GW/m³** in the UO₂ kernel.
* Performed a **3D CFD analysis in Star-CCM+** using polyhedral meshes at three density levels (Coarse/Medium/Fine); validated against a 1D analytical solution (T_max = **1106.35 K**).
* Achieved mesh independence with a maximum relative error of **< 0.0084%** across all mesh levels.
* Extended the model with **temperature-dependent UO₂ thermal conductivity** and quantified the resulting shift in peak temperature.

#### 💧 [Flow Analysis in a PWR Fuel Rod Subchannel](https://github.com/EmreSakarya/pwr-subchannel-flow-cfd)
> **Keywords:** *ANSYS Fluent, Star-CCM+, Turbulent Flow, Darcy-Weisbach, RANS, k-ε, k-ω SST*
* Simulated turbulent coolant flow in a **3.66 m bare-rod PWR subchannel** (D = 9.5 mm, pitch = 12.6 mm, D_h = 11.78 mm) using two independent CFD codes.
* Applied **realizable k-ε** (Star-CCM+) and **k-ω SST** (ANSYS Fluent) turbulence models; compared pressure drop results against the **Darcy-Weisbach / Filonenko** analytical solution.
* Analyzed axial pressure profiles, inlet velocity distributions, and the hydrodynamic entrance region.
* Cross-validated ANSYS Fluent (~873 k cells) and Star-CCM+ (~920 k cells) results — both codes converged to within ~1% of each other and of the analytical reference.

#### ⚛️ [CFD Simulation of a 3×3 PWR Rod Bundle](https://github.com/EmreSakarya/pwr-rod-bundle-cfd)
> **Keywords:** *ANSYS Fluent, Star-CCM+, Conjugate Heat Transfer, k-ω SST, Rod Bundle, PWR*
* Performed a full **conjugate heat transfer CFD simulation** of a 3×3 PWR fuel rod bundle with a non-fuel central guide tube (D_NFR = 12.243 mm), using parameters from a typical 4-loop PWR at 15.51 MPa.
* Applied a **non-uniform radial power distribution** (corner rods: 0.85×, edge rods: 1.00×) with volumetric heat sources up to **2.51×10⁸ W/m³**.
* ANSYS Fluent used 1/4-symmetry (1 subchannel); Star-CCM+ solved the full 3×3 geometry — results cross-validated against analytical correlations (Petukhov friction factor, Dittus-Boelter Nusselt).
* Simultaneously resolved 3D velocity, temperature, and pressure fields within all subchannels.

---

## 📂 Other Engineering Simulations

#### ⚛️ [iPWR Startup Transient Analysis](https://github.com/EmreSakarya/iPWR-Reactor-Simulation-Analysis)
* **Focus:** *IAEA Simulator, Reactivity Management, Transient Analysis*
* Analyzed a controlled power ascension from 0% to 100% using the **IAEA iPWR Simulator**.
* Specifically analyzed the **Start-Up Rate (SUR)** fluctuations and verified safety margins against SCRAM setpoints.

#### ⚛️ [Stochastic Neutron Transport Simulation](https://github.com/EmreSakarya/neutron-transport-monte-carlo)
* **Focus:** *Monte Carlo Method, Statistical Analysis*
* Developed a stochastic algorithm from scratch to model neutron trajectories in slab geometry.
* Validated results against analytical **Exponential Integral** functions and demonstrated $1/\sqrt{N}$ error scaling.

#### ⚡ [Charged Particle Energy Loss Simulation](https://github.com/EmreSakarya/bethe-bloch-analysis)
* **Focus:** *Bethe-Bloch Formula, Adaptive Step-Size*
* Modeled energy deposition of Protons and Alpha particles in water/hydrogen media.
* Simulated **phase transition effects** on stopping power using adaptive numerical integration.

#### 🏭 [Regenerative Rankine Cycle Analysis](https://github.com/EmreSakarya/regenerative-rankine-cycle)
* **Focus:** *Thermodynamics, Reheat, OFWH*
* Analyzed trade-offs between turbine inlet temperature and steam quality using IAPWS-97 standards.
* Integrated logic for thermodynamic efficiency improvements.

---

<div align="center">

*Always open to collaborating on computational physics and engineering simulation projects.*

</div>
