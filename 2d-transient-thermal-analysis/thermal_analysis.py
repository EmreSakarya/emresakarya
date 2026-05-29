"""
2D Transient Thermal Analysis of a Rectangular Plate
=====================================================

Two-stage heat treatment of a 6 cm x 16 cm x 1 cm rectangular plate solved with
the Control Volume Finite Difference Method (CVFDM) and an explicit Euler
time integration scheme:

    Stage 1 (Heating):  localized volumetric heating via spring clamps at the
                        four corner nodes, adiabatic outer boundaries.
    Stage 2 (Cooling):  natural convection cooling inside a liquid-metal bath
                        whose temperature decreases linearly with time.

All thermophysical properties of both the solid plate and the liquid-metal
coolant are temperature dependent and are updated at every node and every
time step (the system is therefore strongly non-linear).

NEM 394 Engineering Project III - Hacettepe University, Nuclear Engineering.
"""

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Solid plate properties (temperature dependent)
# ---------------------------------------------------------------------------
def solid_conductivity(T):
    k0 = 0.0025
    a = 5.0e-6
    b = 5.0e-10
    return k0 + a * T + b * T * T * T


def solid_specific_heat(T):
    return 700.0 + 0.2 * T


def solid_density(T):
    rho0 = 200.0
    beta = 2.0e-6
    return rho0 * (1.0 - beta * (T - 300.0))


def solid_diffusivity(T):
    # alpha = k / (rho * cp)
    return solid_conductivity(T) / (solid_density(T) * solid_specific_heat(T))


# ---------------------------------------------------------------------------
# Liquid metal coolant properties (temperature dependent)
# ---------------------------------------------------------------------------
def liquid_conductivity(T):
    return 124.0 - 0.113 * T + 5.22e-5 * T * T


def liquid_density(T):
    return 971.0 - 0.247 * T


def liquid_expansion(T):
    return 0.247 / (971.0 - 0.247 * T)


def liquid_viscosity(T):
    return 2.94e-4 * np.exp(734.0 / T)


def liquid_specific_heat(T):
    return 1240.0 - 0.27 * T


def liquid_kinematic_viscosity(T):
    return liquid_viscosity(T) / liquid_density(T)


def liquid_diffusivity(T):
    return liquid_conductivity(T) / (liquid_density(T) * liquid_specific_heat(T))


def liquid_prandtl(T):
    return liquid_viscosity(T) * liquid_specific_heat(T) / liquid_conductivity(T)


def bath_temperature(t):
    # The bath cools down linearly with time.
    return 1000.0 - 125.0 * (t / 600.0)


# ---------------------------------------------------------------------------
# Geometry and grid
# ---------------------------------------------------------------------------
Lx = 0.06        # plate width      [m]
Ly = 0.16        # plate length     [m]
Lz = 0.01        # plate thickness  [m]

dx = 0.03        # node spacing in x [m]
dy = 0.04        # node spacing in y [m]

nx = 3           # number of nodes in x  (Lx/dx + 1)
ny = 5           # number of nodes in y  (Ly/dy + 1)

g = 9.81         # gravity [m/s^2]

# Characteristic lengths used in the convection correlations.
# Horizontal faces use the standard As/P definition.
Lc_horizontal = (Lx * Ly) / (2.0 * (Lx + Ly))
Lc_vertical = Lz


# Boundary nodes only own half of a cell.
def cv_width_x(i):
    if i == 0 or i == nx - 1:
        return dx / 2.0
    return dx


def cv_width_y(j):
    if j == 0 or j == ny - 1:
        return dy / 2.0
    return dy


# Heat generation: only the four corner nodes have a source.
qgen = np.zeros((nx, ny))
qgen[0, 0] = 32000.0           # bottom-left clamp
qgen[nx - 1, 0] = 35000.0      # bottom-right clamp
qgen[0, ny - 1] = 28000.0      # top-left clamp
qgen[nx - 1, ny - 1] = 37000.0 # top-right clamp


# Classify the nodes into three groups so they can be plotted separately later.
corner_nodes = []
side_nodes = []
interior_nodes = []
for i in range(nx):
    for j in range(ny):
        is_corner = (i == 0 or i == nx - 1) and (j == 0 or j == ny - 1)
        is_edge = (i == 0 or i == nx - 1 or j == 0 or j == ny - 1)
        if is_corner:
            corner_nodes.append((i, j))
        elif is_edge:
            side_nodes.append((i, j))
        else:
            interior_nodes.append((i, j))


# ---------------------------------------------------------------------------
# PART A - Stability criterion
# The most restrictive case is at the highest temperature (1000 K) where the
# diffusivity is largest (cubic term in k).
# ---------------------------------------------------------------------------
alpha_at_1000 = solid_diffusivity(1000.0)
dt_max = 1.0 / (2.0 * alpha_at_1000 * (1.0 / (dx * dx) + 1.0 / (dy * dy)))

print("---------- PART A: Stability ----------")
print("alpha(300 K)  =", solid_diffusivity(300.0), "m^2/s")
print("alpha(1000 K) =", alpha_at_1000, "m^2/s")
print("dt_max        =", round(dt_max, 3), "s")

# Use 90% of the limit for the main heating run.
dt_heat = 0.90 * dt_max
print("dt_heat used  =", round(dt_heat, 3), "s")


# ---------------------------------------------------------------------------
# PART B - Heating simulation
# Boundaries are adiabatic, so a missing neighbour simply contributes no flux.
# ---------------------------------------------------------------------------
print("\n---------- PART B: Heating ----------")

T = np.full((nx, ny), 300.0)  # initial temperature everywhere

time_now = 0.0
heat_time_list = [0.0]
heat_T_history = [T.copy()]

# Keep stepping until even the coldest node has reached 1000 K.
keep_going = True
while keep_going:
    T_new = T.copy()
    for i in range(nx):
        for j in range(ny):
            T_here = T[i, j]
            k_here = solid_conductivity(T_here)

            wx = cv_width_x(i)
            wy = cv_width_y(j)
            volume = wx * wy * Lz

            # conduction from each existing neighbour
            Q_cond = 0.0
            if i < nx - 1:  # +x
                k_face = 0.5 * (k_here + solid_conductivity(T[i + 1, j]))
                Q_cond += k_face * (wy * Lz) * (T[i + 1, j] - T_here) / dx
            if i > 0:        # -x
                k_face = 0.5 * (k_here + solid_conductivity(T[i - 1, j]))
                Q_cond += k_face * (wy * Lz) * (T[i - 1, j] - T_here) / dx
            if j < ny - 1:  # +y
                k_face = 0.5 * (k_here + solid_conductivity(T[i, j + 1]))
                Q_cond += k_face * (wx * Lz) * (T[i, j + 1] - T_here) / dy
            if j > 0:        # -y
                k_face = 0.5 * (k_here + solid_conductivity(T[i, j - 1]))
                Q_cond += k_face * (wx * Lz) * (T[i, j - 1] - T_here) / dy

            # internal heat generation
            Q_source = qgen[i, j] * volume

            # explicit update
            heat_capacity = solid_density(T_here) * solid_specific_heat(T_here) * volume
            T_new[i, j] = T_here + dt_heat * (Q_cond + Q_source) / heat_capacity

    T = T_new
    time_now += dt_heat
    heat_time_list.append(time_now)
    heat_T_history.append(T.copy())

    if T.min() >= 1000.0:
        keep_going = False

heat_time = np.array(heat_time_list)
heat_T_history = np.array(heat_T_history)
T_after_heating = T.copy()

# locate the slowest and the hottest node
slow_i, slow_j = np.unravel_index(np.argmin(T), T.shape)
hot_i, hot_j = np.unravel_index(np.argmax(T), T.shape)

print("Heating finished at t =", round(time_now, 1), "s   (",
      round(time_now / 3600.0, 2), "hours )")
print("Slowest node =", (slow_i, slow_j), "  T =", round(T[slow_i, slow_j], 2), "K")
print("Hottest node =", (hot_i, hot_j), "  T =", round(T[hot_i, hot_j], 2), "K")


# ---------------------------------------------------------------------------
# PART C - Effect of the time step
# Repeat the heating run for several dt values and record run time / steps.
# ---------------------------------------------------------------------------
print("\n---------- PART C: dt comparison ----------")

dt_try_list = [0.99 * dt_max, 0.90 * dt_max, 0.50 * dt_max, 5.0, 0.1]
partc_dt = []
partc_time = []
partc_steps = []

for dt_try in dt_try_list:
    Tc = np.full((nx, ny), 300.0)
    tt = 0.0
    nstep = 0
    diverged = False

    busy = True
    while busy:
        Tc_new = Tc.copy()
        for i in range(nx):
            for j in range(ny):
                T_here = Tc[i, j]
                k_here = solid_conductivity(T_here)
                wx = cv_width_x(i)
                wy = cv_width_y(j)
                volume = wx * wy * Lz
                Q_cond = 0.0
                if i < nx - 1:
                    k_face = 0.5 * (k_here + solid_conductivity(Tc[i + 1, j]))
                    Q_cond += k_face * (wy * Lz) * (Tc[i + 1, j] - T_here) / dx
                if i > 0:
                    k_face = 0.5 * (k_here + solid_conductivity(Tc[i - 1, j]))
                    Q_cond += k_face * (wy * Lz) * (Tc[i - 1, j] - T_here) / dx
                if j < ny - 1:
                    k_face = 0.5 * (k_here + solid_conductivity(Tc[i, j + 1]))
                    Q_cond += k_face * (wx * Lz) * (Tc[i, j + 1] - T_here) / dy
                if j > 0:
                    k_face = 0.5 * (k_here + solid_conductivity(Tc[i, j - 1]))
                    Q_cond += k_face * (wx * Lz) * (Tc[i, j - 1] - T_here) / dy
                Q_source = qgen[i, j] * volume
                heat_capacity = solid_density(T_here) * solid_specific_heat(T_here) * volume
                Tc_new[i, j] = T_here + dt_try * (Q_cond + Q_source) / heat_capacity

        Tc = Tc_new
        tt += dt_try
        nstep += 1

        # blow-up check
        if (not np.all(np.isfinite(Tc))) or (Tc.max() > 1.0e6):
            diverged = True
            busy = False

        # stop condition
        if Tc.min() >= 1000.0:
            busy = False

        # hard safety cap so the loop never runs forever
        if nstep > 5000000:
            busy = False

    if diverged:
        print("dt =", round(dt_try, 3), "s    -> DIVERGED")
    else:
        print("dt =", round(dt_try, 3), "s    -> t =", round(tt, 1), "s ,", nstep, "steps")
        partc_dt.append(dt_try)
        partc_time.append(tt)
        partc_steps.append(nstep)


# ---------------------------------------------------------------------------
# PART D - Cooling in the liquid metal bath
# No heat generation. Every face loses heat by natural convection to the bath.
# The convection coefficient h is computed locally from the Rayleigh number.
# ---------------------------------------------------------------------------
print("\n---------- PART D: Cooling ----------")


def convection_coefficient(T_surface, T_bath, Lc, orientation):
    """Natural-convection coefficient for one face, three explicit orientations."""
    delta_T = abs(T_surface - T_bath)
    if delta_T < 1.0e-3:
        return 0.0

    T_film = 0.5 * (T_surface + T_bath)
    if T_film < 50.0:
        T_film = 50.0

    beta = liquid_expansion(T_film)
    nu = liquid_kinematic_viscosity(T_film)
    alpha = liquid_diffusivity(T_film)
    k = liquid_conductivity(T_film)
    Pr = liquid_prandtl(T_film)

    Ra = g * beta * delta_T * Lc * Lc * Lc / (nu * alpha)
    if Ra <= 0.0:
        return 0.0

    if orientation == "top":
        if Ra < 1.0e7:
            Nu = 0.54 * Ra ** 0.25
        else:
            Nu = 0.15 * Ra ** (1.0 / 3.0)
    elif orientation == "bottom":
        Nu = 0.27 * Ra ** 0.25
    else:
        # vertical side faces, Churchill-Chu
        denom = (1.0 + (0.492 / Pr) ** (9.0 / 16.0)) ** (4.0 / 9.0)
        Nu = 0.68 + 0.670 * Ra ** 0.25 / denom

    return Nu * k / Lc


def cooling_stable_dt(T_field, t_value):
    """Recompute the allowed step from the current field (convection dominates)."""
    Tb = bath_temperature(t_value)
    worst = 0.0
    for i in range(nx):
        for j in range(ny):
            T_here = T_field[i, j]
            wx = cv_width_x(i)
            wy = cv_width_y(j)
            volume = wx * wy * Lz
            rcpV = solid_density(T_here) * solid_specific_heat(T_here) * volume

            face_area = wx * wy
            h_top = convection_coefficient(T_here, Tb, Lc_horizontal, "top")
            h_bot = convection_coefficient(T_here, Tb, Lc_horizontal, "bottom")
            hA = (h_top + h_bot) * face_area

            if i == 0 or i == nx - 1:
                hv = convection_coefficient(T_here, Tb, Lc_vertical, "vert")
                hA += hv * wy * Lz
            if j == 0 or j == ny - 1:
                hv = convection_coefficient(T_here, Tb, Lc_vertical, "vert")
                hA += hv * wx * Lz

            cond_part = 2.0 * solid_diffusivity(T_here) * (1.0 / (dx * dx) + 1.0 / (dy * dy))
            value = hA / rcpV + cond_part
            if value > worst:
                worst = value

    return 1.0 / worst if worst > 0.0 else 1.0


# start cooling from the heated field
T = T_after_heating.copy()

dt_cool = 0.40 * cooling_stable_dt(T, 0.0)
print("initial cooling dt =", round(dt_cool, 4), "s")

cool_time_list = [0.0]
cool_T_history = [T.copy()]
htop_history = [np.zeros((nx, ny))]
bath_history = [bath_temperature(0.0)]

t_cool = 0.0
step_count = 0
save_every = 25

cooling = True
while cooling:
    # refresh the stable step every once in a while
    if step_count % 20 == 0:
        candidate = 0.40 * cooling_stable_dt(T, t_cool)
        if candidate < dt_cool:
            dt_cool = candidate

    Tb = bath_temperature(t_cool)
    T_new = T.copy()
    htop_now = np.zeros((nx, ny))

    for i in range(nx):
        for j in range(ny):
            T_here = T[i, j]
            k_here = solid_conductivity(T_here)
            wx = cv_width_x(i)
            wy = cv_width_y(j)
            volume = wx * wy * Lz

            # conduction from neighbours (same as heating)
            Q_cond = 0.0
            if i < nx - 1:
                k_face = 0.5 * (k_here + solid_conductivity(T[i + 1, j]))
                Q_cond += k_face * (wy * Lz) * (T[i + 1, j] - T_here) / dx
            if i > 0:
                k_face = 0.5 * (k_here + solid_conductivity(T[i - 1, j]))
                Q_cond += k_face * (wy * Lz) * (T[i - 1, j] - T_here) / dx
            if j < ny - 1:
                k_face = 0.5 * (k_here + solid_conductivity(T[i, j + 1]))
                Q_cond += k_face * (wx * Lz) * (T[i, j + 1] - T_here) / dy
            if j > 0:
                k_face = 0.5 * (k_here + solid_conductivity(T[i, j - 1]))
                Q_cond += k_face * (wx * Lz) * (T[i, j - 1] - T_here) / dy

            # convection on top and bottom faces
            face_area = wx * wy
            h_top = convection_coefficient(T_here, Tb, Lc_horizontal, "top")
            h_bot = convection_coefficient(T_here, Tb, Lc_horizontal, "bottom")
            hA = (h_top + h_bot) * face_area

            # convection on outer vertical faces (boundary nodes only)
            if i == 0 or i == nx - 1:
                hv = convection_coefficient(T_here, Tb, Lc_vertical, "vert")
                hA += hv * wy * Lz
            if j == 0 or j == ny - 1:
                hv = convection_coefficient(T_here, Tb, Lc_vertical, "vert")
                hA += hv * wx * Lz

            Q_conv = -hA * (T_here - Tb)

            heat_capacity = solid_density(T_here) * solid_specific_heat(T_here) * volume
            T_new[i, j] = T_here + dt_cool * (Q_cond + Q_conv) / heat_capacity
            htop_now[i, j] = h_top

    T = T_new
    t_cool += dt_cool
    step_count += 1

    if step_count % save_every == 0:
        cool_time_list.append(t_cool)
        cool_T_history.append(T.copy())
        htop_history.append(htop_now.copy())
        bath_history.append(Tb)

    # stop once the hottest node is below 600 K
    if T.max() <= 600.0:
        cooling = False

cool_time = np.array(cool_time_list)
cool_T_history = np.array(cool_T_history)
htop_history = np.array(htop_history)
bath_history = np.array(bath_history)

print("Cooling finished at t =", round(t_cool, 2), "s   (",
      round(t_cool / 60.0, 2), "min )")
print("final T_max =", round(T.max(), 2), "K  , bath =",
      round(bath_temperature(t_cool), 2), "K")


# ---------------------------------------------------------------------------
# FIGURES
# ---------------------------------------------------------------------------

# ---- Figure 1: solid properties ----
T_axis = np.linspace(300.0, 1100.0, 300)
k_curve = solid_conductivity(T_axis)
rho_curve = solid_density(T_axis)
cp_curve = solid_specific_heat(T_axis)

fig1, ax1 = plt.subplots(1, 3, figsize=(13, 4))
ax1[0].plot(T_axis, k_curve, color="red")
ax1[0].set_xlabel("T [K]"); ax1[0].set_ylabel("k_s [W/m K]")
ax1[0].set_title("Conductivity"); ax1[0].grid(alpha=0.3)
ax1[1].plot(T_axis, rho_curve, color="blue")
ax1[1].set_xlabel("T [K]"); ax1[1].set_ylabel("rho_s [kg/m3]")
ax1[1].set_title("Density"); ax1[1].grid(alpha=0.3)
ax1[2].plot(T_axis, cp_curve, color="green")
ax1[2].set_xlabel("T [K]"); ax1[2].set_ylabel("cp_s [J/kg K]")
ax1[2].set_title("Specific Heat"); ax1[2].grid(alpha=0.3)
fig1.suptitle("Figure 1 - Solid Plate Properties")
fig1.tight_layout()
fig1.savefig("fig01_solid_props.png", dpi=130)


# ---- Figures 2,3,4: heating temperature for each group ----
def plot_heating_group(node_list, title_text, file_name):
    fig, ax = plt.subplots(figsize=(9, 5))
    for (i, j) in node_list:
        if qgen[i, j] > 0.0:
            label_text = "Node(%d,%d)  q=%d kW/m3" % (i, j, int(qgen[i, j] / 1000.0))
        else:
            label_text = "Node(%d,%d)" % (i, j)
        ax.plot(heat_time, heat_T_history[:, i, j], linewidth=1.8, label=label_text)
    ax.axhline(1000.0, color="black", linestyle="--", linewidth=1.2, label="1000 K target")
    ax.set_xlabel("time [s]"); ax.set_ylabel("T [K]")
    ax.set_title(title_text); ax.grid(alpha=0.3); ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(file_name, dpi=130)


plot_heating_group(corner_nodes, "Figure 2 - Heating: Edge (Corner) Nodes", "fig02_heat_edge.png")
plot_heating_group(side_nodes, "Figure 3 - Heating: Side Nodes", "fig03_heat_side.png")
plot_heating_group(interior_nodes, "Figure 4 - Heating: Interior Nodes", "fig04_heat_interior.png")


# ---- Figure 5: dt comparison (Part C) ----
fig5, ax5 = plt.subplots(1, 2, figsize=(11, 4.5))
ax5[0].plot(partc_dt, partc_time, marker="o")
ax5[0].set_xlabel("dt [s]"); ax5[0].set_ylabel("required time [s]")
ax5[0].set_title("Time to reach 1000 K"); ax5[0].grid(alpha=0.3)
ax5[1].plot(partc_dt, partc_steps, marker="s", color="orange")
ax5[1].set_yscale("log")
ax5[1].set_xlabel("dt [s]"); ax5[1].set_ylabel("number of steps")
ax5[1].set_title("Computational cost"); ax5[1].grid(alpha=0.3)
fig5.suptitle("Figure 5 - Effect of Time Step (Part C)")
fig5.tight_layout()
fig5.savefig("fig05_dt_comparison.png", dpi=130)


# ---- Figure 6: liquid properties + Rayleigh ----
Tl_axis = np.linspace(300.0, 1100.0, 300)
kl_c = liquid_conductivity(Tl_axis)
rhol_c = liquid_density(Tl_axis)
cpl_c = liquid_specific_heat(Tl_axis)
mul_c = liquid_viscosity(Tl_axis)
betal_c = liquid_expansion(Tl_axis)
Ra_c = np.where(
    Tl_axis > 330.0,
    g * liquid_expansion(Tl_axis) * 30.0 * Lc_horizontal ** 3
    / (liquid_kinematic_viscosity(Tl_axis) * liquid_diffusivity(Tl_axis)),
    0.0,
)

fig6, ax6 = plt.subplots(2, 3, figsize=(14, 8))
ax6[0, 0].plot(Tl_axis, kl_c, color="red"); ax6[0, 0].set_title("Conductivity"); ax6[0, 0].set_ylabel("k_l [W/m K]")
ax6[0, 1].plot(Tl_axis, rhol_c, color="blue"); ax6[0, 1].set_title("Density"); ax6[0, 1].set_ylabel("rho_l [kg/m3]")
ax6[0, 2].plot(Tl_axis, cpl_c, color="green"); ax6[0, 2].set_title("Specific Heat"); ax6[0, 2].set_ylabel("cp_l [J/kg K]")
ax6[1, 0].plot(Tl_axis, mul_c, color="purple"); ax6[1, 0].set_title("Viscosity"); ax6[1, 0].set_ylabel("mu_l [Pa s]"); ax6[1, 0].set_yscale("log")
ax6[1, 1].plot(Tl_axis, betal_c, color="brown"); ax6[1, 1].set_title("Expansion"); ax6[1, 1].set_ylabel("beta_l [1/K]")
ax6[1, 2].plot(Tl_axis, Ra_c, color="black"); ax6[1, 2].set_title("Rayleigh (dT=30K)"); ax6[1, 2].set_ylabel("Ra"); ax6[1, 2].set_yscale("log")
for r in range(2):
    for c in range(3):
        ax6[r, c].set_xlabel("T [K]"); ax6[r, c].grid(alpha=0.3)
fig6.suptitle("Figure 6 - Liquid Metal Properties and Rayleigh Number")
fig6.tight_layout()
fig6.savefig("fig06_liquid_props.png", dpi=130)


# ---- Figure 7: cooling temperature, all groups in one plot ----
fig7, ax7 = plt.subplots(figsize=(10, 6))
ax7.plot(cool_time, bath_history, color="black", linestyle="--", linewidth=2.0, label="T_bath")
for (i, j) in corner_nodes:
    ax7.plot(cool_time, cool_T_history[:, i, j], color="red", linewidth=1.4)
for (i, j) in side_nodes:
    ax7.plot(cool_time, cool_T_history[:, i, j], color="blue", linewidth=1.2)
for (i, j) in interior_nodes:
    ax7.plot(cool_time, cool_T_history[:, i, j], color="green", linewidth=1.4)
ax7.axhline(600.0, color="red", linestyle=":", linewidth=1.5)
ax7.set_xlabel("time [s]"); ax7.set_ylabel("T [K]")
ax7.set_title("Figure 7 - Cooling: All Node Groups (overlap)")
ax7.grid(alpha=0.3)
ax7.legend()
fig7.tight_layout()
fig7.savefig("fig07_cool_all.png", dpi=130)


# ---- Figure 8: convection coefficient h(t) for edge and side nodes ----
n_h = len(htop_history)
fig8, ax8 = plt.subplots(1, 2, figsize=(14, 5))
for (i, j) in corner_nodes:
    ax8[0].plot(cool_time[0:n_h], htop_history[:, i, j], linewidth=1.8,
                label="Node(%d,%d)" % (i, j))
ax8[0].set_xlabel("time [s]"); ax8[0].set_ylabel("h_top [W/m2 K]")
ax8[0].set_title("(a) Edge (Corner) Nodes"); ax8[0].grid(alpha=0.3); ax8[0].legend(fontsize=9)
for (i, j) in side_nodes:
    ax8[1].plot(cool_time[0:n_h], htop_history[:, i, j], linewidth=1.8,
                label="Node(%d,%d)" % (i, j))
ax8[1].set_xlabel("time [s]"); ax8[1].set_ylabel("h_top [W/m2 K]")
ax8[1].set_title("(b) Side Nodes"); ax8[1].grid(alpha=0.3); ax8[1].legend(fontsize=9, ncol=2)
fig8.suptitle("Figure 8 - Convection Coefficient h(t) During Cooling")
fig8.tight_layout()
fig8.savefig("fig08_h_combined.png", dpi=130)

print("\nAll figures saved.")

if __name__ == "__main__":
    # plt.show()  # uncomment to display the figures interactively
    pass
