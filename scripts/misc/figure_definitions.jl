# Generic colors
pacco_color = :grey20
background_color = (:lightgrey, 1.0)

laskar2004color = :black            # Insolation
lisiecki2005color = :steelblue      # d18O
hodell2023color = :darkorange       # d18O

luthi2008color = :slateblue1     # C
yamamoto2022color = :darkorange2 # C
berends2021color = :palegreen3  # C, SLE
spratt2016color = :orange1       # SLE
bintanja2008color = :dodgerblue1   # SLE, T
barker2011color = :red3             # T

# Variables names
time_label = L"Time (kyr BP)$\,$"
ins_label = L"$I$ (W⋅m$^{-2}$)"
periods_label = L"$P$ (kyr)"
d18O_label = L"δ$^{18}$O (‰) $\,$"
temp_label = L"$T - T_\mathrm{ref}$ (K)"
carbon_label = L"$C$ (ppm)"
icethick_label = L"$H$ (m)"
z_label = L"$z$ (m)"
icevol_label = L"$V_\mathrm{ice}$ (m SLE)"
sed_label = L"$H_\mathrm{sed}$ (m)"
elev_label = L"Elevation (m)$\,$"
mb_label = L"$\dot{m}$ (m⋅yr$^{-1}$)"
snowfall_label = L"$\dot{s}$ (m⋅yr$^{-1}$)"
ablation_label = L"$-\dot{a}$ (m⋅yr$^{-1}$)"
v_label = L"$v$ (m⋅yr$^{-1}$)"

insdot_label = L"$dI/dt$ (W⋅m$^{-2}$⋅yr$^{-1}$)"
icethickdot_label = L"$dH/dt$ (m⋅yr$^{-1}$)"
tempdot_label = L"$dT/dt$ (K⋅yr$^{-1}$)"
icevoldot_label = L"$dV_\mathrm{ice}/dt$ (m SLE⋅yr$^{-1}$)"
mbdot_label = L"$d \dot{m}/dt$ (m⋅yr$^{-2}$)"
vdot_label = L"$v$ (m⋅yr$^{-2}$)"

# Ticks and positions
ticks_periods = [20, 40, 100]

xticks_time_short = Int.(-1.5e3:200:200)    # kyr
xticks_time_long = Int.(-3e3:300:300)       # kyr
xticks_time_800 = Int.(-800:100:0)       # kyr

yticks_insol = [450, 500, 550]
yticks_d18O = [3.5, 4.0, 4.5]
yticks_temp = [-20, -10, 0]
yticks_carbon, yticks_carbon2 = [200, 250, 300], [100, 200, 300, 400]
yticks_ice = [0, 1000, 2000]
yticks_z = [0, 1000, 2000]
yticks_vol = [-120, -60, 0]
yticks_sed = [0, 15, 30]
yticks_snowfall = [0.0, 0.1, 0.2, 0.3, 0.4]
yticks_ablation = [-3, -2, -1, 0]
yticks_mb = [-2, -3, -1, 0, 1]
yticks_v = [0, 100, 200, 300]

yticks_icedot = [-1.0, -0.5, 0.0]
yticks_tempdot = [-0.002, 0.0, 0.002, 0.004, 0.006]
yticks_voldot = [-0.02, 0.0, 0.02]
yticks_insoldot = []
yticks_mbdot = []
yticks_vdot = []

vlines_mpt = [-1.25e3, -0.7e3]

# Axis limits
xlims_time_short = (-1.5e3, 0.0)            # kyr
xlims_time_long = (-3.0e3, 0.0)             # kyr
xlims_time_800 = (-800, 0.0)             # kyr

ylims_ins = (420, 600)
ylims_periods = (0, 155)
ylims_temp = (-22, 15)
ylims_ice, ylims_ice2 = (-100, 2200), (-100, 2500)
ylims_z = (-100, 2200)
ylims_carbon, ylims_carbon2 = (150, 350), (50, 500)
ylims_vol = (-150, 80)
ylims_sed = (-0.5, 40)
ylims_elev = (-300, 2000)
ylims_mb = (-3, 1)
ylims_v = (-10, 350)

# Panel labeling
a_label = L"(a)$\,$"
b_label = L"(b)$\,$"
c_label = L"(c)$\,$"
d_label = L"(d)$\,$"
e_label = L"(e)$\,$"
f_label = L"(f)$\,$"
g_label = L"(g)$\,$"
h_label = L"(h)$\,$"
i_label = L"(i)$\,$"
j_label = L"(j)$\,$"
k_label = L"(k)$\,$"
l_label = L"(l)$\,$"
m_label = L"(m)$\,$"
n_label = L"(n)$\,$"