import numpy as np
import matplotlib.pyplot as plt
import magpylib as magpy

from scipy.spatial.transform import Rotation as R

from atomsmltr.environment import PlaneWaveLaserBeam
from atomsmltr.environment import GaussianLaserBeam
from atomsmltr.environment import ConstantForce
from atomsmltr.atoms import Ytterbium
from atomsmltr.simulation import Configuration, RK4
from atomsmltr.environment import (
    Vertical,
    Horizontal,
    CircularLeft,
    CircularRight,
    Linear,
    Vector,
)
from atomsmltr.environment.fields.magnetic.magpylib import MagpylibWrapper

#------------------------------------------ Atom setup ------------------------------------------#

atom = Ytterbium()
# get transition, to help setting up lasers
main = atom.trans["main"] # 399 nm transition to 1P1
intercombination = atom.trans["intercombination"] # 556 nm transition to 3P1

#------------------------------------------ Green push beam setup ------------------------------------------#

laser_1 = GaussianLaserBeam(
    wavelength=intercombination.wavelength,
    waist=2e-3,
    waist_position=(0, 0, -0.2),
    direction=(0, 0, 1),
)

laser_1.set_power_from_I(intercombination.Isat) # set power to reach Isat
laser_1.polarization = Vertical()
laser_1.tag = "las1"

#------------------------------------------ 2D MOT B-field setup ------------------------------------------#

Br = 1.42  # Tesla
mu0 = 4*np.pi*1e-7
M = Br / mu0

magnet_1 = magpy.magnet.Cuboid(
    position = (-5.3*0.01, -4*0.01, 0), dimension = (2*0.01, 2*0.01, 5.5*0.01), magnetization = (0, 0, M))
magnet_2 = magpy.magnet.Cuboid(
    position = (-5.3*0.01, 4*0.01, 0), dimension = (2*0.01, 2*0.01, 5.5*0.01), magnetization = (0, 0, M))
magnet_3 = magpy.magnet.Cuboid(
    position = (5.3*0.01, -4*0.01, 0), dimension = (2*0.01, 2*0.01, 5.5*0.01), magnetization = (0, 0, -M))
magnet_4 = magpy.magnet.Cuboid(
    position = (5.3*0.01, 4*0.01, 0), dimension = (2*0.01, 2*0.01, 5.5*0.01), magnetization = (0, 0, -M))

magnets = magpy.Collection(magnet_1, magnet_2, magnet_3, magnet_4)

# Rotate so the extraction axis points UP (+z)
rot = R.from_euler('x', 90, degrees=True)
magnets.rotate(rot, anchor=(0, 0, 0))

# Shift the whole 2D MOT assembly 20 cm below the 3D MOT chamber center
magnets.move((0, 0, -0.2))

# wrap it up
mag_field_2D_MOT = MagpylibWrapper(magnets)
mag_field_2D_MOT.tag = "2D MOT coils"

#magnets.show(backend='plotly')

#------------------------------------------ 3D MOT B-field setup ------------------------------------------#

def mot3d_coil(I=35, config="AH", plot_obj=False):

    if config=="H":
        curr_up = -I
        curr_down = -I
    else:
        if config=="AH":
            curr_up = I
            curr_down = -I
        else:
            raise ValueError("Invalid 3D MOT Coil Configuration")

    # Turns
    T = 3
    # Windings
    W = 2

    # Spacing between windings
    s = 1.37 * 0.001 # in mm

    # Diameter
    # d - Inner diameter
    d = 66.167 * 0.001 # in mm

    # Error
    e = 0 # in mm

    # Position of coils
    z_1 = 17.78 * 0.001 + e
    z_2 = -z_1

    coil = magpy.Collection()

    # Note: Unlike the cavity coil function, the 3D MOT function starts populating turns from the bottom of the top coil, and the top of the bottom coil.

    for i in range(0, T):
        for n in range(W):

            # Upper Coils
            winding1 = magpy.current.Circle(
                current = curr_up,
                diameter = d + (2*n + 1) * s,
                position = (0, 0, z_1 + i*s),
            )

            coil.add(winding1)

            # Lower Coils
            winding2 = magpy.current.Circle(
                current = curr_down,
                diameter = d + (2*n + 1) * s,
                position = (0, 0, z_2 - i*s),
            )

            coil.add(winding2)

    if plot_obj:
        coil.show(backend='plotly')

    return coil

motcoil = mot3d_coil(config="AH", plot_obj=False)

# ROTATION: Rotate the 3D MOT coils if needed
# For example, rotating 90 degrees around the X-axis 
# (This would flip the coil symmetry axis from the Z-axis to the Y-axis)
rot_3d_mot = R.from_euler('x', 90, degrees=True)
motcoil.rotate(rot_3d_mot, anchor=(0, 0, 0))

# wrap it up
mag_field_3D_MOT = MagpylibWrapper(motcoil)
mag_field_3D_MOT.tag = "3D MOT coils"

#------------------------------------------ Cavity coil B-field setup ------------------------------------------#

# def cavity_coil(I=30, config="H", plot_obj=False):

#     if config=="H":
#         curr_up = -I
#         curr_down = -I
#     else:
#         if config=="AH":
#             curr_up = -I
#             curr_down = I
#         else:
#             raise ValueError("Invalid 3D MOT Coil Configuration")

#     # Turns
#     T = 9
#     # Windings
#     W = 8

#     # Spacing between windings
#     s = 5.588 * 0.001 # in mm

#     # Diameter
#     # d - Inner diameter
#     d = 273.05 * 0.001 # in mm

#     # Error
#     e = 0 # in mm

#     # Position of coils
#     z_1 = 92.6225 * 0.001 + e
#     z_2 = -z_1

#     coil = magpy.Collection()

#     # Note: Unlike the 3D MOT coil function, the cavity coil function starts populating turns from the middle of both coils.

#     for i in range(0, 2*T, 2):
#         for n in range(W):

#             # Upper Coils
#             winding1 = magpy.current.Circle(
#                 current = curr_up,
#                 diameter = d + (2*n + 1) * s,
#                 position = (0.2483, 0, z_1 + (s)*((i - (T-1))/2)),
#             )

#             coil.add(winding1)

#             # Lower Coils
#             winding2 = magpy.current.Circle(
#                 current = curr_down,
#                 diameter = d + (2*n + 1) * s,
#                 position = (0.2483, 0, z_2 + (s)*((i - (T-1))/2)),
#             )

#             coil.add(winding2)

#     if plot_obj:
#         coil.show(backend='plotly')

#     return coil

# cavitycoil = cavity_coil(config="H", plot_obj=False)

# # wrap it up
# mag_field_cavity = MagpylibWrapper(cavitycoil)
# mag_field_cavity.tag = "Cavity coils"

#------------------------------------------ 3D MOT lasers setup ------------------------------------------#

# Create the 6 beams for a 3D MOT
# 1 & 2: Vertical Pair (+z, -z)
las_z_up = GaussianLaserBeam(wavelength=intercombination.wavelength, waist=6e-3, waist_position=(0, 0, 0), direction=(0,0,1))
las_z_down = GaussianLaserBeam(wavelength=intercombination.wavelength, waist=6e-3, waist_position=(0, 0, 0), direction=(0,0,-1))

# 3 & 4: Horizontal Pair 1 (+x, -x)
las_x_plus = GaussianLaserBeam(wavelength=intercombination.wavelength, waist=6e-3, waist_position=(0, 0, 0), direction=(1,0,0))
las_x_minus = GaussianLaserBeam(wavelength=intercombination.wavelength, waist=6e-3, waist_position=(0, 0, 0), direction=(-1,0,0))

# 5 & 6: Horizontal Pair 2 (+y, -y)
las_y_plus = GaussianLaserBeam(wavelength=intercombination.wavelength, waist=6e-3, waist_position=(0, 0, 0), direction=(0,1,0))
las_y_minus = GaussianLaserBeam(wavelength=intercombination.wavelength, waist=6e-3, waist_position=(0, 0, 0), direction=(0,-1,0))

# Set power (can increase if we want power broadening)
las_z_up.set_power_from_I(5 * intercombination.Isat)
las_z_down.set_power_from_I(5 * intercombination.Isat)
las_x_plus.set_power_from_I(5 * intercombination.Isat)
las_x_minus.set_power_from_I(5 * intercombination.Isat)
las_y_plus.set_power_from_I(5 * intercombination.Isat)
las_y_minus.set_power_from_I(5 * intercombination.Isat)

# Assign Polarizations (Must be Circular and match the anti-Helmholtz field signs!)
las_z_up.polarization = CircularRight()
las_z_down.polarization = CircularRight()
las_x_plus.polarization = CircularRight()
las_x_minus.polarization = CircularRight()
las_y_plus.polarization = CircularLeft()
las_y_minus.polarization = CircularLeft()

# Tag them for coupling identification
las_z_up.tag = "z_up"; las_z_down.tag = "z_down"
las_x_plus.tag = "x_plus"; las_x_minus.tag = "x_minus"
las_y_plus.tag = "y_plus"; las_y_minus.tag = "y_minus"

#------------------------------------------ Gravity ------------------------------------------#

m = Ytterbium().mass  # kg
g = 9.81  # m/s^2
grav_force = (0, 0, -m*g)
gravity = ConstantForce(field_value=grav_force, tag="gravity")

#------------------------------------------ config ------------------------------------------#

config = Configuration()
config.atom = atom
config += laser_1, gravity, mag_field_3D_MOT, mag_field_2D_MOT
config += las_z_up, las_z_down, las_x_plus, las_x_minus, las_y_plus, las_y_minus

# # Set a red-detuning relative to the intercombination transition
detuning_push = 0
config.add_atomlight_coupling("las1", "intercombination", detuning_push) # Arguments: laser = "las1", transition = intercombination", detuning = detuning_push

# Tie all 6 lasers to the intercombination transition simultaneously
detuning_3D_mot = -intercombination.Gamma/2 # red-detuning for trapping
for tag in ["z_up", "z_down", "x_plus", "x_minus", "y_plus", "y_minus"]:
    config.add_atomlight_coupling(tag, "intercombination", detuning_3D_mot)

# - simulation
sim = RK4(config=config)
t = np.linspace(0, 2, 6000) # timesteps for integration
u0 = (0, 0, -0.2, 0, 0, 0) # atom starts with vz=100m/s
res = sim.integrate(u0, t)

# plot
fix, axes = plt.subplots(2, 2, tight_layout=True)

mu_B = 9.274e-24
hbar = 1.055e-34
gJ = 1  # depends on state!
mJ = -1  # choose a sublevel

positions = res.y[:3].T
B_total = mag_field_2D_MOT.get_value(positions) + mag_field_3D_MOT.get_value(positions)
B_mag = np.linalg.norm(B_total, axis=1)

zeeman = (mu_B * gJ * mJ / hbar) * B_mag
doppler = -(2 * np.pi / intercombination.wavelength) * res.y[5]

axes[0, 0].plot(res.t * 1e3, (doppler+zeeman)/(2*np.pi*1e6))
axes[0, 0].set_ylabel("Detuning (MHz)")
axes[0, 1].plot(res.t * 1e3, B_mag * 1e4)
axes[0, 1].set_ylabel("Total B-field Magnitude (G)")

axes[1, 0].plot(res.t * 1e3, res.y[2])
axes[1, 0].set_ylabel("z␣(m)")
axes[1, 1].plot(res.t * 1e3, res.y[5])
axes[1, 1].set_ylabel("vz␣(m/s)")

plt.show()

#------------------------------------------ Coil + magnets visualization ------------------------------------------#

# combined = magpy.Collection()
# combined.add(cavitycoil)
# combined.add(motcoil)
# combined.add(magnets)
# combined.show(backend= 'plotly')

#------------------------------------------ Visualization to check magnet geometry and placement ------------------------------------------#
# magnets_final = magpy.Collection(magnets, coil)
# magnets_final.show(backend='plotly')

# zline = np.linspace(0, 0.2, 1000)  # ±10 cm
# line = np.array([(0, 0, z) for z in zline])

# B_line = magnets.getB(line)   # shape (N, 3)
# Bx_line, By_line, Bz_line = np.moveaxis(B_line, -1, 0)
# Bmag = np.sqrt(Bx_line**2 + By_line**2 + Bz_line**2) * 1e4  # T→G

# plt.figure()
# # plt.plot(xline * 100, Bmag)  # m→cm
# plt.plot(zline * 100, Bx_line*1e4)
# plt.plot(zline * 100, By_line*1e4)
# plt.plot(zline * 100, Bz_line*1e4)
# plt.xlabel("x (cm)")
# plt.ylabel("B (G)")
# plt.title("B along physical x-axis (y = 0, z = 0)")
# plt.grid()
# plt.show()

# Appendix 1

# axes[0, 0].plot(res.t * 1e3, res.y[1])
# axes[0, 0].set_ylabel("y␣(m)")
# axes[0, 1].plot(res.t * 1e3, res.y[4])
# axes[0, 1].set_ylabel("vy␣(m/s)")

# Appendix 2

# diff = res.y[1] - 0.002
# sign_changes = np.where(np.diff(np.sign(diff)))[0]
# print(res.t[sign_changes]*1e3)

# for ax in axes.flat:
#     ax.axvline(x=res.t[sign_changes]*1e3, color='r', linestyle='--')
#     ax.set(xlabel="t␣(ms)")

# y = np.interp(6e-3, res.t, res.y[1])
# vy = np.interp(6e-3, res.t, res.y[4])
# z = np.interp(6e-3, res.t, res.y[2])
# vz = np.interp(6e-3, res.t, res.y[5])

# print(y, vy, z, vz)
