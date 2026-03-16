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

# - setup atom
atom = Ytterbium()
# get transition, to help setting up lasers
main = atom.trans["main"] # 399 nm transition to 1P1
intercombination = atom.trans["intercombination"] # 556 nm transition to 3P1

# - setup laser
laser_1 = GaussianLaserBeam(
    wavelength=intercombination.wavelength,
    waist=2e-3,
    waist_position=(0, 0, 0),
    direction=(0, -np.sin(0), np.cos(0)),
)
laser_1.set_power_from_I(intercombination.Isat) # set power to reach Isat
laser_1.polarization = Vertical()
laser_1.tag = "las1"

# - create the magpylib object (in this case 2D MOT magnets)
Br = 1.42  # Tesla
mu0 = 4*np.pi*1e-7
M = Br / mu0
magnet_1 = magpy.magnet.Cuboid(
    position = (-5.3*0.01, 0, -4*0.01), dimension = (2*0.01, 5.5*0.01, 2*0.01), magnetization = (0, -M, 0))
magnet_2 = magpy.magnet.Cuboid(
    position = (-5.3*0.01, 0, 4*0.01), dimension = (2*0.01, 5.5*0.01, 2*0.01), magnetization = (0, -M, 0))
magnet_3 = magpy.magnet.Cuboid(
    position = (5.3*0.01, 0, -4*0.01), dimension = (2*0.01, 5.5*0.01, 2*0.01), magnetization = (0, M, 0))
magnet_4 = magpy.magnet.Cuboid(
    position = (5.3*0.01, 0, 4*0.01), dimension = (2*0.01, 5.5*0.01, 2*0.01), magnetization = (0, M, 0))

magnets = magpy.Collection(magnet_1, magnet_2, magnet_3, magnet_4)

# wrap it up
mag_field_2D_MOT = MagpylibWrapper(magnets)
mag_field_2D_MOT.tag = "2D MOT coils"

# - create the magpylib object (in this case 3D MOT coils)
def mot3d_coil(I=30, config="H", plot_obj=False):

    if config=="H":
        curr_up = -I
        curr_down = -I
    else:
        if config=="AH":
            curr_up = -I
            curr_down = I
        else:
            raise ValueError("Invalid 3D MOT Coil Configuration")

    # Turns
    T = 9
    # Windings
    W = 8

    # Spacing between windings
    s = 5.588e-3 # in mm
    
    # Width of coil
    L = T * (6e-3) # in mm

    # Diameter
    # d - Inner diameter
    d = 273.05e-3 # in mm

    # D - Outer diameter
    D = d + 2 * W * s # in mm
    D = 370.84e-3 # in mm

    # Error
    e = 0 # in mm

    # Position of coils
    y_1 = -92.6225e-3 + e
    y_2 = -y_1

    coil = magpy.Collection()

    for i in range(0, 2*T, 2):
        for n in range(W):

            # Upper Coils
            winding1 = magpy.current.Circle(
                current = curr_up,
                diameter = d + (2*n + 1) * s,
                position = (0, (y_1 + (s)*((i - (T-1))/2)), 1),
                orientation=R.from_euler('x', 90, degrees=True)
            )

            coil.add(winding1)

            # Lower Coils
            winding2 = magpy.current.Circle(
                current = curr_down,
                diameter = d + (2*n + 1) * s,
                position = (0, (y_2 + (s)*((i - (T-1))/2)), 1),
                orientation=R.from_euler('x', 90, degrees=True)
            )

            coil.add(winding2)

    if plot_obj:
        coil.show(backend='plotly')

    return coil

coil = mot3d_coil(config="AH", plot_obj=False)

# wrap it up
mag_field_3D_MOT = MagpylibWrapper(coil)
mag_field_3D_MOT.tag = "3D MOT coils"

# Visualization to check magnet geometry and placement
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

# let's add gravity, pointing along +y
m = Ytterbium().mass  # kg
g = 9.81  # m/s^2
grav_force = (0,  m * g, 0)
gravity = ConstantForce(field_value=grav_force, tag="gravity")

# - config
config = Configuration()
config.atom = atom
config += laser_1, gravity, mag_field_2D_MOT
config.add_atomlight_coupling("las1", "intercombination", 0) # Arguments: laser = "las1", transition = intercombination", detuning = - or +200 * intercombination.Gamma
# - simulation
sim = RK4(config=config)
t = np.linspace(0, 0.25, 6000) # timesteps for integration
u0 = (0, 0, 0, 0, 0, 0) # atom starts with vz=100m/s
res = sim.integrate(u0, t)

# plot
fix, axes = plt.subplots(3, 2, tight_layout=True)

axes[0, 0].plot(res.t * 1e3, res.y[0])
axes[0, 0].set_ylabel("x␣(m)")
axes[0, 1].plot(res.t * 1e3, res.y[3])
axes[0, 1].set_ylabel("vx␣(m/s)")

axes[1, 0].plot(res.t * 1e3, res.y[1])
axes[1, 0].set_ylabel("y␣(m)")
axes[1, 1].plot(res.t * 1e3, res.y[4])
axes[1, 1].set_ylabel("vy␣(m/s)")

axes[2, 0].plot(res.t * 1e3, res.y[2])
axes[2, 0].set_ylabel("z␣(m)")
axes[2, 1].plot(res.t * 1e3, res.y[5])
axes[2, 1].set_ylabel("vz␣(m/s)")

diff = res.y[1] - 0.002
sign_changes = np.where(np.diff(np.sign(diff)))[0]
print(res.t[sign_changes]*1e3)

for ax in axes.flat:
    ax.axvline(x=res.t[sign_changes]*1e3, color='r', linestyle='--')
    ax.set(xlabel="t␣(ms)")

y = np.interp(6e-3, res.t, res.y[1])
vy = np.interp(6e-3, res.t, res.y[4])
z = np.interp(6e-3, res.t, res.y[2])
vz = np.interp(6e-3, res.t, res.y[5])

print(y, vy, z, vz)

plt.show()

