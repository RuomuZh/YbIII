import numpy as np
import matplotlib.pyplot as plt

# from src.functions.math import func as mathf
# from src.functions.optics import func as optf
# from src.constants import const

import magpylib as magpy
# from magpylib.current import Loop

#------------------------------------------ 2D MOT B-field setup ------------------------------------------#

# Br = 1.42  # Tesla
# mu0 = 4*np.pi*1e-7
# M = Br / mu0

# magnet_1 = magpy.magnet.Cuboid(
#     position = (-5.3*0.01, -20*0.01, -4*0.01), dimension = (2*0.01, 5.5*0.01, 2*0.01), magnetization = (0, -M, 0))
# magnet_2 = magpy.magnet.Cuboid(
#     position = (-5.3*0.01, -20*0.01, 4*0.01), dimension = (2*0.01, 5.5*0.01, 2*0.01), magnetization = (0, -M, 0))
# magnet_3 = magpy.magnet.Cuboid(
#     position = (5.3*0.01, -20*0.01, -4*0.01), dimension = (2*0.01, 5.5*0.01, 2*0.01), magnetization = (0, M, 0))
# magnet_4 = magpy.magnet.Cuboid(
#     position = (5.3*0.01, -20*0.01, 4*0.01), dimension = (2*0.01, 5.5*0.01, 2*0.01), magnetization = (0, M, 0))

# magnets = magpy.Collection(magnet_1, magnet_2, magnet_3, magnet_4)

#------------------------------------------ 3D MOT B-field setup ------------------------------------------#

def mot3d_coil(I=35, config="AH", plot_obj=False):

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
    T = 3
    # Windings
    W = 2

    # Spacing between windings
    s = 1.37 # in mm
    
    # Width of coil
    L = T * (2) # in mm

    # Diameter
    # d - Inner diameter
    d = 66.167 # in mm

    # D - Outer diameter
    D = d + 2 * W * s # in mm
    #D = 370.84 # in mm

    # Error
    e = 0 # in mm

    # Position of coils
    z_1 = 17.78 + e
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

#------------------------------------------ Cavity coil B-field setup ------------------------------------------#

def cavity_coil(I=30, config="H", plot_obj=False):

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
    s = 5.588 # in mm
    
    # Width of coil
    L = T * (6) # in mm

    # Diameter
    # d - Inner diameter
    d = 273.05 # in mm

    # D - Outer diameter
    D = d + 2 * W * s # in mm
    D = 370.84 # in mm

    # Error
    e = 0 # in mm

    # Position of coils
    z_1 = 92.6225 + e
    z_2 = -z_1

    coil = magpy.Collection()

    # Note: Unlike the 3D MOT coil function, the cavity coil function starts populating turns from the middle of both coils.

    for i in range(0, 2*T, 2):
        for n in range(W):

            # Upper Coils
            winding1 = magpy.current.Circle(
                current = curr_up,
                diameter = d + (2*n + 1) * s,
                position = (248.3, 0, z_1 + (s)*((i - (T-1))/2)),
            )

            coil.add(winding1)

            # Lower Coils
            winding2 = magpy.current.Circle(
                current = curr_down,
                diameter = d + (2*n + 1) * s,
                position = (248.3, 0, z_2 + (s)*((i - (T-1))/2)),
            )

            coil.add(winding2)

    if plot_obj:
        coil.show(backend='plotly')

    return coil

#------------------------------------------ Coil objects ------------------------------------------#

cavitycoil = cavity_coil(config="H", plot_obj=False)
motcoil = mot3d_coil(config="AH", plot_obj=False)

#------------------------------------------ Coil + magnets visualization ------------------------------------------#

# combined = magpy.Collection()
# combined.add(cavitycoil)
# combined.add(motcoil)
# #combined.add(magnets)
# combined.show(backend= 'plotly')


#------------------------------------------ Simulating B-fields along particular axes ------------------------------------------#

xline = np.linspace(-300, 300, 1000)  # ±1 m
line = np.array([(x, 0, 0) for x in xline])

B_line = motcoil.getB(line) + cavitycoil.getB(line)  # shape (N, 3)
Bx_line, By_line, Bz_line = np.moveaxis(B_line, -1, 0)
Bmag = np.sqrt(Bx_line**2 + By_line**2 + Bz_line**2) * 1e4  # T→G

Bx_G = Bx_line*1e7 # Since in the function, we calculate B-fields scaled up by 1000 (m instead of mm) and in the units of T, we multiply 10^4 * 10^3 = 10^7
By_G = By_line*1e7
Bz_G = Bz_line*1e7
x_cm = xline/10

dBx_dx = np.gradient(Bx_G, x_cm)
dBy_dx = np.gradient(By_G, x_cm)
dBz_dx = np.gradient(Bz_G, x_cm)

#------------------------------------------ Plot B-fields ------------------------------------------#

plt.figure()
# plt.plot(xline * 100, Bmag)  # m→cm
plt.plot(x_cm, Bx_G, 'red', label=r'$B_x$')
plt.plot(x_cm, By_G, 'blue', label=r'$B_y$')
plt.plot(x_cm, Bz_G, 'green', label=r'$B_z$')
plt.xlabel("x (cm)")
plt.ylabel("B (G)")
plt.title("B Field along the conveyor axis")
plt.grid()

plt.legend(loc='best', fontsize=10)
plt.show()

#------------------------------------------ Plot B-field gradient ------------------------------------------#

# plt.figure()
# # plt.plot(xline * 100, Bmag)  # m→cm
# plt.plot(x_cm, dBx_dx, 'red', label=r'$dB_x/dx$')
# plt.plot(x_cm, dBy_dx, 'blue', label=r'$dB_y/dx$')
# plt.plot(x_cm, dBz_dx, 'green', label=r'$dB_z$/dx')
# plt.xlabel("x (cm)")
# plt.ylabel("dB/dx (G/cm)")
# plt.title("B Field Gradients along the conveyor axis")
# plt.grid()

# plt.legend(loc='best', fontsize=10)
# plt.show()


#------------------------------------------ Streamplot stuff ------------------------------------------#

# import matplotlib.pyplot as plt

# def plot_plane(plane="xy", min=[-300, -300], max=[300, 300], res=100, title="Magnetic Field Streamplot (XY plane)", save_plot=False, fname="mot3d_ah_i_200a", dir="results\\"):

#     fig, axs = plt.subplots(1, 1, figsize=(26,10))

#     # create grid
#     ax1_ts = np.linspace(min[0], max[0], res)
#     ax2_ts = np.linspace(min[0], max[1], res)
#     comp = [0, 1]

#     if plane=="xy":
#         grid = np.array([[(x, y, 0) for x in ax1_ts] for y in ax2_ts])
#         comp = [0, 1]
#     else:
#         if plane=="yz":
#             grid = np.array([[(0, y, z) for y in ax1_ts] for z in ax2_ts])
#             comp = [1, 2]
#         else:
#             if plane=="xz":
#                 grid = np.array([[(x, 0, z) for x in ax1_ts] for z in ax2_ts])
#                 comp = [0, 2]
#             else:
#                 raise ValueError("Invalid Plane")

#     # compute and plot field of coil
#     B = magpy.getB(coil, grid) *1e2*1e4 # mT to G
#     Bamp = np.linalg.norm(B, axis=2)
#     Bamp /= np.amax(Bamp)

#     sp = axs.streamplot(
#         grid[:,:,comp[0]], grid[:,:,comp[1]], B[:,:,comp[0]], B[:,:,comp[1]],
#         density=2,
#         color=np.where(np.linalg.norm(B, axis=2) > 100, 1000, 10*np.linalg.norm(B, axis=2)),
#         linewidth=np.sqrt(Bamp)*3,
#         cmap='jet',
#     )

#     # figure styling
#     axs.set(
#         title='Magnetic field of coils [G]',
#         xlabel=plane[0] + '-position [mm]',
#         ylabel=plane[1] + '-position [mm]',
#         aspect=1,
#     )

#     plt.colorbar(sp.lines, ax=axs, label='[G]')

#     if save_plot:
#         plt.savefig(dir + fname + '.png', bbox_inches='tight')

#     plt.tight_layout()
#     plt.show()


# plot_plane(plane="xy", fname="mot3d_ah_i_200a_xy", save_plot=False)
