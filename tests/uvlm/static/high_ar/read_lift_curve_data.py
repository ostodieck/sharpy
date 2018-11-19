"""
Read Lift Curve Test Data
"""

import numpy as np
import os
import sharpy.utils.h5utils as h5
import sharpy.utils.algebra as algebra
import matplotlib.pyplot as plt


def comp_total_forces(data, moment_reference, AR):
    """Compute the total aerodynamic forces coefficients from SHARPy data
    """

    # Initialise force and moment vector
    ftot = np.zeros((3,))
    mtot = np.zeros((3,))

    rho = data.settings['StaticUvlm']['rho']
    u_inf = data.settings['StaticUvlm']['velocity_field_input']['u_inf']
    span = np.max(data.structure.timestep_info[0].pos[:, 1])*2
    chord = span/AR

    aeroforces = data.aero.timestep_info[0].forces

    for surface in range(len(aeroforces)):

        _, chord_cells, span_cells = aeroforces[surface].shape

        # Sum forces
        for axis_dir in range(3):
            ftot[axis_dir] += aeroforces[surface][axis_dir, :, :].sum()

        # Calculate moments
        # Aerogrid  coordinates
        zeta = data.aero.timestep_info[0].zeta[surface]
        # Calculate moments with respect to reference location
        for cell_x in range(chord_cells):
            for cell_y in range(span_cells):
                moment_arm = zeta[:, cell_x, cell_y] - moment_reference
                mtot += np.cross(moment_arm, aeroforces[surface][:3, cell_x, cell_y])

    cftot = ftot/0.5/rho/u_inf**2/span/chord
    cmtot = mtot/0.5/rho/u_inf**2/span/chord**2

    return cftot, cmtot

M = 10
alpha = 4
u_inf = 50
AR_range = np.arange(10, 100, 10)

iter = 0

cftot = np.zeros((len(AR_range), 3))
cmtot = np.zeros((len(AR_range), 3))
CL_alpha = np.zeros(len(AR_range))

for AR in AR_range:

    # Case Name
    case_name = "flat_plate_K%.4d_a%.4d_uinf%.4d_AR%.2d" %(M*AR*10, alpha*100, u_inf, AR)
    route = os.path.abspath('.') + '/cases/'
    dir_name = route+case_name+'/'

    data = h5.readh5(dir_name+'flat_plate.data.h5').data

    cftot[iter,:], cmtot[iter, :] = comp_total_forces(data, [0,0,0], AR)

    CL_alpha[iter] = cftot[iter,2]/alpha/np.pi*180

    iter += 1

fig = plt.figure(figsize=(12,10))

plt.plot(AR_range, CL_alpha,
         lw=2,
         color='k')

AR_range_fine = np.linspace(np.min(AR_range), np.max(AR_range), 50)
AR_theory = np.pi*2*AR_range_fine/(AR_range_fine+2)

plt.plot(AR_range_fine, AR_theory,
         lw=2,
         ls='--')

plt.xlabel('Aspect Ratio, AR')
plt.ylabel('Lift Curve Slope, $C_{L_alpha}')
plt.grid()

plt.show()
