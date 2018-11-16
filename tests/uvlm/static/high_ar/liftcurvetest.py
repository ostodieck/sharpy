"""
Test the lift curve slope of a high aspect ratio flat plate
"""

import numpy as np
import cases.templates.flying_wings as flyingwings
import os
import sharpy.utils.algebra as algebra
import sharpy.sharpy_main

# Flight conditions
u_inf = 50
alpha = 4
rho = 1.225
q_inf = 0.5*rho*u_inf**2
wing_orientation = algebra.euler2quat(np.array([0, alpha*np.pi/180, 0]))


# AR_range = np.arange(10,100,10)
AR_range = [10,20]

CL = np.zeros(len(AR_range))
CL_alpha = np.zeros(len(AR_range))
K = np.zeros(len(AR_range))
iter = 0

for AR in AR_range:

    # Wing properties
    chord = 1
    span = AR*chord

    S_ref = span*chord


    # Number of bound panels
    m = 10
    ar_panel = 1

    ar_ratio = AR // ar_panel

    M, N, Mstarfact = m, ar_ratio*m, 15
    K[iter] = M*N
    Nsurfaces = 2

    # Case name
    case_name = "flat_plate_K%.4d_a%.4d_uinf%.4d_AR%.2d" %(K[iter], alpha*100, u_inf, AR)
    route = os.path.abspath('.') + '/cases/'
    dir_name = route+case_name+'/'
    os.system('rm -rf dir_name')
    os.system('mkdir -p %s' %dir_name)

    # Create wing
    wing = flyingwings.FlyingWing(M = M,
                                  N=N,
                                  Mstar_fact=Mstarfact,
                                  alpha=alpha,
                                  u_inf=u_inf,
                                  main_chord=chord,
                                  aspect_ratio=AR,
                                  b_ref=span,
                                  rho=rho,
                                  n_surfaces=Nsurfaces,
                                  route = dir_name,
                                  case_name='flat_plate')
    # wing.root_airfoil_P = 0
    # wing.root_airfoil_M = 0
    # wing.tip_airfoil_P = 0
    # wing.tip_airfoil_M = 0

    wing.clean_test_files()
    wing.update_derived_params()
    wing.generate_fem_file()
    wing.generate_aero_file()

    ### solution flow
    wing.set_default_config_dict()
    wing.config['SHARPy']['flow'] = [
        'BeamLoader', 'AerogridLoader', 'StaticUvlm', 'BeamPlot',
        'AerogridPlot', 'SaveData', 'AeroForcesCalculator']
    wing.config['SaveData'] = {'folder': dir_name}
    wing.config['AerogridLoader']['freestream_dir'] = [1., 0., 0.]
    # wing.config['AerogridLoader']['mstar'] = 1
    wing.config['SHARPy']['write_screen'] = False
    wing.config['BeamLoader']['orientation'] = wing_orientation
    wing.config['StaticUvlm']['velocity_field_input']['u_inf_direction'] = [1., 0, 0]
    wing.config['StaticUvlm']['horseshoe'] = 'off'
    wing.config['StaticCoupled']['aero_solver_settings']['velocity_field_input']['u_inf_direction'] = [1., 0., 0.]
    if M < 5:
        wing.config['StaticCoupled']['tolerance'] = 1e-7
    wing.config.write()

    ### solve
    data = sharpy.sharpy_main.main(['path_to_solver_useless',
                                    dir_name + 'flat_plate' + '.solver.txt'])
    fz_aero = 0
    fx_aero = 0

    for surf in range(Nsurfaces):
        fz_aero += data.aero.timestep_info[0].forces[0][2, :, :].sum()
        fx_aero += data.aero.timestep_info[0].forces[0][0, :, :].sum()

    CL[iter] = fz_aero/q_inf/S_ref
    CL_alpha[iter] = CL[iter]/alpha/np.pi*180

    print("AR = %.2d" %AR)
    print("C_L = %.3f" %CL[iter])
    print("C_L_alpha = %.3f" %CL_alpha[iter])
    iter += 1