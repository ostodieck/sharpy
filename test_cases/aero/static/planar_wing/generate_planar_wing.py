import h5py as h5
import numpy as np
import configparser

aspect_ratio = 10.0
chord = 1.0
length = aspect_ratio/chord

def clean_test_files(route, case_name):
    fem_file_name = route + '/' + case_name + '.fem.h5'
    if os.path.isfile(fem_file_name):
        os.remove(fem_file_name)

    aero_file_name = route + '/' + case_name + '.aero.h5'
    if os.path.isfile(aero_file_name):
        os.remove(aero_file_name)

    solver_file_name = route + '/' + case_name + '.solver.txt'
    if os.path.isfile(solver_file_name):
        os.remove(solver_file_name)

    flightcon_file_name = route + '/' + case_name + '.flightcon.txt'
    if os.path.isfile(flightcon_file_name):
        os.remove(flightcon_file_name)


def generate_files(route, case_name, num_elem=10, num_node_elem=3):
    clean_test_files(route, case_name)
    num_node, coordinates = generate_fem_file(route,
                                              case_name,
                                              num_elem,
                                              num_node_elem
                                              )
    generate_aero_file(route, case_name, num_elem, num_node, coordinates)
    generate_solver_file(route, case_name)
    generate_flightcon_file(route, case_name)


def generate_fem_file(route, case_name, num_elem, num_node_elem=3):

    num_node = (num_node_elem - 1)*num_elem + 1
    # import pdb; pdb.set_trace()
    angle = 90*np.pi/180.0
    dihedral = 0*np.pi/180.0
    n_elem_beam = int(num_elem/2)
    n_node_beam = int((num_node - 1)/2)

    x = np.zeros((num_node,))
    y = np.zeros((num_node,))
    z = np.zeros((num_node,))
    x[0:n_node_beam+1] = (np.linspace(0, length, n_node_beam+1))*np.cos(angle)
    y[0:n_node_beam+1] = (np.linspace(0, length, n_node_beam+1))*np.sin(angle)*np.cos(dihedral)
    z[0:n_node_beam+1] = (np.linspace(0, length, n_node_beam+1))*np.sin(dihedral)
    x[n_node_beam+1:] = ((np.linspace(0, -length, n_node_beam+1))*np.cos(angle))[1:]
    y[n_node_beam+1:] = ((np.linspace(0, -length, n_node_beam+1))*np.sin(angle)*np.cos(dihedral))[1:]
    z[n_node_beam+1:] = ((np.linspace(0, -length, n_node_beam+1))*np.sin(dihedral))[1:]

    structural_twist = np.zeros_like(x)

    frame_of_reference_delta = np.zeros((num_node, 3))
    for inode in range(num_node):
        # frame_of_reference_delta[inode, :] = [0, 1, 0]
        if inode < n_node_beam + 1:
            frame_of_reference_delta[inode, :] = [-np.sin(angle), np.cos(angle), 0]
        else:
            frame_of_reference_delta[inode, :] = [np.sin(angle), -np.cos(angle), 0]

    scale = 1

    x *= scale
    y *= scale
    z *= scale

    conn = np.zeros((num_elem, num_node_elem), dtype=int)
    for ielem in range(num_elem):
        conn[ielem, :] = (np.ones((3,)) * ielem * (num_node_elem - 1)
                          + [0, 2, 1])
        if ielem == n_elem_beam:
            conn[ielem, 0] = 0

    # stiffness array
    # import pdb; pdb.set_trace()
    num_stiffness = 1
    ea = 4.8e8
    ga = 3.231e8
    gj = 1.0e6
    ei = 9.346e6
    base_stiffness = np.diag([ea, ga, ga, gj, ei, ei])
    stiffness = np.zeros((num_stiffness, 6, 6))
    # import pdb; pdb.set_trace()
    for i in range(num_stiffness):
        stiffness[i, :, :] = base_stiffness

    # element stiffness
    elem_stiffness = np.zeros((num_elem,), dtype=int)

    # mass array
    num_mass = 1
    m_bar = 100
    j = 10
    base_mass = np.diag([m_bar, m_bar, m_bar, j, j, j])
    mass = np.zeros((num_mass, 6, 6))
    for i in range(num_mass):
        mass[i, :, :] = base_mass
    # element masses
    elem_mass = np.zeros((num_elem,), dtype=int)

    # bocos
    boundary_conditions = np.zeros((num_node, 1), dtype=int)
    boundary_conditions[0] = 1
    boundary_conditions[n_node_beam + 1] = -1
    boundary_conditions[-1] = -1

    # beam number
    beam_number = np.zeros((num_elem, 1), dtype=int)
    beam_number[n_elem_beam:] = 1

    # new app forces scheme (only follower)
    n_app_forces = 0
    node_app_forces = np.array([])
    app_forces = np.zeros((n_app_forces, 6))
    # app_forces[0, :] = [0, 0, 0, 0, 0, -11744.5275328e3]

    # lumped masses input
    n_lumped_mass = 0
    lumped_mass_nodes = np.array([], dtype=int)
    lumped_mass = np.zeros((n_lumped_mass, ))
    # lumped_mass[0] = 600e3/9.81
    lumped_mass_inertia = np.zeros((n_lumped_mass, 3, 3))
    lumped_mass_position = np.zeros((n_lumped_mass, 3))

    with h5.File(route + '/' + case_name + '.fem.h5', 'a') as h5file:
        coordinates = h5file.create_dataset('coordinates', data = np.column_stack((x, y, z)))
        conectivities = h5file.create_dataset('connectivities', data = conn)
        num_nodes_elem_handle = h5file.create_dataset(
            'num_node_elem', data = num_node_elem)
        num_nodes_handle = h5file.create_dataset(
            'num_node', data = num_node)
        num_elem_handle = h5file.create_dataset(
            'num_elem', data = num_elem)
        stiffness_db_handle = h5file.create_dataset(
            'stiffness_db', data = stiffness)
        stiffness_handle = h5file.create_dataset(
            'elem_stiffness', data = elem_stiffness)
        mass_db_handle = h5file.create_dataset(
            'mass_db', data = mass)
        mass_handle = h5file.create_dataset(
            'elem_mass', data = elem_mass)
        frame_of_reference_delta_handle = h5file.create_dataset(
            'frame_of_reference_delta', data=frame_of_reference_delta)
        structural_twist_handle = h5file.create_dataset(
            'structural_twist', data=structural_twist)
        bocos_handle = h5file.create_dataset(
            'boundary_conditions', data=boundary_conditions)
        beam_handle = h5file.create_dataset(
            'beam_number', data=beam_number)
        app_forces_handle = h5file.create_dataset(
            'app_forces', data=app_forces)
        node_app_forces_handle = h5file.create_dataset(
            'node_app_forces', data=node_app_forces)
        lumped_mass_nodes_handle = h5file.create_dataset(
            'lumped_mass_nodes', data=lumped_mass_nodes)
        lumped_mass_handle = h5file.create_dataset(
            'lumped_mass', data=lumped_mass)
        lumped_mass_inertia_handle = h5file.create_dataset(
            'lumped_mass_inertia', data=lumped_mass_inertia)
        lumped_mass_position_handle = h5file.create_dataset(
            'lumped_mass_position', data=lumped_mass_position)
    return num_node, coordinates


def generate_aero_file(route, case_name, num_elem, num_node, coordinates):
    # example airfoil
    num_node = (3 - 1)*num_elem + 1
    # import pdb; pdb.set_trace()
    n_elem_beam = int(num_elem/2)
    n_node_beam = int((num_node - 1)/2)

    naca_x, naca_y = generate_naca_camber(P=0, M=0)
    # airfoil distribution
    airfoil_distribution = []
    for i in range(num_node):
        airfoil_distribution.append(0)

    surface_distribution = np.zeros((num_elem,), dtype=int)
    surface_distribution[n_elem_beam:] = 1

    surface_m = np.zeros((2,), dtype=int)
    surface_m[0] = 8
    surface_m[1] = 8

    m_distribution = 'uniform'

    aero_node = np.ones(num_node, dtype=bool)

    # twist distribution
    twist = np.linspace(0, 0, num_node)*np.pi/180

    # chord distribution
    chord_dist = chord*np.ones((num_node,), dtype=float)

    # elastic axis distribution
    elastic_axis = 0.5*np.ones((num_node,))

    # import pdb; pdb.set_trace()
    with h5.File(route + '/' + case_name + '.aero.h5', 'a') as h5file:
        airfoils_group = h5file.create_group('airfoils')
        # add one airfoil
        naca_airfoil = airfoils_group.create_dataset('0', data=np.column_stack((naca_x, naca_y)))

        # chord
        chord_input = h5file.create_dataset('chord', data = chord_dist)
        dim_attr = chord_input .attrs['units'] = 'm'

        # twist
        twist_input = h5file.create_dataset('twist', data = twist)
        dim_attr = twist_input.attrs['units'] = 'rad'

        # airfoil distribution
        airfoil_distribution_input = h5file.create_dataset('airfoil_distribution', data=airfoil_distribution)

        surface_distribution_input = h5file.create_dataset('surface_distribution', data=surface_distribution)
        surface_m_input = h5file.create_dataset('surface_m', data=surface_m)
        m_distribution_input = h5file.create_dataset('m_distribution', data=m_distribution.encode('ascii', 'ignore'))

        aero_node_input = h5file.create_dataset('aero_node', data=aero_node)
        elastic_axis_input = h5file.create_dataset('elastic_axis', data=elastic_axis)


def generate_naca_camber(M=0, P=0):
    m = M/100
    p = P/10

    def naca(x, m, p):
        if x < 1e-6:
            return 0.0
        elif x < p:
            return m/(p*p)*(2*p*x - x*x)
        elif x > p and x < 1+1e-6:
            return m/((1-p)*(1-p))*(1 - 2*p + 2*p*x - x*x)

    # import pdb; pdb.set_trace()
    x_vec = np.linspace(0, 1, 1000)
    y_vec = np.array([naca(x, m, p) for x in x_vec])

    mat = np.column_stack((x_vec, y_vec))
    # np.savetxt(route + '/naca.csv', mat, delimiter=',')

    return x_vec, y_vec


def generate_solver_file(route, case_name):
    file_name = route + '/' + case_name + '.solver.txt'
    config = configparser.ConfigParser()
    config['SHARPy'] = {'case': 'planar_wing',
                        'route': './tests/aero/static/planar_wing',
                        'flow': 'StaticUvlm, BeamPlot, AeroGridPlot, AeroForcesSteadyCalculator',
                        'plot': 'on'}
    config['StaticUvlm'] = {'print_info': 'on',
                            'M_distribution': 'uniform',
                            'Mstar': 1,
                            'rollup': 'off',
                            'aligned_grid': 'on',
                            'prescribed_wake': 'off'}
    config['AeroGridPlot'] = {'route': './tests/aero/static/planar_wing/output',
                              'on_screen': 'off'}
    config['BeamPlot'] = {'route': './tests/aero/static/planar_wing/output',
                              'on_screen': 'off'}
    config['AeroForcesSteadyCalculator'] = {}

    with open(file_name, 'w') as configfile:
        config.write(configfile)


def generate_flightcon_file(route, case_name):
    file_name = route + '/' + case_name + '.flightcon.txt'
    config = configparser.ConfigParser()
    config['FlightCon'] = {'u_inf': 1.0,
                           'alpha': 25.0,
                           'beta': 0.0,
                           'rho_inf': 1.225,
                           'c_ref': 1.0,
                           'b_ref': 1.0}

    with open(file_name, 'w') as configfile:
        config.write(configfile)

if __name__ == '__main__':
    import os
    dir_path = os.path.dirname(os.path.realpath(__file__))
    generate_files(dir_path + '/', 'planar_wing', 10, 3)
    print('The test case has been successfully generated!!!')