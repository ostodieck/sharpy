"""@package PyBeam.Solver.NonlinearStatic
@brief      Nonlinear static solvers.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       25/10/2012
@pre        None
@warning    None

@modified   Alfonso del Carre
"""

# import DerivedTypes
# import BeamIO
# import beam.utils.derivedtypes as derivedtypes
# import XbeamLib
# import BeamInit
import ctypes as ct

import numpy as np

import sharpy.beam.utils.beamlib as beamlib
from sharpy.presharpy.utils.settings import str2bool
from sharpy.utils.solver_interface import solver, BaseSolver


@solver
class NonLinearStatic(BaseSolver):
    solver_id = 'NonLinearStatic'
    solver_type = 'structural'

    def __init__(self):
        pass

    def initialise(self, data):
        self.data = data
        self.settings = data.settings[self.solver_id]
        self.convert_settings()
        data.beam.generate_aux_information()

    def run(self):
        print('Running non linear static solver...')
        beamlib.cbeam3_solv_nlnstatic(self.data.beam, self.settings)
        self.data.beam.update()
        print('...Finished')
        return self.data

    def convert_settings(self):
        self.settings['print_info'] = ct.c_bool(str2bool(self.settings['print_info']))
        self.settings['out_b_frame'] = ct.c_bool(str2bool(self.settings['out_b_frame']))
        self.settings['out_a_frame'] = ct.c_bool(str2bool(self.settings['out_a_frame']))
        self.settings['elem_proj'] = ct.c_int(int(self.settings['elem_proj']))
        self.settings['max_iterations'] = ct.c_int(int(self.settings['max_iterations']))
        self.settings['num_load_steps'] = ct.c_int(int(self.settings['num_load_steps']))
        self.settings['delta_curved'] = ct.c_double(float(self.settings['delta_curved']))
        self.settings['min_delta'] = ct.c_double(float(self.settings['min_delta']))
        self.settings['newmark_damp'] = ct.c_double(float(self.settings['newmark_damp']))
        self.settings['gravity_on'] = ct.c_bool(str2bool(self.settings['gravity_on']))
        self.settings['gravity'] = ct.c_double(float(self.settings['gravity']))
        self.settings['gravity_dir'] = np.fromstring(self.settings['gravity_dir'], sep=',', dtype=ct.c_double)
