import ctypes as ct
import numpy as np

from sharpy.utils.solver_interface import solver, BaseSolver, solver_from_string
import sharpy.utils.settings as settings

import sharpy.structure.utils.xbeamlib as xbeamlib
import sharpy.utils.multibody as mb
import sharpy.structure.utils.lagrangeconstraints as lagrangeconstraints
import sharpy.utils.exceptions as exc


_BaseStructural = solver_from_string('_BaseStructural')


@solver
class NonLinearDynamicMultibody(_BaseStructural):
    """
    Nonlinear dynamic multibody

    Nonlinear dynamic step solver for multibody structures.

    """
    solver_id = 'NonLinearDynamicMultibody'
    solver_classification = 'structural'

    settings_types = _BaseStructural.settings_types.copy()
    settings_default = _BaseStructural.settings_default.copy()
    settings_description = _BaseStructural.settings_description.copy()

    settings_table = settings.SettingsTable()
    __doc__ += settings_table.generate(settings_types, settings_default, settings_description)

    def __init__(self):
        self.data = None
        self.settings = None

        # Total number of unknowns in the Multybody sistem
        self.sys_size = None

        # Total number of equations associated to the Lagrange multipliers
        self.lc_list = None
        self.num_LM_eq = None

        self.gamma = None
        self.beta = None

    def initialise(self, data, custom_settings=None):

        self.data = data
        if custom_settings is None:
            self.settings = data.settings[self.solver_id]
        else:
            self.settings = custom_settings
        settings.to_custom_types(self.settings, self.settings_types, self.settings_default)

        # load info from dyn dictionary
        self.data.structure.add_unsteady_information(
            self.data.structure.dyn_dict, self.settings['num_steps'].value)

        # Define Newmark constants
        self.gamma = 0.5 + self.settings['newmark_damp'].value
        self.beta = 0.25*(self.gamma + 0.5)*(self.gamma + 0.5)

        # Define the number of equations
        self.lc_list = lagrangeconstraints.initialize_constraints(self.data.structure.ini_mb_dict)
        self.num_LM_eq = lagrangeconstraints.define_num_LM_eq(self.lc_list)

        # Define the number of dofs
        self.define_sys_size()

    def add_step(self):
        self.data.structure.next_step()

    def next_step(self):
        pass

    def define_sys_size(self):

        MBdict = self.data.structure.ini_mb_dict
        self.sys_size = self.data.structure.num_dof.value

        for ibody in range(self.data.structure.num_bodies):
            if (MBdict['body_%02d' % ibody]['FoR_movement'] == 'free'):
                self.sys_size += 10

    def assembly_MB_eq_system(self, MB_beam, MB_tstep, ts, dt, Lambda, Lambda_dot, MBdict):
        self.lc_list = lagrangeconstraints.initialize_constraints(MBdict)
        self.num_LM_eq = lagrangeconstraints.define_num_LM_eq(self.lc_list)

        MB_M = np.zeros((self.sys_size+self.num_LM_eq, self.sys_size+self.num_LM_eq), dtype=ct.c_double, order='F')
        MB_C = np.zeros((self.sys_size+self.num_LM_eq, self.sys_size+self.num_LM_eq), dtype=ct.c_double, order='F')
        MB_K = np.zeros((self.sys_size+self.num_LM_eq, self.sys_size+self.num_LM_eq), dtype=ct.c_double, order='F')
        MB_Asys = np.zeros((self.sys_size+self.num_LM_eq, self.sys_size+self.num_LM_eq), dtype=ct.c_double, order='F')
        MB_Q = np.zeros((self.sys_size+self.num_LM_eq,), dtype=ct.c_double, order='F')
        #ipdb.set_trace()
        first_dof = 0
        last_dof = 0
        # Loop through the different bodies
        for ibody in range(len(MB_beam)):

            # Initialize matrices
            M = None
            C = None
            K = None
            Q = None

            # Generate the matrices for each body
            if MB_beam[ibody].FoR_movement == 'prescribed':
                last_dof = first_dof + MB_beam[ibody].num_dof.value
                M, C, K, Q = xbeamlib.cbeam3_asbly_dynamic(MB_beam[ibody], MB_tstep[ibody], self.settings)

            elif MB_beam[ibody].FoR_movement == 'free':
                last_dof = first_dof + MB_beam[ibody].num_dof.value + 10
                M, C, K, Q = xbeamlib.xbeam3_asbly_dynamic(MB_beam[ibody], MB_tstep[ibody], self.settings)


            ############### Assembly into the global matrices
            # Flexible and RBM contribution to Asys
            MB_M[first_dof:last_dof, first_dof:last_dof] = M.astype(dtype=ct.c_double, copy=True, order='F')
            MB_C[first_dof:last_dof, first_dof:last_dof] = C.astype(dtype=ct.c_double, copy=True, order='F')
            MB_K[first_dof:last_dof, first_dof:last_dof] = K.astype(dtype=ct.c_double, copy=True, order='F')

            #Q
            MB_Q[first_dof:last_dof] = Q

            first_dof = last_dof

        # Define the number of equations
        # Generate matrices associated to Lagrange multipliers
        LM_C, LM_K, LM_Q = lagrangeconstraints.generate_lagrange_matrix(
            self.lc_list,
            MB_beam,
            MB_tstep,
            ts,
            self.num_LM_eq,
            self.sys_size,
            dt,
            Lambda,
            Lambda_dot,
            "dynamic")

        # Include the matrices associated to Lagrange Multipliers
        MB_C += LM_C
        MB_K += LM_K
        MB_Q += LM_Q

        MB_Asys = MB_K + MB_C*self.gamma/(self.beta*dt) + MB_M/(self.beta*dt*dt)

        return MB_Asys, MB_Q

    def integrate_position(self, MB_beam, MB_tstep, dt):
        vel = np.zeros((6,),)
        acc = np.zeros((6,),)
        for ibody in range(0, len(MB_tstep)):
            # I think this is the right way to do it, but to make it match the rest I change it temporally
            if True:
                # MB_tstep[ibody].mb_quat[ibody,:] =  algebra.quaternion_product(MB_tstep[ibody].quat, MB_tstep[ibody].mb_quat[ibody,:])
                acc[0:3] = (0.5-self.beta)*np.dot(MB_beam[ibody].timestep_info.cga(),MB_beam[ibody].timestep_info.for_acc[0:3])+self.beta*np.dot(MB_tstep[ibody].cga(),MB_tstep[ibody].for_acc[0:3])
                vel[0:3] = np.dot(MB_beam[ibody].timestep_info.cga(),MB_beam[ibody].timestep_info.for_vel[0:3])
                MB_tstep[ibody].for_pos[0:3] += dt*(vel[0:3] + dt*acc[0:3])
            else:
                MB_tstep[ibody].for_pos[0:3] += dt*np.dot(MB_tstep[ibody].cga(),MB_tstep[ibody].for_vel[0:3])

        # Use next line for double pendulum (fix position of the second FoR)
        # MB_tstep[ibody].for_pos[0:3] = np.dot(algebra.quat2rotation(MB_tstep[0].quat), MB_tstep[0].pos[-1,:])
        # print("tip final pos: ", np.dot(algebra.quat2rotation(MB_tstep[0].quat), MB_tstep[0].pos[-1,:]))
        # print("FoR final pos: ", MB_tstep[ibody].for_pos[0:3])
        # print("pause")

    def extract_resultants(self):
        # TODO: code
        pass

    def compute_forces_constraints(self, MB_beam, MB_tstep, ts, dt, Lambda, Lambda_dot):
        try:
            self.lc_list[0]
        except IndexError:
            return

        # TODO the output of this routine is wrong. check at some point.
        LM_C, LM_K, LM_Q = lagrangeconstraints.generate_lagrange_matrix(self.lc_list, MB_beam, MB_tstep, ts, self.num_LM_eq, self.sys_size, dt, Lambda, Lambda_dot, "dynamic")
        F = -np.dot(LM_C[:, -self.num_LM_eq:], Lambda_dot) - np.dot(LM_K[:, -self.num_LM_eq:], Lambda)

        first_dof = 0
        for ibody in range(len(MB_beam)):
            # Forces associated to nodes
            body_numdof = MB_beam[ibody].num_dof.value
            body_freenodes = np.sum(MB_beam[ibody].vdof > -1)
            last_dof = first_dof + body_numdof
            MB_tstep[ibody].forces_constraints_nodes[(MB_beam[ibody].vdof > -1), :] = F[first_dof:last_dof].reshape(body_freenodes, 6, order='C')

            # Forces associated to the frame of reference
            if MB_beam[ibody].FoR_movement == 'free':
                # TODO: How are the forces in the quaternion equation interpreted?
                MB_tstep[ibody].forces_constraints_FoR[ibody, :] = F[last_dof:last_dof+10]
                last_dof += 10

            first_dof = last_dof
            # print(MB_tstep[ibody].forces_constraints_nodes)
        # TODO: right now, these forces are only used as an output, they are not read when the multibody is splitted

    def run(self, structural_step=None, dt=None):
        if structural_step is None:
            structural_step = self.data.structure.timestep_info[-1]

        if structural_step.mb_dict is not None:
            MBdict = structural_step.mb_dict
        else:
            MBdict = self.data.structure.ini_mb_dict

        if dt is None:
            dt = self.settings['dt'].value
        else:
            self.settings['dt'] = ct.c_float(dt)

        self.lc_list = lagrangeconstraints.initialize_constraints(MBdict)
        self.num_LM_eq = lagrangeconstraints.define_num_LM_eq(self.lc_list)

        # TODO: only working for constant forces
        MB_beam, MB_tstep = mb.split_multibody(
            self.data.structure,
            structural_step,
            MBdict,
            self.data.ts)

        # Lagrange multipliers parameters
        num_LM_eq = self.num_LM_eq
        Lambda = np.zeros((num_LM_eq,), dtype=ct.c_double, order='F')
        Lambda_dot = np.zeros((num_LM_eq,), dtype=ct.c_double, order='F')

        # Initialize
        q = np.zeros((self.sys_size + num_LM_eq,), dtype=ct.c_double, order='F')
        dqdt = np.zeros((self.sys_size + num_LM_eq,), dtype=ct.c_double, order='F')
        dqddt = np.zeros((self.sys_size + num_LM_eq,), dtype=ct.c_double, order='F')

        # Predictor step
        mb.disp2state(MB_beam, MB_tstep, q, dqdt, dqddt)

        q += dt*dqdt + (0.5 - self.beta)*dt*dt*dqddt
        dqdt += (1.0 - self.gamma)*dt*dqddt
        dqddt = np.zeros((self.sys_size + num_LM_eq,), dtype=ct.c_double, order='F')
        if not num_LM_eq == 0:
            Lambda = q[-num_LM_eq:].astype(dtype=ct.c_double, copy=True, order='F')
            Lambda_dot = dqdt[-num_LM_eq:].astype(dtype=ct.c_double, copy=True, order='F')
        else:
            Lambda = 0
            Lambda_dot = 0

        # Newmark-beta iterations
        old_Dq = 1.0
        LM_old_Dq = 1.0

        converged = False
        for iteration in range(self.settings['max_iterations'].value):
            # Check if the maximum of iterations has been reached
            if iteration == self.settings['max_iterations'].value - 1:
                error = ('Solver did not converge in %d iterations.\n res = %e \n LM_res = %e' % 
                        (iteration, res, LM_res))
                raise exc.NotConvergedSolver(error)

            # Update positions and velocities
            mb.state2disp(q, dqdt, dqddt, MB_beam, MB_tstep)
            MB_Asys, MB_Q = self.assembly_MB_eq_system(MB_beam,
                                                       MB_tstep,
                                                       self.data.ts,
                                                       dt,
                                                       Lambda,
                                                       Lambda_dot,
                                                       MBdict)

            # Compute the correction
            # ADC next line not necessary
            # Dq = np.zeros((self.sys_size+num_LM_eq,), dtype=ct.c_double, order='F')
            # MB_Asys_balanced, T = scipy.linalg.matrix_balance(MB_Asys)
            # invT = np.matrix(T).I
            # MB_Q_balanced = np.dot(invT, MB_Q).T

            Dq = np.linalg.solve(MB_Asys, -MB_Q)
            # least squares solver
            # Dq = np.linalg.lstsq(np.dot(MB_Asys_balanced, invT), -MB_Q_balanced, rcond=None)[0]

            # Evaluate convergence
            if iteration:
                res = np.max(np.abs(Dq[0:self.sys_size]))/old_Dq
                if np.isnan(res):
                    raise exc.NotConvergedSolver('Multibody res = NaN')
                if num_LM_eq:
                    LM_res = np.max(np.abs(Dq[self.sys_size:self.sys_size+num_LM_eq]))/LM_old_Dq
                else:
                    LM_res = 0.0
                if (res < self.settings['min_delta'].value) and (LM_res < self.settings['min_delta'].value):
                    converged = True

            # Compute variables from previous values and increments
            # TODO:decide If I want other way of updating lambda
            # this for least sq
            # q[:, np.newaxis] += Dq
            # dqdt[:, np.newaxis] += self.gamma/(self.beta*dt)*Dq
            # dqddt[:, np.newaxis] += 1.0/(self.beta*dt*dt)*Dq

            # this for direct solver
            q += Dq
            dqdt += self.gamma/(self.beta*dt)*Dq
            dqddt += 1.0/(self.beta*dt*dt)*Dq

            if not num_LM_eq == 0:
                Lambda = q[-num_LM_eq:].astype(dtype=ct.c_double, copy=True, order='F')
                Lambda_dot = dqdt[-num_LM_eq:].astype(dtype=ct.c_double, copy=True, order='F')
            else:
                Lambda = 0
                Lambda_dot = 0

            if converged:
                break

            if not iteration:
                old_Dq = np.max(np.abs(Dq[0:self.sys_size]))
                if old_Dq < 1.0:
                    old_Dq = 1.0
                if num_LM_eq:
                    LM_old_Dq = np.max(np.abs(Dq[self.sys_size:self.sys_size+num_LM_eq]))
                else:
                    LM_old_Dq = 1.0

        mb.state2disp(q, dqdt, dqddt, MB_beam, MB_tstep)
        # end: comment time stepping

        # End of Newmark-beta iterations
        self.integrate_position(MB_beam, MB_tstep, dt)
        # lagrangeconstraints.postprocess(self.lc_list, MB_beam, MB_tstep, MBdict, "dynamic")
        lagrangeconstraints.postprocess(self.lc_list, MB_beam, MB_tstep, "dynamic")
        self.compute_forces_constraints(MB_beam, MB_tstep, self.data.ts, dt, Lambda, Lambda_dot)
        if self.settings['gravity_on']:
            for ibody in range(len(MB_beam)):
                xbeamlib.cbeam3_correct_gravity_forces(MB_beam[ibody], MB_tstep[ibody], self.settings)
        mb.merge_multibody(MB_tstep, MB_beam, self.data.structure, structural_step, MBdict, dt)

        # structural_step.q[:] = q[:self.sys_size].copy()
        # structural_step.dqdt[:] = dqdt[:self.sys_size].copy()
        # structural_step.dqddt[:] = dqddt[:self.sys_size].copy()

        return self.data

