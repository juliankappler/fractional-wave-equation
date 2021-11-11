#!/usr/bin/env python

import numpy as np
from scipy.optimize import root
from scipy.special import gamma
import math
import time


class integrator_base:
	'''
	This is the base class that contains the shared codebase for both the
	linear and nonlinear integrators
	'''
	def __init__(self,parameters):
		'''
		Set simulation parameters, which are provided
		as dictionary in the argument parameters
		'''
		#
		#######################################
		# Spatial and temporal discretization #
		#######################################
		self.L = parameters['L']
		self.T = parameters['T']
		#
		self.Nx = parameters['Nx']
		self.Nt = parameters['Nt']
		#
		self.dx = self.L/(self.Nx + 1)
		self.dt = self.T/self.Nt
		#
		self.x_array = np.arange(self.Nx+2,dtype=float)*self.dx
		self.t_array = np.arange(self.Nt+1,dtype=float)*self.dt
		#
		#########################
		# Fractional derivative #
		#########################
		self.alpha = parameters['alpha']
		#
		#####################
		# Initial condition #
		#####################
		self.current_state = np.zeros(self.Nx+2,dtype=float)
		self.initial_velocity = np.zeros(self.Nx+2,dtype=float)

		# some constants
		self.Lambda = gamma(3-self.alpha)*self.dt**self.alpha
		self.omega0 = self.omega(0)
		#

		try:
			self.verbose = parameters['verbose']
		except KeyError:
			self.verbose = True

	def set_parameters(self,parameters):
		'''
		Change parameters of an existing instance of this class
		'''
		try:
			self.L = parameters['L']
		except KeyError:
			pass
		#
		try:
			self.T = parameters['T']
		except KeyError:
			pass
		#
		try:
			self.Nx = parameters['Nx']
		except KeyError:
			pass
		#
		try:
			self.Nt = parameters['Nt']
		except KeyError:
			pass
		#
		self.dx = self.L/(self.Nx + 1)
		self.dt = self.T/self.Nt
		#
		self.x_array = np.arange(self.Nx+2,dtype=float)*self.dx
		self.t_array = np.arange(self.Nt+1,dtype=float)*self.dt
		#
		try:
			self.alpha = parameters['alpha']
		except KeyError:
			pass
		#
		self.Lambda = gamma(3-self.alpha)*self.dt**self.alpha
		self.omega0 = self.omega(0)
		#
		try:
			self.verbose = parameters['verbose']
		except KeyError:
			self.verbose = True

	def get_parameters(self,print_parameters=False):
		'''
		return parameters of an existing instance of this class
		'''
		output_dictionary = {'L':self.L,
             'Nx':self.Nx,
             'T':self.T,
             'Nt':self.Nt,
             'alpha':self.alpha,
			 'verbose':self.verbose}
		if print_parameters:
			print("Parameters set for this instance:")
			print("L       = {0}".format(self.L))
			print("Nx      = {0}".format(self.Nx))
			print("T       = {0}".format(self.T))
			print("Nt      = {0}".format(self.Nt))
			print("alpha   = {0}".format(self.alpha))
			print("verbose = {0}".format(self.verbose))
		return output_dictionary

	def print_time_remaining(self,step,end='\r'):
		elapsed_time = time.time() - self.system_time_at_start_of_simulation
		m_elapsed, s_elapsed = divmod(elapsed_time, 60)
		h_elapsed, m_elapsed = divmod(m_elapsed, 60)
		remaining_time = (self.Nt/step -1)*elapsed_time
		m_remaining, s_remaining = divmod(remaining_time, 60)
		h_remaining, m_remaining = divmod(m_remaining, 60)
		print("Running simulation. Progress: {0}%, elapsed time: {1:d}:{2:02d}:{3:02d}, remaining time: {4:d}:{5:02d}:{6:02d}\t\t\t".format(int(step/self.Nt*100.),
																														int(np.round(h_elapsed)),int(np.round(m_elapsed)),int(np.round(s_elapsed)),
																														int(np.round(h_remaining)),int(np.round(m_remaining)),int(np.round(s_remaining))),
			  end=end)




	def omega(self,j):
		return (j+1)**(2-self.alpha)-j**(2-self.alpha)

	def fractional_sum_over_past(self):
		n = len(self.state_history)-2 # because state_history has a -dt and
							# a zero time entry
		return_sum = np.zeros(self.Nx+2,dtype=float)
		# add -1 term (note that b_{-1}^n = -omega(n))
		return_sum += -self.omega(n)*self.state_history[0]
		# add 0 term (note that b_0^n = 2*omega(n)-omega(n-1))
		if 1 <= n:
			return_sum += (2*self.omega(n)-self.omega(n-1)) \
								* self.state_history[1]
		else:
			return_sum += 2*self.omega(n)*self.state_history[1]
		# add terms for 1 <= j <= n-1
		for i,e in enumerate(self.state_history[:-1]):
			if 2 <= i: # the first two we already summed
				j = i-1 # i = 2 is the j=1 term.
				return_sum += (-self.omega(n-j+1) \
							    + 2*self.omega(n-j) \
								- self.omega(n-j-1) ) * e
		# add last term (note that omega(0) = 1)
		if 1 <= n: # not needed in the first timestep
			return_sum += (-self.omega(1)+2)*self.state_history[-1]
		return return_sum

	'''
	def fractional_sum_over_past(self,step):
		# at integration step = n, the array self.state_history
		# contains n + 1 nonzero entries:
		# self.state_history[0] = U at time t = -dt
		# self.state_history[1] = U at time t = 0
		# . . .
		# self.state_history[n] = U at time t = (n-1)*dt

		n = len(self.state_history)-2 # because state_history has a -dt and
							# a zero time entry
		return_sum = np.zeros(self.Nx+2,dtype=float)
		# add -1 term (note that b_{-1}^n = -omega(n))
		return_sum += -self.omega(n)*self.state_history[0]
		# add 0 term (note that b_0^n = 2*omega(n)-omega(n-1))
		if 1 <= n:
			return_sum += (2*self.omega(n)-self.omega(n-1)) \
								* self.state_history[1]
		else:
			return_sum += 2*self.omega(n)*self.state_history[1]
		# add terms for 1 <= j <= n-1
		for i,e in enumerate(self.state_history[:-1]):
			if 2 <= i: # the first two we already summed
				j = i-1 # i = 2 is the j=1 term.
				return_sum += (-self.omega(n-j+1) \
							    + 2*self.omega(n-j) \
								- self.omega(n-j-1) ) * e
		# add last term (note that omega(0) = 1)
		if 1 <= n: # not needed in the first timestep
			return_sum += (-self.omega(1)+2)*self.state_history[-1]
		return return_sum
	''';


	def return_results(self):
		output_dictionary = {'t':self.t_array,
							'x':self.x_array,
							'y':np.array(self.state_history[1:]),
							'L':self.L,'T':self.T,
							'dx':self.dx,'Nx':self.Nx,
							'dt':self.dt,'Nt':self.Nt,
							'alpha':self.alpha,
							}
		return output_dictionary



class nonlinear(integrator_base):
	#

	def root_of_this_function_defines_next_configuration(self,
									current_state_candidate,
									fractional_sum_over_past_without_bc,
									D,current_U0,current_UL):
		#
		current_state_candidate_with_boundary_conditions = np.zeros(self.Nx+2,dtype=float)
		current_state_candidate_with_boundary_conditions[0] = current_U0
		current_state_candidate_with_boundary_conditions[-1] = current_UL
		current_state_candidate_with_boundary_conditions[1:-1] = current_state_candidate
		#
		Psi = current_state_candidate_with_boundary_conditions
		dPsi = Psi[2:]-Psi[:-2]   # dPsi[i] = Psi[i+1]-Psi[i], i.e. dPsi[0] = Psi[x_1] - LeftBoundaryCondition
		# use the array with the derivatives to calculate
		# D(\partial_x u) for each of the N+1 intervals
		dArray = D(dPsi/(2.*self.dx))/self.dx**2
		return_array = self.omega0*current_state_candidate \
					+ self.Lambda*dArray*(2*current_state_candidate \
						-current_state_candidate_with_boundary_conditions[:-2] \
						-current_state_candidate_with_boundary_conditions[2:]) \
					- fractional_sum_over_past_without_bc
		return return_array



	def simulate(self,D,U0,UL):
		# Since the fractional derivative is nonlocal in time, we need to
		# store all past states. These are saved in the following list:
		self.state_history = []
		self.state_history.append(self.current_state.copy() \
									- self.dt*self.initial_velocity)
		self.state_history.append(self.current_state.copy())
		# Note that state_history[0] is at time -dt,
		# while state_history[1] is the initial condition

		self.system_time_at_start_of_simulation = time.time()

		for step,next_time in enumerate(self.t_array[1:]):
			if self.verbose:
				if ( step % int(self.Nt/100) == 0 ) and step > 0:
					self.print_time_remaining(step)
			self.integrate(step,next_time,D,U0,UL)
		#
		if self.verbose:
			self.print_time_remaining(step+1,end='\n')
		return self.return_results()

	def integrate(self,step,next_time,D,U0,UL):
		# Note that our solver takes the *next* configuration to
		# calculate the new field configuration. That means at step
		# 0, we already use the boundary condition at step 1. Thats
		# why in integrate(), the boundary condition at step+1 is used.
		previous_state = self.state_history[-1].copy()
		fractional_sum_over_past_without_bc = self.fractional_sum_over_past()[1:-1].copy()
		current_root = np.zeros(self.Nx+2,dtype=float)
		current_root[0] = U0( next_time )
		current_root[-1] = UL( next_time )
		current_root[1:-1] = root(fun=self.root_of_this_function_defines_next_configuration,
								x0 = previous_state[1:-1],
								args=(fractional_sum_over_past_without_bc,
										D,
										current_root[0],
										current_root[-1])).x
		self.state_history.append(current_root)
		return 0

'''
class linear(integrator):
	#
	## THIS IS FOR LINEAR MATRIX INVERSION MODE
	#W = CreateW(CurConf)
	#Prop = sp.linalg.inv(K+dt*W)*K
	#Prop = sp.linalg.inv(K+Lambda*W)*K

	#	CurConf = Prop.dot(LastConf)
	#	History.append(Prop.dot(SumOverPast()))



	def root_of_this_function_defines_next_configuration(self,
									current_state_candidate,
									fractional_sum_over_past_without_bc):
		#
		current_state_candidate_with_boundary_conditions = np.zeros(N+2,dtype=float)
		current_state_candidate_with_boundary_conditions[0] = boundary_condition_function(step+1)
		current_state_candidate_with_boundary_conditions[1:-1] = current_state_candidate
		#
		Psi = current_state_candidate_with_boundary_conditions
		dPsi = Psi[2:]-Psi[:-2]   # dPsi[i] = Psi[i+1]-Psi[i], i.e. dPsi[0] = Psi[x_1] - LeftBoundaryCondition
		# use the array with the derivatives to calculate
		# D(\partial_x u) for each of the N+1 intervals
		dArray = vD(dPsi/(2.*dx))/dx**2
		return_array = omega0*current_state_candidate \
					+ Lambda*dArray*(2*current_state_candidate \
						-current_state_candidate_with_boundary_conditions[:-2] \
						-current_state_candidate_with_boundary_conditions[2:]) \
					- fractional_sum_over_past_without_bc
		return return_array

	def integrate(self,step):
		# Note that our solver takes the *next* configuration to
		# calculate the new field configuration. That means at step
		# 0, we already use the boundary condition at step 1. Thats
		# why in integrate(), the boundary condition at step+1 is used.
		previous_state = self.state_history[-1].copy()
		fractional_sum_over_past_without_bc = self.fractional_sum_over_past()[1:-1].copy()
		current_root = np.zeros(self.Nx+2,dtype=float)
		current_root[0] = boundary_condition_function( (step+1)*self.dt )
		current_root[1:-1] = root(fun=self.root_of_this_function_defines_next_configuration,
								x0 = previous_state[1:-1],
								args=(fractional_sum_over_past_without_bc)).x
		self.state_history.append(current_root)
		return 0

''';
