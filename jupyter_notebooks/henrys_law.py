import numpy as np 


class HenrysLaw(object):

	def __init__(self, h_cp, t_dep, ref_temp, name=None):
		self.name = name
		self._ref_h_cp = h_cp
		self._t_dep = t_dep
		self.ref_temp = ref_temp

	def h_cp(self, temp=None):
		"""Calculate H_CP at defined temperature.

		Should have units of mol/m^3/Pa if following Sander paper. 
		"""
		_t = temp or self.ref_temp
		temp_term = (1.0/_t) - (1.0/self.ref_temp)
		exp_term = np.exp(self._t_dep * temp_term)
		return self._ref_h_cp * exp_term

	def molar_conc(self, p, temp=None):
		"""Calculate molar concentration from partial pressure p and temp."""
		_h_cp = self.h_cp(temp)
		conc = p*_h_cp  # mol/m^3

		# 1 L = (dm)^3 = (0.1 m)^3 = 1e-3 m^3
		return 1e-3 * conc



_REF_T = 298.15
O2_SOL = HenrysLaw(1.2e-5, 1700, _REF_T, name='O2')
CO2_SOL = HenrysLaw(3.3e-4, 2400, _REF_T, name='CO2')
