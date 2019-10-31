import numpy as np

from henrys_law import O2_SOL, CO2_SOL



class RuBisCOParams(object):
    def __init__(self, kcat_C, K_C, K_O, S, name=None, mw=65e3):
        """
        Args:
            kcat_C: units uM/s.
            K_C: units uM.
            K_O: units uM.
            S: unitless.
            name: a name for this Rubisco.
            mw: mass per active site in g/mol.
            
        Specificity factor S is defined as
            S = kcat_C * K_O / (kcat_O * K_C)
        """
        self.kcat_C = kcat_C
        self.K_C = K_C
        self.K_O = K_O
        self.S = S
        self.kcat_O = kcat_C * K_O / (K_C * S)
        self.name = name
        self.mw = mw

    def vs(self, co2, o2):
        """Concentrations in micromolar. Returns per-active site rates."""
        vc = ((self.kcat_C) /
              (1.0 + (self.K_C / co2) + (self.K_C * o2 / (self.K_O * co2))))
        vo = ((self.kcat_O) /
              (1.0 + (self.K_O / o2) + (self.K_O * co2 / (self.K_C * o2))))
        return vc, vo

    def net_vc(self, co2, o2):
        """Returns the net carboxylation rate assuming c2 photorespiration.

        Args:
            co2: uM co2 concentration.
            o2: uM o2 concentration.

        Returns: net carboxylation rate of this rubisco.
        """
        vc, vo = self.vs(co2, o2)
        return vc - vo/2
    
    def __str__(self):
        return '<Rubisco kcat_C=%.2f kcat_O=%.2f K_C=%.2f K_O=%.2f S=%.1f>' % (
            self.kcat_C, self.kcat_O, self.K_C, self.K_O, self.S)


# Occhilanni 2015
PCC7942 = RuBisCOParams(kcat_C=14.4, K_C=172, K_O=585, S=43, name='PCC7942')
# Savir 2010
RUBRUM = RuBisCOParams(kcat_C=7.3, K_C=80, K_O=406, S=12.3, name='Rubrum', mw=52e3)
SPINACH = RuBisCOParams(kcat_C=3.7, K_C=14, K_O=480, S=80, name='Spinach')
MAIZE = RuBisCOParams(kcat_C=4.4, K_C=34, K_O=810, S=78, name='Maize')
TOBACCO = RuBisCOParams(kcat_C=3.4, K_C=10.7, K_O=295, S=82, name='Tobacco')
_R = 8.314e-3 # Universal gas constant kJ/mol/K


class ArrheniusParam(object):
    """A parameter with Arrhenius T-dependence. T in Kelvin"""
    def __init__(self, ref_value, ref_temp, e_act):
        self.ref_value = ref_value
        self.ref_temp = ref_temp
        self.e_act = e_act

    def value(self, temp=None):
        my_temp = temp or self.ref_temp
        exp_term = np.exp((-self.e_act/(_R*my_temp)) * (self.ref_temp - my_temp) / self.ref_temp)
        return self.ref_value*exp_term

    def __str__(self):
        return '<ArrheniusParam ref_val=%.2f ref_temp=%.2f e_act=%.2f>' % (self.ref_value, self.ref_temp, self.e_act)


class ArrheniusRubisco(object):
    def __init__(self, kcat_C, K_C, K_O, S, ref_temp=None, name=None, mw=65e3):
        """
        Args:
            kcat_C: units uM/s. ArrheniusParam
            K_C: units uM. ArrheniusParam
            K_O: units uM. ArrheniusParam
            S: unitless. ArrheniusParam
            ref_temp: temperature in Kelvin. default is 298.15 K. 
            name: a name for this Rubisco.
            mw: mass per active site in g/mol.
            
        Specificity factor S is defined as
            S = kcat_C * K_O / (kcat_O * K_C)
        """
        self.kcat_C = kcat_C
        self.K_C = K_C
        self.K_O = K_O
        self.S = S
        self.ref_temp = ref_temp or 298.15 
        self.name = name
        self.mw = mw

    def concrete_temp(self, temp):
        kcat_C = self.kcat_C.value(temp)
        K_C = self.K_C.value(temp)
        K_O = self.K_O.value(temp)
        S = self.S.value(temp)

        # Need to convert K_C and K_O into umolar from Pa. 
        K_C_umolar = CO2_SOL.molar_conc(K_C, temp) * 1e6
        K_O_umolar = O2_SOL.molar_conc(K_O, temp) * 1e6

        # Need to convert S to soluble - S_C/O has K_O in numerate and K_C in denominator
        sol_ratio =  O2_SOL.h_cp(temp) / CO2_SOL.h_cp(temp)
        S_sol = S * sol_ratio

        return RuBisCOParams(kcat_C, K_C_umolar, K_O_umolar, S_sol, name=self.name, mw=self.mw)

    def vs(self, co2, o2, temp=None):
        """Concentrations in micromolar. Returns per-active site rates."""
        temp = temp or self.ref_temp
        rps = self.concrete_temp(temp)
        return rps.vs(co2, o2)

    def net_vc(self, co2, o2, temp=None):
        """Returns the net carboxylation rate assuming c2 photorespiration.

        Args:
            co2: uM co2 concentration.
            o2: uM o2 concentration.

        Returns: net carboxylation rate of this rubisco.
        """
        temp = temp or self.ref_temp
        rps = self.concrete_temp(temp)
        return rps.net_vc(co2, o2)

    def __str__(self):
        return '<ArrheniusRubisco kcat_C=%s K_C=%s K_O=%s S=%s>' % (
            self.kcat_C, self.K_C, self.K_O, self.S)


_REF_T = 298.15
ATHALIANA_T = ArrheniusRubisco(
    kcat_C=ArrheniusParam(3.1, _REF_T, 59.6),
    K_C=ArrheniusParam(36, _REF_T, 63.0),
    K_O=ArrheniusParam(23100, _REF_T, 16.9),
    S=ArrheniusParam(2003, _REF_T, -28.7),
    ref_temp=_REF_T,
    name='A. thaliana')