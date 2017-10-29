"""
Ion object. Holds all methods properties of a CHIANTI ion.
"""

import numpy as np
from scipy.interpolate import splrep, splev, interp1d
import astropy.units as u
import astropy.constants as const
import periodictable

from .base import IonBase


class Ion(IonBase):
    """
    Ion class

    Parameters
    ----------
    ion_name: `str`
    temperature: ~astropy.Quantity

    The ion object is the fundamental unit of the fiasco library. An Ion object contains
    all of the properties and methods needed to access important information about each ion
    from the CHIANTI database.
    """
    
    @u.quantity_input
    def __init__(self, ion_name, temperature: u.K, *args, **kwargs):
        super().__init__(ion_name, *args, **kwargs)
        self.temperature = temperature
        # Use selected datasets
        if super().ip is not None:
            self._ip = super().ip[kwargs.get('ip_filename','chianti')]*const.h.cgs*const.c.cgs
        else:
            self._ip = None
    
    @property
    def ip(self):
        return self._ip
    
    @ip.setter
    def ip(self,value):
        self._ip = value
        
    @staticmethod
    def burgess_tully_descale(x, y, energy_ratio, c, scaling_type):
        """
        Convert scaled Burgess-Tully parameters to physical quantities. For more details see
        [1]_.

        References
        ----------
        .. [1] Burgess, A. and Tully, J. A., 1992, A&A, `254, 436 <http://adsabs.harvard.edu/abs/1992A%26A...254..436B>`_ 

        """
        nots = splrep(x, y, s=0)
        if scaling_type == 1:
            x_new = 1.0 - np.log(c)/np.log(energy_ratio + c)
            upsilon = splev(x_new, nots, der=0)*np.log(energy_ratio + np.e)
        elif scaling_type == 2:
            x_new = energy_ratio/(energy_ratio + c)
            upsilon = splev(x_new, nots, der=0)
        elif scaling_type == 3:
            x_new = energy_ratio/(energy_ratio + c)
            upsilon = splev(x_new, nots, der=0)/(energy_ratio + 1.0)
        elif scaling_type == 4:
            x_new = 1.0 - np.log(c)/np.log(energy_ratio + c)
            upsilon = splev(x_new, nots, der=0)*np.log(energy_ratio + c)
        elif scaling_type == 5:
            # dielectronic
            x_new = energy_ratio/(energy_ratio + c)
            upsilon = splev(x_new,nots,der=0)/energy_ratio
        elif scaling_type == 6:
            # protons
            x_new = energy_ratio/(energy_ratio + c)
            upsilon = 10**splev(x_new, nots, der=0)
        else:
            raise ValueError('Unrecognized BT92 scaling option.')

        return upsilon
    
    def direct_ionization_rate(self):
        """
        Calculate direct ionization rate in cm3/s
        
        Needs an equation reference or explanation
        """
        if self.diparams is None:
            return np.zeros(self.temperature.shape)*u.cm**3/u.s
        xgl,wgl = np.polynomial.laguerre.laggauss(12)
        kBT = const.k_B.cgs*self.temperature
        energy = np.outer(xgl,kBT)*kBT.unit + self.ip
        cross_section = self.direct_ionization_cross_section(energy)
        term1 = np.sqrt(8./np.pi/const.m_e.cgs)*np.sqrt(kBT)*np.exp(-self.ip/kBT)
        term2 = ((wgl*xgl)[:,np.newaxis]*cross_section).sum(axis=0)
        term3 = (wgl[:,np.newaxis]*cross_section).sum(axis=0)*self.ip/kBT
        
        return term1*(term2 + term3)
    
    @u.quantity_input
    def direct_ionization_cross_section(self, energy: u.erg):
        """
        Calculate direct ionization cross-section.

        The cross-sections are calculated one of two ways:

        - Using the method of [1]_ for hydrogenic and He-like ions
        - Using the scaled cross-sections of [2]_ for all other ions

        References
        ----------
        .. [1] Fontes, C. J., et al., 1999, Phys. Rev. A., `59 1329 <https://journals.aps.org/pra/abstract/10.1103/PhysRevA.59.1329>`_
        .. [2] Dere, K. P., 2007, A&A, `466, 771 <http://adsabs.harvard.edu/abs/2007A%26A...466..771D>`_
        """
        is_hydrogenic = (self.atomic_number - self.charge_state == 1) and (self.atomic_number >= 6)
        is_he_like = (self.atomic_number - self.charge_state == 2) and (self.atomic_number >= 10)
        
        if is_hydrogenic or is_he_like:
            # Fontes cross sections
            U = energy/self.ip
            A = 1.13
            B = 1 if is_hydrogenic else 2
            F = 1 if self.atomic_number < 20 else (140 + (self.atomic_number/20)**3.2)/141
            if self.atomic_number >= 16:
                c = -0.28394
                d = 1.95270
                C = 0.20594
                if self.atomic_number > 20:
                    C += ((self.atomic_number - 20)/50.5)**1.11
                D = 3.70590
            else:
                c = -0.80414
                d = 2.32431
                C = 0.14424
                D = 3.82652
            Qrp = 1./U*(A*np.log(U) + D*(1. - 1./U)**2 + C*U*(1. - 1./U)**4 + (c/U + d/U**2)*(1. - 1./U))
            
            return B*(np.pi*const.a0.cgs**2)*F*Qrp/(self.ip.to(u.Ry).value**2)
                
        else:
            # Cross-sections from diparams file
            cross_section_total = np.zeros(energy.shape)
            if self.diparams is None:
                return cross_section_total*u.cm**2
            for ip,bt_c,bt_e,bt_cross_section in zip(self.diparams['ip'],self.diparams['bt_c'],self.diparams['bt_e'],
                                                     self.diparams['bt_cross_section']):
                U = energy/(ip.to(u.erg))
                scaled_energy = 1. - np.log(bt_c)/np.log(U - 1. + bt_c)
                f_interp = interp1d(bt_e.value,bt_cross_section.value,kind='cubic',fill_value='extrapolate')
                scaled_cross_section = f_interp(scaled_energy.value)*bt_cross_section.unit
                # Only nonzero at energies above the ionization potential
                scaled_cross_section *= (U.value > 1.)
                cross_section = scaled_cross_section*(np.log(U) + 1.)/U/(ip**2)
                if not hasattr(cross_section_total,'unit'):
                    cross_section_total = cross_section_total*cross_section.unit
                cross_section_total += cross_section
                
            return cross_section_total
            
    def excitation_autoionization_rate(self):
        """
        Calculate ionization rate due to excitation autoionization
        """
        if self.easplups is None:
            return np.zeros(self.temperature.shape)*u.cm**3/u.s
        # Collision constant
        c = (const.h.cgs**2)/((2.*np.pi*const.m_e.cgs)**(1.5)*np.sqrt(const.k_B.cgs))
        kBTE = u.Quantity([(const.k_B.cgs*self.temperature)/(de.to(u.erg)) for de in self.easplups['delta_energy']])
        # Descale upsilon
        shape = self.easplups['bt_upsilon'].shape
        xs = np.tile(np.linspace(0,1,shape[1]),shape[0]).reshape(shape)
        args = [xs,self.easplups['bt_upsilon'].value,kBTE.value,self.easplups['bt_c'].value,self.easplups['bt_type']]
        upsilon = u.Quantity(list(map(self.burgess_tully_descale,*args)))
        # Rate coefficient
        rate = c*upsilon*np.exp(-1/kBTE)/np.sqrt(self.temperature[np.newaxis,:])
        
        return rate.sum(axis=0)
    
    def ionization_rate(self):
        """
        Total ionization rate.

        Includes contributions from both direct ionization and excitation-autoionization

        See Also
        --------
        direct_ionization_rate
        excitation_autoionization_rate
        """
        return self.direct_ionization_rate() + self.excitation_autoionization_rate()
    
    def radiative_recombination_rate(self):
        """
        Radiative recombination rate

        The recombination rate due to interaction with the ambient radiation field
        is calculated using a set of fit parameters using one of two methods:

        - Method of [1]_, (show expression)
        - Method of [2]_, (show expression)

        References
        ----------
        .. [1] Badnell, N. R., 2006, APJS, `167 334 <http://adsabs.harvard.edu/abs/2006ApJS..167..334B>`_
        .. [2] Shull, J. M., Van Steenberg, M., 1982, `48 95 <http://adsabs.harvard.edu/abs/1982ApJS...48...95S>`_
        """
        if self.rrparams is None:
            return np.zeros(self.temperature.shape)*u.cm**3/u.s
        if self.rrparams['fit_type'][0] == 1 or self.rrparams['fit_type'][0] == 2:
            A = self.rrparams['A_fit']
            B = self.rrparams['B_fit']
            if self.rrparams['fit_type'] == 2:
                B = B + self.rrparams['C_fit']*np.exp(-self.rrparams['T2_fit']/self.temperature)
            T0 = self.rrparams['T0_fit']
            T1 = self.rrparams['T1_fit']
            
            return A/(np.sqrt(self.temperature/T0)*(1 + np.sqrt(self.temperature/T0))**(1. - B)
                      *(1. + np.sqrt(self.temperature/T1))**(1. + B))
        elif self.rrparams['fit_type'][0] == 3:
            return self.rrparams['A_fit']*(self.temperature/(1e4*u.K))**(-self.rrparams['eta_fit'])
        else:
            raise ValueError('Unrecognized fit type {}'.format(self.rrparams['fit_type']))
    
    def dielectronic_recombination_rate(self):
        if self.drparams is None:
            return np.zeros(self.temperature.shape)*u.cm**3/u.s
        if self.drparams['fit_type'][0] == 1:
            E_over_T = np.outer(self.drparams['E_fit'],1./self.temperature)*(self.drparams['E_fit'].unit/self.temperature.unit)
            return self.temperature**(-1.5)*(self.drparams['c_fit'][:,np.newaxis]*np.exp(-E_over_T)).sum(axis=0)
        elif self.drparams['fit_type'][0] == 2:
            A = self.drparams['A_fit']
            B = self.drparams['B_fit']
            T0 = self.drparams['T0_fit']
            T1 = self.drparams['T1_fit']
            return A*self.temperature**(-1.5)*np.exp(-T0/self.temperature)*(1. + B*np.exp(-T1/self.temperature))
        else:
            raise ValueError('Unrecognized fit type {}'.format(self.drparams['fit_type']))
    
    def recombination_rate(self):
        """
        Total recombination rate.

        Includes contributions from dielectronic recombination and radiative recombination.

        See Also
        --------
        radiative_recombination_rate
        dielectronic_recombination_rate
        """
        return self.radiative_recombination_rate() + self.dielectronic_recombination_rate()
