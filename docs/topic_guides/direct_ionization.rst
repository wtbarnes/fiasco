.. _fiasco-topic-guide-direct-ionization-rate:

The direct ionization rate and cross-section
============================================

The ionization rate due to the collisions with free electrons can be written as the integral of the
velocity-weighted collisional cross-section over the Maxwell-Boltzmann distribution.
Following Section 3.5.1 of :cite:t:`del_zanna_solar_2018`, this can be written as,

.. math::

    C^I = \sqrt{\frac{8}{\pi m_e}}(k_BT)^{-3/2}\int_I^{\infty}\mathrm{d}E\,E\sigma_I(E)\exp{\left(-\frac{E}{k_BT}\right)}

where :math:`E` is the energy of the incident electron, :math:`I` is the ionization energy of the initially bound electron,
and :math:`\sigma_I` is the ionization cross-section.

Making the substitution :math:`x=(E-I)/k_BT`, the above integral can be rewritten as,

.. math::

    \begin{aligned}
        C^I = \sqrt{\frac{8k_BT}{\pi m_e}}\exp{\left(-\frac{I}{k_BT}\right)}&\left(\int_0^{\infty}\mathrm{d}x\,x\sigma_{I}(k_BTx+I)e^{-x} \right. \\
                                                                            &\left. + \frac{I}{k_BT}\int_0^{\infty}\mathrm{d}x\,\sigma_{I}(k_BTx+I)e^{-x}\right).
    \end{aligned}

Each of these integrals is of the form such that they can be evaluated using Gauss-Laguerre quadrature,

.. math::

    \int_0^\infty\mathrm{d}x e^{-x}f(x) \approx \sum_{i=1}^n w_if(x_i),

where :math:`x_i` is the :math:`i`-th root of the Laguerre polynomial and :math:`w_i` are weights.
:math:`x_i` and :math:`w_i` can be computed using `numpy.polynomial.laguerre.laggauss`.
Typically, using a degree of :math:`n=12` is sufficient for this approximation.

.. note::

    There is a typo in the expression for the ionization rate integral in Eq. 32 of :cite:t:`del_zanna_solar_2018`.

The direction ionization cross-section, :math:`\sigma_I`, is computed according to the method of :cite:t:`dere_ionization_2007`
which employs a scaling similar to that used by :cite:t:`burgess_analysis_1992`.
Rearranging Eq. 3 of :cite:t:`dere_ionization_2007`,

.. math::

    \sigma_I = \frac{\Sigma (\log{u} + 1)}{uI^2}

where :math:`u=E/I` is the energy of the incident electron scaled by ionization potential and :math:`\Sigma` is the scaled
cross-section which is defined over,

.. math::

    U = 1 - \frac{\log{f}}{\log{u - 1 + f}}

where :math:`f` is a fitting parameter. :math:`U,f,\Sigma` are all stored in the CHIANTI database such that :math:`\sigma_I`
can be computed for a given :math:`E`. These scaled cross-section data are then interpolated to a given energy array.

The total rate is the summation of :math:`C^I` over all electronic configurations of a given ion for which there is a defined
cross-section.
At a minimum, this includes the outer-shell electron though contributions from inner-shell electrons are also included for some ions.
Sections 3.3 and 3.4 of :cite:t:`young_chianti_2025` provide more details on the calculation of the direct ionization
cross-section and rate and how this is done in CHIANTI.
