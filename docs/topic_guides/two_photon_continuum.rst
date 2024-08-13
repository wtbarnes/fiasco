.. _fiasco-topic-guide-two-photon-continuum:

Calculating the two-photon continuum emission
=============================================

This topic guide provides a detailed description of the two-photon continuum emission calculation
as implemented in :meth:`fiasco.Ion.two_photon`.

In hydrogen-like ions, the transition :math:`2S_{1/2} \rightarrow 1S_{1/2} + h\nu` cannot occur
as an electric dipole transition, but only as a much slower magnetic dipole transition.
The dominant transition then becomes :math:`2S_{1/2} \rightarrow 1S_{1/2} + h\nu_{1} + h\nu_{2}`.

In helium-like ions, the transition from :math:`1s2s ^{1}S_{0} \rightarrow 1s^{2}\ ^{1}S_{0} + h\nu`
is forbidden under quantum selection rules since :math:`\Delta J = 0`.
Similarly, the dominant transition becomes
:math:`1s2s ^{1}S_{0} \rightarrow 1s^{2}\ ^{1}S_{0} + h\nu_{1} + h\nu_{2}`.

In both cases, the energy of the two photons emitted equals the energy difference of the two levels.
As a consequence, no photons can be emitted beneath the rest wavelength for a given transition.
See the introduction of :cite:t:`drake_spontaneous_1986` for a concise description of the process.

The emission is given by

.. math::

    C_{2p}(\lambda, T, n_{e}) = hc \frac{n_{j}(X^{+m}) A_{ji} \lambda_{0} \psi(\frac{\lambda_{0}}{\lambda})}{\psi_{\text{norm}}\lambda^{3}}

where :math:`\lambda_{0}` is rest wavelength of the (disallowed) transition,
:math:`A_{ji}` is the Einstein spontaneous emission coefficient,
:math:`\psi` is so-called spectral distribution function, given approximately by

.. math::

    \psi(y) \approx 2.623 \sqrt{\cos{\Big(\pi(y-\frac{1}{2})\Big)}}

according to :cite:p:`gronenschild_calculated_1978` and :math:`\psi_{\text{norm}}` is a normalization factor such that

.. math::

    \frac{1}{\psi_{\text{norm}}} \int_{0}^{1} \psi(y) dy = 2

for hydrogen-like ions and :math:`1` for helium-like ions.
Finally, :math:`n_{j}(X^{+m})` is the density of ions with charge state :math:`m` of element :math:`X`
in excited state :math:`j` and is given by

.. math::

    n_{j}(X^{+m}) = \frac{n_{j}(X^{+m})}{n(X^{+m})} \frac{n(X^{+m})}{n(X)} \frac{n(X)}{n_{H}} \frac{n_{H}}{n_{e}} n_{e}.
