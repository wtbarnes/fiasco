.. _fiasco-topic-guide-freebound-gaunt-factor:

The implementation of the total free-bound Gaunt factor in fiasco
=================================================================

The calculation of the wavelength-averaged total free-bound Gaunt factor, :math:`G_{fb}`,
as suggested by :cite:t:`mewe_calculated_1986` and used in CHIANTI to compute the total free-bound
radiative losses (:meth:`fiasco.Ion.free_bound_radiative_loss`) makes an approximation to an infinite
sum that can be improved upon.
Specifically, Equation 14 of :cite:t:`mewe_calculated_1986` has a simple analytic solution.
They make the approximation,

.. math::

    f_{1}(Z, n, n_{0} ) = \sum_{1}^{\infty} n^{-3} - \sum_{1}^{n_{0}} n^{-3} = \zeta(3) - \sum_{1}^{n_{0}} n^{-3} \approx 0.21 n_{0}^{-1.5}

where :math:`\zeta(x)` is the Riemann zeta function.
However, the second sum is analytic,

.. math::

    \sum_{1}^{n_{0}} n^{-3} = \zeta(3) + \frac{1}{2}\psi^{(2)}(n_{0}+1)

where :math:`\psi^{n}(x)` is the :math:`n`-th derivative of the digamma function (or polygamma function).
As such, we can write the full solution of Equation 14 as,

.. math::

    f_{1}(Z, n, n_{0}) = \zeta(3) - \sum_{1}^{n_{0}} n^{-3} = - \frac{1}{2}\psi^{(2)}(n_{0}+1)

The final expression is therefore simplified and more accurate than the approximation used by :cite:t:`mewe_calculated_1986`.

Note also that, unlike in :cite:t:`mewe_calculated_1986`, the calculation of the total free-bound Gaunt factor,
used by :meth:`~fiasco.Ion.free_bound_radiative_loss`, does not include the abundances and ionization equilibria.
