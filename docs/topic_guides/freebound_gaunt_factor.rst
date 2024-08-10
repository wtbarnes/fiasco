.. _fiasco-topic-guide-freebound-gaunt-factor:

The implementation of the total free-bound Gaunt factor in fiasco
==================================

The calculation of the total free-bound Gaunt factor, used by the free-bound radiative loss function, does not include the
abundances and ionization equilibria, unlike in :cite:t:`mewe_freebound_1986`.

Additionally, the calculation suggested by :cite:t:`mewe_freebound_1986` and used in CHIANTI makes an approximation to an infinite
sum that can be improved upon.  Specifically, Equation 14 of :cite:t:`mewe_freebound_1986` has a simple
analytic solution.  They approximate
.. math::
    f_{1}(Z, n, n_{0} ) = \sum_{1}^{\infty} n^{-3} - \sum_{1}^{n_{0}} n^{-3} = \zeta(3) - \sum_{1}^{n_{0}} n^{-3} \approx 0.21 n_{0}^{-1.5}

where :math: `\zeta(x)` is the Riemann zeta function.

However, the second sum is analytic, :math: `\sum_{1}^{n_{0}} n^{-3} = \zeta(3) + \frac{1}{2}\psi^{(2)}(n_{0}+1)`
where :math: `\psi^{n}(x)` is the n-th derivative of the digamma function (a.k.a. the polygamma function).
So, we can write the full solution as:
.. math::
    f_{1}(Z, n, n_{0}) = \zeta(3) - \sum_{1}^{n_{0}} n^{-3} = - \frac{1}{2}\psi^{(2)}(n_{0}+1)

The final expression is therefore simplified and more accurate than :cite:t:`mewe_freebound_1986`.
