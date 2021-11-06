### Numerical calculation of 2D Larkin-Imry-Ma state.

- Vector fields `L` and `A`, `|L|=1`;
- `A` is chaotic;
- `L` interacts with `A`: `E = Uin * (L*A)^2`;
- `L` has gradient energy: `E = Ugr * ((dL/dx)^2 + (dL/dy)^2)`.

One can start from some initial distribution of `L` and look for the
equilibrium distribution of `L`. `Ugr` is a control parameter.

There are a few non-trivial effects which can be seen here:

1. Larkin-Imry-Ma state: if gradient energy is strong then `L` is
uniform, if gradient energy is weak, `L` is chaotic (following `A` field).
Between these limits there is an intermediate state where `L` form some
structures of size `R` which can be estimated in the following way:

By rotating a uniform area of size `R` (with `N \propto R^2` impurities)
you can decrease interaction energy by `\propto Uin/sqrt(N) \propto Uin/R` and
increase gradient energy by `\propto Ugr/R^2`. There is an
optimal `R` where these two energies are equal: `R \propto Ugr/Uin`.

2. You can choose one of two different interactions: `E = Uin * (L*A)^2`
or `E = Uin * (L*A)`. If you decrease gradient energy and then increase
it again, in the first case you will return to the original distribution
of `L` (uniform or with vortices), in the second case you will have lots
of random vortices. The reason is that in the first case interaction with
`A` will never rotate `L` more then 90 degrees and information about
original orientation will survive.

### Examples:

picture1: start from uniform L, gradient energy Ugr goes up-down-up:
https://slazav.xyz/tmp/lim1u.gif

picture2: start from L with a vortex, Ugr goes up-down-up. The vortex is stable:
https://slazav.xyz/tmp/lim1v.gif

picture3: start from L with a gradient of phase:
https://slazav.xyz/tmp/lim2s.gif

picture4: start from L with a few vortices:
https://slazav.xyz/tmp/lim3s.gif

