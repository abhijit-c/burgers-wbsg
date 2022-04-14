# Burgers-WBSG
An implementation of a well-balanced stochastic Galerkin method for an inviscid
Burgers equation with random inputs. Methodology based off of the paper by Jin,
Xiu, and Zhu[^1].

## Timetable
| Dates         | Task                      |
|---------------|---------------------------|
| 04/11 - 04/14 | Read and understand paper |
| 04/15 - 04/21 | Implement                 |
| 04/22 - 04/24 | Build Presentation        |

# Abhi's Notes on the Paper

- Looking at page 6, it seems like they're used a forward Euler scheme for the
  discretization of both the flux term and the source term. I remember Dr.
  Chertock mentioning we could reuse our central upwind... I suppose she
  probably meant our forward upwind code?

- (3.14) is the main equation we need to implement. Along with precomputing the
  expectations $E[\Phi_k \Phi_m \Phi_n]$; store those in some 3D Matrix?

- The interface method they built to do this WB-SG seems to be pretty custom to
  the given PDE (atleast general to the flux). In particular, we need the form
  $b'(x) q(u)$ in order to do the manipulation 2.7 - 2.8. Of course, if we just
  stick to the paper examples, this isn't a big deal.

- Regarding (3.6), the $\Phi_M(z)$ is probably a typo, should be $\Phi_m(z)$,

- When doing the numerical example (4.2), what do they choose for the right
  boundary condition? The left one is directly specified. Ahh, since they're
  using Forward Euler, they don't need a ghost cell on the right, hence they can
  just run the method to update the right endpoints,

- In section 4.1, t get the true mean and standard deviation, I guess we take
  u = 2- b, substitute in the b given, and use typical expectation formula and
  integrate them over the domain? I suppose that's why they specified by so that
  it was supported only over a small part of the domain.

## References

[^1]: Jin, S., Xiu, D. & Zhu, X. A Well-Balanced Stochastic Galerkin Method for
  Scalar Hyperbolic Balance Laws with Random Inputs. J Sci Comput 67, 1198–1218
  (2016). https://doi.org/10.1007/s10915-015-0124-2
