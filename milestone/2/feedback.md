# High-level description feedback

Cool project idea.

Another method used to relax the N^2 N-body problem are cell lists
(https://en.wikipedia.org/wiki/Cell_lists).  The moment conserving interpolation
kernel I was talking about with Kevin is a B-spline based kernel derived in
Monaghan 1985
(https://www.sciencedirect.com/science/article/pii/0021999185900063?via%3Dihub),
the kernel is callem M4' and often used in particle-to-mesh and mesh-to-particle
operators.

When you compare against python be sure to benchmark the python version on the
compute cluster as well to be fair.  Looking forward to your results.
