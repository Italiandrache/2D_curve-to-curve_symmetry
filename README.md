# 2D_curve-to-curve_symmetry
This Python code determines the symmetry of a given parametric curve with respect to another parametric differentiable curve in a 2D plane. The latter curve is treated as the axis of reflection.

I originally conceived this idea in 2015. The code is relatively slow because it employs sympy to perform the necessary mathematical computations. Therefore, I recommend utilizing multiprocessing for improved performance.

KNOWN ISSUES:
There may be instances where the code fails to function correctly, for example when the intersections between the perpendicular line to the tangent at a given point on the mirror curve and the curve to be mirrored are infinite.
