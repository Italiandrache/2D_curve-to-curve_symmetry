# 2D_curve-to-curve_symmetry
This Python script determines the symmetric image of a given parametric curve with respect to another parametric differentiable curve in a 2D Cartesian plane. The latter curve is treated as the axis of reflection. The mirror curve equations are defined by ```xMirror_t``` and ```yMirror_t``` and his name by ```mirrorName```; for the to-be-mirrored curve we have ```xToBeMirrored_q```, ```yToBeMirrored_q``` and ```toBeMirroredName``` instead.

I originally conceived this idea in 2015. The code is relatively slow because it employs sympy to perform the necessary mathematical computations. Therefore, I recommend utilizing multiprocessing for improved performance. You can define the number of concurrently active processes via the ```nProcesses``` value. For optimum performance don't exceed the number of physical cores of your CPU (x2 if SMT (AMD) or Hyper-Threading (Intel) is enabled).

The code operates by performing the following main steps:

 - Evaluating the derivative of the mirror curve at each point*.
 - Determining the intersections between each perpendicular line* to the tangents of the mirror curve and the to-be-mirrored curve**.
 - Flipping the vectors that connect the intersection points with the point on the mirror curve from where the perpendicular line was cast***.

*Considering that, in general, it is impossible to compute the derivative for an infinite number of points within a finite time, the code provides the capability to define intervals of the parameter that generates the curve and adjust the point density within those intervals. You can define these values as tuples inside ```tRangeValuesList```. Beware to define ```tRangePlot``` as well, as the parameter interval used to just plot the mirror curve. The ability to selectively define the ```tRange``` parameter intervals and adjust their point density provides an efficient approach to performing the computations required for achieving a precise result. This flexibility avoids unnecessarily lengthy computation times by focusing only on the calculations necessary to attain the desired level of accuracy in different areas of the curve.

**The to-be-mirrored curve might be defined by any two parametric expressions, but the existence of an exact solution set for a system consisting of an affine equation and a nonlinear non-affine equation relies on the nature of the equations and their properties. In certain cases, numerical methods or approximation techniques are necessary to find solutions. However, a system of affine equations always possesses a set of analytic solutions; that's why in the code the to-be-mirrored curve is approximated by a set of line segments. Additionally, by increasing the point density in the range of the to-be-mirrored curve parameter, we can approach the original curve arbitrarily closely, except maybe for unusual curves like the Weierstrass function. You can define the range and the point density of the to-be-mirrored curve parameter as tuples inside ```qRangeValueList```.

***KNOWN ISSUES: There may be instances where the code fails to function correctly, for example when the intersections between the perpendicular line to the tangent at a given point on the mirror curve and the curve to be mirrored are an infinite amount; and there may be some other edge-cases I forgot to consider.

Please note that the steps mentioned above describe the basic procedure of the code, but the actual implementation involves additional considerations and optimizations.

Finally, the code generates a plot showcasing the two parametric curves along with the symmetric curve (beware of defining appropriate x and y boundaries of the plot via ```plt.xlim(xMin, xMax)``` and ```plt.ylim(yMin, yMax)```). This plot is saved in the folder where the script was executed from. Additionally, a .csv file is generated containing the coordinates of the points of the mirrored curve, allowing you to easily export the data if needed.
