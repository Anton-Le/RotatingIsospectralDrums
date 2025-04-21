# Additional tests of functionality. #       {#tests}


As remarked upon in the thesis (Appendix B) the integrated tests demonstrating
the correct functionality of the software are presented here. 
Additionally unit tests which serve as a check of basic functions are listed as well.


## The Poisson Problem ## 

The basic integration test is the classic Poisson problem on the unit square
\f[ -\Delta \varphi (x,y) = f(x,y)\\
(x,y) \in [0,1]^2\;. \f]

Assuming correct implementation of the mesh input and refinement 
the Poisson problem is solved for a given inhomogeneity \f$ f(x,y) \f$ numerically
and the solution is compared to the known analytic solution \f$ u(x,y) \f$.

The tests were developed by first picking an appropriate function \f$ u \in C^\infty \f$,
and then acting on it with the appropriate differential operator to obtain \f$ f \f$.
The result was then encoded in the test.

The verification of correctness is done by computing the \f$ L^2 \f$ or, when possible,
the \f$ H^1 \f$ norms of the errors of numerical solutions. Additionally the \f$ L^2 \f$ norm
of the inhomogeneity is computed. Using these quantities the behaviour of the 
method for decreasing mesh width can be studied and compared to theoretical predictions.

Each test has been implemented in such a way that if it is compiled with #P1BASIS ,
#P2BASIS or without a preprocessor directive (defaults to #P3BASIS) a mesh-refinement
study is run when the resulting binary is invoked.
The resulting errors are printed on screen and may be collected by hand for visualization purposes.

\attention Not all mesh width' are used for each order of the ansatz functions.
This is due to the fact that for coarse meshes and low orders the resulting matrices
consist of only one entry.


### Polynomial Solution ###     

The functions \f$ u,f \f$ for the first test were taken from \cite Kressner2010 .
They are \f[
u(x,y) = x^3(x-1)y(y-1)\\
f(x,y) = (12 x^2-6x)y(y-1) + 2 x^3(x-1)\;.
\f]

This pair is suitable as a first integrated test, for the polynomials and their products
can be integrated exactly by the chosen quadrature rule, which eliminates the
quadrature error from the list of potential problems. 

This test has been implemented in \ref QualityCheck.cpp .

In the image below the \f$ L^2,H^1\f$ errors are shown as functions of the mesh width. 
It can be clearly seen that the errors display the expected behaviour.
This shows that the basic implementation of the FEM is solid.


![The behaviour of errors w.r.t. the mesh width for different orders of ansatz polynomials.](/home/solid/secondary/lebedev/Dokumente/Master_Thesis/Plots/tests/Comparisons/PoissonPolynomialTest.png)

### Trigonometric Solution ###  

The trigonometric solution refers to the choice of the eigenfunctions of the
square as the solutions. Unlike the polynomial of the preceding section the 
Taylor expansion of the trigonometric functions is infinite. This guarantees
that the functions can not be integrated exactly by the Gauss quadrature.
Additionally these functions are oscillatory, which serves as a good test for
a static approximation of eigenmodes of a drum.

The functions were chosen to be
\f[
u(x,y) = \sin(\pi x)\sin(2\pi y)\\
f(x,y) = 5\pi^2 \sin(\pi x)\sin(2\pi y)\;.
\f]

This test is implemented in \ref TrigQualityCheck.cpp .

From the following image one can infer that the FEM implementation is suitable even for
oscillatory functions and thus for the consideration of eigenmodes.

![The behaviour of errors w.r.t. the mesh width for different orders of ansatz polynomials.](/home/solid/secondary/lebedev/Dokumente/Master_Thesis/Plots/tests/Comparisons/PoissonTrigonometricTest.png)

### Poisson with transport and trigonometric solution ###

The last test uses a different differential operator. In this case the following 
equation is to be solved:
\f[ \Delta u + \alpha (-y\partial_x + x\partial_y)u = f\;. \f]

This equation is a step into the direction of eq (3.20) from the thesis in so far
that it already utilizes the transport term \f$ -y\partial_x + x\partial_y \f$.

The functions used for this test are
\f[
u(x,y) = x^3(1-x)y(y-1)\\
f(x,y) = x(6y(1-y) + \alpha x^4(2y-1) + x^3(2+\alpha - 2\alpha y) + 3x(y-1)y(\alpha y +4)
+x^2(-2+4y^2\alpha -4y^3\alpha) )\;.
\f]

This test is implemented in \ref TransportCheck.cpp . The parameter value is chosen to be quite large: \f$ \alpha = 10 \f$.

The final results show that the current implementation of simple Lagrange finite elements
is, at least in principle, suitable for the treatment of equations with a discrete version of a rotational vector field \f$\partial_\varphi\f$.

![The behaviour of errors w.r.t. the mesh width for different orders of ansatz polynomials.](/home/solid/secondary/lebedev/Dokumente/Master_Thesis/Plots/tests/Comparisons/PoissonTransportTest.png)

## Miscellaneous Tests. ##

The miscellaneous tests are comprised mostly of basic checks of functionality,
i.e. unit tests. Contrary to the integration tests they do not check the interplay between 
a lot of functions.
These unit tests for \f$ \mathcal{P}_{2,3} \f$ finite elements have been implemented
in \ref P2_MeshTest.cpp and \ref P3_MeshTest.cpp .

In both cases reference finite element matrices are hard-coded into the routines.
The matrices generated by the functions MassMatrixSetup() , StiffnessMatrixSetup() 
and DampingMatrixSetup() are then compared to the reference matrices. The admissible error
is set to \f$ \epsilon=10^{-10} \f$.

Additionally the mesh refinement and the insertion of additional points for higher-order FE
is tested by acting on a refined or coarse mesh with ref2d() , p2basepoints() or p3basepoints() .