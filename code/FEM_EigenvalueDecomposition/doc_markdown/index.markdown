# Rotating Isospectral Drums - The Practical Parts # {#mainpage}

This is the documentation of the source code, the additional supporting files and
the quality checks of the software used in the master's thesis.

\remark This documentation contains mathematical formulae which are written in LaTeX
and rendered with MathJax. Thus an internet connection is required for a correct display.

## The Software ##

### Description ###
The present code implements a finite element approach to the discretization of
partial differential equations on arbitrary domains.

The code implements triangular Lagrange finite elements of orders 1-3 in two dimensions.
It allows for a simple automatisation of parameter studies through the use of a parameter
file \ref parameters.txt as well as of parameters which can be passed at launch time.

Provided a mesh, an offset a set of angular velocities and a spatial shift (which may vanish)
of the axis of rotation the program determines first and foremost the eigenvalues (or frequencies)
of the stationary domain. Afterwards the same number of eigenvalues with smallest magnitude
is determined for each angular velocity provided.
In such a way one can follow the evolution of the eigenvalues of eq. (3.20) with increasing
\f$ \omega \f$.

### Compilation ###

To compile the software a C++ compiler with C++11 support is required.
Furthermore the boost program options library should be present on the system.

Libraries for linear algebra [Eigen] an eigenvalue decomposition
[ARPACK++] are provided in the **lib** directory of the medium.
The latter library will require an installation. Eigen is a template library and may be used
out of the box. Optimization options may be added to taste, though `-O3` compiler option is recommended.

The order of the ansatz polynomials can be choosen by defining P1BASIS or P2BASIS at compile-time.
If the definition is absent a cubic basis is used.

### Execution ###

The program utilizes the boost program_options library to parse command-line arguments.
The set of arguments recognized by the program is encoded in prog_opt_init().
The hard-coded values therein are the last fall-back. The next level of user interaction
is the file \ref parameters.txt , where user-provided values for the parameters are provided.
The file has to consist of key=value pairs with each pair being placed in a separate row.
The hard-coded parameter values as well as the values provided in \ref parameters.txt
are superceded by those passed to the program on the command line in the form `-param <value>`.

The current iteration of the software expects a file (\ref AngularFrequencies.dat ) with a set of angular frequencies (one per row) for which
the eigenvalues of a domain have to be determined. The name of the file should be provided as the `-omega`
parameter and the file should be located in the same directory where the program is executed.

The program accepts the following parameters
- *N* (unsigned int) number of refinement steps.
 + default: 3
- *CF* (string) name of the file with vertex coordinates.
 + default: Coordinates.dat
- *EF* (string) name of the file with triangle<->vertex assignments.
 + default: Tri.dat
- *BFD* (string) name of the file with boundary element coordinates of the Dirichlet boundary.
 + default: Boundary_Dirichlet.dat
- *BFN* (string) name of the file with boundary element coordinates of the Neumann boundary.
 + default: Boundary_Neumann.dat
- *n* (double) index of refraction.
 + default: 1.0
- *omega* (string) name of the file containing a list of angular frequencies to be simulated.
 + default: AngularFrequencies.dat
- *Neig* (unsigned int) number of eigenvalues to compute.
 + default: 10
- *tol* (double) relative error of eigenvalue computations (ARPACK).
 + default: 1e-8
- *recomp* {0,1} flag indicating whether the system matrices should be recomputed.
 + default: 1
- *FP* (string) prefix for matrix storage files.
 + default: NONE
- *xShft* (double) x-shift of the geometry.
 + default: 0
- *yShft* (double) y-shift of the geometry.
 + default: 0



## Output ##

When compiled the code will produce multiple output files when run:

- Matrix-Market files (*_stiffness.mm, *_mass.mm, *_damping.mm). The system matrices are
stored to disc using the matrix-market format. If the \a recompute flag is set to 0
the matrices are not computed anew but are instead read from the files. An exception to this
is the full mass matrix, for it requires an update for each angular frequency \f$ \omega \f$.
\attention These files tend to be very large. Care should be taken to ensure sufficient storage
space on the hard  drive. The saveMatrix() functions are responsible for writing the matrices to disc.

- output.dat A data file containing the angular velocities and the N smallest eigenvalues for each velocity in
the following format
| \f$\omega_0\f$ | \f$\omega_1\f$ |  ...  |  \f$\omega_m\f$ |
| :--------------: | :--------------: | :----: | :--------------: |
| \f$\lambda_1\f$  | \f$\lambda_1\f$  |  ...  |  \f$\lambda_1\f$ |
| ...  |  ...  |  ...  |  ... |
| \f$\lambda_N\f$  |  \f$\lambda_N\f$  |  ...  | : \f$\lambda_N\f$ |

The output is handled by saveToFile().

- Eigenvalues_*.dat Files containing the values of the eigenfunctions computed in case linear FE are chosen by defining P1BASIS 
at compile time.

- Plot_Elements*.dat and Plot_Vertices*.dat Files containing the definitions of the mesh
for which the eigenmodes were written. Again available only for linear FE.


## Numerical Experiments ##

Numerical experiments can be run manually, as long as all the necessary parameters are provided.
To facilitate mesh-refinement studies a Python wrapper was written to run a set of simulations en-bloc.

### Execution ###

The wrapper script is \ref NumExpRun.py . It allows one to specify a list of indices
of refraction `-ior <A>,<B>,<C>`, the levels of mesh refinement `-nref 2,3,4` the number of \f$ \omega \f$-values to be spaced logarithmically in the interval
\f$ [0,3\cdot 10^8]\f$ as `-samples 100` as well as the directory with reference parameters `-refdir <path>` and the name of the executable file located therein `-runfile`.
Additionally a working directory can be specified with `-workdir <absolute path>` and a domain prefix ("ID1" or "ID2") can be provided with `-domain`.

If the necessary parameters are provided (or the default values are used) the script will launch numerical simulations for appropriate pairings of the different parameters.
Specifically it will run a mesh-refinement study for each index of refraction using the number of \f$ \omega \f$-samples  provided.
The parameter file \ref parameters.txt as well as the mesh definitions and angular frequencies (\ref AngularFrequencies.dat) are copied from the reference directory
into the directory of the numerical experiment. The latter is created as a subdirectory of the working directory.

At the end of the simulation cycle the directory in which the script was launched will contain multiple subdirectories with names beginning with `Exp_`
containing the computed results for the different indices of refraction.


### Evaluation Guidelines ###

Prior to using the IPython notebooks provided in the **eval** directory for the evaluation of numerical results the raw results obtained in the previous step
have to be preprocessed by executing \ref Evalscript.py in the directory containing `Exp_*` subdirectories. This script will consolidate data
from the pure text files written by my code into compressed *.npz-archives. The latter may then be evaluated with ease with the IPython notebooks
provided with the thesis.

## Addendum ##

Unit and integration tests are described briefly on the additional tests page.
There the additional quality checks of the method are described in some detail.

The description of the parameter file and mesh files is provided on the [supporting files] page.

# License #

The present software is provided under the GNU General Public License 3.0 (GPLv3). The full text of said license
may be found in \ref license.txt .

I hope that academics will, in the near future, stop guarding their code in the same way Gollum guarded
*The Ring*. Whenever results were obtained through simulations, or more generally through numerical means, 
it is the duty of the researcher to publish his or her code to facilitate reproducibility.


[Eigen]: http://eigen.tuxfamily.org "Eigen linear algebra template library"
[ARPACK++]: https://github.com/m-reuter/arpackpp "ARPACK++ - C++ wrappers for ARPACK"
[supporting files]: @ref supportfiles "supporting files"