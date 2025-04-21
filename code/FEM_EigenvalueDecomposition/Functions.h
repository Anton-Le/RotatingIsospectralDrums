/*! \file Functions.h
 *  \brief Function declarations.
 * \copyright GNU General Public License 3.0 (GPLv3)
 * 
 * File:   Functions.h
 * Author: Anton Lebedev
 *
 * Created on 6. October 2015 11:15
 * 
 * Declaration of functions used to acquire and process
 * the data stored in the geometry/geometry2 structs.
 */

#ifndef FUNCTIONS_H
#define	FUNCTIONS_H

#include "Geometry.h"

/*! \brief Function for retrieving options from cmd-line or file.
 * 
 * A function for retrieveing program options using
 * boost::program_options.
 * 
 * \param[in] argc  argument count, standard parameter passed from main.
 * \param[in] argv[]  array of characters containing the command-line arguments.
 *                Standard parameter passed from main.
 * \param[out] par  parameter structure which will contain the retrieved parameters.
 * 
 * The function fills the par struct with the parameters provided by the
 * user in a parameters.txt file or with default values. 
 */
void prog_opt_init(const int argc, const char *const argv[], params& par);

/*! \brief Function for reading the geometry data from files. 
 * 
 * \param[in] VertexFilename - Name of the file containing the vertex coordinates.
 * \param[in] ElementFilename - Name of the file containing the association of the element
 * vertices to the mesh vertices.
 * \param[in] BoundaryFilenameDirichlet - Name of the file containing the association of the
 * boundary element vertices to the mesh vertices for the Dirichlet boundary.
 * \param[in] BoundaryFilenameNeumann - Name of the file containing the association of the 
 * boundary element vertices to the mesh vertices for the Neumann boundary.
 * \return geo - geometry2 structure with problem data.
 *  
 * The function reads the mesh data from the files passed as arguments, stores
 * it in std::vector objects (no manual specification of sizes required) and then
 * transfers the data to the geometry2 object. */
geometry2 readData2D(const std::string VertexFilename, const std::string ElementFilename, const std::string BoundaryFilenameDirichlet,\
        const std::string BoundaryFilenameNeumann);

/*! \brief Function for reading the list of angular frequencies.
 * 
 * A function which retrieves the list of all angular frequencies from a file
 * and stores it into a vector.
 * 
 * \param[in] omegaFile  Name of the file containing the list.
 * \return omegaList Vector of doubles containing the angular frequencies.
 */
std::vector<dtype> readOmega(const std::string omegaFile);

/*! \brief Output eigenvalues to a file.
 * 
 * A function which takes an array of eigenvalues and a vector of angular frequencies
 * and writes them to output.dat
 * 
 * \param[in] Neig number of computed eigenvalues.
 * \param[in] omega  vector of angular frequencies.
 * \param[in] eigval rectangular array of complex eigenvalues. Size #angular frequencies x 2*Neig.
 */
void saveToFile(const uint Neig, const std::vector<dtype> &omega, cdtype **eigval);

/*! \brief Mesh refinement function.
 * 
 * \param[in,out] Pars - structure containing mesh data.
 * \param[in] Nref - number of refinement steps (>=0).
 * 
 * The function takes a coarse mesh stored in \a Pars ,
 * refines it \a Nref times and overwrites the original mesh
 * with the refined version.
 * Note that due to the use of Euler's formula for planar graphs the
 * function only works for meshes of simply connected domains of the
 * same topology class as a disc.
 */
void ref2d(geometry2& Pars, const uint Nref);

/*! \brief Vertex insertion for quadratic elements.
 *
 * \param[in,out] Pars - structure containing mesh data.
 * 
 * The function takes a **final**, i.e., fully refined, mesh stored in \a Pars
 * and overwrites it with a version where additional nodes for \f$ \mathcal{P}_2 \f$
 * finite elements have been added. The nodes are the midpoints of the edges.
 * 
 * \sa p3basepoints, Phi_2, dPhi_2
 */
void p2basepoints(geometry2& Pars);

/*! \brief Vertex insertion for cubic elements.
 *
 * \param Pars - structure containing mesh data.
 * \return Pars - structure containing an expanded mesh
 * with edge midpoints added.
 * 
 * The function takes a fully refined mesh stored in \a Pars
 * and overwrites it with a version with additional nodes for \f$ \mathcal{P}_2 \f$
 * added to each triangle.
 * 
 * \sa p2basepoints, Phi_3, dPhi_3
 */
void p3basepoints(geometry2& Pars);

//functions for the problem

/*! \brief Basis function on the standard simplex.
 * 
 * \param[in] i number of the base function (0,1,2).
 * \param[in] u x-value on the standard simplex.
 * \param[in] v y-value on the standard simplex.
 * \return Phi_i(u,v) - function value.
 * 
 * The function computes the value of a basis function \f$ \phi_i \f$ on the standard simplex
 * in \f$ \mathbb{R}^2 \f$ for linear ansatz functions.
 */
dtype Phi_1(const uint i, const dtype& u, const dtype& v);

/*! \brief Gradient of Phi_i.
 * 
 * \param[in] i index of the base function (0,1,2).
 * \param[in] u,v The coordinates on the std-simplex in 2D.
 * \param[out] x vector to contain the gradient.
  * 
 * The function takes the coordinates \a u,v on a standard simplex and
 * computes the gradient of the \a i-th basis function \f$ \phi \f$ on it.
 * The gradient is returned in \a x , which is overwritten on exit.
 */
void dPhi_1(const uint i, Eigen::Vector2d& x, const dtype& u, const dtype& v);


/*! \brief Basis function on the standard simplex.
 * 
 * \param[in] i number of the base function (0..5).
 * \param[in] u,v x,y-value on the standard simplex.
 * \return Function value: Phi_i(u,v).
 * 
 * The function computes the value of the \a i-th basis function 
 * at the coordinates \a u,v in the standard simplex on \f$ \mathbb{R}^2 \f$.
 * The function uses quadratic Largange polynomials.
 */
dtype Phi_2(const uint i, const dtype& u, const dtype& v);

/*! \brief Gradient of Phi_i.
 * 
 * \param[in] i index of the base function (0..5).
 * \param[in,out] x vector to contain the gradient.
 * 
 * The function computes the gradient of the \a i-th base function Phi_i(u,v)
 * on the standard simplex. The position at which the gradient is
 * computed is provided using \a u , \a v and the vector \a x is overwritten
 * with the gradient of \a Phi_2 on exit.
 * 
 * \sa Phi_2, dPhi_1, dPhi_3
 */
void dPhi_2(const uint i, Eigen::Vector2d& x, const dtype& u, const dtype& v);

/*! \brief Basis function on the standard simplex.
 * 
 * \param[in] i number of the base function (0..9).
 * \param[in] u,v x,y-value on the standard simplex.
 * \return Phi_i(u,v) - function value.
 * 
 * The function computes the value of the \a i-th basis function
 * on the standard simplex in 2D at the point \a (u,v) therein.
 * The function implements cubic ansatz polynomials.
 * 
 * \sa Phi_2, dPhi_1, dPhi_3
 */
dtype Phi_3(const uint i, const dtype& u, const dtype& v);

/*! \brief Gradient of Phi_i.
 * 
 * \param[in] i index of the base function (0..9).
 * \param[out] x vector to contain the gradient.
 * 
 * Computes the value of the gradient of \a i-th cubic basis function at the
 * point \a (u,v) and stores it in the vector \a x , which is overwritten on exit.
 */
void dPhi_3(const uint i, Eigen::Vector2d& x, const dtype& u, const dtype& v);

/*! \brief Compute the jacobian matrix of the transform
 * 
 * \param[in] el - vertex indices of the target triangle.
 * \param[in] vertex - array of vertices of the mesh.
 * \param[out] J - \f$ 2\times 2\f$ Matrix.
 * 
 * The function computes the Jacobi matrix \a J of the transformation \f$ \chi \f$
 * from the 2D standard simplex to the triangle \a el defined by the vertices
 * provided in \a vertex . The matrix \a J is overwritten on exit
 */
void Jacobian(const Tri& el, const coord<>* vertex, Eigen::Matrix2d& J);

/*! \brief Compute the inverse TRANSPOSED jacobian matrix of the transform.
 * 
 * \param[in] el - vertex indices of the target triangle.
 * \param[in] vertex - array of vertices of the mesh.
 * \param[out] J - \f$2\times 2\f$ Matrix.
 * 
 * Computes the inverse **transposed** Jacobi matrix of the mapping from
 * the standard simplex to the target triangle. The function requires the
 * indices of the triangle vertices and the array containing the coordinates
 * of the mesh vertices. The matrix \a J is overwritten in the routine.
 */
void InvJacobian(const Tri& el, const coord<>* vertex, Eigen::Matrix2d& J);

/*!\brief Transformation function.
 * 
 * The function transforms a given pair of coordinate values on the standard simplex
 * to a pair of coordinates in the target triangle.
 * Executes an affine transformation. The functions Jacobian and InvJacobian
 * are 'derivatives' of this transform.
 * 
 * \param[in] el - vertex indices of the target triangle.
 * \param[in] vertex - array of vertices of the mesh.
 * \param[in,out] x - vector containing the coordinates in the simplex.
 * 
 * The function takes a list of vertices \a vertex and the indices of the
 * element \a el and computes the image of the point stored in \a x
 * under the transformation \f$ \chi: \hat{T}\rightarrow T \f$.
 * The vector \a x is overwritten with the coordinates in the real triangle \a el.
 */
void g(const Tri& el, const coord<>* vertex, Eigen::Vector2d& x);

/*! \brief 1D basis function on [-1,1].
 * 
 * \attention The function implements a base function in 1 dimension and is not suitable for
 * higher-dimensional calculations. It is necessary when computing boundary integrals in 2D, but
 * is valid only for \f$ \mathcal{P}_1 \f$ finite elements.
 * 
 * The function implements the base function \Phi_i \in P1 on the reference element I=[-1,1].
 * 
 * \param i - {0,1} the index of the function of the basis. 0 corresponds to the
 * left node (x = -1) whilst 1 corresponds to the right node (x = 1).
 * \param x - local coordinate on [-1,1].
 * \return Phi_i(x) - the value of the basis function i at x.
 */
dtype Phi1D(const uint &i, const dtype& x);

/*! \brief Derivative of the 1D basis function.
 * 
 * \attention The function is valid in 1D or for boundary integrals in 2D, but there
 * only for \f$ \mathcal{P}_1 \f$ finite elements.
 * 
 * The function computes the derivative of the basis function Phi_i
 * on I=[-1,1]. Which is either -1/2 or 1/2.
 * 
 * \param i - {0,1} the index of the function of the basis. 0 corresponds to the
 * left node (x = -1) whilst 1 corresponds to the right node (x = 1).
 * \return +-0.5  - the derivative of Phi_i(x), depending only on i. 
 */
dtype dPhi1D(const uint &i);

//interval mapping
/*! \brief Mapping from I=[-1,1] to \f$ \mathbb{R}^2 \f$.
 * 
 * \param[in] i Index (global) of the starting (left) point of the line.
 * \param[in] j Index (global) of the ending (right) point of the line.
 * \param[in] s Local coordinate on I=[-1,1].
 * \param[in] *x Pointer to the array containing all the vertices of the triangulation.
 * \param[out] y Eigen vector to which the result of the mapping is stored.
 * \return y 2D vector with y = 1/2 * (x[j] - x[i])*s + 1/2 *(x[j] + x[i]). 
 * 
 * The function computes the mapping from I=[-1,1] to R^2 for a parametrized
 * edge starting at \a x[i] and ending at \a x[j].
 * On exit the global coordinates (x,y) of the point with the edge-coordinate \a s
 * are stored in the first and second component of \a y.
 * 
 * The mapping is:
 * \f[
  y = \frac{1}{2}\left( x[j] -x[i] \right) s + \frac{1}{2}\left( x[j]+ x[i] \right).
 * \f]
 */
void g1D(const uint &i, const uint &j, const dtype &s, const coord<>* x, Eigen::Vector2d& y);

/*! \brief Jacobian of the 1D mapping g1D.
 * 
 * \param[in] i Index (global) of the starting (left) point of the line.
 * \param[in] j Index (global) of the ending (right) point of the line.
 * \param[in] x Pointer to the array containing all the vertices of the triangulation.
 * \return 0.5||x[j]-x[i]|| The Jacobian of the affine transformation.
 * 
  * The function computes the Jacobian (measure transform) of the mapping from I=[-1,1]
 * to a line defined by x[i], x[j].
 */
dtype dg1D(const uint &i, const uint &j, const coord<>* x);


//------------------------------------------------------------------------------
/*! \brief Determine number of non-zero entries per row.
 * 
 * \param[in] geo Structure of type geometry2 which contains the geometry data.
 * \param[out] NnzPerRow Number of non-zero entries per row of the matrix, for each row.
 * 
 * The function determines the number of non-vanishing entries for each
 * row of a matrix (M, A, V) for one component.
 * 
 * From NnzPerRow a sparse (CRS) matrix can be set-up for the problem by
 * concatenating two of these vectors to account for real and imaginary parts.
 */
void DetermineSparsityPattern(const geometry2 &geo, Eigen::VectorXi &NnzPerRow);

/*! \brief Set-up of the mass matrix.
 * 
 * \param[out] M *Uncompressed* and initialized matrix in CSR/CRS format.
 * \param geo Structure of type geometry2 containing the geometry of the domain.
 * \param[in] GD A nx3 array containing in its rows the nodes and weights of the gauss
 * quadrature.
 * \param[in] QuadOrder Number of nodes of the quadrature (#rows of GD).
 * \param[in] freeNodes Number of free (non-Dirichlet) nodes in the mesh.
 * \param[in] NnzPerRow 
 * \parblock 
 * Eigen integer vector of length \a freeNodes.
 * A vector containing the estimate for the number of non-zeros per row
 * of the matrix.
 * \endparblock
 * 
 * The function fills a given matrix \a M with entries of a mass matrix
 * of a FEM. The routine requires a Gauss-type quadrature in 2D provided
 * by a list \a GD of weights and nodes along with the number \a QuadOrder of these.
 * The number of \a freeNodes is used for the resizing of the matrix and 
 * \a NnzPerRow contains the number of non-zero entries per row which allows
 * a more efficient filling of the matrix by avoiding constant resize-and-copy operations.
 * The matrix \a M is resized and overwritten.
 *  
 */
void MassMatrixSetup(CSRMat& M, const geometry2& geo, const dtype GD[][3], const uint QuadOrder,
        const uint& freeNodes, const Eigen::VectorXi& NnzPerRow);

/*! \brief Set-up of the stiffness matrix.
 * 
 * \param[in,out] A **Uncompressed** and initialized matrix in CSR/CRS format.
 * \param[in] geo Structure of type geometry2 containing the geometry data of the domain.
 * \param[in] GD A nx3 array containing in it's rows the nodes and weights of the gauss
 * quadrature.
 * \param[in] QuadOrder Number of nodes of the quadrature (#rows of GD).
 * \param[in] freeNodes Number of free (not bound by Dirichlet conditions) nodes
 * of the mesh.
 * \param[in] NnzPerRow Eigen integer vector of length \a freeNodes .
 * 
 * 
 * The function fills a given matrix \a A with entries of a FEM stiffness matrix 
 * where nodes and base functions are taken from the \a geo structure.
 * The parameters \a NnzPerRow and \a freeNodes serve the purpose of accelerating
 * the filling of the matrix by providing a precomputed estimate of the size of
 * the memory needed for each row of the matrix. This allows one
 * to avoid frequent reallocations and memcopy.
 */
void StiffnessMatrixSetup(CSRMat& A, const geometry2& geo, const dtype GD[][3], const uint QuadOrder, 
        const uint& freeNodes, const Eigen::VectorXi& NnzPerRow);

/*! \brief Set-up of the damping (or transport) matrix.
 * 
 * \param[out] C **Uncompressed** and initialized matrix in CRS/CSR format.
 * \param[in] geo Struct of type geometry2 containing the geometry data of the domain.
 * \param[in] GD a nx3 array containing in it's rows the nodes and weights of the gauss
 * quadrature.
 * \param[in] QuadOrder Number of nodes of the quadrature (# of rows of GD).
 * \param[in] freeNodes Number of free (not bound by Dirichlet conditions) nodes
 * of the mesh.
 * \param[in] NnzPerRow Eigen integer vector of length \a freeNodes .
 * 
 * 
 * The function fills the provided matrix \a C with entries of a FEM gyroscopic matrix.
 * The matrix is a discrete version of an azimuthal unit vector field  \f$\partial_\varphi \f$
 * in the real plane. It implements eq. (6.12c) from the thesis.
 * 
 * Though Dirichlet boundary conditions are not explicitly assumed in the function, its
 * form is nevertheless indicative of this assumption seeing as 
 * terms implementing necessary boundary integrals are missing.
 * 
 * The order of the ansatz functions and the number of associated nodes per triangle
 * is given globally using PkELEMENTS preprocessor directive.
 * 
 * The parameters \a QuadOrder and \a GD are used to pass the data necessary for
 * an implementation of a numerical quadrature and \a NnzPerRow and \a freeNodes
 * are used for resizing and efficient memory handling.
 * 
 * On exit the matrix \a C is overwritten with new entries.
 */
void DampingMatrixSetup(CSRMat& C, const geometry2& geo, const dtype GD[][3], const uint QuadOrder, 
        const uint& freeNodes, const Eigen::VectorXi& NnzPerRow);

/*! \brief Set-up of the full mass matrix.
 * 
 * \param[out] M **Uncompressed** and initialized matrix.
 * \param[in] geo Struct of type geometry2 containing the geometry data of the domain.
 * \param[in] GD A nx3 array containing in it's rows the nodes and weights of the gauss
 * quadrature.
 * \param[in] QuadOrder Number of nodes of the quadrature (# of rows of GD).
 * \param[in] freeNodes Number of free (not bound by Dirichlet conditions) nodes
 * of the mesh.
 * \param[in] NnzPerRow Eigen integer vector of length = number of free nodes in the grid.
 * \param[in] omega Angular frequency of the domain.
 * \param[in] IOR Index of refraction of the medium the domain is comosed of.
 * 
 * 
 * The function fills a matrix \a M in CSR/CSC format with entries corresponding to a discrete
 * version of the full mass matrix (\f$ \gamma\cdot\left( n^2-\left(\frac{r\omega}{c}\right)^2\right) \f$).
 * For the mass matrix no assumptions as to the type of the boundary conditions are necessary.
 * The function implements eq. (6.12d) of the accompanying thesis.
 * 
 * Again, \a QuadOrder and \a GD are used for the implementation of a numerical quadrature and
 * \a freeNodes and \a NnzPerRow for memory handling.
 * 
 * The angular frequency parameter \a omega and the index of refraction \a IOR are necessary
 * to compute the coefficient which scales the whole matrix.
 * 
 * The entries depend, of course, on the order of the ansatz polynomials \f$ \mathcal{P}_k \f$.
 *  */
void MassMatrixSetup_Full(CSRMat& M, const geometry2& geo, const dtype GD[][3], const uint QuadOrder,
                          const uint& freeNodes, const Eigen::VectorXi& NnzPerRow, const dtype& omega, const dtype& IOR);

/*! \brief Set-up of the auxiliary stiffness matrix.
 * 
 * \param[out] A **Uncompressed** and initialized matrix.
 * \param[in] geo Struct of type geometry2 containing the geometry data of the domain.
 * \param[in] GD a nx3 array containing in it's rows the nodes and weights of the gauss
 * quadrature.
 * \param[in] QuadOrder Number of nodes of the quadrature (# of rows of GD).
 * \param[in] freeNodes Number of free (not bound by Dirichlet conditions) nodes
 * of the mesh.
 * \param[in] NnzPerRow Eigen integer vector of length = number of free nodes in the grid.
 * 
 * 
 * The function fills a matrix \a A with entries corresponding to a discretization of
 * \f$ \partial_\varphi^2 \f$ on a FE mesh provided in \a geo .
 * 
 * The matrix is then scaled by \f$ \left(\frac{\omega}{c}\right)^2 \f$ and subtracted from
 * the original stiffness matrix (which is a discrete version of the Laplacian), hence the prefix "auxiliary".
 * This function is the implementation of eq. (6.12b) of the accompanying thesis.
 * 
 * The significance of other parameters is as for 
 * \sa MassMatrixSetup, StiffnessMatrixSetup, DampingMatrixSetup
 */
void StiffnessMatrixSetup_Auxiliary(CSRMat& A, const geometry2& geo, const dtype GD[][3], const uint QuadOrder,
                                    const uint& freeNodes, const Eigen::VectorXi& NnzPerRow);

/*! \brief Function for shifting the mesh in space.
 * 
 * \param[in,out] geo Struct of type geometry2 containing the geometry data to be modified.
 * \param[in] x,y The shift vector.
 * 
* The function shifts all vertices of the mesh stored in \a geo
 * by \a x in the x-direction and \a y in the y-direction.
 * The mesh \a geo is overwritten.
 * 
 * \attention This function should be called before any other operation, save for
 * refinement and node insertion, is performed with the mesh. Specifically it should
 * be called before any matrix is set-up using the mesh.
 */
void shift(geometry2 &geo, const dtype& x, const dtype& y);

/*! \brief Function which writes eigenvectors to disc for a particular angular frequency.
 * 
 * \param[in] Fileprefix Prefix for the name file name.
 * \param[in] omega Angular frequency (file suffix).
 * \param[in] Neig Number of eigenvectors to be written.
 * \param[in] geo Structure contining the mesh.
 * \param[in] EV Array of dimension Neig x freeNodes containing the computed eigenvectors.
 * \param[in] writeMesh Flag indicating whether the mesh should be written as well.
 * \param[in] GD a nx3 array containing in it's rows the nodes and weights of the gauss
 * quadrature.
 * \param[in] QuadOrder Number of nodes of the quadrature (# of rows of GD).
 * 
 * The function stores computed eigenvalues on disc in a simple text file
 * with the name *Eigenvalues_<Fileprefix>_omega_<omega>.dat*.
 * If the **writeMesh** flag is set the mesh is saved to
 * *Plot_Vertices_<Fileprefix>.dat* and *Plot_Elements_<Fileprefix>.dat*.
 * 
 */
void saveEV(const std::string Fileprefix, const dtype omega, const uint Neig, const geometry2& geo, cdtype** const EV,
        const bool writeMesh, const uint freeNodes, const dtype GD[][3], const uint QuadOrder);
#endif	/* FUNCTIONS_H */

