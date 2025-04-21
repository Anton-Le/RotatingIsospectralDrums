/*! \file Misc.h
 * \brief Auxiliary structures.
 *  
 * File:   Misc.h
 * Author: Anton Lebedev
 *
 * Created on 6. October 2015 11:15
 * 
 * Declaration of auxiliary structures.
 * Global includes and shorthand typedefs.
 */

#ifndef MISC_H
#define	MISC_H

#include <string>
#include <complex>
#define EIGEN_VECTORIZE_SSE_4_2
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>


/*! \brief A macro for computing a^2, taken from Numerical Recipes for C++. */
#define SQU(a) ((a)*(a))

/*! \fn const uint triVert() 
 * 
 * Number of vertices in a triangle.
 */
#ifdef P2BASIS
    const uint triVert(6);
#elif P3BASIS
    const uint triVert(10);
#else
    const uint triVert(3);
#endif
    
//!  short-hand declaration of double to be able to switch for float.
typedef double dtype;
//! short-hand declaration of complex numbers
typedef std::complex<double> cdtype;
//! short-hand notation for a Eigen vector of variable size.
typedef Eigen::VectorXd Vec;
//! short-hand notation for a Eigen vector of variable size.
typedef Eigen::VectorXcd cVec;
//! short-hand notation for an Eigen sparse matrix in row-major (C) format
typedef Eigen::SparseMatrix<dtype> CSRMat;
//! short-hand notation for an Eigen sparse matrix in row-major (C) format
typedef Eigen::SparseMatrix<cdtype> cCSRMat;
//! short-hand notation for an Eigen dense matrix of variable size
typedef Eigen::MatrixXd DensMat;
//! short-hand notation for an Eigen dense matrix of variable size
typedef Eigen::MatrixXcd cDensMat;

/*! \brief A structure to keep program parameters in.
 * The structure is used store the runtime parameters
 * read using boost::program_options from command-line
 * or a file */
struct params {
    uint Nref;                          //!< number of refinement iterations for the grid.
    std::string CoordinateFilename;     //!< name of the file containing the coordinates.
    std::string ElementsFilename;       //!< name of the file conatining the element<-> coordinate associations.
    std::string BoundaryFilenameD;       //!< name of the file conatining Dirichlet boundary <-> coordinate associations.
    std::string BoundaryFilenameN;       //!< name of the file conatining Neumann boundary <-> coordinate associations.
    dtype n;                            //!< index of refraction.
    std::string omega;                        //!< angular velocity of the cavity.
    uint Neig;                          //!< number of eigenvalues to compute.
    dtype tol;                          //!< relative tolerance of eigenvalue computations.
    uint recompute;                     //!< flag indicating whether the system matrices should be recomputed.
    std::string MatPrefix;              //!< prefix used for the matrix storage files.
    dtype xShft;                        //!< x-shift value to permit shifting the geometry
    dtype yShft;                        //!< y-shift value to permit shifting the geometry
};


/*! \brief Speed of light. */
const dtype c(299792458);

#endif	/* MISC_H */

