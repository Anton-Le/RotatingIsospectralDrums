/*! \file TransportCheck.cpp
 * \author Anton Lebedev
 * \copyright GNU General Public License 3.0 (GPLv3)
 * \brief Integration test with an augmented Poisson equation.
 * 
 * File:   TransportCheck.cpp
 * Author: lebedev
 *
 * Created on 22.07.2016, 12:30:55
 */

#include <cstdlib>
#include <iostream>
#include "Geometry.h"
#include "Functions.h"
#include <Eigen/Cholesky>
#include <Eigen/Sparse>
#include "GaussQuadParam.h"

//base functions from Functions.cpp
extern dtype(*Phi)(const uint, const dtype&, const dtype&);
extern void (*dPhi)(const uint, Eigen::Vector2d&, const dtype&, const dtype&);

const dtype Gauss[][3] = {
    {0.3333333333333333, 0.3333333333333333, 0.0378610912003147 / 2},
    {0.7974269853530872, 0.1012865073234563, 0.0376204254131829 / 2},
    {0.1012865073234563, 0.7974269853530872, 0.0376204254131829 / 2},
    {0.1012865073234563, 0.1012865073234563, 0.0376204254131829 / 2},
    {0.0597158717897698, 0.4701420641051151, 0.0783573522441174 / 2},
    {0.4701420641051151, 0.0597158717897698, 0.0783573522441174 / 2},
    {0.4701420641051151, 0.4701420641051151, 0.0783573522441174 / 2},
    {0.5357953464498992, 0.2321023267750504, 0.1162714796569659 / 2},
    {0.2321023267750504, 0.5357953464498992, 0.1162714796569659 / 2},
    {0.2321023267750504, 0.2321023267750504, 0.1162714796569659 / 2},
    {0.9410382782311209, 0.0294808608844396, 0.0134442673751655 / 2},
    {0.0294808608844396, 0.9410382782311209, 0.0134442673751655 / 2},
    {0.0294808608844396, 0.0294808608844396, 0.0134442673751655 / 2},
    {0.7384168123405100, 0.2321023267750504, 0.0375097224552317 / 2},
    {0.7384168123405100, 0.0294808608844396, 0.0375097224552317 / 2},
    {0.2321023267750504, 0.7384168123405100, 0.0375097224552317 / 2},
    {0.2321023267750504, 0.0294808608844396, 0.0375097224552317 / 2},
    {0.0294808608844396, 0.7384168123405100, 0.0375097224552317 / 2},
    {0.0294808608844396, 0.2321023267750504, 0.0375097224552317 / 2},
};

const uint Order(19);
const dtype alpha(10);

using namespace std;

inline dtype f(const dtype& x, const dtype& y);
inline dtype u(const dtype& x, const dtype& y);
void du(Eigen::Vector2d& v, const dtype& x, const dtype& y);

/*! \brief A wrapper for basis generation.
 * 
 * The function is a wrapper for NULL, p2basis, p3basis. It
 * inserts the appropriate basis nodes into a mesh depending on
 * the degree of the Lagrange ansatz polynomials chosen at compile time.
 * \param geo [in,out] Mesh to be modified.
 *  */
void pKbasis(geometry2 & geo);

/*! \brief Computes the approx. L^2 deviation of a numerical solution from an analytic one.
 * 
 * Function which computes an approximation of the L^2 error $||u-u_{h}||_{L^2}$
 *  of a numerical solution using the analytic solution and Gauss quadrature
 * to approximate the integrals.
 * NOTE! It is vital that the quadrature is of sufficiently high order due to
 * the occurence of products of polynomials.
 * \param vals [in] Vector containing the computed values of the approximation.
 * \param geo [in] Struct of type geometry2 containing the mesh of the domain.
 * \param GD [in] Array of nodes and weights of Gauss quadrature.
 * \param QuadOrder [in] Number of nodes of the quadrature
 * \return Err L^2 error approximation.
 *  */
dtype L2FuncErr(const Vec& vals, const geometry2& geo, const dtype GD[][3], const uint QuadOrder);

/*! \brief Computes the approx. L^2 deviation of the gradients of a numerical solution from an analytic one.
 * 
 * Function which computes an approximation of the L^2 error $||\nabla u- \nabla u_{h}||_{L^2}$
 *  of a numerical solution using the analytic solution and Gauss quadrature
 * to approximate the integrals.
 * NOTE! It is vital that the quadrature is of sufficiently high order due to
 * the occurence of products of polynomials.
 * \param vals [in] Vector containing the computed values of the approximation.
 * \param geo [in] Struct of type geometry2 containing the mesh of the domain.
 * \param GD [in] Array of nodes and weights of Gauss quadrature.
 * \param QuadOrder [in] Number of nodes of the quadrature
 * \return Err L^2 error approximation.
 *  */
dtype L2DerivErr(const Vec& vals, const geometry2& geo, const dtype GD[][3], const uint QuadOrder);


dtype H1FuncErr(const Vec& vals, const geometry2& geo, const dtype GD[][3], const uint QuadOrder);


dtype L2NormRHS(const geometry2& geo, const dtype GD[][3], const uint QuadOrder);

void EvaluateU(const geometry2& geo, Vec& u);

void setUpRHS(const geometry2& geo, const dtype GD[][3], const uint QuadOrder, Vec& RHS);

//------------------------------------------------------------------------------

int main(int argc, char** argv) {
    std::string CoordFname("Koordinaten_ref.dat"), TriFname("Elemente_ref.dat"), BdFname("Boundary_ref.dat");

    params RefPar;
    //parameter input
    prog_opt_init(argc, argv, RefPar);

    //hard-coded parameters
    RefPar.BoundaryFilenameD = "Boundary_refRect.dat";
    RefPar.BoundaryFilenameN = "";
    RefPar.ElementsFilename = "Elements_refRect.dat";
    RefPar.CoordinateFilename = "Coord_refRect.dat";
    RefPar.tol = 1e-10;

    //maximum refinement order
    const uint maxRef(7);

    //additional constants
    dtype meshSize(0);
    dtype L2f(0), L2Err(0), H1Err(0);

    //read and store non-refined geometry
    geometry2 geo = readData2D(RefPar.CoordinateFilename, RefPar.ElementsFilename, RefPar.BoundaryFilenameD, RefPar.BoundaryFilenameN);
    geometry2 refGeo;
    refGeo = geo;


    std::cout << "%SUITE_STARTING% QualityCheck" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    //iterate over the refinements
    Eigen::SparseLU<CSRMat> solverLLT;
    for (uint ref(2); ref < maxRef; ++ref) {
        //1. refine mesh
        RefPar.Nref = ref;
        geo = refGeo;
        ref2d(geo, RefPar.Nref);
        meshSize = 1.0 / std::pow(2, RefPar.Nref);
        //2. insert appropriate basis nodes
        pKbasis(geo);
        //3. create the stiffness matrix
        uint freeNodes = geo.updateFIdx();
        Eigen::VectorXi NnzPerRow(freeNodes);
        DetermineSparsityPattern(geo, NnzPerRow);
        CSRMat A(freeNodes, freeNodes);
        A.reserve(NnzPerRow);
        StiffnessMatrixSetup(A, geo, Gauss, Order, freeNodes, NnzPerRow);
        A.prune(RefPar.tol, 1e-2);
        A.makeCompressed();
        CSRMat C(freeNodes, freeNodes);
        C.reserve(NnzPerRow);
        DampingMatrixSetup(C, geo, Gauss, Order, freeNodes, NnzPerRow);
        C.prune(RefPar.tol, 1e-3);
        C.makeCompressed();
        
        CSRMat B = -A+alpha*C;
        //4. create the RHS F.
        Vec F(freeNodes);
        setUpRHS(geo, Gauss, Order, F);
        //5. initialize the CG solver
        solverLLT.compute(B);
        //6. Solve Ax=F.
        Vec x(freeNodes);
        x = solverLLT.solve(F);

        //7. Compute the L^2 norm of the inhomogeneity f
        L2f = L2NormRHS(geo, Gauss, Order);
        //8. Using x from step 6 compute the L^2 error of the numerical
        // solution and its gradient as well as the H^1 error.
        L2Err = L2FuncErr(x, geo, Gauss, Order);
        H1Err = H1FuncErr(x, geo, Gauss, Order);
        //9. print
        std::cout << "Refinement : " << ref << std::endl;
        std::cout << "h :\t" << meshSize << std::endl;
        std::cout << "L2 :\t" << L2Err << std::endl;
        std::cout << "H1 :\t" << H1Err << std::endl;
#ifdef P1BASIS
        //if P1 basis is defined compute the L^infty error
        Vec uAnalytic(freeNodes);
        EvaluateU(geo, uAnalytic);
        //compute the deviation
        Eigen::ArrayXd error = uAnalytic.array() - x.array();
        dtype LinftyErr = error.abs().maxCoeff();
        std::cout << "L^inf:\t" << LinftyErr << std::endl;
#endif
        std::cout << "f :\t" << L2f << std::endl;
        std::cout << "\n";
    }


    std::cout << "%SUITE_FINISHED% time=0" << std::endl;
    return 0;
}

dtype f(const dtype& x, const dtype& y) {
    //test function for the analytic solution
    static dtype f = 0;
    f = 0;
    f = x*(6*y*(1.0-y)+alpha*pow(x,4)*(2*y-1.) + pow(x,3)*(2+alpha-2*alpha*y) + 3*x*(y-1.)*y*(alpha*y+4.)
            +x*x*(-2.+4*y*y*alpha-4*pow(y,3)*alpha));
    return f;
}

dtype u(const dtype& x, const dtype& y) {
    //analytic solution
    static dtype u = 0;
    u = 0;
    u = x * x * x * (x - 1.0) * y * (y - 1.0);
    return u;
}

void du(Eigen::Vector2d& v, const dtype& x, const dtype& y) {
    v.setZero();
    v(0) = y * (y - 1.)*(4 * x * x * x - 3 * x * x);
    v(1) = x * x * x * (x - 1.)*(2 * y - 1.);
}

//definition of the functions

void pKbasis(geometry2 & geo) {
    //call appropriate functions for P2, P3 elements and do nothing for P1
#ifdef P2BASIS
    p2basepoints(geo);
#elif P3BASIS
    p3basepoints(geo);
#else
    ;
#endif

}

dtype L2NormRHS(const geometry2& geo, const dtype GD[][3], const uint QuadOrder) {
    dtype NrmSq(0), integral(0), jac(0);

    Eigen::Vector2d tmp;
    Eigen::Matrix2d J;
    //Sum over all elements
    for (uint k(0); k < geo.Nelem; ++k) {
        Jacobian(geo.element[k], geo.vertices, J);
        jac = fabs(J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1));

        //compute the integral over the element
        integral = 0;
        for (uint s(0); s < QuadOrder; ++s) {
            //map the quadrature nodes from the std simplex to the triangle
            tmp(0) = GD[s][0];
            tmp(1) = GD[s][1];
            g(geo.element[k], geo.vertices, tmp);
            integral += GD[s][2] * SQU(f(tmp(0), tmp(1)));
        }
        integral *= jac;
        //add to the total
        NrmSq += integral;
    }
    //return the ||f||_L^2 = sqrt( ||f||^2 )
    return std::sqrt(NrmSq);
}

void setUpRHS(const geometry2& geo, const dtype GD[][3], const uint QuadOrder, Vec& RHS) {

    dtype integral(0), jac(0);
    uint FuncIdx1 = 0;

    Eigen::Vector2d tmp;
    Eigen::Matrix2d J;

    for (uint k(0); k < geo.Nelem; ++k) { //iteration over the elements(subintervals)
        //compute the volume of the element
        Jacobian(geo.element[k], geo.vertices, J);
        //abosulte value of the determinant
        jac = fabs(J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1));

        for (uint i(0); i < triVert; ++i) { //iteration over the boundaries

            //determine function index
            FuncIdx1 = geo.fIdx[ geo.element[k].vert[i] ];
            //INHOMOGENOUS DIRICHLET ADDITION
            if (geo.isInDirichletBoundary[ geo.element[k].vert[i] ] == true) {
                continue; //skip the rest
            }

            //integrate
            integral = 0.0;

            for (uint j(0); j < QuadOrder; ++j) {
                //map to the appropriate element
                tmp(0) = GD[j][0];
                tmp(1) = GD[j][1];
                g(geo.element[k], geo.vertices, tmp);
                integral += GD[j][2] * f(tmp(0), tmp(1)) * Phi(i, GD[j][0], GD[j][1]);
            }
            integral *= jac; //transform the measure (constant)
            RHS(FuncIdx1) += integral; //we skip first node

        }
    }
}

void EvaluateU(const geometry2& geo, Vec& uVec){

    uint FuncIdx(0);
    for(uint j(0); j<geo.Nvert; ++j){
        //skip dirichlet boundary
        if(geo.isInDirichletBoundary[j]==true)
            continue;
        FuncIdx = geo.fIdx[ j ];
        uVec(FuncIdx) = u(geo.vertices[j].x[0], geo.vertices[j].x[1]);
    }
    
}

dtype L2FuncErr(const Vec& vals, const geometry2& geo, const dtype GD[][3], const uint QuadOrder) {
    dtype NrmSq(0), integral(0), jac(0), arg(0);
    uint FuncIdx1 = 0;

    Eigen::Vector2d tmp;
    Eigen::Matrix2d J;
    //Sum over all elements
    for (uint k(0); k < geo.Nelem; ++k) {
        Jacobian(geo.element[k], geo.vertices, J);
        jac = fabs(J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1));

        //compute the integral over the element
        integral = 0;
        for (uint s(0); s < QuadOrder; ++s) {
            //map the quadrature nodes from the std simplex to the triangle
            tmp(0) = GD[s][0];
            tmp(1) = GD[s][1];
            g(geo.element[k], geo.vertices, tmp);

            //sum over the basis functions on the std-simplex
            arg = 0;
            for (uint j(0); j < triVert; ++j) {
                //skip dirichlet boundary
                if (geo.isInDirichletBoundary[ geo.element[k].vert[j] ] == true)
                    continue;
                FuncIdx1 = geo.fIdx[ geo.element[k].vert[j] ];
                arg += vals(FuncIdx1) * Phi(j, GD[s][0], GD[s][1]);
            }
            //switch the sign and add the analytic value
            arg -= u(tmp(0), tmp(1));


            integral += GD[s][2] * SQU(arg);
        }
        integral *= jac;
        //add to the total
        NrmSq += integral;
    }
    //return the ||f||_L^2 = sqrt( ||f||^2 )
    return std::sqrt(NrmSq);
}

dtype L2DerivErr(const Vec& vals, const geometry2& geo, const dtype GD[][3], const uint QuadOrder) {
    dtype NrmSq(0), integral(0), jac(0);
    uint FuncIdx1 = 0;

    Eigen::Vector2d tmp, dPhi1, dPhi2, arg, tmp2;
    Eigen::Matrix2d J, Jinv;
    //Sum over all elements
    for (uint k(0); k < geo.Nelem; ++k) {
        Jacobian(geo.element[k], geo.vertices, J);
        jac = fabs(J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1));
        InvJacobian(geo.element[k], geo.vertices, Jinv);

        //compute the integral over the element
        integral = 0;
        for (uint s(0); s < QuadOrder; ++s) {
            //map the quadrature nodes from the std simplex to the triangle
            tmp(0) = GD[s][0];
            tmp(1) = GD[s][1];
            g(geo.element[k], geo.vertices, tmp);

            //sum over the gradients of the basis functions
            arg.setZero();
            for (uint j(0); j < triVert; ++j) {
                //skip dirichlet boundary
                if (geo.isInDirichletBoundary[ geo.element[k].vert[j] ] == true)
                    continue;
                tmp2(0) = GD[s][0];
                tmp2(1) = GD[s][1];
                dPhi(j, tmp2, GD[s][0], GD[s][1]);
                dPhi1 = Jinv*tmp2;

                FuncIdx1 = geo.fIdx[ geo.element[k].vert[j] ];
                arg += vals(FuncIdx1) * dPhi1;
            }
            arg *= -1;
            du(tmp2, tmp(0), tmp(1));
            //dPhi2 = Jinv*tmp2;
            dPhi2 = tmp2;
            arg += dPhi2;

            integral += GD[s][2] * SQU(arg.norm());
        }
        integral *= jac;
        //add to the total
        NrmSq += integral;
    }
    //return the ||f||_L^2 = sqrt( ||f||^2 )
    return std::sqrt(NrmSq);
}

dtype H1FuncErr(const Vec& vals, const geometry2& geo, const dtype GD[][3], const uint QuadOrder) {
    dtype Err(0), DerivErr(0);
    Err = L2FuncErr(vals, geo, GD, QuadOrder);
    DerivErr = L2DerivErr(vals, geo, GD, QuadOrder);
    return std::sqrt(SQU(Err) + SQU(DerivErr));
}