/*
 * File:   main.cpp
 * Author: Anton Lebedev
 *
 * Created on 6. October 2015 11:15
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "Geometry.h"
#include "Functions.h"
#include "GaussQuadParam.h"

#include "../lib/Eigen_3.2.6/unsupported/Eigen/SparseExtra"
#include <Eigen/UmfPackSupport>
#include "arcomp.h"
#include "arrgsym.h"
#include "arrsnsym.h"

using namespace std;



//------------------------------------------------------------------------------

int main(int argc, char** argv) {

    params RefPar;
    //parameter input
    prog_opt_init(argc, argv, RefPar);

    //read-in
    geometry2 geo = readData2D(RefPar.CoordinateFilename, RefPar.ElementsFilename, RefPar.BoundaryFilenameD, RefPar.BoundaryFilenameN);

    //pull the flag for recomputation of the system matrices
    bool recompute(RefPar.recompute);

    //shift if required
    if ( (RefPar.xShft != 0) && (RefPar.yShft != 0) )
        shift(geo, RefPar.xShft, RefPar.yShft);

    //refine
    ref2d(geo, RefPar.Nref);
    //add midpoints for P_k FE.
#ifdef P2BASIS
    p2basepoints(geo);
#elif P3BASIS
    p3basepoints(geo);
#else
    ;
#endif

    //update function indices - compute number of free nodes
    uint freeNodes = geo.updateFIdx();


    //vector to store # of non-zeros per row of the stiffness matrix A
    Eigen::VectorXi NnzPerRow(freeNodes);
    DetermineSparsityPattern(geo, NnzPerRow);

    //--------------------------------------------------------------------------
    //matrix storage filename handling
    std::string Afile, Mfile, Cfile, tmpString;
    std::string InterpolationOrder;
    std::stringstream converter;

    //convert the # of refinement steps into a string
    converter << RefPar.Nref;
    converter >> tmpString;
    //define strings indicating the order of the interpolation polynomial
#ifdef P2BASIS
    InterpolationOrder = "P2";
#elif P3BASIS
    InterpolationOrder = "P3";
#else
    InterpolationOrder = "P1";
#endif
    Afile = RefPar.MatPrefix + "_" + InterpolationOrder + "_" + tmpString + "_stiffness.mm";
    Mfile = RefPar.MatPrefix + "_" + InterpolationOrder + "_" + tmpString + "_mass.mm";
    Cfile = RefPar.MatPrefix + "_" + InterpolationOrder + "_" + tmpString + "_damping.mm";

    std::ifstream InputFile;
    
    Eigen::initParallel();

    //------------------------------------------------------------------------------
    //---------------- Matrix creation ---------------------------------------------
    std::cout << "Creating matrix A" << std::endl;
    //create stiffness matrix
    CSRMat A(freeNodes, freeNodes);
    A.reserve(NnzPerRow);
    //iterator variable
    uint k(0);

    InputFile.open(Afile.c_str());
    if (!recompute && InputFile.good()) {
        //file exists and can be opened, now close it again
        InputFile.close();
        //read matrix
        Eigen::loadMarket<CSRMat>(A, Afile);
    }
    else {
        std::cout << "Filling matrix A" << std::endl;
        StiffnessMatrixSetup(A, geo, GD, QuadOrder, freeNodes, NnzPerRow);

        //store matrix
        Eigen::saveMarket<CSRMat>(A, Afile);
    }
    //create the mass matrix
    std::cout << "Creating matrix M" << std::endl;
    CSRMat M(freeNodes, freeNodes);
    M.reserve(NnzPerRow);


    InputFile.open(Mfile.c_str());
    if (!recompute && InputFile.good()) {
        //file exists and can be opened, now close it again
        InputFile.close();
        //read matrix
        Eigen::loadMarket<CSRMat>(M, Mfile);
    }
    else {
        std::cout << "Filling matrix M" << std::endl;
        MassMatrixSetup(M, geo, GD, QuadOrder, freeNodes, NnzPerRow);

        //store matrix
        Eigen::saveMarket<CSRMat>(M, Mfile);
    }

    //prune vanishing entries
    A.prune(RefPar.tol, 1e-3);
    M.prune(RefPar.tol, 1e-3);
    //compress the global sparse matrices
    A.finalize();
    M.finalize();
    A.makeCompressed();
    M.makeCompressed();
    
    //scale mass matrix by IOR^2
    M *= SQU(RefPar.n);

    //Use ARPACK to determine the eigenvalues
    dtype tol(RefPar.tol); //tolerance
    dtype sigma(0); //shift for later use
    uint Neig = RefPar.Neig;
    uint precision = static_cast<uint> (-std::log10(RefPar.tol)) + 1;

    //use ARPACK++
    Eigen::SparseMatrix<dtype, Eigen::ColMajor> Atilde = -A;
    Eigen::VectorXd w(freeNodes), z(freeNodes);
    Eigen::UmfPackLU< Eigen::SparseMatrix<dtype, Eigen::ColMajor> > Asolver;
    Eigen::SimplicialLDLT<CSRMat> Asolver2;
    Asolver.compute(Atilde);
    
    ARrcSymGenEig<dtype> prob('S', Atilde.cols(), Neig, 0.0);

    while (!prob.ArnoldiBasisFound()) {

        // Calling ARPACK FORTRAN code. Almost all work needed to
        // find an Arnoldi basis is performed by TakeStep.

        prob.TakeStep();

        switch (prob.GetIdo()) {
            case -1:

                // Performing w <- OP*B*v for the first time.
                // This product must be performed only if GetIdo is equal to
                // -1. GetVector supplies a pointer to the input vector, v,
                // and PutVector a pointer to the output vector, w.
                w = Eigen::Matrix<dtype, Eigen::Dynamic, 1>::Map(prob.GetVector(), w.rows());
                z = M*w;
                w = Asolver.solve(z);
                Eigen::Matrix<dtype, Eigen::Dynamic, 1>::Map(prob.PutVector(), w.rows()) = w;
                break;

            case 1:

                // Performing w <- OP*B*v when Bv is available.
                // This product must be performed whenever GetIdo is equal to
                // 1. GetProd supplies a pointer to the previously calculated
                // product Bv and PutVector a pointer to the output vector w.
                w = Eigen::Matrix<dtype, Eigen::Dynamic, 1>::Map(prob.GetProd(), w.rows());
                z = Asolver.solve(w);
                Eigen::Matrix<dtype, Eigen::Dynamic, 1>::Map(prob.PutVector(), w.rows()) = z;
                break;

            case 2:
                w = Eigen::Matrix<dtype, Eigen::Dynamic, 1>::Map(prob.GetVector(), w.rows());
                z = M*w;
                Eigen::Matrix<dtype, Eigen::Dynamic, 1>::Map(prob.PutVector(), w.rows()) = z;
        }
    }
    prob.FindEigenvalues();

    int nconvRe = prob.ConvergedEigenvalues();
    std::cout << "Found " << nconvRe << " eigenvalues: " << std::endl;
    for (int i(0); i < nconvRe; ++i) {
        std::cout << "  lambda[" << (i + 1) << "]: " << std::setprecision(precision) << prob.Eigenvalue(i)\
            << "\t" << std::setprecision(precision) << std::pow(std::abs(prob.Eigenvalue(i)), 2) << std::endl;

    }

    //create the transport matrix
    std::cout << "Creating matrix C" << std::endl;
    CSRMat C(freeNodes, freeNodes);
    C.reserve(NnzPerRow);

    InputFile.open(Cfile.c_str());
    if (!recompute && InputFile.good()) {
        //file exists and can be opened, now close it again
        InputFile.close();
        //read matrix
        Eigen::loadMarket<CSRMat>(C, Cfile);
    }
    else {
        std::cout << "Filling matrix C" << std::endl;
        DampingMatrixSetup(C, geo, GD, QuadOrder, freeNodes, NnzPerRow);
        Eigen::saveMarket<CSRMat>(C, Cfile);
    }
    //prune vanishing entries
    C.prune(RefPar.tol, 1e-3);
    //compress
    C.finalize();
    C.makeCompressed();


    //create the global matrices for the linearized eigenvalue problem
    CSRMat L(2 * freeNodes, 2 * freeNodes); //twice the size system matrix
    CSRMat R(2 * freeNodes, 2 * freeNodes); //RHS matrix
    //vector with total number of nonzeros
    Eigen::VectorXi NnzPerRow_R(2 * freeNodes);
    Eigen::VectorXi NnzPerRow_L(2 * freeNodes);


    std::cout << "Setting up the index vectors" << std::endl;
    //collect the number of non zeros per row
    Eigen::VectorXi NnzPerRowA(A.rows());
    Eigen::VectorXi NnzPerRowM(A.rows());
    Eigen::VectorXi NnzPerRowC(A.rows());

    std::cout << "Setting up the index vectors" << std::endl;
    //set up the vectors
    for (uint j(0); j < freeNodes; ++j) {
        NnzPerRow_R(j) = NnzPerRow_R(j + freeNodes) = NnzPerRow(j);
        NnzPerRow_L(j) = NnzPerRow(j) + NnzPerRow(j);
        NnzPerRow_L(j + freeNodes) = NnzPerRow(j);
    }

    dtype omega(0);
    dtype alpha(0);

    //set the value using user input
    alpha = -2 * omega / c;

    //read angular frequencies from file

    std::vector<dtype> OmegaList = readOmega(RefPar.omega);
    //reserve space
    R.reserve(NnzPerRow_R);
    L.reserve(NnzPerRow_R);


    //set matrices from triplets
    uint LNnz(0), RNnz(0);
    LNnz = C.nonZeros() + 2 * A.nonZeros();
    RNnz = A.nonZeros() + M.nonZeros();
    std::vector<Eigen::Triplet<dtype> > Ltriplets;
    std::vector<Eigen::Triplet<dtype> > Rtriplets;

    //create all required matrices
    Eigen::SparseMatrix<dtype, Eigen::ColMajor> Acsc(2 * freeNodes, 2 * freeNodes);
    Eigen::SparseMatrix<dtype, Eigen::ColMajor> Mcsc(2 * freeNodes, 2 * freeNodes);
    //define the matrices necessary for the Hamiltonian decomposition
    //due to the indended use of UMFPACK and ARPACK we will need them in CSC format
    Eigen::SparseMatrix<dtype, Eigen::ColMajor> B1(2 * freeNodes, 2 * freeNodes);
    Eigen::SparseMatrix<dtype, Eigen::ColMajor> B2(2 * freeNodes, 2 * freeNodes);
    Eigen::SparseMatrix<dtype, Eigen::ColMajor> D(2 * freeNodes, 2 * freeNodes);
    //vectors
    //support vectors
    z.resize(2*freeNodes);
    w.resize(2*freeNodes);
    Vec x(freeNodes), y(freeNodes), x2(freeNodes), y2(freeNodes);
    //Eigen::UmfPackLU< Eigen::SparseMatrix<dtype, Eigen::ColMajor> > Dsolver;
#ifdef P1BASIS
    //create the array to store eigenvectors
    cdtype **EV = new cdtype*[freeNodes];
    for(uint i(0); i< freeNodes; ++i)
        EV[i] = new cdtype[2*Neig];
    //zero it for good measure
    for(uint i(0); i< freeNodes; ++i)
        for(uint j(0); j<2*Neig; ++j)
            EV[i][j] = 0;
    //filename of the file storing the eigenvectors
    std::string EigVecFile;
    
    //flag for mesh writing
    bool writeMesh(true);
#endif

    //array to store the eigenvalues
    cdtype **eigvalList = new cdtype*[2 * Neig];
    for (uint j(0); j < 2 * Neig; ++j)
        eigvalList[j] = new cdtype[OmegaList.size()];

    std::cout << "Iterating over " << OmegaList.size() << " angular frequencies" << std::endl;
    
    //for the linearized equation we have to factorize A only once, then use it
//    Asolver2.analyzePattern(A);
//    Asolver2.factorize(A);
    
    uint eigvalIter(0);
    //iterate over the omega values
    for (std::vector<dtype>::iterator omegaIt = OmegaList.begin(); omegaIt != OmegaList.end(); omegaIt++) {
        std::cout << "Angular frequency " << *omegaIt << std::endl;
        alpha = -2.0 * (*omegaIt) / c;
        
        //set up auxiliary stiffness matrix 
        StiffnessMatrixSetup_Auxiliary(Atilde, geo, GD, QuadOrder, freeNodes, NnzPerRow);
        Atilde.prune(RefPar.tol, 1e-3);
        Atilde.makeCompressed();
        //set up the full mass matrix
        MassMatrixSetup_Full(M, geo, GD, QuadOrder, freeNodes, NnzPerRow, *omegaIt, RefPar.n);
        M.prune(RefPar.tol, 1e-3);
        M.makeCompressed();
        
        Atilde *= -SQU(*omegaIt/c);
        Atilde += A;
        //update the solver!
        if(eigvalIter == 0)
            Asolver2.analyzePattern(Atilde);
        
        Asolver2.factorize(Atilde);
        if(Asolver2.info() != Eigen::Success){
            std::cerr << " Failed to factorize the augmented stiffness matrix using Cholesky factorization..." <<std::endl;
            return 1;
        }

        Ltriplets.reserve(LNnz);
        Rtriplets.reserve(RNnz);
        //fill the triplets
        std::cout << "Filling triplets " << std::endl;
        uint s(0);
        //RHS
        //Stiffness part
        for (k = 0; k < Atilde.outerSize(); ++k)
            for (CSRMat::InnerIterator it(Atilde, k); it; ++it)
                Rtriplets.push_back(Eigen::Triplet<dtype>(it.row(), it.col() + freeNodes, -1.0 * it.value()));

        //mass part
        for (k = 0; k < M.outerSize(); ++k)
            for (CSRMat::InnerIterator it(M, k); it; ++it)
                Rtriplets.push_back(Eigen::Triplet<dtype>(it.row() + freeNodes, it.col(), it.value()));


        s = 0;
        //LHS
        //transport part
        for (k = 0; k < C.outerSize(); ++k)
            for (CSRMat::InnerIterator it(C, k); it; ++it)
                Ltriplets.push_back(Eigen::Triplet<dtype>(it.row(), it.col() + freeNodes, alpha * it.value()));
        //mass part
        for (k = 0; k < M.outerSize(); ++k)
            for (CSRMat::InnerIterator it(M, k); it; ++it) {
                Ltriplets.push_back(Eigen::Triplet<dtype>(it.row(), it.col(), it.value()));
                Ltriplets.push_back(Eigen::Triplet<dtype>(it.row() + freeNodes, it.col() + freeNodes, it.value()));
            }

        //shrink triplets to size
        Rtriplets.shrink_to_fit();
        Ltriplets.shrink_to_fit();
        //set matrices from triplets
        std::cout << "Filling matrices" << std::endl;
        R.setFromTriplets(Rtriplets.begin(), Rtriplets.end());
        L.setFromTriplets(Ltriplets.begin(), Ltriplets.end());

        std::cout << "Done setting up matrices " << std::endl;
        L.finalize();
        R.finalize();
        L.makeCompressed();
        R.makeCompressed();


        //clear triplets
        Ltriplets.clear();
        Rtriplets.clear();

        //SA views of matrix A and M
        Acsc = R;
        Mcsc = L;
        Acsc.makeCompressed();
        Mcsc.makeCompressed();

        /* -------------------------------------------------------------------------------------------------
         Implementation of the structure-preserving shift-and-invert method using the callback interface of 
         * ARPACK and the shift proposed in Tisseur and Meerbergen "The Quadratic Eigenvalue Problem" (SIAM 2001)
         ---------------------------------------------------------------------------------------------------*/

        //now set up the matrices B1 and B2
        Ltriplets.reserve(LNnz);
        Rtriplets.reserve(LNnz);

        std::cout << "Setting up auxiliary matrices for the Hamiltonian problem" << std::endl;
        std::cout << "Filling triplets" << std::endl;
        //RHS + LHS due to symmetry
        s = 0;
        //mass part
        for (k = 0; k < M.outerSize(); ++k)
            for (CSRMat::InnerIterator it(M, k); it; ++it) {
                //B1
                Ltriplets.push_back(Eigen::Triplet<dtype>(it.row() + freeNodes, it.col() + freeNodes, it.value()));
                //B2
                Rtriplets.push_back(Eigen::Triplet<dtype>(it.row(), it.col(), it.value()));
                ++s;
            }

        //transport part
        for (k = 0; k < C.outerSize(); ++k)
            for (CSRMat::InnerIterator it(C, k); it; ++it) {
                //B1
                Ltriplets.push_back(Eigen::Triplet<dtype>(it.row(), it.col() + freeNodes, 0.5 * alpha * it.value()));
                //B2
                Rtriplets.push_back(Eigen::Triplet<dtype>(it.row(), it.col() + freeNodes, 0.5 * alpha * it.value()));
                ++s;
            }


        //fill in the identity matrices
        for (k = 0; k < freeNodes; ++k) {
            //B1
            Ltriplets.push_back(Eigen::Triplet<dtype>(k, k, 1.0));
            //B2
            Rtriplets.push_back(Eigen::Triplet<dtype>(k + freeNodes, k + freeNodes, 1.0));
            ++s;
        }
        //we've overestimated the number of non-vanishing elements for
        // the B1 and B2 matrices, shrink to fit
        Ltriplets.shrink_to_fit();
        Rtriplets.shrink_to_fit();
        //set matrices
        std::cout << "Setting up matrices" << std::endl;
        B1.setFromTriplets(Ltriplets.begin(), Ltriplets.end());
        B2.setFromTriplets(Rtriplets.begin(), Rtriplets.end());

        //clear the triplets
        Ltriplets.clear();
        Rtriplets.clear();

        //finalize matrices
        B1.finalize();
        B2.finalize();
        D.finalize();
        B1.makeCompressed();
        B2.makeCompressed();



        //initialize ARPACK callback interface
        std::cout << "Initializing ARPACK callback interface" << std::endl;
        ARrcNonSymStdEig<dtype> HamiltonianProb(Acsc.cols(), 2 * Neig, sigma);
        HamiltonianProb.ChangeTol(tol);

        //search for eigenvalues
        std::cout << "Constructing subspaces..." << std::endl;
        while (!HamiltonianProb.ArnoldiBasisFound()) {
            //take a step
            HamiltonianProb.TakeStep();

            if ((HamiltonianProb.GetIdo() == 1) || (HamiltonianProb.GetIdo() == -1)) {
                //Apply the  shift operator
                //Application of the shift operator following [MW] in a simplified form
                //fetch the vector and multiply it with B1
                x.setZero();
                y.setZero();
                w = B1 * Eigen::Matrix<dtype, Eigen::Dynamic, 1>::Map(HamiltonianProb.GetVector(), w.rows());
                //split the vector in x and y
                for(k=0; k<freeNodes; ++k){
                    x(k) = w(k);
                    y(k) = w(k+freeNodes);
                }
                //step1. x_2 <- A^-1 *x, x_2<- -x_2
                x2 = Asolver2.solve(x);
                x2 *= -1;
                if(Asolver2.info() != Eigen::Success){
                    std::cerr << " Failed to solve the linear system using Cholesky factorization of the augmented stiffness matrix..." << std::endl;
                    return 1;
                }
                //step 2. y_2 = alpha*C*x2;
                y2 = C*x2;
                y2 *= alpha;
                //step 3. y <- y+y_2
                y += y2;
                //step 4. y_2 = A^-1 *y, y_2 <- -y_2;
                y2 = Asolver2.solve(y);
                y2 *= -1;
                if(Asolver2.info() != Eigen::Success){
                    std::cerr << " Failed to solve the linear system using Cholesky factorization of the augmented stiffness matrix..." << std::endl;
                    return 1;
                }
                //join x_2, y_2 into w
                for(k=0; k<freeNodes; ++k){
                   z(k) = x2(k);
                   z(k+freeNodes) = y2(k);
                }
                //multiply w by B_2 and store into the ARPACK designated area
                Eigen::Matrix<dtype, Eigen::Dynamic, 1>::Map(HamiltonianProb.PutVector(), z.rows()) = B2*z;
            }
        }
        std::cout << "Searching for eigenvalues..." << std::endl;
        //HamiltonianProb.FindEigenvalues();

        int nconvEV = HamiltonianProb.FindEigenvectors();
        int nconv = HamiltonianProb.ConvergedEigenvalues();
        std::cout << "Found " << nconv << " eigenvalues: " << std::endl;
        for (int i(0); i < nconv; ++i) {
            std::cout << "  lambda[" << (i + 1) << "]: " << std::setprecision(precision) << HamiltonianProb.Eigenvalue(i)\
            << "\t" << std::setprecision(precision) << std::pow(std::abs(HamiltonianProb.Eigenvalue(i)), 2) << std::endl;
            eigvalList[i][eigvalIter] = HamiltonianProb.Eigenvalue(i);
        }
#ifdef P1BASIS       
        //determine the eigenvectors every 5th iteration
        if(eigvalIter% 5 == 0){
            
            //int nconvEV = HamiltonianProb.FindEigenvectors();
            std::cout << "Found " << nconvEV << " eigenvectors, storing..." <<std::endl;
            for(int i(0); i < nconvEV; ++i){
                for(uint j(0); j<freeNodes; j++){
                    //the unscaled eigenvectors of the QEVP are stored 
                    //in the second part of the eigenvectors of the linearized problem
                    EV[j][i] = HamiltonianProb.Eigenvector(i, j+freeNodes);
                }
            }
            

            std::cout << "Writing eigenvectors for omega " << *omegaIt << std::endl;
            //EigVecFile = "Eigenmodes_omega_" + tmpString + ".dat";
            //std::cout<< "Filename "<< EigVecFile <<std::endl;
            //write as many eigenvectors as have been determined
            //here: nconv
            saveEV(RefPar.MatPrefix, *omegaIt, nconvEV, geo, EV, writeMesh, freeNodes, GD, QuadOrder);
            //set the writeMesh flag to false after first iteration
            writeMesh = false;
        }
#endif
        eigvalIter++;
    }

    std::cout << "Shift: " << sigma << std::endl;
    //save
    saveToFile(2 * Neig, OmegaList, eigvalList);
    //clean up
#ifdef P1BASIS
    for(uint i(0); i< freeNodes; ++i)
        delete[] EV[i];

    delete[] EV;
#endif
    for (uint j(0); j < 2 * Neig; ++j)
        delete[] eigvalList[j];
    delete[] eigvalList;

    return 0;
}
