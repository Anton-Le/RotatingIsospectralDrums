/*! \file P2_MeshTest.cpp
 * \author Anton Lebedev
 * \copyright GNU General Public License 3.0 (GPLv3)
 * \brief Unit tests for the quadratic basis.
 * 
 * File:   P2_MeshTest.cpp
 * Author: lebedev
 *
 * Created on 18.01.2016, 14:18:23
 */

#include <stdlib.h>
#include <iostream>
#include "Geometry.h"
#include "Functions.h"
#include "GaussQuadParam.h"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
/*
 * Simple C++ Test Suite
 */
//const uint triVert(6);

void createReferenceMesh(geometry2& refgeo) {
    /*!\brief Create the reference mesh for P2 elements.
     * 
     * Computed by hand on [7.4] p.3.
     * No. Elements : 2
     * No. Boundary elements: 8
     * No. Vertices: 9
     * Coordinates:
     * (0,0), (1,0), (1,1), (0,1), (1,1/2), (1/2,1/2), (1/2,0), (1/2,1), (0,1/2)
     * Elements:
     * 0: 1 2 0 4 5 6
     * 1: 2 3 0 7 8 5
     * Boundary:
     * (0,6), (1,4), (2,7), (3,8), (6,1), (4,2), (7,3), (8,0)
     *  */
    refgeo.resize(9, 2, 0, 8);

    //set the coordinates
    //0
    refgeo.vertices[0].x[0] = 0;
    refgeo.vertices[0].x[1] = 0;
    //1
    refgeo.vertices[1].x[0] = 1;
    refgeo.vertices[1].x[1] = 0;
    //2
    refgeo.vertices[2].x[0] = 1;
    refgeo.vertices[2].x[1] = 1;
    //3
    refgeo.vertices[3].x[0] = 0;
    refgeo.vertices[3].x[1] = 1;
    //4
    refgeo.vertices[4].x[0] = 1;
    refgeo.vertices[4].x[1] = 0.5;
    //5
    refgeo.vertices[5].x[0] = 0.5;
    refgeo.vertices[5].x[1] = 0.5;
    //6
    refgeo.vertices[6].x[0] = 0.5;
    refgeo.vertices[6].x[1] = 0;
    //7
    refgeo.vertices[7].x[0] = 0.5;
    refgeo.vertices[7].x[1] = 1;
    //8
    refgeo.vertices[8].x[0] = 0;
    refgeo.vertices[8].x[1] = 0.5;

    //set the elements
    //0
    refgeo.element[0].vert[0] = 1;
    refgeo.element[0].vert[1] = 2;
    refgeo.element[0].vert[2] = 0;
    refgeo.element[0].vert[3] = 4;
    refgeo.element[0].vert[4] = 5;
    refgeo.element[0].vert[5] = 6;
    //1
    refgeo.element[1].vert[0] = 2;
    refgeo.element[1].vert[1] = 3;
    refgeo.element[1].vert[2] = 0;
    refgeo.element[1].vert[3] = 7;
    refgeo.element[1].vert[4] = 8;
    refgeo.element[1].vert[5] = 5;

    //set the boundary
    //0
    refgeo.boundaryNeumann[0].vert[0] = 0;
    refgeo.boundaryNeumann[0].vert[1] = 6;
    //1
    refgeo.boundaryNeumann[1].vert[0] = 1;
    refgeo.boundaryNeumann[1].vert[1] = 4;
    //2
    refgeo.boundaryNeumann[2].vert[0] = 2;
    refgeo.boundaryNeumann[2].vert[1] = 7;
    //3
    refgeo.boundaryNeumann[3].vert[0] = 3;
    refgeo.boundaryNeumann[3].vert[1] = 8;
    //4
    refgeo.boundaryNeumann[4].vert[0] = 6;
    refgeo.boundaryNeumann[4].vert[1] = 1;
    //5
    refgeo.boundaryNeumann[5].vert[0] = 4;
    refgeo.boundaryNeumann[5].vert[1] = 2;
    //6
    refgeo.boundaryNeumann[6].vert[0] = 7;
    refgeo.boundaryNeumann[6].vert[1] = 3;
    //7
    refgeo.boundaryNeumann[7].vert[0] = 8;
    refgeo.boundaryNeumann[7].vert[1] = 0;
}

void testCopy(const geometry2& in, geometry2& out) {
    out = in;

    for (uint j(0); j < in.Nvert; ++j) {
        //compare vertices
        if ((out.vertices[j].x[0] != in.vertices[j].x[0]) || (out.vertices[j].x[1] != in.vertices[j].x[1])) {
            std::cout << "%TEST_FAILED% time=0 testname=testCopy (P2_MeshTest) message=Vertices not identical" << std::endl;
            return;
        }
        //compare function indices
        if (out.fIdx[j] != in.fIdx[j]) {
            std::cout << "%TEST_FAILED% time=0 testname=testCopy (P2_MeshTest) message=function indices differ" << std::endl;
            return;
        }
        //compare entries of dirichlet boundary
        if (out.isInDirichletBoundary[j] != in.isInDirichletBoundary[j]) {
            std::cout << "%TEST_FAILED% time=0 testname=testCopy (P2_MeshTest) message=Entries in the Dirichlet bd. masking list differ." << std::endl;
            return;
        }
    }
    //compare Elements
    for (uint j(0); j < in.Nelem; ++j) {
        for (uint k(0); k < triVert; ++k) {
            if (out.element[j].vert[k] != in.element[j].vert[k]) {
                std::cout << "%TEST_FAILED% time=0 testname=testCopy (P2_MeshTest) message=Elements differ" << std::endl;
                return;
            }
        }
    }
    //compare boundary elements
    for (uint j(0); j < in.NBdElemD; ++j) {
        if ((out.boundaryDirichlet[j].vert[0] != in.boundaryDirichlet[j].vert[0]) || (out.boundaryDirichlet[j].vert[1] != in.boundaryDirichlet[j].vert[1])) {
            std::cout << "%TEST_FAILED% time=0 testname=testCopy (P2_MeshTest) message=Dirichlet boundary elements differ" << std::endl;
            return;
        }
    }
    for (uint j(0); j < in.NBdElemN; ++j) {
        if ((out.boundaryNeumann[j].vert[0] != in.boundaryNeumann[j].vert[0]) || (out.boundaryNeumann[j].vert[1] != in.boundaryNeumann[j].vert[1])) {
            std::cout << "%TEST_FAILED% time=0 testname=testCopy (P2_MeshTest) message=Neumann boundary elements differ" << std::endl;
            return;
        }
    }

}

void testP2Refinement(geometry2& refgeo, geometry2& p2geo, params& Pars) {
    //refine the reference to a degree Nref+1
    uint referenceRefinement = Pars.Nref + 1;

    ref2d(refgeo, referenceRefinement);

    //refine the P2 mesh
    ref2d(p2geo, Pars.Nref);
    //add P_2 coordinates
    p2basepoints(p2geo);

    //compare the number of vertices
    if (p2geo.Nvert != refgeo.Nvert) {
        std::cout << "%TEST_FAILED% time=0 testname=testP2Refinement (P2_MeshTest) message=Number of vertices differs!" << std::endl;
        return;
    }
    //both meshes have the same number of vertices.
    //now check if the vertices themselves are identical
    for (uint j(0); j < refgeo.Nvert; ++j) {
        if ((p2geo.vertices[j].x[0] != refgeo.vertices[j].x[0]) || (p2geo.vertices[j].x[1] != refgeo.vertices[j].x[1])) {
            std::cout << "%TEST_FAILED% time=0 testname=testP2Refinement (P2_MeshTest) message=Vertices not identical" << std::endl;
            return;
        }
    }
    //check if the number of free nodes differs
    uint freeNodesSimplical(0), freeNodesP2(0);
    freeNodesSimplical = refgeo.updateFIdx();
    freeNodesP2 = p2geo.updateFIdx();
    if (freeNodesP2 != freeNodesSimplical) {
        std::cout << "%TEST_FAILED% time=0 testname=testP2Refinement (P2_MeshTest) message=Number of non-bound nodes differs" << std::endl;
        return;
    }
}

void testP2Refinement2(geometry2& refgeo, geometry2& geo) {
    std::cout << "Refining mesh" << std::endl;
    //add P_2 coordinates
    p2basepoints(geo);
    //compare the number of vertices
    if (geo.Nvert != refgeo.Nvert) {
        std::cout << "%TEST_FAILED% time=0 testname=testP2Refinement (P2_MeshTest) message=Number of vertices differs!" << std::endl;
        return;
    }
    //both meshes have the same number of vertices.
    //now check if the vertices themselves are identical
    for (uint j(0); j < refgeo.Nvert; ++j) {
        if ((geo.vertices[j].x[0] != refgeo.vertices[j].x[0]) || (geo.vertices[j].x[1] != refgeo.vertices[j].x[1])) {
            std::cout << "%TEST_FAILED% time=0 testname=testP2Refinement (P2_MeshTest) message=Vertices not identical" << std::endl;
            return;
        }
    }

    //compare the number of boundary elements
    if (geo.NBdElemD != refgeo.NBdElemD) {
        std::cout << "%TEST_FAILED% time=0 testname=testP2Refinement (P2_MeshTest) message=Number of Dirichlet boundary elements differs!" << std::endl;
        return;
    }
    //compare the number of boundary elements
    if (geo.NBdElemN != refgeo.NBdElemN) {
        std::cout << "%TEST_FAILED% time=0 testname=testP2Refinement (P2_MeshTest) message=Number of Neumann boundary elements differs!" << std::endl;
        return;
    }

    //check if the Neumann boundary differs
    for (uint j(0); j < refgeo.NBdElemN; ++j) {
        if ((geo.boundaryNeumann[j].vert[0] != refgeo.boundaryNeumann[j].vert[0]) || (geo.boundaryNeumann[j].vert[1] != refgeo.boundaryNeumann[j].vert[1])) {
            std::cout << "%TEST_FAILED% time=0 testname=testP2Refinement (P2_MeshTest) message=Boundary not identical" << std::endl;
            return;
        }
    }
    //check if the number of free nodes differs
    uint freeNodesSimplical(0), freeNodesP2(0);
    freeNodesSimplical = refgeo.updateFIdx();
    freeNodesP2 = geo.updateFIdx();
    if (freeNodesP2 != freeNodesSimplical) {
        std::cout << "%TEST_FAILED% time=0 testname=testP2Refinement (P2_MeshTest) message=Number of non-bound nodes differs" << std::endl;
        return;
    }
}

void createReferenceMassMatrix(DensMat& M);

void createReferenceStiffnessMatrix(DensMat& A);

void createReferenceDampingMatrix(DensMat& C);

void createReferenceMassMatrixRefined(DensMat& M);

void testMassMatrix(const DensMat& Mref, const CSRMat& M);

void testDampingMatrix(const DensMat& Mref, const CSRMat& M);

void testStiffnessMatrix(const DensMat& Aref, const CSRMat& A);

int main(int argc, char** argv) {
    std::string CoordFname("Koordinaten_ref.dat"), TriFname("Elemente_ref.dat"), BdFname("Boundary_ref.dat");

    params RefPar;
    //parameter input
    prog_opt_init(argc, argv, RefPar);
    //hard-coding the test parameters
    RefPar.BoundaryFilenameN = "Boundary_refRect.dat";
    RefPar.BoundaryFilenameD = "";
    RefPar.ElementsFilename = "Elements_refRect.dat";
    RefPar.CoordinateFilename = "Coord_refRect.dat";
    RefPar.Nref = 1;
    RefPar.n = 0;


    //read-in
    geometry2 geo = readData2D(RefPar.CoordinateFilename, RefPar.ElementsFilename, RefPar.BoundaryFilenameD, RefPar.BoundaryFilenameN);
    geometry2 geoP2;
    //create the reference P2 mesh
    geometry2 refgeoP2;
    createReferenceMesh(refgeoP2);

    //store original geometry, geometry2 struct has no copy-constructor!
    geometry2 geoOrig;
    geoOrig = geo;


    std::cout << "%SUITE_STARTING% P2_MeshTest" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% testCopy (P2_MeshTest)" << std::endl;
    testCopy(geo, geoP2);
    std::cout << "%TEST_FINISHED% time=0 testCopy (P2_MeshTest)" << std::endl;

    //test the refinement 
    std::cout << "%TEST_STARTED% testP2Refinement (P2_MeshTest)" << std::endl;
    testP2Refinement(geo, geoP2, RefPar);
    std::cout << "%TEST_FINISHED% time=0 testP2Refinement (P2_MeshTest)" << std::endl;

    geoP2 = geoOrig;

    std::cout << "%TEST_STARTED% testP2Refinement (P2_MeshTest)" << std::endl;
    testP2Refinement2(refgeoP2, geoP2);
    std::cout << "%TEST_FINISHED% time=0 testP2Refinement (P2_MeshTest)" << std::endl;

    //create reference matrix
    DensMat Mref(9, 9);
    createReferenceMassMatrix(Mref);
    //create the mass matrix
    uint freeNodes = geoP2.updateFIdx();
    Eigen::VectorXi NnzPerRow(freeNodes);
    DetermineSparsityPattern(geoP2, NnzPerRow);
    
    CSRMat M(freeNodes, freeNodes);
    M.reserve(NnzPerRow);
    MassMatrixSetup(M, geoP2, GD, QuadOrder, freeNodes, NnzPerRow);
    //test
    std::cout << "%TEST_STARTED% testMassMatrix (P2_MeshTest)" << std::endl;
    testMassMatrix(Mref, M);
    std::cout << "%TEST_FINISHED% time=0 testMassMatrix (P2_MeshTest)" << std::endl;

    //damping matrix test
    DensMat Cref(9,9);
    createReferenceDampingMatrix(Cref);
    CSRMat C(freeNodes, freeNodes);
    C.reserve(NnzPerRow);
    DampingMatrixSetup(C, geoP2, GD, QuadOrder, freeNodes, NnzPerRow);
    //test
    std::cout << "%TEST_STARTED% testDampingMatrix (P2_MeshTest)" << std::endl;
    testDampingMatrix(Cref, C);
    std::cout << "%TEST_FINISHED% time=0 testDampingMatrix (P2_MeshTest)" << std::endl;

    DensMat Aref(9,9);
    createReferenceStiffnessMatrix(Aref);
    CSRMat A(freeNodes, freeNodes);
    A.reserve(NnzPerRow);
    StiffnessMatrixSetup(A, geoP2, GD, QuadOrder, freeNodes, NnzPerRow);
    std::cout << "%TEST_STARTED% testStiffnessMatrix (P2_MeshTest)" << std::endl;
    testStiffnessMatrix(Aref, A);
    std::cout << "%TEST_FINISHED% time=0 testStiffnessMatrix (P2_MeshTest)" << std::endl;

    //test P2 basis on a refined mesh
    
    //reset the mesh
    geo = geoOrig;
    //refine
    ref2d(geo, RefPar.Nref);
    p2basepoints(geo);
    //reset the sparse matrix
    freeNodes = geo.updateFIdx();
    NnzPerRow.resize(freeNodes);
    DetermineSparsityPattern(geo, NnzPerRow);
    M.resize(freeNodes, freeNodes);
    M.reserve(NnzPerRow);
    //set up matrix
    MassMatrixSetup(M, geo, GD, QuadOrder, freeNodes, NnzPerRow);
    
    //set up the reference
    createReferenceMassMatrixRefined(Mref);
    std::cout << "%TEST_STARTED% testMassMatrixRefined (P2_MeshTest)" << std::endl;
    testMassMatrix(Mref, M);
    std::cout << "%TEST_FINISHED% time=0 testMassMatrixRefined (P2_MeshTest)" << std::endl;
    
    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

// definition of functions

void createReferenceMassMatrix(DensMat& M) {
    M.resize(9, 9);
    M.setZero();
    M << 1. / 30, -1. / 360, -1. / 180, -1. / 360, -1. / 90, 0, 0, -1. / 90, 0,
            -1. / 360, 1. / 60, -1. / 360, 0, 0, -1. / 90, 0, 0, 0,
            -1. / 180, -1. / 360, 1. / 30, -1. / 360, 0, 0, -1. / 90, 0, -1. / 90,
            -1. / 360, 0, -1. / 360, 1. / 60, 0, -1. / 90, 0, 0, 0,
            -1. / 90, 0, 0, 0, 4. / 45, 2. / 45, 2. / 45, 0, 0,
            0, -1. / 90, 0, -1. / 90, 2. / 45, 8. / 45, 2. / 45, 2. / 45, 2. / 45,
            0, 0, -1. / 90, 0, 2. / 45, 2. / 45, 4. / 45, 0, 0,
            -1. / 90, 0, 0, 0, 0, 2. / 45, 0, 4. / 45, 2. / 45,
            0, 0, -1. / 90, 0, 0, 2. / 45, 0, 2. / 45, 4. / 45;
}

void createReferenceStiffnessMatrix(DensMat& A) {
    A.resize(9, 9);
    A.setZero();
    A << 1, 1. / 6, 0, 1. / 6, 0, 0, -2. / 3, 0, -2. / 3,
            1. / 6, 1, 1. / 6, 0, -2. / 3, 0, -2. / 3, 0, 0,
            0, 1. / 6, 1, 1. / 6, -2. / 3, 0, 0, -2. / 3, 0,
            1. / 6, 0, 1. / 6, 1, 0, 0, 0, -2. / 3, -2. / 3,
            0, -2./3, -2./3, 0, 8./3, -4./3, 0, 0, 0,
            0, 0, 0, 0, -4. / 3, 16. / 3, -4. / 3, -4. / 3, -4. / 3,
            -2. / 3, -2. / 3, 0, 0, 0, -4. / 3, 8. / 3, 0, 0,
            0, 0, -2. / 3, -2. / 3, 0, -4. / 3, 0, 8. / 3, 0,
            -2. / 3, 0, 0, -2. / 3, 0, -4. / 3, 0, 0, 8. / 3;
}

void createReferenceDampingMatrix(DensMat& C) {
    C.resize(9, 9);
    C.setZero();
    C << 0, 7. / 360, 0, -7. / 360, 1. / 45, 0, -1. / 90, -1. / 45, 1. / 90,
            -7. / 360, 0, -11. / 360, 0, 1. / 9, -1. / 30, 1. / 45, 0, 0,
            0, 11. / 360, 0, -17. / 360, -7. / 45, 0, 1. / 90, 7. / 45, -1. / 90,
            7. / 360, 0, 17. / 360, 0, 0, 1. / 30, 0, -1. / 9, -1. / 45,
            -1. / 45, -1. / 9, 7. / 45, 0, 0, 2. / 9, -2. / 45, 0, 0,
            0, 1. / 30, 0, -1. / 30, -2. / 9, 0, -8. / 45, 2. / 9, 8. / 45,
            1. / 90, -1. / 45, -1. / 90, 0, 2. / 45, 8. / 45, 0, 0, 0,
            1. / 45, 0, -7. / 45, 1. / 9, 0, -2. / 9, 0, 0, 2. / 45,
            -1. / 90, 0, 1. / 90, 1. / 45, 0, -8. / 45, 0, -2. / 45, 0;
}

//testing functions

void testMassMatrix(const DensMat& Mref, const CSRMat& M) {

    //check whether the dimensions are identical
    if ((Mref.rows() != M.rows()) || (Mref.cols() != M.cols())) {
        std::cout << " %TEST_FAILED% time=0 testname=testMassMatrix (P2_MeshTest) message=Matrix dimensions differ!" << std::endl;
        return;
    }
    //if matrix dimensions agree we compute the difference between the two
    DensMat Mdens(Mref.rows(), Mref.cols());
    Mdens = M.toDense();
    Mdens -= Mref;

    //check the deviation
    dtype tol(1e-10);
    if ((std::fabs(Mdens.maxCoeff()) > tol) || (std::fabs(Mdens.minCoeff()) > tol)) {
        std::cout << " %TEST_FAILED% time=0 testname=testMassMatrix (P2_MeshTest) message=Matrix elements differ!" << std::endl;
        std::cout << "Reference matrix \n" << Mref << std::endl;
        std::cout << "Computed matrix \n" << M.toDense() << std::endl;
        std::cout << "Difference \n" << Mdens << std::endl;
        return;
    }

}

void testDampingMatrix(const DensMat& Cref, const CSRMat& C) {
    //check whether the dimensions are identical
    if ((Cref.rows() != C.rows()) || (Cref.cols() != C.cols())) {
        std::cout << " %TEST_FAILED% time=0 testname=testDampingMatrix (P2_MeshTest) message=Matrix dimensions differ!" << std::endl;
        return;
    }
    //if matrix dimensions agree we compute the difference between the two
    DensMat Cdens(Cref.rows(), Cref.cols());
    Cdens = C.toDense();
    Cdens -= Cref;

    //check the deviation
    dtype tol(1e-10);
    if ((std::fabs(Cdens.maxCoeff()) > tol) || (std::fabs(Cdens.minCoeff()) > tol)) {
        std::cout << " %TEST_FAILED% time=0 testname=testDampingMatrix (P2_MeshTest) message=Matrix elements differ!" << std::endl;
        std::cout << "Reference matrix \n" << Cref << std::endl;
        std::cout << "Computed matrix \n" << C.toDense() << std::endl;
        std::cout << "Difference \n" << Cdens << std::endl;
        return;
    }
}

void testStiffnessMatrix(const DensMat& Aref, const CSRMat& A) {

    //check whether the dimensions are identical
    if ((Aref.rows() != A.rows()) || (Aref.cols() != A.cols())) {
        std::cout << " %TEST_FAILED% time=0 testname=testStiffnessMatrix (P2_MeshTest) message=Matrix dimensions differ!" << std::endl;
        return;
    }
    //if matrix dimensions agree we compute the difference between the two
    DensMat Adens(Aref.rows(), Aref.cols());
    Adens = A.toDense();
    Adens -= Aref;

    //check the deviation
    dtype tol(1e-10);
    if ((std::fabs(Adens.maxCoeff()) > tol) || (std::fabs(Adens.minCoeff()) > tol)) {
        std::cout << " %TEST_FAILED% time=0 testname=testStiffnessMatrix (P2_MeshTest) message=Matrix elements differ!" << std::endl;
        std::cout << "Reference matrix \n" << Aref << std::endl;
        std::cout << "Computed matrix \n" << A.toDense() << std::endl;
        std::cout << "Difference \n" << Adens << std::endl;
        return;
    }
}

void createReferenceMassMatrixRefined(DensMat& M){
    M.resize(25,25);
    M.setZero();
    M << 1./120,0,0,0,0,-1./720,-1./1440,0,-1./1440,0,-1./360,0,0,-1./360,0,0,0,0,0,0,0,0,0,0,0,
0,1./240,0,0,-1./1440,0,-1./1440,0,0,0,0,-1./360,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,1./120,0,-1./1440,-1./720,0,-1./1440,0,-1./360,0,0,0,0,-1./360,0,0,0,0,0,0,0,0,0,0,
0,0,0,1./240,0,0,0,-1./1440,-1./1440,0,0,0,-1./360,0,0,0,0,0,0,0,0,0,0,0,0,
0,-1./1440,-1./1440,0,1./80,-1./720,-1./720,0,0,0,-1./360,0,0,0,0,0,-1./360,0,-1./360,0,0,0,0,0,0,
-1./720,0,-1./720,0,-1./720,1./40,-1./720,-1./720,-1./720,0,0,-1./360,-1./360,0,0,0,0,-1./360,0,0,-1./360,-1./360,0,0,-1./360,
-1./1440,-1./1440,0,0,-1./720,-1./720,1./80,0,0,-1./360,0,0,0,0,0,-1./360,0,0,0,-1./360,0,0,0,0,0,
0,0,-1./1440,-1./1440,0,-1./720,0,1./80,-1./720,0,0,0,0,-1./360,0,0,0,0,-1./360,0,0,0,0,-1./360,0,
-1./1440,0,0,-1./1440,0,-1./720,0,-1./720,1./80,0,0,0,0,0,-1./360,0,0,0,0,-1./360,0,0,-1./360,0,0,
0,0,-1./360,0,0,0,-1./360,0,0,2./45,1./90,1./90,0,0,0,0,0,1./90,1./90,0,0,0,0,0,0,
-1./360,0,0,0,-1./360,0,0,0,0,1./90,2./45,1./90,0,0,0,0,0,0,0,1./90,1./90,0,0,0,0,
0,-1./360,0,0,0,-1./360,0,0,0,1./90,1./90,2./45,0,0,0,1./90,1./90,0,0,0,0,0,0,0,0,
0,0,0,-1./360,0,-1./360,0,0,0,0,0,0,2./45,1./90,1./90,0,0,0,0,0,0,0,1./90,1./90,0,
-1./360,0,0,0,0,0,0,-1./360,0,0,0,0,1./90,2./45,1./90,0,0,0,0,1./90,0,0,0,0,1./90,
0,0,-1./360,0,0,0,0,0,-1./360,0,0,0,1./90,1./90,2./45,0,0,0,1./90,0,0,1./90,0,0,0,
0,0,0,0,0,0,-1./360,0,0,0,0,1./90,0,0,0,1./45,1./90,0,0,0,0,0,0,0,0,
0,0,0,0,-1./360,0,0,0,0,0,0,1./90,0,0,0,1./90,1./45,0,0,0,0,0,0,0,0,
0,0,0,0,0,-1./360,0,0,0,1./90,0,0,0,0,0,0,0,1./45,1./90,0,0,0,0,0,0,
0,0,0,0,-1./360,0,0,-1./360,0,1./90,0,0,0,0,1./90,0,0,1./90,2./45,0,0,1./90,0,0,0,
0,0,0,0,0,0,-1./360,0,-1./360,0,1./90,0,0,1./90,0,0,0,0,0,2./45,1./90,0,0,0,1./90,
0,0,0,0,0,-1./360,0,0,0,0,1./90,0,0,0,0,0,0,0,0,1./90,1./45,0,0,0,0,
0,0,0,0,0,-1./360,0,0,0,0,0,0,0,0,1./90,0,0,0,1./90,0,0,1./45,0,0,0,
0,0,0,0,0,0,0,0,-1./360,0,0,0,1./90,0,0,0,0,0,0,0,0,0,1./45,1./90,0,
0,0,0,0,0,0,0,-1./360,0,0,0,0,1./90,0,0,0,0,0,0,0,0,0,1./90,1./45,0,
0,0,0,0,0,-1./360,0,0,0,0,0,0,0,1./90,0,0,0,0,0,1./90,0,0,0,0,1./45;
}