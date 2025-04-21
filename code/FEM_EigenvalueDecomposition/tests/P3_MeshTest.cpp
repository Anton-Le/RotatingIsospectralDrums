/*! \file P3_MeshTest.cpp
 * \author Anton Lebedev
 * \copyright GNU General Public License 3.0 (GPLv3)
 * \brief Unit tests for the cubic basis.
 * 
 * File:   P3_MeshTest.cpp
 * Author: aquinox
 *
 * Created on 31.05.2016, 20:55:27
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
//10-point triangle
//const uint triVert(10);


void createReferenceMesh(geometry2& refgeo){
    /*!\brief Create the reference mesh for P2 elements.
     * 
     * Computed by hand on [7.4] p.10.
     * No. Elements : 2
     * No. Boundary elements: 12
     * No. Vertices: 16
     * Coordinates:
     * (0,0), (1,0), (1,1), (0,1), (1,1/3), (1,2/3), (1/3,1/3), (2/3,2/3), (1/3,0), (2/3,0), (2/3,1), (1/,3,1), (0,1/3), (0,2/3)
     * (2/3,1/3) ,(1/3,2/3)
     * Elements:
     * 0: 1 2 0 4 5 6 7 8 9 14
     * 1: 2 3 0 10 11 12 13 7 6 15
     * Boundary:
     * (0,8) (1,4) (2,10) (3,13) (8,9) (9,1) (4,5) (5,2) (10,11) (11,3) (13,12) (12,0)
     *  */
    refgeo.resize(16,2,0,12);
    
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
    refgeo.vertices[4].x[1] = 1./3.;
    //5
    refgeo.vertices[5].x[0] = 1;
    refgeo.vertices[5].x[1] = 2./3.;
    //6
    refgeo.vertices[6].x[0] = 1./3.;
    refgeo.vertices[6].x[1] = 1./3.;
    //7
    refgeo.vertices[7].x[0] = 2./3.;
    refgeo.vertices[7].x[1] = 2./3.;
    //8
    refgeo.vertices[8].x[0] = 1./3.;
    refgeo.vertices[8].x[1] = 0;
    //9
    refgeo.vertices[9].x[0] = 2./3.;
    refgeo.vertices[9].x[1] = 0;
    //10
    refgeo.vertices[10].x[0] = 2./3.;
    refgeo.vertices[10].x[1] = 1;
    //11
    refgeo.vertices[11].x[0] = 1./3.;
    refgeo.vertices[11].x[1] = 1;
    //12
    refgeo.vertices[12].x[0] = 0;
    refgeo.vertices[12].x[1] = 1./3.;
    //13
    refgeo.vertices[13].x[0] = 0;
    refgeo.vertices[13].x[1] = 2./3.;
    //14
    refgeo.vertices[14].x[0] = 2./3.;
    refgeo.vertices[14].x[1] = 1./3.;
    //15
    refgeo.vertices[15].x[0] = 1./3.;
    refgeo.vertices[15].x[1] = 2./3.;
    
    //set the elements
    //0
    refgeo.element[0].vert[0] = 1;
    refgeo.element[0].vert[1] = 2;
    refgeo.element[0].vert[2] = 0;
    refgeo.element[0].vert[3] = 4;
    refgeo.element[0].vert[4] = 5;
    refgeo.element[0].vert[5] = 6;
    refgeo.element[0].vert[6] = 7;
    refgeo.element[0].vert[7] = 8;
    refgeo.element[0].vert[8] = 9;
    refgeo.element[0].vert[9] = 14;
    //1
    refgeo.element[1].vert[0] = 2;
    refgeo.element[1].vert[1] = 3;
    refgeo.element[1].vert[2] = 0;
    refgeo.element[1].vert[3] = 10;
    refgeo.element[1].vert[4] = 11;
    refgeo.element[1].vert[5] = 12;
    refgeo.element[1].vert[6] = 13;
    refgeo.element[1].vert[7] = 7;
    refgeo.element[1].vert[8] = 6;
    refgeo.element[1].vert[9] = 15;
    
    //set the boundary
    //0
    refgeo.boundaryNeumann[0].vert[0] = 0;
    refgeo.boundaryNeumann[0].vert[1] = 8;
    //1
    refgeo.boundaryNeumann[1].vert[0] = 1;
    refgeo.boundaryNeumann[1].vert[1] = 4;
    //2
    refgeo.boundaryNeumann[2].vert[0] = 2;
    refgeo.boundaryNeumann[2].vert[1] = 10;
    //3
    refgeo.boundaryNeumann[3].vert[0] = 3;
    refgeo.boundaryNeumann[3].vert[1] = 13;
    //4
    refgeo.boundaryNeumann[4].vert[0] = 8;
    refgeo.boundaryNeumann[4].vert[1] = 9;
    //5
    refgeo.boundaryNeumann[5].vert[0] = 9;
    refgeo.boundaryNeumann[5].vert[1] = 1;
    //6
    refgeo.boundaryNeumann[6].vert[0] = 4;
    refgeo.boundaryNeumann[6].vert[1] = 5;
    //7
    refgeo.boundaryNeumann[7].vert[0] = 5;
    refgeo.boundaryNeumann[7].vert[1] = 2;
    //8
    refgeo.boundaryNeumann[8].vert[0] = 10;
    refgeo.boundaryNeumann[8].vert[1] = 11;
    //9
    refgeo.boundaryNeumann[9].vert[0] = 11;
    refgeo.boundaryNeumann[9].vert[1] = 3;
    //10
    refgeo.boundaryNeumann[10].vert[0] = 13;
    refgeo.boundaryNeumann[10].vert[1] = 12;
    //11
    refgeo.boundaryNeumann[11].vert[0] = 12;
    refgeo.boundaryNeumann[11].vert[1] = 0;
}

void testP3base(const geometry2& refgeo, geometry2& p3geo, params& Pars){
    /* Necessary checks:
     * 
     * 1. Number of vertices.
     * 2. Number of boundary edges.
     * 3. Number of elements. (should be unchanged)
     * 4. Correctness of coordinates.
     * 5. Correctness of edge definitions.
     * 6. Correctness of element definitions.
     */ 
     //1. check number of vertices
    if(p3geo.Nvert != refgeo.Nvert){
	std::cout << "%TEST_FAILED% time=0 testname=testP3base (P3_MeshTest) message=Number of vertices erroneous." << std::endl;
	return;
      }
    else if((p3geo.NBdElemN != refgeo.NBdElemN) || (p3geo.NBdElemD != refgeo.NBdElemD)){
	std::cout << "%TEST_FAILED% time=0 testname=testP3base (P3_MeshTest) message=Number of boundary elements erroneous." << std::endl;
        return;
      }
    
    //check elements
    for(uint j(0); j<p3geo.Nelem; ++j)
        for(uint i(0); i<10; ++i)
            if( p3geo.element[j].vert[i] != refgeo.element[j].vert[i]){
                std::cout << "%TEST_FAILED% time=0 testname=testP3base (P3_MeshTest) message=Elements differ." << std::endl;
                return;
            }
    for( uint j(0); j<p3geo.NBdElemD; ++j)
        if( ( p3geo.boundaryDirichlet[j].vert[0] != refgeo.boundaryDirichlet[j].vert[0] ) || ( p3geo.boundaryDirichlet[j].vert[1] != refgeo.boundaryDirichlet[j].vert[1] ) ){
            std::cout << "%TEST_FAILED% time=0 testname=testP3base (P3_MeshTest) message=Dirichlet boundary elements differ!" << std::endl;
                return;
        }
    for( uint j(0); j<p3geo.NBdElemN; ++j)
        if( ( p3geo.boundaryNeumann[j].vert[0] != refgeo.boundaryNeumann[j].vert[0] ) || ( p3geo.boundaryNeumann[j].vert[1] != refgeo.boundaryNeumann[j].vert[1] ) ){
            std::cout << "%TEST_FAILED% time=0 testname=testP3base (P3_MeshTest) message=Neumann boundary elements differ!" << std::endl;
                return;
        }
    
        //check the coordinates
    for(uint j(0); j<p3geo.Nvert; ++j)
        if( (p3geo.vertices[j].x[0] != refgeo.vertices[j].x[0]) || (p3geo.vertices[j].x[1] != refgeo.vertices[j].x[1]) ){
            std::cout << "%TEST_FAILED% time=0 testname=testP3base (P3_MeshTest) message=Vertices differ." << std::endl;
            std::cerr << " x \t x_ref \t y \t y_ref"<<std::endl;
            for( uint s(0); s<p3geo.Nvert; ++s){
                std::cerr<<p3geo.vertices[s].x[0] << "\t " << refgeo.vertices[s].x[0];
                std::cerr<< "\t "<<p3geo.vertices[s].x[1] << "\t " << refgeo.vertices[s].x[1]<<std::endl;
            }
            return;
        }
}

void createReferenceMassMatrix(DensMat& M);

void createReferenceStiffnessMatrix(DensMat& A);

void createReferenceDampingMatrix(DensMat& C);

void createReferenceMassMatrixAlgOrder(DensMat& M);

void createReferenceStiffnessMatrixAlgOrder(DensMat& A);

void createReferenceDampingMatrixAlgOrder(DensMat& C);

void createReferenceDampingMatrixAlgOrderShifted(DensMat& C);

void testMassMatrix(const DensMat& Mref, const CSRMat& M);

void testDampingMatrix(const DensMat& Mref, const CSRMat& M);

void testStiffnessMatrix(const DensMat& Aref, const CSRMat& A);

int main(int argc, char** argv) {
    std::string CoordFname("Koordinaten_ref.dat"), TriFname("Elemente_ref.dat"), BdFname("Boundary_ref.dat");

    params RefPar;
    //parameter input
    prog_opt_init(argc, argv, RefPar);
    
    //hard-coded parameters
    RefPar.BoundaryFilenameN="Boundary_refRect.dat";
    RefPar.BoundaryFilenameD = "";
    RefPar.ElementsFilename="Elements_refRect.dat";
    RefPar.CoordinateFilename="Coord_refRect.dat";
    RefPar.Nref = 0;


    //read-in
    geometry2 geo = readData2D(RefPar.CoordinateFilename, RefPar.ElementsFilename, RefPar.BoundaryFilenameD, RefPar.BoundaryFilenameN);
    geometry2 geoOrig;
    geometry2 refgeoP3;
    geoOrig  = geo;
    //add P_3 coordinates
    p3basepoints(geo);
    
    
    std::cout << "%SUITE_STARTING% P3_MeshTest" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;
    
    //create reference mesh
    createReferenceMesh(refgeoP3);
    
    //test the refinement 
    std::cout << "%TEST_STARTED% testP3base (P3_MeshTest)" << std::endl;
    testP3base(refgeoP3, geo, RefPar);
    std::cout << "%TEST_FINISHED% time=0 testP3base (P3_MeshTest)\n" << std::endl;
    
    //test the matrices against a natural assignment of the global basis functions
    // to local ones
    std::cout << "###################################################################\n"<<std::endl;
    std::cout << "%TEST_STARTED% Natural node ordering matrices.\n" << std::endl;
    std::cout << "###################################################################\n"<<std::endl;
    
    std::cout << "%TEST_STARTED% testMassMatrix (P3_MeshTest)" << std::endl;
    DensMat Mref(16,16);
    createReferenceMassMatrix(Mref);
    
    //determine matrix characteristics
    uint freeNodes = geo.updateFIdx();
    Eigen::VectorXi NnzPerRow(freeNodes);
    DetermineSparsityPattern(geo, NnzPerRow);
    //create the mass matrix
    CSRMat M(freeNodes, freeNodes);
    M.reserve(NnzPerRow);
    MassMatrixSetup(M, geo, GD, QuadOrder, freeNodes, NnzPerRow);
    testMassMatrix(Mref, M);
    std::cout << "%TEST_FINISHED% time=0 testMassMatrix (P3_MeshTest)" << std::endl;
    
    
    std::cout << "%TEST_STARTED% testStiffnessMatrix (P3_MeshTest)" << std::endl;
    DensMat Aref(16,16);
    createReferenceStiffnessMatrix(Aref);
    
    CSRMat A(freeNodes, freeNodes);
    A.reserve(NnzPerRow);
    StiffnessMatrixSetup(A, geo, GD, QuadOrder, freeNodes, NnzPerRow);
    testStiffnessMatrix(Aref, A);
    std::cout << "%TEST_FINISHED% time=0 testStiffnessMatrix (P3_MeshTest)" << std::endl;
    
    std::cout << "%TEST_STARTED% testDampingMatrix (P3_MeshTest)" << std::endl;
    DensMat Cref(16,16);
    createReferenceDampingMatrix(Cref);
    CSRMat C(freeNodes, freeNodes);
    C.reserve(NnzPerRow);
    DampingMatrixSetup(C, geo, GD, QuadOrder, freeNodes, NnzPerRow);
    testDampingMatrix(Cref, C);
    std::cout << "%TEST_FINISHED% time=0 testDampingMatrix (P3_MeshTest)" << std::endl;
     std::cout << "###################################################################\n"<<std::endl;
    std::cout << "%TEST_FINISHED% Natural node ordering matrices.\n" << std::endl;
     std::cout << "###################################################################\n"<<std::endl;
    
    //test the matrices against an algorithmic assignment of global basis functions
    //to local ones.
    std::cout << "###################################################################\n"<<std::endl;
    std::cout << "%TEST_STARTED% Algorithmic node ordering matrices.\n" << std::endl;
    std::cout << "###################################################################\n"<<std::endl;
    
    std::cout << "%TEST_STARTED% testMassMatrix (P3_MeshTest)" << std::endl;
    createReferenceMassMatrixAlgOrder(Mref);
    testMassMatrix(Mref, M);
    std::cout << "%TEST_FINISHED% time=0 testMassMatrix (P3_MeshTest)" << std::endl;
    
    std::cout << "%TEST_STARTED% testStiffnessMatrix (P3_MeshTest)" << std::endl;
    createReferenceStiffnessMatrixAlgOrder(Aref);
    testStiffnessMatrix(Aref, A);
    std::cout << "%TEST_FINISHED% time=0 testStiffnessMatrix (P3_MeshTest)" << std::endl;
    
    std::cout << "%TEST_STARTED% testDampingMatrix (P3_MeshTest)" << std::endl;
    createReferenceDampingMatrixAlgOrder(Cref);
    testDampingMatrix(Cref, C);
    std::cout << "%TEST_FINISHED% time=0 testDampingMatrix (P3_MeshTest)" << std::endl;
    
    std::cout << "%TEST_STARTED% testDampingMatrixShifted (P3_MeshTest)" << std::endl;
    createReferenceDampingMatrixAlgOrderShifted(Cref);
    
    //shift the geometry
    dtype shiftX(5.0), shiftY(5.0);
    shift(geo, shiftX, shiftY);
    //recreate the matrix
    C.setZero();
    DampingMatrixSetup(C, geo, GD, QuadOrder, freeNodes, NnzPerRow);
    //test
    testDampingMatrix(Cref, C);
    std::cout << "%TEST_FINISHED% time=0 testDampingMatrixShifted (P3_MeshTest)" << std::endl;
    
    std::cout << "###################################################################\n"<<std::endl;
    std::cout << "%TEST_FINISHED% Algorithmic node ordering matrices.\n" << std::endl;
    std::cout << "###################################################################\n"<<std::endl;
        
    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

//functions for creating the reference matrix with naturally ordered basis functions
void createReferenceMassMatrix(DensMat& M){
    M.resize(16,16);
    M.setZero();
    M << 19./1680,11./13440,11./6720,11./13440,9./4480,9./4480,3./1120,0,3./2240,0,9./4480,9./4480,3./2240,0,3./1120,3./1120,
11./13440,19./3360,11./13440,0,3./2240,0,9./4480,9./4480,0,3./2240,0,0,0,0,3./1120,0,
11./6720,11./13440,19./1680,11./13440,0,3./2240,0,3./1120,9./4480,9./4480,3./2240,0,9./4480,9./4480,3./1120,3./1120,
11./13440,0,11./13440,19./3360,0,0,9./4480,9./4480,0,0,0,3./2240,0,3./2240,0,3./1120,
9./4480,3./2240,0,0,9./224,-9./640,-9./2240,-9./896,-9./896,9./448,0,0,0,0,27./2240,0,
9./4480,0,3./2240,0,-9./640,9./224,-9./896,9./448,-9./2240,-9./896,0,0,0,0,27./2240,0,
3./1120,9./4480,0,9./4480,-9./2240,-9./896,9./112,-9./320,9./448,-9./896,-9./896,-9./2240,9./448,-9./896,27./2240,27./2240,
0,9./4480,3./1120,9./4480,-9./896,9./448,-9./320,9./112,-9./896,-9./2240,9./448,-9./896,-9./896,-9./2240,27./2240,27./2240,
3./2240,0,9./4480,0,-9./896,-9./2240,9./448,-9./896,9./224,-9./640,0,0,0,0,27./2240,0,
0,3./2240,9./4480,0,9./448,-9./896,-9./896,-9./2240,-9./640,9./224,0,0,0,0,27./2240,0,
9./4480,0,3./2240,0,0,0,-9./896,9./448,0,0,9./224,-9./640,-9./2240,-9./896,0,27./2240,
9./4480,0,0,3./2240,0,0,-9./2240,-9./896,0,0,-9./640,9./224,-9./896,9./448,0,27./2240,
3./2240,0,9./4480,0,0,0,9./448,-9./896,0,0,-9./2240,-9./896,9./224,-9./640,0,27./2240,
0,0,9./4480,3./2240,0,0,-9./896,-9./2240,0,0,-9./896,9./448,-9./640,9./224,0,27./2240,
3./1120,3./1120,3./1120,0,27./2240,27./2240,27./2240,27./2240,27./2240,27./2240,0,0,0,0,81./560,0,
3./1120,0,3./1120,3./1120,0,0,27./2240,27./2240,0,0,27./2240,27./2240,27./2240,27./2240,0,81./560;
}

void createReferenceStiffnessMatrix(DensMat& A){
    A.resize(16,16);
    A.setZero();
    A << 17./20,-7./80,0,-7./80,-3./80,-3./80,3./40,3./40,-27./40,27./80,-3./80,-3./80,-27./40,27./80,0,0,
-7./80,17./20,-7./80,0,-51./80,3./8,-3./40,-3./40,3./8,-51./80,0,0,0,0,0,0,
0,-7./80,17./20,-7./80,27./80,-27./40,3./40,3./40,-3./80,-3./80,-27./40,27./80,-3./80,-3./80,0,0,
-7./80,0,-7./80,17./20,0,0,-3./40,-3./40,0,0,3./8,-51./80,3./8,-51./80,0,0,
-3./80,-51./80,27./80,0,27./8,-27./16,27./80,27./80,0,0,0,0,0,0,-81./40,0,
-3./80,3./8,-27./40,0,-27./16,27./8,27./80,-27./16,0,0,0,0,0,0,0,0,
3./40,-3./40,3./40,-3./40,27./80,27./80,27./4,-27./20,-27./16,27./80,27./80,27./80,-27./16,27./80,-81./40,-81./40,
3./40,-3./40,3./40,-3./40,27./80,-27./16,-27./20,27./4,27./80,27./80,-27./16,27./80,27./80,27./80,-81./40,-81./40,
-27./40,3./8,-3./80,0,0,0,-27./16,27./80,27./8,-27./16,0,0,0,0,0,0,
27./80,-51./80,-3./80,0,0,0,27./80,27./80,-27./16,27./8,0,0,0,0,-81./40,0,
-3./80,0,-27./40,3./8,0,0,27./80,-27./16,0,0,27./8,-27./16,0,0,0,0,
-3./80,0,27./80,-51./80,0,0,27./80,27./80,0,0,-27./16,27./8,0,0,0,-81./40,
-27./40,0,-3./80,3./8,0,0,-27./16,27./80,0,0,0,0,27./8,-27./16,0,0,
27./80,0,-3./80,-51./80,0,0,27./80,27./80,0,0,0,0,-27./16,27./8,0,-81./40,
0,0,0,0,-81./40,0,-81./40,-81./40,0,-81./40,0,0,0,0,81./10,0,
0,0,0,0,0,0,-81./40,-81./40,0,0,0,-81./40,0,-81./40,0,81./10;
}

void createReferenceDampingMatrix(DensMat& C){
    C.resize(16,16);
    C.setZero();
    C << 0,-29./4480,0,29./4480,-9./2240,-51./4480,0,0,-3./448,-3./4480,51./4480,9./2240,3./448,3./4480,3./320,-3./320,
29./4480,0,47./4480,0,141./2240,-99./2240,27./4480,93./4480,-9./2240,3./2240,0,0,0,0,-9./560,0,
0,-47./4480,0,31./1920,51./896,-51./448,0,0,-9./4480,-3./320,51./448,-51./896,9./4480,3./320,3./448,-3./448,
-29./4480,0,-31./1920,0,0,0,-27./4480,-93./4480,0,0,99./2240,-141./2240,9./2240,-3./2240,0,9./560,
9./2240,-141./2240,-51./896,0,0,621./4480,-117./4480,-9./112,-9./4480,-9./160,0,0,0,0,81./448,0,
51./4480,99./2240,51./448,0,-621./4480,0,-279./4480,117./560,9./2240,99./4480,0,0,0,0,-81./2240,0,
0,-27./4480,0,27./4480,117./4480,279./4480,0,0,-9./112,9./224,-279./4480,-117./4480,9./112,-9./224,-27./280,27./280,
0,-93./4480,0,93./4480,9./112,-117./560,0,0,153./4480,207./4480,117./560,-9./112,-153./4480,-207./4480,-27./112,27./112,
3./448,9./2240,9./4480,0,9./4480,-9./2240,9./112,-153./4480,0,-27./4480,0,0,0,0,-27./2240,0,
3./4480,-3./2240,3./320,0,9./160,-99./4480,-9./224,-207./4480,27./4480,0,0,0,0,0,459./2240,0,
-51./4480,0,-51./448,-99./2240,0,0,279./4480,-117./560,0,0,0,513./4480,-9./2240,-99./4480,0,81./2240,
-9./2240,0,51./896,141./2240,0,0,117./4480,9./112,0,0,-513./4480,0,9./4480,9./160,0,-81./448,
-3./448,0,-9./4480,-9./2240,0,0,-9./112,153./4480,0,0,9./2240,-9./4480,0,27./4480,0,27./2240,
-3./4480,0,-3./320,3./2240,0,0,9./224,207./4480,0,0,99./4480,-9./160,-27./4480,0,0,-459./2240,
-3./320,9./560,-3./448,0,-81./448,81./2240,27./280,27./112,27./2240,-459./2240,0,0,0,0,0,0,
3./320,0,3./448,-9./560,0,0,-27./280,-27./112,0,0,-81./2240,81./448,-27./2240,459./2240,0,0;
}

//functions for creating the reference matrix with algorithmically ordered basis functions

void createReferenceMassMatrixAlgOrder(DensMat& M){
    M.resize(16,16);
    M.setZero();
    M << 19./1680,11./13440,11./6720,11./13440,9./4480,9./4480,0,3./1120,3./2240,0,9./4480,9./4480,0,3./2240,3./1120,3./1120,
11./13440,19./3360,11./13440,0,3./2240,0,9./4480,9./4480,0,3./2240,0,0,0,0,3./1120,0,
11./6720,11./13440,19./1680,11./13440,0,3./2240,3./1120,0,9./4480,9./4480,3./2240,0,9./4480,9./4480,3./1120,3./1120,
11./13440,0,11./13440,19./3360,0,0,9./4480,9./4480,0,0,0,3./2240,3./2240,0,0,3./1120,
9./4480,3./2240,0,0,9./224,-9./640,-9./896,-9./2240,-9./896,9./448,0,0,0,0,27./2240,0,
9./4480,0,3./2240,0,-9./640,9./224,9./448,-9./896,-9./2240,-9./896,0,0,0,0,27./2240,0,
0,9./4480,3./1120,9./4480,-9./896,9./448,9./112,-9./320,-9./896,-9./2240,9./448,-9./896,-9./2240,-9./896,27./2240,27./2240,
3./1120,9./4480,0,9./4480,-9./2240,-9./896,-9./320,9./112,9./448,-9./896,-9./896,-9./2240,-9./896,9./448,27./2240,27./2240,
3./2240,0,9./4480,0,-9./896,-9./2240,-9./896,9./448,9./224,-9./640,0,0,0,0,27./2240,0,
0,3./2240,9./4480,0,9./448,-9./896,-9./2240,-9./896,-9./640,9./224,0,0,0,0,27./2240,0,
9./4480,0,3./2240,0,0,0,9./448,-9./896,0,0,9./224,-9./640,-9./896,-9./2240,0,27./2240,
9./4480,0,0,3./2240,0,0,-9./896,-9./2240,0,0,-9./640,9./224,9./448,-9./896,0,27./2240,
0,0,9./4480,3./2240,0,0,-9./2240,-9./896,0,0,-9./896,9./448,9./224,-9./640,0,27./2240,
3./2240,0,9./4480,0,0,0,-9./896,9./448,0,0,-9./2240,-9./896,-9./640,9./224,0,27./2240,
3./1120,3./1120,3./1120,0,27./2240,27./2240,27./2240,27./2240,27./2240,27./2240,0,0,0,0,81./560,0,
3./1120,0,3./1120,3./1120,0,0,27./2240,27./2240,0,0,27./2240,27./2240,27./2240,27./2240,0,81./560;
}

void createReferenceStiffnessMatrixAlgOrder(DensMat& A){
    A.resize(16,16);
    A.setZero();
    A << 17./20,-7./80,0,-7./80,-3./80,-3./80,3./40,3./40,-27./40,27./80,-3./80,-3./80,27./80,-27./40,0,0,
-7./80,17./20,-7./80,0,-51./80,3./8,-3./40,-3./40,3./8,-51./80,0,0,0,0,0,0,
0,-7./80,17./20,-7./80,27./80,-27./40,3./40,3./40,-3./80,-3./80,-27./40,27./80,-3./80,-3./80,0,0,
-7./80,0,-7./80,17./20,0,0,-3./40,-3./40,0,0,3./8,-51./80,-51./80,3./8,0,0,
-3./80,-51./80,27./80,0,27./8,-27./16,27./80,27./80,0,0,0,0,0,0,-81./40,0,
-3./80,3./8,-27./40,0,-27./16,27./8,-27./16,27./80,0,0,0,0,0,0,0,0,
3./40,-3./40,3./40,-3./40,27./80,-27./16,27./4,-27./20,27./80,27./80,-27./16,27./80,27./80,27./80,-81./40,-81./40,
3./40,-3./40,3./40,-3./40,27./80,27./80,-27./20,27./4,-27./16,27./80,27./80,27./80,27./80,-27./16,-81./40,-81./40,
-27./40,3./8,-3./80,0,0,0,27./80,-27./16,27./8,-27./16,0,0,0,0,0,0,
27./80,-51./80,-3./80,0,0,0,27./80,27./80,-27./16,27./8,0,0,0,0,-81./40,0,
-3./80,0,-27./40,3./8,0,0,-27./16,27./80,0,0,27./8,-27./16,0,0,0,0,
-3./80,0,27./80,-51./80,0,0,27./80,27./80,0,0,-27./16,27./8,0,0,0,-81./40,
27./80,0,-3./80,-51./80,0,0,27./80,27./80,0,0,0,0,27./8,-27./16,0,-81./40,
-27./40,0,-3./80,3./8,0,0,27./80,-27./16,0,0,0,0,-27./16,27./8,0,0,
0,0,0,0,-81./40,0,-81./40,-81./40,0,-81./40,0,0,0,0,81./10,0,
0,0,0,0,0,0,-81./40,-81./40,0,0,0,-81./40,-81./40,0,0,81./10;
}

void createReferenceDampingMatrixAlgOrder(DensMat& C){
    C.resize(16,16);
    C.setZero();
    C <<0,-29./4480,0,29./4480,-9./2240,-51./4480,0,0,-3./448,-3./4480,51./4480,9./2240,3./4480,3./448,3./320,-3./320,
29./4480,0,47./4480,0,141./2240,-99./2240,93./4480,27./4480,-9./2240,3./2240,0,0,0,0,-9./560,0,
0,-47./4480,0,31./1920,51./896,-51./448,0,0,-9./4480,-3./320,51./448,-51./896,3./320,9./4480,3./448,-3./448,
-29./4480,0,-31./1920,0,0,0,-93./4480,-27./4480,0,0,99./2240,-141./2240,-3./2240,9./2240,0,9./560,
9./2240,-141./2240,-51./896,0,0,621./4480,-9./112,-117./4480,-9./4480,-9./160,0,0,0,0,81./448,0,
51./4480,99./2240,51./448,0,-621./4480,0,117./560,-279./4480,9./2240,99./4480,0,0,0,0,-81./2240,0,
0,-93./4480,0,93./4480,9./112,-117./560,0,0,153./4480,207./4480,117./560,-9./112,-207./4480,-153./4480,-27./112,27./112,
0,-27./4480,0,27./4480,117./4480,279./4480,0,0,-9./112,9./224,-279./4480,-117./4480,-9./224,9./112,-27./280,27./280,
3./448,9./2240,9./4480,0,9./4480,-9./2240,-153./4480,9./112,0,-27./4480,0,0,0,0,-27./2240,0,
3./4480,-3./2240,3./320,0,9./160,-99./4480,-207./4480,-9./224,27./4480,0,0,0,0,0,459./2240,0,
-51./4480,0,-51./448,-99./2240,0,0,-117./560,279./4480,0,0,0,513./4480,-99./4480,-9./2240,0,81./2240,
-9./2240,0,51./896,141./2240,0,0,9./112,117./4480,0,0,-513./4480,0,9./160,9./4480,0,-81./448,
-3./4480,0,-3./320,3./2240,0,0,207./4480,9./224,0,0,99./4480,-9./160,0,-27./896,0,-459./2240,
-3./448,0,-9./4480,-9./2240,0,0,153./4480,-9./112,0,0,9./2240,-9./4480,27./896,0,0,27./2240,
-3./320,9./560,-3./448,0,-81./448,81./2240,27./112,27./280,27./2240,-459./2240,0,0,0,0,0,0,
3./320,0,3./448,-9./560,0,0,-27./112,-27./280,0,0,-81./2240,81./448,459./2240,-27./2240,0,0;
}

void createReferenceDampingMatrixAlgOrderShifted(DensMat& C){
    C.resize(16,16);
    C.setZero();
    C << 0,-1607./13440,0,1607./13440,-159./2240,-351./4480,0,0,-39./64,1257./4480,351./4480,159./2240,-1257./4480,39./64,201./2240,-201./2240,
1607./13440,0,901./13440,0,123./320,-639./2240,99./640,627./4480,-549./2240,723./2240,0,0,0,0,-99./560,0,
0,-901./13440,0,579./4480,303./896,-321./448,0,0,-309./4480,-171./2240,321./448,-303./896,171./2240,309./4480,39./448,-39./448,
-1607./13440,0,-579./4480,0,0,0,-99./640,-627./4480,0,0,639./2240,-123./320,-723./2240,549./2240,0,99./560,
159./2240,-123./320,-303./896,0,0,4401./4480,-153./224,-1737./4480,-549./4480,-9./160,0,0,0,0,135./64,0,
351./4480,639./2240,321./448,0,-4401./4480,0,927./560,-2439./4480,9./2240,639./4480,0,0,0,0,-621./2240,0,
0,-99./640,0,99./640,153./224,-927./560,0,0,2313./4480,261./640,927./560,-153./224,-261./640,-2313./4480,-27./14,27./14,
0,-627./4480,0,627./4480,1737./4480,2439./4480,0,0,-171./112,9./14,-2439./4480,-1737./4480,-9./14,171./112,-999./560,999./560,
39./64,549./2240,309./4480,0,549./4480,-9./2240,-2313./4480,171./112,0,-2727./4480,0,0,0,0,-81./320,0,
-1257./4480,-723./2240,171./2240,0,9./160,-639./4480,-261./640,-9./14,2727./4480,0,0,0,0,0,4779./2240,0,
-351./4480,0,-321./448,-639./2240,0,0,-927./560,2439./4480,0,0,0,459./640,-639./4480,-9./2240,0,621./2240,
-159./2240,0,303./896,123./320,0,0,153./224,1737./4480,0,0,-459./640,0,9./160,549./4480,0,-135./64,
1257./4480,0,-171./2240,723./2240,0,0,261./640,9./14,0,0,639./4480,-9./160,0,-783./896,0,-4779./2240,
-39./64,0,-309./4480,-549./2240,0,0,2313./4480,-171./112,0,0,9./2240,-549./4480,783./896,0,0,81./320,
-201./2240,99./560,-39./448,0,-135./64,621./2240,27./14,999./560,81./320,-4779./2240,0,0,0,0,0,0,
201./2240,0,39./448,-99./560,0,0,-27./14,-999./560,0,0,-621./2240,135./64,4779./2240,-81./320,0,0;
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