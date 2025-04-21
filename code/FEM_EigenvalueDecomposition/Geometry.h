/* \file Geometry.h
 * \brief File defining geometric structures.
 * 
 * \version 1.0
 * \date 10.09.2016
 * \copyright GNU General Public License 3.0 (GPLv3).
 * 
 * File:   Geometry.h.
 * Author: Anton Lebedev
 *
 * Created on 6. October 2015 11:15
 * 
 * Declaration of basic structures to hold and process meshes
 * in 1D (line elements) and 2D (triangle meshes).
 */

#ifndef GEOMETRY_H
#define	GEOMETRY_H
#include "Misc.h"

/*! \brief A template for coordinates.
 * 
 * Template struct for coordinates in \f$ \mathbb{R}^n\f$.
 * 
 * \tparam T Datatype of the coordinates.
 * \tparam n Dimension of the coordinate space.
 * 
 * No member functions defined explicitly,
 * compiler will handle that automatically. 
 */
template <typename T = dtype, unsigned n = 2 >
struct coord {
    T x[n];     //!< coordinate vector.
};

//template for element indices

/*! \brief A template for element indices.
 * Template struct for indices of the elements.
 * 
 * \tparam T - datatype of the index.
 * \tparam n - number of vertices.
 */
template <typename T = unsigned, unsigned n = 3 >
struct element {
    T vert[n];  //!< vertex-index array.
};
//! shorthand notation for line elements
typedef element<uint, 2> Line;
//! shorthand notation for triangle elements
typedef element<uint, triVert> Tri;
//! shorthand notation for quadrilateral elements
typedef element<uint, 4> Quad;

/*! \brief Data structure for 1D mesh.
 * 
 * The data structure is used to store
 * vertices, elements and function-to-vertex mappings
 * of a 1D inhomogeneous mesh.
 */
struct geometry {
    uint Nvert; //!<  number of vertices in the mesh (No.elements + 1).
    uint Nelem; //!< number of elements in the mesh.
    dtype x0, xMax;     //!< boundary values.
    dtype *vertices;    //!< array of vertices.
    uint *fIdx;         //!< vertex <-> function mapping array.
    Line *nodes;        //!< array of nodes .

    geometry(uint _Nvert, uint _Nelem) : Nvert(_Nvert), Nelem(_Nelem) {
        vertices = new dtype[Nvert];
        fIdx = new uint[Nvert];

        for (uint j(0); j < Nvert; ++j) {
            vertices[j] = 0.0;
            fIdx[j] = 0;
        }
        nodes = new Line[Nelem];
        for (uint j(0); j < Nelem; ++j) {
            //nodes[j].vert[0] = 0;
            //nodes[j].vert[1] = 0;
        }
        x0 = 0.0;
        xMax = 1.0;
    }

    const geometry& operator=(const geometry& src) {
        this->Nelem = src.Nelem;
        this->Nvert = src.Nvert;
        //clear
        if (vertices != NULL)
            delete[] vertices;
        if (fIdx != NULL)
            delete[] fIdx;
        if (nodes != NULL)
            delete[] nodes;

        //create structures
        this->vertices = new dtype[Nvert];
        this->fIdx = new uint[Nvert];
        this->nodes = new Line[Nelem];

        for (uint j(0); j < Nvert; ++j) {
            this->vertices[j] = src.vertices[j];
            this->fIdx[j] = src.fIdx[j];
        }

        for (uint j(0); j < Nelem; ++j) {
            this->nodes[j] = src.nodes[j];
        }

        this->x0 = src.x0;
        this->xMax = src.xMax;
        return *this;
    }

    void resize(const uint _Nvert, const uint _Nelem) {
        if (vertices != NULL)
            delete[] vertices;
        if (fIdx != NULL)
            delete[] fIdx;
        if (nodes != NULL)
            delete[] nodes;
        this->Nvert = _Nvert;
        this->Nelem = _Nelem;
        vertices = new dtype[Nvert];
        fIdx = new uint[Nvert];
        nodes = new Line[Nelem];
    }

    ~geometry() {
        //destructor - handles memory 
        if (vertices != NULL)
            delete[] vertices;
        if (fIdx != NULL)
            delete[] fIdx;
        if (nodes != NULL)
            delete[] nodes;
    }

};

/*! \brief Data structure for a 2D mesh.
 * 
 * The data structure is used to store
 * data about the mesh and the values of a function
 * on the vertices of the a 2D inhomogeneous mesh.
 */
struct geometry2 {
    uint Nvert; //!<  number of vertices in the mesh.
    uint Nelem; //!< number of elements in the mesh.
    uint NBdElemD;       //!< number of elements composing the Dirichlet domain boundary.
    uint NBdElemN;       //!< number of elements composing the Neumann domain boundary.
    coord<> *vertices;  //!< array of vertices.
    std::complex<dtype> *fVal;        //!< array of function values on the vertices of the mesh.
    uint *fIdx;         //!< vertex <-> function mapping array.
    bool *isInDirichletBoundary; //!< inidcator array for vertices in the Dirichlet boundary.
    Tri *element;     //!< array of elements.
    Line *boundaryDirichlet;     //!< array of boundary elements for the Dirichlet boundary.
    Line *boundaryNeumann;      //!< array of boundary elements for the Neumann boundary.

    /*! \brief A constructor.
     * 
     * The default constructor creates an entity of the struct and
     * allocates the required amount of storage space for the arrays.
     * 
     * \param _Nvert = 10 - number of vertices.
     * \param _Nelem = 10 - number of elements (triangles).
     * \param _BdElemD = 2 - number of boundary elements (lines) of the Dirichlet boundary.
     * \param _BdElemN = 3 - number of boundary elements (lines) of the Neumann boundary.
     */
    geometry2(uint _Nvert = 10, uint _Nelem = 10, uint _BdElemD = 2, uint _BdElemN = 3);
    /*! \brief Destructor.
     * 
     * Handles the release of the memory allocated for the arrays of the structure.
     */
    ~geometry2();
    /*! \brief Assignment operator.
     * 
     * Redefined assignemnt operator. Handles allocation and
     * de-allocation of memory for the dynamic arrays and the 
     * appropriate copying of data between the source and the current object.
     * Handles self-assignment preemptively.
     * \param src - const. reference to the source object.
     * \return &(*this) - reference to the current object.
     */
    const geometry2& operator=(const geometry2& src);
    /*! \brief Resize function.
     * 
     * The function resizes and zeroes the arrays of the current object.
     * All prior information stored in the object is LOST!
     * 
     * \param _Nvert - new number of vertices.
     * \param _Nelem - new number of elements (triangles).
     * \param _NBdElemD - new number of boundary elements (lines) of the Dirichlet boundary.
     * \param _NBdElemN - new number of boundary elements (lines) of the Neumann boundary.
     */
    void resize(const uint _Nvert, const uint _Nelem, const uint _NBdElem, const uint _NBdElemN);
    /*! \brief Update function for the vertex<-> function association. 
     * 
     * The function updates the indices of the base functions in the FE space
     * for each vertex.
     * \return Number of free (not bound) vertices. */
    uint updateFIdx();

};


#endif	/* GEOMETRY_H */

