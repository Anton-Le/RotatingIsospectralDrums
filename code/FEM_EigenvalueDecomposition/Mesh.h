/* 
 * File:   Mesh.h
 * Author: lebedev
 *
 * Created on 21. März 2014, 14:34
 */

#ifndef MESH_H
#define	MESH_H

#include "Misc.h"
#include "Geometry.h"
#include <mgl2/mgl.h>
#include <mgl2/qt.h>
#include <string>

/*! \brief Mesh structure designated for 3D plots.
 * 
 * This data structure is used to store a version of the grid(mesh)
 * fit for plotting with MathGL.
 */

struct MGLTriMesh : public mglDraw {
    uint Ntri;  //!<    Number of elements (triangles) in the mesh.
    uint Nvert; //!<    Number of vertices in the mesh.
    mglPoint *vertices; //!<    Array of vertices in 3-space.
    Tri *Triangle;    //!<    Array of nodes.
    std::string LineStyle;      //!<    Variable storing the line-style of the mesh plot.
    std::string VertexStyle;    //!<    Variable storing the vertex style of the mesh plot.

    /*! \brief A dummy constructor.
     * 
     * \param _Nvert - number of vertices in the mesh.
     * \param _Nelem - number of elements(triangles) in the mesh.
     * 
     * \return *this - the current object initialized.
     * 
     * The constructor reserves space for _Nvert vertices and _Nelem
     * elements. It furthermore sets the style of the mesh lines to be
     * dashed grey and the style of the vertices to be a red (small) 'x'. */
    MGLTriMesh(const uint _Nvert, const uint _Nelem);
    /*! \brief The destructor.
     *
     * The default destructor frees the space reserved by the constructor for
     * the arrays of elements and vertices. 
     */
    ~MGLTriMesh();
    /*! \brief Function for updating the z-coordinates.
     * 
     * \param fVal - an array of size Nvert containing the
     * values of the wave function at each vertex.
     * 
     * \return max - the maximum value of |psi|^2 on the grid usable i.e.,
     * for scaling purposes.
     * 
     * The function updates the z-coordinates of the vertices to correspond
     * to the values of |psi|^2(x,y). 
     */
    dtype updateFval(const std::complex<dtype>* fVal);

    /*! \brief Redefined MGL Draw function.
     * 
     * \param gr - pointer to the MGLGraph object which will be used for
     * drawing.
     * 
     * The function draws the mesh essentially by iterating over the elements
     * (triangles) and drawing each triangle - without caring for doubly drawn
     * edges. */
    int Draw(mglGraph* gr);
};

#endif	/* MESH_H */

