#include "Mesh.h"
#include "Geometry.h"
#include <iostream>

//-------------struct

MGLTriMesh::MGLTriMesh(const uint _Nvert, const uint _Nelem) : Ntri(_Nelem), Nvert(_Nvert) {
    this->VertexStyle = "*r";
    this->LineStyle = "-x";

    //reserve space
    this->vertices = new mglPoint[Nvert];
    this->Triangle = new Tri[Ntri];
}

MGLTriMesh::~MGLTriMesh() {
    //clean-up
    if (vertices != NULL)
        delete[] vertices;
    if (Triangle != NULL)
        delete[] Triangle;

}

dtype MGLTriMesh::updateFval(const std::complex<dtype>* fVal) {
    dtype max(0);
    for(uint j(0); j<Nvert; ++j)
        vertices[j].z = std::norm(fVal[j]);
    for(uint j(0); j<Nvert; ++j)
        if(vertices[j].z > max)
            max = vertices[j].z;
    return max;
}

int MGLTriMesh::Draw(mglGraph* gr) {
//gr->Title("Mesh Plot");
    gr->Rotate(50, +10);
    gr->Box("k", false);
    gr->Grid("xyz", "B--");
    //gr->SetRotatedText(false);
    gr->SetFontDef("r:iC");
    gr->SetFontSize(4);
    //set the ranges


    //plot
    for (uint j(0); j < Ntri; ++j) {
        //for each triangle draw three lines
        //x1->x2
        gr->Line(vertices[ Triangle[j].vert[0] ], vertices[ Triangle[j].vert[1] ], this->LineStyle.c_str());
        //x2->x3
        gr->Line(vertices[ Triangle[j].vert[1] ], vertices[ Triangle[j].vert[2] ], this->LineStyle.c_str());
        //x3->x1
        gr->Line(vertices[ Triangle[j].vert[2] ], vertices[ Triangle[j].vert[0] ], this->LineStyle.c_str());
        //vertices
        gr->Mark(vertices[Triangle[j].vert[0] ], this->VertexStyle.c_str());
        gr->Mark(vertices[Triangle[j].vert[1] ], this->VertexStyle.c_str());
        gr->Mark(vertices[Triangle[j].vert[2] ], this->VertexStyle.c_str());
    }
    gr->Label('x', " x ", 1);
    gr->Label('y', " y ", 1);
    gr->Label('z', "|\\psi (x,y)|^2");
    //ticks
    gr->Adjust();


    return 0;
}