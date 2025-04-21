#include "Geometry.h"

//helper function for skipping boundary nodes

inline bool bd_VertCheck(const Line& el, const uint& j) {
    return (el.vert[0] == j || el.vert[1] == j);
}

const geometry2& geometry2::operator=(const geometry2& src) {
    if (this == &src)
        return *this;

    this->Nelem = src.Nelem;
    this->Nvert = src.Nvert;
    this->NBdElemD = src.NBdElemD;
    this->NBdElemN = src.NBdElemN;
    //clear
    if (vertices != NULL)
        delete[] vertices;
    if (fIdx != NULL)
        delete[] fIdx;
    if (isInDirichletBoundary != NULL)
        delete[] isInDirichletBoundary;
    if (element != NULL)
        delete[] element;
    if (fVal != NULL)
        delete[] fVal;
    if (boundaryDirichlet != NULL)
        delete[] boundaryDirichlet;
    if (boundaryNeumann != NULL)
        delete[] boundaryNeumann;

    //create structures
    this->vertices = new coord<>[Nvert];
    this->fIdx = new uint[Nvert];
    this->isInDirichletBoundary = new bool[Nvert];
    this->fVal = new std::complex<dtype>[Nvert];
        this->element = new Tri[Nelem];
    this->boundaryDirichlet = new Line[NBdElemD];
    this->boundaryNeumann = new Line[NBdElemN];

    for (uint j(0); j < Nvert; ++j) {
        this->vertices[j] = src.vertices[j];
        this->fIdx[j] = src.fIdx[j];
        this->isInDirichletBoundary[j] = src.isInDirichletBoundary[j];
        this->fVal[j] = src.fVal[j];
    }

        for (uint j(0); j < Nelem; ++j) {
            this->element[j] = src.element[j];
        }
    for (uint j(0); j < NBdElemD; ++j)
        this->boundaryDirichlet[j] = src.boundaryDirichlet[j];

    for (uint j(0); j < NBdElemN; ++j)
        this->boundaryNeumann[j] = src.boundaryNeumann[j];

    return *this;
}

geometry2::geometry2(uint _Nvert, uint _Nelem, uint _BdElemD, uint _BdElemN) : Nvert(_Nvert), Nelem(_Nelem), NBdElemD(_BdElemD), NBdElemN(_BdElemN) {
    vertices = new coord<>[Nvert];
    fIdx = new uint[Nvert];
    isInDirichletBoundary = new bool[Nvert];
    fVal = new std::complex<dtype>[Nvert];
    boundaryDirichlet = new Line[NBdElemD];
    boundaryNeumann = new Line[NBdElemN];
        element = new Tri[Nelem];

    for (uint j(0); j < Nvert; ++j) {
        fIdx[j] = 0;
        fVal[j] = 0.0;
        isInDirichletBoundary[j] = false;
        vertices[j].x[0] = vertices[j].x[1] = 0.0;
    }

    for (uint j(0); j < NBdElemD; ++j)
        boundaryDirichlet[j].vert[0] = boundaryDirichlet[j].vert[1] = 0;
    for (uint j(0); j < NBdElemN; ++j)
        boundaryNeumann[j].vert[0] = boundaryNeumann[j].vert[1] = 0;
}

geometry2::~geometry2() {
    //destructor - handles memory 
    if (vertices != NULL)
        delete[] vertices;
    if (fIdx != NULL)
        delete[] fIdx;
    if (element != NULL)
        delete[] element;
    if (boundaryDirichlet != NULL)
        delete[] boundaryDirichlet;
    if (boundaryNeumann != NULL)
        delete[] boundaryNeumann;
    if (fVal != NULL)
        delete[] fVal;
    if (isInDirichletBoundary != NULL)
        delete[] isInDirichletBoundary;
}

void geometry2::resize(const uint _Nvert, const uint _Nelem, const uint _NBdElemD, const uint _NBdElemN) {
    if (vertices != NULL)
        delete[] vertices;
    if (fIdx != NULL)
        delete[] fIdx;
    if (element != NULL)
        delete[] element;
    if (boundaryDirichlet != NULL)
        delete[] boundaryDirichlet;
    if (boundaryNeumann != NULL)
        delete[] boundaryNeumann;
    if (fVal != NULL)
        delete[] fVal;
    if (isInDirichletBoundary != NULL)
        delete[] isInDirichletBoundary;

    this->Nvert = _Nvert;
    this->Nelem = _Nelem;
    this->NBdElemD = _NBdElemD;
    this->NBdElemN = _NBdElemN;
    vertices = new coord<>[Nvert];
    fIdx = new uint[Nvert];
    fVal = new std::complex<dtype>[Nvert];
        element = new Tri[Nelem];
    boundaryDirichlet = new Line[NBdElemD];
    boundaryNeumann = new Line[NBdElemN];
    isInDirichletBoundary = new bool[Nvert];

    //fill with blanks
    for (uint j(0); j<this->Nvert; ++j) {
        vertices[j].x[0] = vertices[j].x[1] = 0.0;
        fIdx[j] = 0;
        fVal[j] = 0.0;
        isInDirichletBoundary[j] = false;
    }
    for (uint j(0); j<this->Nelem; ++j)
            element[j].vert[0] = element[j].vert[1] = element[j].vert[2] = 0;
    for (uint j(0); j<this->NBdElemD; ++j)
        boundaryDirichlet[j].vert[0] = boundaryDirichlet[j].vert[1] = 0;
    for (uint j(0); j < NBdElemN; ++j)
        boundaryNeumann[j].vert[0] = boundaryNeumann[j].vert[1] = 0;
}

uint geometry2::updateFIdx() {
    //assign function numbers skipping boundaries
    uint funcIdx(0);
    std::vector<Line> bdElements(this->NBdElemD);

    //copy data
    for (uint j(0); j<this->NBdElemD; ++j) {
        bdElements[j] = this->boundaryDirichlet[j];
    }


    for (uint j(0); j < this->Nvert; ++j) {
        //clear
        this->fIdx[j] = 0;
        this->fVal[j] = 0.0;
        this->isInDirichletBoundary[j] = false;
        if (std::any_of(bdElements.begin(), bdElements.end(), std::bind(bd_VertCheck, std::placeholders::_1, j))) {
            isInDirichletBoundary[j] = true;
            continue;
        }

        this->fIdx[j] = funcIdx;
        funcIdx++;
    }
    return funcIdx;
}
