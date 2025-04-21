#include "Functions.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <boost/program_options.hpp>
#include <vector>
#include <algorithm>
#include <utility>
#include <map>
#include <limits>
#include <sstream>

//pointers to the chosen functions
#ifdef P1BASIS
dtype (*Phi)(const uint, const dtype&, const dtype&) = &Phi_1;
void (*dPhi)(const uint, Eigen::Vector2d&, const dtype&, const dtype&) = &dPhi_1;
#elif P2BASIS
dtype (*Phi)(const uint, const dtype&, const dtype&) = &Phi_2;
void (*dPhi)(const uint, Eigen::Vector2d&, const dtype&, const dtype&) = &dPhi_2;
#else
dtype (*Phi)(const uint, const dtype&, const dtype&) = &Phi_3;
void (*dPhi)(const uint, Eigen::Vector2d&, const dtype&, const dtype&) = &dPhi_3;
#endif

void prog_opt_init(const int argc, const char *const argv[], params &par) {
    using namespace boost::program_options;

    //parameter description
    options_description options("Options for the 2D elliptic FEM solver.");
    options.add_options()("help,h", "display this help message")
            ("N", value<uint>(&(par.Nref))->default_value(3), "(unsigned int) number of refinement steps")
            ("CF", value<std::string>(&(par.CoordinateFilename))->default_value("Coordinates.dat"), "(string) name of the file with vertex coordinates")
            ("EF", value<std::string>(&(par.ElementsFilename))->default_value("Tri.dat"), "(string) name of the file with triangle<->vertex assignments")
            ("BFD", value<std::string>(&(par.BoundaryFilenameD))->default_value("Boundary_Dirichlet.dat"), "(string) name of the file with boundary element coordinates"
            " of the Dirichlet boundary")
            ("BFN", value<std::string>(&(par.BoundaryFilenameN))->default_value("Boundary_Neumann.dat"), "(string) name of the file with boundary element coordinates"
            " of the Neumann boundary")
            ("n", value<dtype>(&(par.n))->default_value(1.0), "(double) index of refraction")
            ("omega", value<std::string>(&(par.omega))->default_value("AngularFrequencies.dat"), "(string) name of the file containing a list of angular frequencies to be simulated.")
            ("Neig", value<uint>(&(par.Neig))->default_value(10), "(unsigned int) number of eigenvalues to compute")
            ("tol", value<dtype>(&(par.tol))->default_value(1e-8), "(double) relative error of eigenvalue computations (ARPACK)")
            ("recomp", value<uint>(&(par.recompute))->default_value(1), "{0,1} flag indicating whether the system matrices should be recomputed")
            ("FP", value<std::string>(&(par.MatPrefix))->default_value(""), "(string) prefix for matrix storage files.")
            ("xShft", value<dtype>(&(par.xShft))->default_value(0), "(double) x-shift of the geometry")
            ("yShft", value<dtype>(&(par.yShft))->default_value(0), "(double) y-shift of the geometry");

    variables_map varmap;

    try {
        // Parse command-line options
        store(parse_command_line(argc, argv, options), varmap);
        // Parse config file (if any))
        std::ifstream ifs("parameters.txt");
        if (ifs.good()) store(parse_config_file(ifs, options), varmap);
        //close config file
        ifs.close();
    } catch (std::exception& ex) {
        std::cout << ex.what() << std::endl << options << std::endl;
        exit(1);
    }

    // check whether help was requested
    if (varmap.count("help") > 0) {
        std::cout << options << std::endl;
        exit(0);
    }

    // save parameter values in the 'par' variable
    try {
        notify(varmap);
    } catch (std::exception& ex) {
        std::cout << ex.what() << std::endl << options << std::endl;
        exit(1);
    }
    //check whether the user requested a matrix load without specifying a file prefix
    //assert( (par.MatPrefix.empty() || par.recompute==0) );
}

//2D Mesh read-in

inline bool bd_VertCheck(const Line& el, const uint& j) {
    return (el.vert[0] == j || el.vert[1] == j);
}

geometry2 readData2D(const std::string VertexFilename, const std::string ElementFilename, const std::string BoundaryFilenameDirichlet,\
        const std::string BoundaryFilenameNeumann) {
    std::ifstream file;

    typedef coord<> coordinate;
    coordinate vertex;
    std::vector< coordinate > vertices;

    //read vertices
    file.open(VertexFilename.c_str());
    while (file >> vertex.x[0] >> vertex.x[1]) {
        vertices.push_back(vertex);
    }
    file.close();

    std::vector<Tri> elements;
    Tri el;
    //elements
    file.open(ElementFilename.c_str());
    while (file >> el.vert[0] >> el.vert[1] >> el.vert[2]) {
        elements.push_back(el);
    }
    file.close();

    //Dirichlet boundary
    Line bdEl;
    std::vector<Line> bdElementsD;

    file.open(BoundaryFilenameDirichlet.c_str());
    while (file >> bdEl.vert[0] >> bdEl.vert[1]) {
        bdElementsD.push_back(bdEl);
    }
    file.close();

    //Neumann boundary
    std::vector<Line> bdElementsN;

    file.open(BoundaryFilenameNeumann.c_str());
    while (file >> bdEl.vert[0] >> bdEl.vert[1]) {
        bdElementsN.push_back(bdEl);
    }
    file.close();

    //fit
    vertices.shrink_to_fit();
    elements.shrink_to_fit();
    bdElementsD.shrink_to_fit();
    bdElementsN.shrink_to_fit();

    //create our array
    geometry2 data(vertices.size(), elements.size(), bdElementsD.size(), bdElementsN.size());
    //copy data
    for (uint j(0); j < data.Nvert; ++j)
        data.vertices[j] = vertices[j];

    //store function Numbers skipping boundaries
    for (uint j(0); j < data.Nelem; ++j) {
        data.element[j] = elements[j];
    }
    for (uint j(0); j < data.NBdElemD; ++j)
        data.boundaryDirichlet[j] = bdElementsD[j];
    for (uint j(0); j < data.NBdElemN; ++j)
        data.boundaryNeumann[j] = bdElementsN[j];

    //assign function numbers skipping Dirichlet boundaries
    uint funcIdx(0);
    for (uint j(0); j < data.Nvert; ++j) {
        if (std::any_of(bdElementsD.begin(), bdElementsD.end(), std::bind(bd_VertCheck, std::placeholders::_1, j)))
            continue;

        data.fIdx[j] = funcIdx;
        funcIdx++;
    }
    return data;
}

std::vector<dtype> readOmega(const std::string omegaFile) {
    std::ifstream file;

    file.open(omegaFile.c_str());
    dtype omega(0);
    std::vector<dtype> omegaList;
    while (file >> omega) {
        omegaList.push_back(omega);
    }
    file.close();
    omegaList.shrink_to_fit();
    return omegaList;
}

void saveToFile(const uint Neig, const std::vector<dtype> &omega, cdtype **eigval) {
    std::ofstream file;

    file.open("output.dat");
    for (uint j(0); j < omega.size(); ++j)
        file << std::setprecision(DBL_DIG) << omega[j] << '\t';
    file << std::endl;
    for (uint n(0); n < Neig; ++n) {
        for (uint j(0); j < omega.size(); ++j) {
            file << std::setprecision(DBL_DIG) << eigval[n][j].real() << '+' << std::setprecision(DBL_DIG) << eigval[n][j].imag() << "i\t";
        }
        file << std::endl;
    }
    file.close();
}

//mesh refinement

using namespace std;
using namespace Eigen;


//helper function for dumping the result

void red_refine(dtype** Koordinaten,
                int** Elemente,
                int** Randdaten_Dirichlet,
                int** Randdaten_Neumann,
                geometry2& Pars) {

    // Vordefinieren der Variablen:
    uint imin(0), imax(0);
    uint nElem_neu(0), nKoord_neu(0), nRand_D_neu(0), nRand_N_neu(0);

    // ##############################################################
    // Neues Aufstellen der Elemente:
    // ##############################################################
    nElem_neu = Pars.Nelem;
    nKoord_neu = Pars.Nvert;

    // Vordefinieren von Vektoren/Matrizen
    uint Nachfolger[3]; //vektor der indices des nachfolgenden knotens
    Nachfolger[0] = 1;
    Nachfolger[1] = 2;
    Nachfolger[2] = 0;

    uint Knoten[6]; //vektor der Indices aller knoten eines Dreiecks
    uint KnotenRD[2]; //vektor der Indices aller Knoten eines Randes
    // Vordefinieren von K:

    //C++11 variante
    std::map< std::pair<uint, uint>, uint> K;
    std::pair<uint, uint> HashKey;


    for (uint i = 0; i < Pars.Nelem; ++i) {
        // Uebertragen auf Knoten:
        Knoten[0] = Elemente[i][0];
        Knoten[1] = Elemente[i][1];
        Knoten[2] = Elemente[i][2];

        for (int j = 0; j < 3; ++j) {
            // Bestimmen von imin, imax;
            //orient the edge using the global index of its boundary vertices
            imin = min(Elemente[i][j], Elemente[i][Nachfolger[j]]);
            imax = max(Elemente[i][j], Elemente[i][Nachfolger[j]]);

            HashKey.first = imin;
            HashKey.second = imax;

            if (K.count(HashKey) == 0) { //edge has not yet been subdivided

                // Updaten des KnotenVektors:
                Knoten[3 + j] = nKoord_neu;

                // Erstellen einer neuen Koordinate - Mittelpunkt der Seite
                for (uint k = 0; k < 2; ++k) {
                    //BUG off by 1
                    Koordinaten[nKoord_neu][k] = 0.5 * (Koordinaten[imin][k] + Koordinaten[imax][k]);
                }

                // Eintragen in K
                K.emplace(HashKey, nKoord_neu);

                // Updaten von nKoord_neu
                ++nKoord_neu;

            }
            else {

                Knoten[3 + j] = K.at(HashKey);
            }
        }

        // Updaten der Elemente:
        Elemente[i][0] = Knoten[3];
        Elemente[i][1] = Knoten[4];
        Elemente[i][2] = Knoten[5];

        Elemente[nElem_neu][0] = Knoten[0];
        Elemente[nElem_neu][1] = Knoten[3];
        Elemente[nElem_neu][2] = Knoten[5];
        ++nElem_neu; // Neues Element erstellt

        Elemente[nElem_neu][0] = Knoten[3];
        Elemente[nElem_neu][1] = Knoten[1];
        Elemente[nElem_neu][2] = Knoten[4];
        ++nElem_neu; // Neues Element erstellt

        Elemente[nElem_neu][0] = Knoten[4];
        Elemente[nElem_neu][1] = Knoten[2];
        Elemente[nElem_neu][2] = Knoten[5];
        ++nElem_neu; // Neues Element erstellt
    }


    // Updaten in der Parameter Struktur:
    Pars.Nvert = nKoord_neu;
    Pars.Nelem = nElem_neu;


    // ##############################################################
    // Neues Aufstellen der Dirichlet Randdaten:
    // ##############################################################

    KnotenRD[0] = KnotenRD[1] = 0;
    nRand_D_neu = Pars.NBdElemD;
    if (nRand_D_neu != 0) {
        for (uint i = 0; i < Pars.NBdElemD; ++i) {
            // Bestimmen von imin, imax;
            imin = min(Randdaten_Dirichlet[i][0], Randdaten_Dirichlet[i][1]);
            imax = max(Randdaten_Dirichlet[i][0], Randdaten_Dirichlet[i][1]);

            HashKey.first = imin;
            HashKey.second = imax;

            KnotenRD[0] = K.at(HashKey);
            KnotenRD[1] = Randdaten_Dirichlet[i][1];

            Randdaten_Dirichlet[i][1] = KnotenRD[0];
            Randdaten_Dirichlet[nRand_D_neu][0] = KnotenRD[0];
            Randdaten_Dirichlet[nRand_D_neu][1] = KnotenRD[1];

            // Erhoehen der Anzahl Randdaten:
            ++nRand_D_neu;
        }

        // Updaten in der Parameter Struktur:
        Pars.NBdElemD = nRand_D_neu;
    }

    //    // ##############################################################
    //    // Neues Aufstellen der Neumann Randdaten:
    //    // ##############################################################
    KnotenRD[0] = KnotenRD[1] = 0;
    nRand_N_neu = Pars.NBdElemN;
    if (nRand_N_neu != 0) {
        for (uint i = 0; i < Pars.NBdElemN; ++i) {
            // Bestimmen von imin, imax;

            imin = min(Randdaten_Neumann[i][0], Randdaten_Neumann[i][1]);
            imax = max(Randdaten_Neumann[i][0], Randdaten_Neumann[i][1]);

            HashKey.first = imin;
            HashKey.second = imax;

            KnotenRD[0] = K.at(HashKey);
            KnotenRD[1] = Randdaten_Neumann[i][1];


            Randdaten_Neumann[i][1] = KnotenRD[0];
            Randdaten_Neumann[nRand_N_neu][0] = KnotenRD[0];
            Randdaten_Neumann[nRand_N_neu][1] = KnotenRD[1];

            // Erhoehen der Anzahl Randdaten:
            ++nRand_N_neu;
        }

        // Updaten in der Parameter Struktur:
        Pars.NBdElemN = nRand_N_neu;
    }
}

void ref2d(geometry2& Pars, const uint Nref) {

    /* Vor dem Aufruf der Funktion Gebietsverfeinerung sollte geprueft werden,
     ob Pars.Verfeinerung > 0 erfuellt ist. Falls nicht fuehrt dieses Programm
     nicht zur gewuenschten Verfeinerung.
     */
    if (Nref == 0)
        return;

    // ##############################################################
    // Berechnen, wiegross die verfeinerten Strukturen werden:
    uint vorfaktor, Anz_Elemente, Anz_Knoten, Anz_Kanten, Anz_Rand_D, Anz_Rand_N;
    vorfaktor = pow(2, Nref);

    Anz_Elemente = vorfaktor * vorfaktor * Pars.Nelem;
    Anz_Rand_D = vorfaktor * Pars.NBdElemD;
    Anz_Rand_N = vorfaktor * Pars.NBdElemN; // Falls NeumannRand
    Anz_Kanten = (3 * Anz_Elemente + Anz_Rand_D) / 2;
    Anz_Kanten = (3 * Anz_Elemente + Anz_Rand_D + Anz_Rand_N) / 2; // Falls NeumannRand
    Anz_Knoten = 1 + Anz_Kanten - Anz_Elemente;

    //allocate memory for the refined mesh
    int **Elemente, **Randdaten_D, **Randdaten_N;
    dtype **Koordinaten;

    Elemente = new int*[Anz_Elemente];
    for (uint j(0); j < Anz_Elemente; ++j)
        Elemente[j] = new int[3];

    Randdaten_D = new int*[Anz_Rand_D];
    for (uint j(0); j < Anz_Rand_D; ++j)
        Randdaten_D[j] = new int[2];

    Randdaten_N = new int*[Anz_Rand_N];
    for (uint j(0); j < Anz_Rand_N; ++j)
        Randdaten_N[j] = new int[2];

    Koordinaten = new dtype*[Anz_Knoten];
    for (uint j(0); j < Anz_Knoten; ++j)
        Koordinaten[j] = new dtype[2];

    //blank the allocated arrays
    for (uint j(0); j < Anz_Elemente; ++j)
        Elemente[j][0] = Elemente[j][1] = Elemente[j][2] = 0;
    for (uint j(0); j < Anz_Rand_D; ++j)
        Randdaten_D[j][0] = Randdaten_D[j][1] = 0;
    for (uint j(0); j < Anz_Rand_N; ++j)
        Randdaten_N[j][0] = Randdaten_N[j][1] = 0;
    for (uint j(0); j < Anz_Knoten; ++j)
        Koordinaten[j][0] = Koordinaten[j][1] = 0.0;

    // ##############################################################
    // Uebertragen in Elemente und Koordinaten:
    for (uint j(0); j < Pars.Nvert; ++j) {
        Koordinaten[j][0] = Pars.vertices[j].x[0];
        Koordinaten[j][1] = Pars.vertices[j].x[1];
    }
    for (uint j(0); j < Pars.Nelem; ++j) {
        Elemente[j][0] = Pars.element[j].vert[0];
        Elemente[j][1] = Pars.element[j].vert[1];
        Elemente[j][2] = Pars.element[j].vert[2];
    }
    for (uint j(0); j < Pars.NBdElemD; ++j) {
        Randdaten_D[j][0] = Pars.boundaryDirichlet[j].vert[0];
        Randdaten_D[j][1] = Pars.boundaryDirichlet[j].vert[1];
    }
    for (uint j(0); j < Pars.NBdElemN; ++j) {
        Randdaten_N[j][0] = Pars.boundaryNeumann[j].vert[0];
        Randdaten_N[j][1] = Pars.boundaryNeumann[j].vert[1];
    }

    // ##############################################################
    // Durchfuehren der Verfeinerung:
    for (uint v = 0; v < Nref; ++v) {
        red_refine(Koordinaten,
                   Elemente,
                   Randdaten_D,
                   Randdaten_N,
                   Pars);
    }


    // ##############################################################
    // Speichern:
    /*Ausgeben_KER(Koordinaten, "Coordinates_refined.dat",
            Elemente, "Elements_refined.dat",
            Randdaten_D, "Dirichlet_refined.dat",
            Randdaten_N, "Neumann_refined.dat",
            Pars);
     */

    //feed back
    //copy back to the geometry
    Pars.resize(Pars.Nvert, Pars.Nelem, Pars.NBdElemD, Pars.NBdElemN);

    for (uint j(0); j < Pars.Nvert; ++j) {
        Pars.vertices[j].x[0] = Koordinaten[j][0];
        Pars.vertices[j].x[1] = Koordinaten[j][1];

    }
    for (uint j(0); j < Pars.Nelem; ++j) {
        Pars.element[j].vert[0] = Elemente[j][0];
        Pars.element[j].vert[1] = Elemente[j][1];
        Pars.element[j].vert[2] = Elemente[j][2];
    }
    if (Pars.NBdElemD != 0) {
        for (uint j(0); j < Pars.NBdElemD; ++j) {
            Pars.boundaryDirichlet[j].vert[0] = Randdaten_D[j][0];
            Pars.boundaryDirichlet[j].vert[1] = Randdaten_D[j][1];
        }
    }
    if (Pars.NBdElemN != 0) {
        for (uint j(0); j < Pars.NBdElemN; ++j) {
            Pars.boundaryNeumann[j].vert[0] = Randdaten_N[j][0];
            Pars.boundaryNeumann[j].vert[1] = Randdaten_N[j][1];
        }
    }

    //clean-up

    for (uint j(0); j < Anz_Elemente; ++j)
        delete[] Elemente[j];
    delete[] Elemente;

    for (uint j(0); j < Anz_Rand_D; ++j)
        delete[] Randdaten_D[j];
    delete[] Randdaten_D;

    for (uint j(0); j < Anz_Rand_N; ++j)
        delete[] Randdaten_N[j];
    delete[] Randdaten_N;

    for (uint j(0); j < Anz_Knoten; ++j)
        delete[] Koordinaten[j];
    delete[] Koordinaten;
}

//functions for adding edge centers to the mesh to use for P_2 FE.

void EdgeSplit(dtype** Koordinaten,
               uint** Elemente,
               uint** Randdaten_Dirichlet,
               uint** Randdaten_Neumann,
               geometry2& Pars) {

    // Vordefinieren der Variablen:
    uint imin, imax;
    uint nKoord_neu, nRand_D_neu, nRand_N_neu;

    // ##############################################################
    // Neues Aufstellen der Elemente:
    // ##############################################################
    nKoord_neu = Pars.Nvert;

    // Vordefinieren von Vektoren/Matrizen
    uint Nachfolger[3]; //vektor der indices des nachfolgenden knotens
    Nachfolger[0] = 1;
    Nachfolger[1] = 2;
    Nachfolger[2] = 0;

    uint Knoten[6]; //vektor der Indices aller knoten eines Dreiecks
    uint KnotenRD[2]; //vektor der Indices aller Knoten eines Randes
    // Vordefinieren von K:

    //C++11 variante
    std::map< std::pair<uint, uint>, uint> K;
    std::pair<uint, uint> HashKey;


    for (uint i = 0; i < Pars.Nelem; ++i) {
        // Uebertragen auf Knoten:
        Knoten[0] = Elemente[i][0];
        Knoten[1] = Elemente[i][1];
        Knoten[2] = Elemente[i][2];

        for (int j = 0; j < 3; ++j) {
            // Bestimmen von imin, imax;
            //orient the edge using the global index of its boundary vertices
            imin = min(Elemente[i][j], Elemente[i][Nachfolger[j]]);
            imax = max(Elemente[i][j], Elemente[i][Nachfolger[j]]);

            HashKey.first = imin;
            HashKey.second = imax;

            if (K.count(HashKey) == 0) { //edge has not yet been subdivided

                // Updaten des KnotenVektors:
                Knoten[3 + j] = nKoord_neu;

                // Erstellen einer neuen Koordinate - Mittelpunkt der Seite
                for (int k = 0; k < 2; ++k) {
                    Koordinaten[nKoord_neu][k] = 0.5 * (Koordinaten[imin][k] + Koordinaten[imax][k]);
                }

                // Eintragen in K
                K.emplace(HashKey, nKoord_neu);

                // Updaten von nKoord_neu
                ++nKoord_neu;

            }
            else {

                Knoten[3 + j] = K.at(HashKey);
            }
        }

        // Updaten der Elemente:
        Elemente[i][0] = Knoten[0];
        Elemente[i][1] = Knoten[1];
        Elemente[i][2] = Knoten[2];

        Elemente[i][3] = Knoten[3];
        Elemente[i][4] = Knoten[4];
        Elemente[i][5] = Knoten[5];
    }


    // Updaten in der Parameter Struktur:
    Pars.Nvert = nKoord_neu;


    // ##############################################################
    // Neues Aufstellen der Dirichlet Randdaten:
    // ##############################################################

    KnotenRD[0] = KnotenRD[1] = 0;
    nRand_D_neu = Pars.NBdElemD;
    if (nRand_D_neu != 0) {
        for (uint i = 0; i < Pars.NBdElemD; ++i) {
            // Bestimmen von imin, imax;
            imin = min(Randdaten_Dirichlet[i][0], Randdaten_Dirichlet[i][1]);
            imax = max(Randdaten_Dirichlet[i][0], Randdaten_Dirichlet[i][1]);

            HashKey.first = imin;
            HashKey.second = imax;

            KnotenRD[0] = K.at(HashKey);
            KnotenRD[1] = Randdaten_Dirichlet[i][1];

            Randdaten_Dirichlet[i][1] = KnotenRD[0];
            Randdaten_Dirichlet[nRand_D_neu][0] = KnotenRD[0];
            Randdaten_Dirichlet[nRand_D_neu][1] = KnotenRD[1];

            // Erhoehen der Anzahl Randdaten:
            ++nRand_D_neu;
        }

        // Updaten in der Parameter Struktur:
        Pars.NBdElemD = nRand_D_neu;
    }

    //    // ##############################################################
    //    // Neues Aufstellen der Neumann Randdaten:
    //    // ##############################################################
    KnotenRD[0] = KnotenRD[1] = 0;
    nRand_N_neu = Pars.NBdElemN;
    if (nRand_N_neu != 0) {
        for (uint i = 0; i < Pars.NBdElemN; ++i) {
            // Bestimmen von imin, imax;

            imin = min(Randdaten_Neumann[i][0], Randdaten_Neumann[i][1]);
            imax = max(Randdaten_Neumann[i][0], Randdaten_Neumann[i][1]);

            HashKey.first = imin;
            HashKey.second = imax;

            KnotenRD[0] = K.at(HashKey);
            KnotenRD[1] = Randdaten_Neumann[i][1];


            Randdaten_Neumann[i][1] = KnotenRD[0];
            Randdaten_Neumann[nRand_N_neu][0] = KnotenRD[0];
            Randdaten_Neumann[nRand_N_neu][1] = KnotenRD[1];

            // Erhoehen der Anzahl Randdaten:
            ++nRand_N_neu;
        }

        // Updaten in der Parameter Struktur:
        Pars.NBdElemN = nRand_N_neu;
    }
}

void p2basepoints(geometry2& Pars) {
    /* Essentially this function is a copy of the ref2d function, except with a
     * fixed "refinement" degree */
    // Berechnen, wiegross die verfeinerten Strukturen werden:
    uint vorfaktor(2), Anz_Elemente, Anz_Knoten, Anz_Kanten, Anz_Rand_D, Anz_Rand_N;

    Anz_Elemente = SQU(vorfaktor) * Pars.Nelem;
    Anz_Rand_D = vorfaktor * Pars.NBdElemD;
    Anz_Rand_N = vorfaktor * Pars.NBdElemN; // Falls NeumannRand
    Anz_Kanten = (3 * Anz_Elemente + Anz_Rand_D) / 2;
    Anz_Kanten = (3 * Anz_Elemente + Anz_Rand_D + Anz_Rand_N) / 2; // Falls NeumannRand
    Anz_Knoten = 1 + Anz_Kanten - Anz_Elemente;

    //reset # of elements to the old value
    Anz_Elemente = Pars.Nelem;
    //allocate memory for the refined mesh
    uint **Elemente, **Randdaten_D, **Randdaten_N;
    dtype **Koordinaten;

    Elemente = new uint*[Anz_Elemente];
    for (uint j(0); j < Anz_Elemente; ++j)
        Elemente[j] = new uint[6];

    Randdaten_D = new uint*[Anz_Rand_D];
    for (uint j(0); j < Anz_Rand_D; ++j)
        Randdaten_D[j] = new uint[2];

    Randdaten_N = new uint*[Anz_Rand_N];
    for (uint j(0); j < Anz_Rand_N; ++j)
        Randdaten_N[j] = new uint[2];

    Koordinaten = new dtype*[Anz_Knoten];
    for (uint j(0); j < Anz_Knoten; ++j)
        Koordinaten[j] = new dtype[2];

    //blank the allocated arrays
    for (uint j(0); j < Anz_Elemente; ++j)
        for (uint k(0); k < 6; ++k)
            Elemente[j][k] = 0;

    for (uint j(0); j < Anz_Rand_D; ++j)
        Randdaten_D[j][0] = Randdaten_D[j][1] = 0;
    for (uint j(0); j < Anz_Rand_N; ++j)
        Randdaten_N[j][0] = Randdaten_N[j][1] = 0;
    for (uint j(0); j < Anz_Knoten; ++j)
        Koordinaten[j][0] = Koordinaten[j][1] = 0.0;

    // ##############################################################
    // Uebertragen in Elemente und Koordinaten:
    for (uint j(0); j < Pars.Nvert; ++j) {
        Koordinaten[j][0] = Pars.vertices[j].x[0];
        Koordinaten[j][1] = Pars.vertices[j].x[1];
    }
    for (uint j(0); j < Pars.Nelem; ++j) {
        Elemente[j][0] = Pars.element[j].vert[0];
        Elemente[j][1] = Pars.element[j].vert[1];
        Elemente[j][2] = Pars.element[j].vert[2];
    }
    for (uint j(0); j < Pars.NBdElemD; ++j) {
        Randdaten_D[j][0] = Pars.boundaryDirichlet[j].vert[0];
        Randdaten_D[j][1] = Pars.boundaryDirichlet[j].vert[1];
    }
    for (uint j(0); j < Pars.NBdElemN; ++j) {
        Randdaten_N[j][0] = Pars.boundaryNeumann[j].vert[0];
        Randdaten_N[j][1] = Pars.boundaryNeumann[j].vert[1];
    }

    // ##############################################################
    // Durchfuehren der Verfeinerung:
    EdgeSplit(Koordinaten, Elemente, Randdaten_D, Randdaten_N, Pars);

    // ##############################################################

    //feed back
    //copy back to the geometry
    //std::cout<< "Koordinaten"<<std::endl;

    Pars.resize(Pars.Nvert, Pars.Nelem, Pars.NBdElemD, Pars.NBdElemN);

    for (uint j(0); j < Pars.Nvert; ++j) {
        Pars.vertices[j].x[0] = Koordinaten[j][0];
        Pars.vertices[j].x[1] = Koordinaten[j][1];
        //std::cout<<j<<": ("<<Koordinaten[j][0]<<","<<Koordinaten[j][1]<<")"<<std::endl;
    }
    for (uint j(0); j < Pars.Nelem; ++j) {
        for (uint k(0); k < triVert; ++k) {
            Pars.element[j].vert[k] = Elemente[j][k];
        }
    }
    if (Pars.NBdElemD != 0) {
        for (uint j(0); j < Pars.NBdElemD; ++j) {
            Pars.boundaryDirichlet[j].vert[0] = Randdaten_D[j][0];
            Pars.boundaryDirichlet[j].vert[1] = Randdaten_D[j][1];
        }
    }
    if (Pars.NBdElemN != 0) {
        for (uint j(0); j < Pars.NBdElemN; ++j) {
            Pars.boundaryNeumann[j].vert[0] = Randdaten_N[j][0];
            Pars.boundaryNeumann[j].vert[1] = Randdaten_N[j][1];
        }
    }

    //clean-up

    for (uint j(0); j < Anz_Elemente; ++j)
        delete[] Elemente[j];
    delete[] Elemente;

    for (uint j(0); j < Anz_Rand_D; ++j)
        delete[] Randdaten_D[j];
    delete[] Randdaten_D;

    for (uint j(0); j < Anz_Rand_N; ++j)
        delete[] Randdaten_N[j];
    delete[] Randdaten_N;

    for (uint j(0); j < Anz_Knoten; ++j)
        delete[] Koordinaten[j];
    delete[] Koordinaten;
}

void Edge2Split(dtype** Koordinaten,
                uint** Elemente,
                uint** Randdaten_Dirichlet,
                uint** Randdaten_Neumann,
                geometry2& Pars) {

    // Vordefinieren der Variablen:
    uint imin, imax;
    uint nKoord_neu, nRand_D_neu, nRand_N_neu;

    // ##############################################################
    // Neues Aufstellen der Elemente:
    // ##############################################################
    nKoord_neu = Pars.Nvert;

    // Vordefinieren von Vektoren/Matrizen
    uint Nachfolger[3]; //vektor der indices des nachfolgenden knotens
    Nachfolger[0] = 1;
    Nachfolger[1] = 2;
    Nachfolger[2] = 0;

    uint Knoten[10]; //vektor der Indices aller knoten eines Dreiecks
    uint KnotenRD[2]; //vektor der Indices aller Knoten eines Randes
    // Vordefinieren von K:

    //C++11 variante
    std::map< std::pair<uint, uint>, std::pair<uint, uint> > K;
    std::pair<uint, uint> HashKey;
    std::pair<uint, uint> EdgeVertices;


    for (uint i = 0; i < Pars.Nelem; ++i) {
        // Uebertragen auf Knoten:
        Knoten[0] = Elemente[i][0];
        Knoten[1] = Elemente[i][1];
        Knoten[2] = Elemente[i][2];

        for (int j = 0; j < 3; ++j) {
            // Bestimmen von imin, imax;
            //orient the edge using the global index of its boundary vertices
            imin = min(Elemente[i][j], Elemente[i][Nachfolger[j]]);
            imax = max(Elemente[i][j], Elemente[i][Nachfolger[j]]);

            HashKey.first = imin;
            HashKey.second = imax;

            if (K.count(HashKey) == 0) { //edge has not yet been subdivided

                // create new coordinates 1/3|2/3 of the edge
                for (int k = 0; k < 2; ++k) {
                    Koordinaten[nKoord_neu][k] = (1.0 / 3.0) * (2.0 * Koordinaten[imin][k] + Koordinaten[imax][k]);
                    Koordinaten[nKoord_neu + 1][k] = (1.0 / 3.0) * (Koordinaten[imin][k] + 2.0 * Koordinaten[imax][k]);
                }


                EdgeVertices.first = nKoord_neu;
                EdgeVertices.second = nKoord_neu + 1;

                // Store into the bookkeeping array K
                K.emplace(HashKey, EdgeVertices);

                // Updaten des KnotenVektors:
                Knoten[3 + 2 * j] = nKoord_neu++;
                //offset 6
                Knoten[3 + 2 * j + 1] = nKoord_neu++;


            }
            else {
                //swap the order to account for an intrinsic ordering of the nodes
                //and an assymatry of the functions on the edges.
                EdgeVertices = K.at(HashKey);
                Knoten[3 + 2 * j] = EdgeVertices.second;
                //offset 6
                Knoten[3 + 2 * j + 1] = EdgeVertices.first;
            }
        }

        // Updaten der Elemente:
        Elemente[i][0] = Knoten[0];
        Elemente[i][1] = Knoten[1];
        Elemente[i][2] = Knoten[2];

        Elemente[i][3] = Knoten[3];
        Elemente[i][4] = Knoten[4];
        Elemente[i][5] = Knoten[5];

        Elemente[i][6] = Knoten[6];
        Elemente[i][7] = Knoten[7];
        Elemente[i][8] = Knoten[8];
    }

    //add the coordinates at the centre of each triangle
    for (uint i(0); i < Pars.Nelem; ++i) {
        for (uint k = 0; k < 2; ++k) {
            Koordinaten[nKoord_neu][k] = (1.0 / 3.0) * (Koordinaten[ Elemente[i][0] ][k] +\
 Koordinaten[ Elemente[i][1] ][k] + Koordinaten[ Elemente[i][2] ][k]);
        }
        //add the final vertex to the triangle
        Elemente[i][9] = nKoord_neu;
        ++nKoord_neu;
    }

    // Updaten in der Parameter Struktur:
    Pars.Nvert = nKoord_neu;


    // ##############################################################
    // Neues Aufstellen der Dirichlet Randdaten:
    // ##############################################################

    KnotenRD[0] = KnotenRD[1] = 0;
    nRand_D_neu = Pars.NBdElemD;
    if (nRand_D_neu != 0) {
        for (uint i = 0; i < Pars.NBdElemD; ++i) {
            // Bestimmen von imin, imax;
            imin = min(Randdaten_Dirichlet[i][0], Randdaten_Dirichlet[i][1]);
            imax = max(Randdaten_Dirichlet[i][0], Randdaten_Dirichlet[i][1]);

            HashKey.first = imin;
            HashKey.second = imax;
            EdgeVertices = K.at(HashKey);

            if (imin == Randdaten_Dirichlet[i][0]) {
                //the orientation of the boundary edge definition 
                //coincides with the convention made for the refinement
                KnotenRD[0] = EdgeVertices.first;
                KnotenRD[1] = Randdaten_Dirichlet[i][1];

                Randdaten_Dirichlet[i][1] = KnotenRD[0];
                //2nd third of the edge
                Randdaten_Dirichlet[nRand_D_neu][0] = KnotenRD[0];
                Randdaten_Dirichlet[nRand_D_neu][1] = EdgeVertices.second;
                //increase # of boundary edges
                ++nRand_D_neu;
                //3rd third of the edge
                Randdaten_Dirichlet[nRand_D_neu][0] = EdgeVertices.second;
                Randdaten_Dirichlet[nRand_D_neu][1] = KnotenRD[1]; //the old right vertex

                // Erhoehen der Anzahl Randdaten:
                ++nRand_D_neu;
            }
            else {
                //the orientation of the boundary edge is opposite to the
                //convention made for refinement
                KnotenRD[0] = Randdaten_Dirichlet[i][1];
                KnotenRD[1] = EdgeVertices.second;

                Randdaten_Dirichlet[i][1] = KnotenRD[1];
                //2nd third of the edge
                Randdaten_Dirichlet[nRand_D_neu][0] = KnotenRD[1];
                Randdaten_Dirichlet[nRand_D_neu][1] = EdgeVertices.first;
                //increase # of boundary edges
                ++nRand_D_neu;
                //3rd third of the edge
                Randdaten_Dirichlet[nRand_D_neu][0] = EdgeVertices.first;
                Randdaten_Dirichlet[nRand_D_neu][1] = KnotenRD[0]; //the old right vertex

                // Erhoehen der Anzahl Randdaten:
                ++nRand_D_neu;
            }
        }

        // Updaten in der Parameter Struktur:
        Pars.NBdElemD = nRand_D_neu;
    }

    //    // ##############################################################
    //    // Neues Aufstellen der Neumann Randdaten:
    //    // ##############################################################
    KnotenRD[0] = KnotenRD[1] = 0;
    nRand_N_neu = Pars.NBdElemN;
    if (nRand_N_neu != 0) {
        for (uint i = 0; i < Pars.NBdElemN; ++i) {
            // Bestimmen von imin, imax;

            imin = min(Randdaten_Neumann[i][0], Randdaten_Neumann[i][1]);
            imax = max(Randdaten_Neumann[i][0], Randdaten_Neumann[i][1]);

            HashKey.first = imin;
            HashKey.second = imax;

            //retrieve values
            EdgeVertices = K.at(HashKey);
            if (imin == Randdaten_Neumann[i][0]) {
                //appropriate orientation
                //temporary storage
                KnotenRD[0] = EdgeVertices.first;
                KnotenRD[1] = Randdaten_Neumann[i][1];

                //1st third of the edge
                Randdaten_Neumann[i][1] = KnotenRD[0];
                //2nd third of the edge
                Randdaten_Neumann[nRand_N_neu][0] = KnotenRD[0];
                Randdaten_Neumann[nRand_N_neu][1] = EdgeVertices.second;
                //incease # of boundary edges
                ++nRand_N_neu;
                //3rd third of the edge
                Randdaten_Neumann[nRand_N_neu][0] = EdgeVertices.second;
                Randdaten_Neumann[nRand_N_neu][1] = KnotenRD[1];
                // Erhoehen der Anzahl Randdaten:
                ++nRand_N_neu;
            }
            else {
                //flipped orientation
                //temporary storage
                KnotenRD[0] = Randdaten_Neumann[i][1];
                KnotenRD[1] = EdgeVertices.second;

                //1st third of the edge
                Randdaten_Neumann[i][1] = KnotenRD[1];
                //2nd third of the edge
                Randdaten_Neumann[nRand_N_neu][0] = KnotenRD[1];
                Randdaten_Neumann[nRand_N_neu][1] = EdgeVertices.first;
                //incease # of boundary edges
                ++nRand_N_neu;
                //3rd third of the edge
                Randdaten_Neumann[nRand_N_neu][0] = EdgeVertices.first;
                Randdaten_Neumann[nRand_N_neu][1] = KnotenRD[0];
                // Erhoehen der Anzahl Randdaten:
                ++nRand_N_neu;
            }
        }

        // Updaten in der Parameter Struktur:
        Pars.NBdElemN = nRand_N_neu;
    }
}

void p3basepoints(geometry2& Pars) {
    /* Essentially this function is a copy of the ref2d function, except with a
     * fixed "refinement" degree */
    // Berechnen, wiegross die verfeinerten Strukturen werden:
    uint vorfaktor(3), Anz_Elemente, Anz_Knoten, Anz_Kanten, Anz_Rand_D, Anz_Rand_N;

    Anz_Elemente = Pars.Nelem;
    Anz_Rand_D = vorfaktor * Pars.NBdElemD;
    Anz_Rand_N = vorfaktor * Pars.NBdElemN; // Falls NeumannRand
    Anz_Kanten = (9 * Anz_Elemente + Anz_Rand_D) / 2;
    Anz_Kanten = (9 * Anz_Elemente + Anz_Rand_D + Anz_Rand_N) / 2; // Falls NeumannRand
    //2 nodes per edge + 1 node per element interior
    Anz_Knoten = 1 + Anz_Kanten - Anz_Elemente; //Euler formula for simply connected graphs with exterior domain 'removed'
    //Euler's formula is deficient here, as we have interior vertices 1 per element
    Anz_Knoten += Anz_Elemente;


    //reset # of elements to the old value
    Anz_Elemente = Pars.Nelem;
    //allocate memory for the refined mesh
    uint **Elemente, **Randdaten_D, **Randdaten_N;
    dtype **Koordinaten;

    Elemente = new uint*[Anz_Elemente];
    for (uint j(0); j < Anz_Elemente; ++j)
        //10 vertices per element
        Elemente[j] = new uint[10];

    Randdaten_D = new uint*[Anz_Rand_D];
    for (uint j(0); j < Anz_Rand_D; ++j)
        Randdaten_D[j] = new uint[2];

    Randdaten_N = new uint*[Anz_Rand_N];
    for (uint j(0); j < Anz_Rand_N; ++j)
        Randdaten_N[j] = new uint[2];

    Koordinaten = new dtype*[Anz_Knoten];
    for (uint j(0); j < Anz_Knoten; ++j)
        Koordinaten[j] = new dtype[2];

    //blank the allocated arrays
    for (uint j(0); j < Anz_Elemente; ++j)
        for (uint k(0); k < 10; ++k)
            Elemente[j][k] = 0;

    for (uint j(0); j < Anz_Rand_D; ++j)
        Randdaten_D[j][0] = Randdaten_D[j][1] = 0;
    for (uint j(0); j < Anz_Rand_N; ++j)
        Randdaten_N[j][0] = Randdaten_N[j][1] = 0;
    for (uint j(0); j < Anz_Knoten; ++j)
        Koordinaten[j][0] = Koordinaten[j][1] = 0.0;

    // ##############################################################
    // Uebertragen in Elemente und Koordinaten:
    for (uint j(0); j < Pars.Nvert; ++j) {
        Koordinaten[j][0] = Pars.vertices[j].x[0];
        Koordinaten[j][1] = Pars.vertices[j].x[1];
    }
    for (uint j(0); j < Pars.Nelem; ++j) {
        Elemente[j][0] = Pars.element[j].vert[0];
        Elemente[j][1] = Pars.element[j].vert[1];
        Elemente[j][2] = Pars.element[j].vert[2];
    }
    for (uint j(0); j < Pars.NBdElemD; ++j) {
        Randdaten_D[j][0] = Pars.boundaryDirichlet[j].vert[0];
        Randdaten_D[j][1] = Pars.boundaryDirichlet[j].vert[1];
    }
    for (uint j(0); j < Pars.NBdElemN; ++j) {
        Randdaten_N[j][0] = Pars.boundaryNeumann[j].vert[0];
        Randdaten_N[j][1] = Pars.boundaryNeumann[j].vert[1];
    }

    // ##############################################################
    // Durchfuehren der Verfeinerung:
    Edge2Split(Koordinaten, Elemente, Randdaten_D, Randdaten_N, Pars);

    // ##############################################################

    //feed back
    //copy back to the geometry
    //std::cout<< "Koordinaten"<<std::endl;

    Pars.resize(Pars.Nvert, Pars.Nelem, Pars.NBdElemD, Pars.NBdElemN);

    for (uint j(0); j < Pars.Nvert; ++j) {
        Pars.vertices[j].x[0] = Koordinaten[j][0];
        Pars.vertices[j].x[1] = Koordinaten[j][1];
    }
    for (uint j(0); j < Pars.Nelem; ++j) {
        //10 vertices!
        for (uint k(0); k < 10; ++k) {
            Pars.element[j].vert[k] = Elemente[j][k];
        }
    }
    if (Pars.NBdElemD != 0) {
        for (uint j(0); j < Pars.NBdElemD; ++j) {
            Pars.boundaryDirichlet[j].vert[0] = Randdaten_D[j][0];
            Pars.boundaryDirichlet[j].vert[1] = Randdaten_D[j][1];
        }
    }
    if (Pars.NBdElemN != 0) {
        for (uint j(0); j < Pars.NBdElemN; ++j) {
            Pars.boundaryNeumann[j].vert[0] = Randdaten_N[j][0];
            Pars.boundaryNeumann[j].vert[1] = Randdaten_N[j][1];
        }
    }

    //clean-up

    for (uint j(0); j < Anz_Elemente; ++j)
        delete[] Elemente[j];
    delete[] Elemente;

    for (uint j(0); j < Anz_Rand_D; ++j)
        delete[] Randdaten_D[j];
    delete[] Randdaten_D;

    for (uint j(0); j < Anz_Rand_N; ++j)
        delete[] Randdaten_N[j];
    delete[] Randdaten_N;

    for (uint j(0); j < Anz_Knoten; ++j)
        delete[] Koordinaten[j];
    delete[] Koordinaten;
}

//functions for the problem

dtype Phi_1(const uint i, const dtype& u, const dtype& v) {
    dtype res(0.0);
    switch (i) {
        case 0:
            res = 1.0 - (u + v);
            break;
        case 1:
            res = u;
            break;
        case 2:
            res = v;
            break;
        default:
            std::cerr << "Erroneous index in Phi. Is " << i << " should be 0..2" << std::endl;
    }
    return res;
}

void dPhi_1(const uint i, Eigen::Vector2d& x, const dtype& u, const dtype& v) {
    switch (i) {
        case 0:
            x(0) = -1.0;
            x(1) = -1.0;
            break;
        case 1:
            x(0) = 1;
            x(1) = 0;

            break;
        case 2:
            x(0) = 0;
            x(1) = 1;

            break;
        default:
            std::cerr << "Erroneous index in gradPhi. Is " << i << " should be 0,1,2" << std::endl;
    }
}

dtype Phi_2(const uint i, const dtype& u, const dtype& v) {
    dtype res(0.0);
    switch (i) {
        case 0:
            res = (1.0 - u - v)*(1.0 - 2.0 * u - 2.0 * v);
            //res = 1.0 - (u + v);
            break;
        case 1:
            res = (2.0 * u - 1.0) * u;
            //  res = u;
            break;
        case 2:
            res = (2.0 * v - 1.0) * v;
            //res = v;
            break;
        case 3:
            res = 4.0 * (1.0 - u - v) * u;
            break;
        case 4:
            res = 4.0 * u*v;
            break;
        case 5:
            res = 4.0 * (1.0 - u - v) * v;
            break;
        default:
            std::cerr << "Erroneous index in Phi. Is " << i << " should be 0..5" << std::endl;
    }
    return res;
}

void dPhi_2(const uint i, Eigen::Vector2d& x, const dtype& u, const dtype& v) {
    switch (i) {
        case 0:
            x(0) = -3. + 4. * u + 4. * v;
            x(1) = -3. + 4. * u + 4. * v;
            break;
        case 1:
            x(0) = -1. + 4. * u;
            x(1) = 0.0;
            break;
        case 2:
            x(0) = 0.0;
            x(1) = -1. + 4. * v;
            break;
        case 3:
            x(0) = -4. * (-1. + 2. * u + v);
            x(1) = -4. * u;
            break;
        case 4:
            x(0) = 4. * v;
            x(1) = 4. * u;
            break;
        case 5:
            x(0) = -4. * v;
            x(1) = -4. * (-1. + u + 2. * v);
            break;
        default:
            std::cerr << "Erroneous index in gradPhi. Is " << i << " should be 0,1,2,3,4,5" << std::endl;
    }
}

dtype
Phi_3(const uint i, const dtype& u, const dtype& v) {
    dtype res(0.0);
    switch (i) {
        case 0:
            res = (1.0 - u - v)*(1.0 - 3 * u - 3 * v)*(2.0 - 3 * u - 3 * v);
            break;
        case 1:
            res = (3 * u - 1.0) * u * (3 * u - 2.0);
            break;
        case 2:
            res = (3 * v - 1.0) * v * (3 * v - 2.0);
            break;
        case 3:
            res = 9 * (1.0 - u - v) * u * (2.0 - 3 * u - 3 * v);
            break;
        case 4:
            res = 9 * (1.0 - u - v) * u * (3 * u - 1.0);
            break;
        case 5:
            res = 9 * u * v * (3 * u - 1.0);
            break;
        case 6:
            res = 9 * u * v * (3 * v - 1.0);
            break;
        case 7:
            res = 9 * (1.0 - u - v) * v * (3 * v - 1.0);
            break;
        case 8:
            res = 9 * (1.0 - u - v) * v * (2.0 - 3 * u - 3 * v);
            break;
        case 9:
            res = 54 * u * v * (1.0 - u - v);
            break;
        default:
            std::cerr << "Erroneous index in Phi. Is " << i << " should be 0..9" << std::endl;
    }
    return res / 2.0;
}

void
dPhi_3(const uint i, Eigen::Vector2d& x, const dtype& u, const dtype& v) {
    switch (i) {
        case 0:
            x(0) = -11.0 + 36.0 * u - 27.0 * SQU(u) + 36.0 * v - 54.0 * u * v - 27.0 * SQU(v);
            x(1) = x(0);
            x /= 2;
            break;
        case 1:
            x(0) = 0.5 * (2.0 - 18.0 * u + 27.0 * SQU(u));
            x(1) = 0.0;
            break;
        case 2:
            x(0) = 0.0;
            x(1) = 0.5 * (2.0 - 18.0 * v + 27.0 * SQU(v));
            break;
        case 3:
            x(0) = 2.0 - 10.0 * u + 9.0 * SQU(u) - 5.0 * v + 12.0 * u * v + 3.0 * SQU(v);
            x(1) = -5.0 * u + 6.0 * SQU(u) + 6.0 * u*v;
            x *= 9.0 / 2.0;
            break;
        case 4:
            x(0) = 1.0 - 8.0 * u + 9.0 * SQU(u) - v + 6.0 * u*v;
            x(1) = -u + 3.0 * SQU(u);
            x *= -9.0 / 2.0;
            break;
        case 5:
            x(0) = -v + 6.0 * u*v;
            x(1) = -u + 3.0 * SQU(u);
            x *= 9.0 / 2.0;
            break;
        case 6:
            x(0) = -v + 3 * SQU(v);
            x(1) = -u + 6 * u*v;
            x *= 9.0 / 2.0;
            break;
        case 7:
            x(0) = -v + 3.0 * SQU(v);
            x(1) = 1.0 - u - 8.0 * v + 6.0 * u * v + 9.0 * SQU(v);
            x *= -9.0 / 2.0;
            break;
        case 8:
            x(0) = -5.0 * v + 6.0 * u * v + 6.0 * SQU(v);
            x(1) = 2.0 - 5.0 * u + 3.0 * SQU(u) - 10.0 * v + 12.0 * u * v + 9.0 * SQU(v);
            x *= 9.0 / 2.0;
            break;
        case 9:
            x(0) = -v + 2.0 * u * v + SQU(v);
            x(1) = -u + 2.0 * u * v + SQU(u);
            x *= -27.0;
            break;
        default:
            std::cerr << "Erroneous index in gradPhi. Is " << i << " should be 0..9" << std::endl;
    }
}

void Jacobian(const Tri& el, const coord<>* vertex, Eigen::Matrix2d& J) {
    J(0, 0) = vertex[ el.vert[1] ].x[0] - vertex[ el.vert[0] ].x[0];
    J(0, 1) = vertex[ el.vert[2] ].x[0] - vertex[ el.vert[0] ].x[0];

    J(1, 0) = vertex[ el.vert[1] ].x[1] - vertex[ el.vert[0] ].x[1];
    J(1, 1) = vertex[ el.vert[2] ].x[1] - vertex[ el.vert[0] ].x[1];
}

void InvJacobian(const Tri& el, const coord<>* vertex, Eigen::Matrix2d& J) {
    J(0, 0) = vertex[ el.vert[2] ].x[1] - vertex[ el.vert[0] ].x[1];
    J(1, 0) = vertex[ el.vert[2] ].x[0] - vertex[ el.vert[0] ].x[0];
    J(1, 0) *= -1.0;

    J(0, 1) = vertex[ el.vert[1] ].x[1] - vertex[ el.vert[0] ].x[1];
    J(0, 1) *= -1.0;
    J(1, 1) = vertex[ el.vert[1] ].x[0] - vertex[ el.vert[0] ].x[0];
    dtype detJ = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
    J /= detJ;
}

void g(const Tri& el, const coord<>* vertex, Eigen::Vector2d& x) {
    dtype u(0), v(0);
    u = x(0);
    v = x(1);
    x(0) = vertex[ el.vert[0] ].x[0];
    x(1) = vertex[ el.vert[0] ].x[1];

    x(0) += (vertex[ el.vert[1] ].x[0] - vertex[ el.vert[0] ].x[0]) * u;
    x(0) += (vertex[ el.vert[2] ].x[0] - vertex[ el.vert[0] ].x[0]) * v;

    x(1) += (vertex[ el.vert[1] ].x[1] - vertex[ el.vert[0] ].x[1]) * u;
    x(1) += (vertex[ el.vert[2] ].x[1] - vertex[ el.vert[0] ].x[1]) * v;
}

//1D mappings and base functions

dtype Phi1D(const uint &i, const dtype& x) {
    return ((i == 0) ? (0.5 * (1.0 - x)) : (0.5 * (1.0 + x)));
}

dtype dPhi1D(const uint &i) {
    return ( (i == 0) ? (-0.5) : (0.5));
}

//interval mapping

void g1D(const uint &i, const uint &j, const dtype &s, const coord<>* x, Eigen::Vector2d& y) {

    y(0) = (x[j].x[0] - x[i].x[0]);
    y(1) = (x[j].x[1] - x[i].x[1]);

    y *= s;
    y(0) += (x[j].x[0] + x[i].x[0]);
    y(1) += (x[j].x[1] + x[i].x[1]);
    y *= 0.5;
}

dtype dg1D(const uint &i, const uint &j, const coord<>* x) {

    dtype res = SQU(x[j].x[0] - x[i].x[0]) + SQU(x[j].x[1] - x[i].x[1]);
    return (0.5 * sqrt(res));
}

//------------------------------------------------------------------------------

void DetermineSparsityPattern(const geometry2 &geo, Eigen::VectorXi &NnzPerRow) {
    //pull # of rows
    uint Nrows = NnzPerRow.rows();
    uint FuncIdx(0);
    //clear
    NnzPerRow.setZero();
    //iterate over the Neumann boundary
    if (geo.NBdElemN != 0) {
        for (uint k(0); k < geo.NBdElemN; ++k) {//iterate over the boundary elements (lines)

            //skip the vertex if it is also in the Dirichlet boundary
            //left boundary vertex
            if (geo.isInDirichletBoundary[ geo.boundaryNeumann[k].vert[0] ] != true) {
                FuncIdx = geo.fIdx[ geo.boundaryNeumann[k].vert[0] ];
                NnzPerRow(FuncIdx) += 1;
            }
            //right boundary vertex
            if (geo.isInDirichletBoundary[ geo.boundaryNeumann[k].vert[1] ] != true) {
                FuncIdx = geo.fIdx[ geo.boundaryNeumann[k].vert[1] ];
                NnzPerRow(FuncIdx) += 1;
            }

        }
    }
    // scale by 2
    NnzPerRow /= 2;
    //iterate over all elements (triangles)
    for (uint k(0); k < geo.Nelem; ++k) {
        for (uint i(0); i < triVert; ++i) {//iterate over  the vertices in a boundary.
            //P1 elements -> 3 vertices in an element, P2 elements -> 6 vertices in an element
            FuncIdx = geo.fIdx[ geo.element[k].vert[i] ];
            //skip Dirichlet vertices
            if (geo.isInDirichletBoundary[ geo.element[k].vert[i] ] == true)
                continue;

            //else increment counter
            /* We'll have at most as many non-vanishing overlap integrals
             * per element as there are nodes in that triangle. 1 removed for
             * the function itself to not grossly overestimate the number of
             * non-vanishing elements.
             *  */
            NnzPerRow(FuncIdx) += triVert - 1;
        }
    }
    //the final increment - the function itself
    for (uint j(0); j < Nrows; ++j)
        NnzPerRow(j) += 1;

}

void MassMatrixSetup(CSRMat& M, const geometry2& geo, const dtype GD[][3], const uint QuadOrder,
                     const uint& freeNodes, const Eigen::VectorXi& NnzPerRow) {

    M.setZero();
    if (M.isCompressed() == true)
        M.uncompress();
#pragma omp parallel shared(M)
    {
        dtype jac = 0.0;
        dtype integral = 0.0;
        uint FuncIdx1 = 0;
        uint FuncIdx2 = 0;
        uint k(0);

        Eigen::Vector2d tmp;
        Eigen::Matrix2d J;
        CSRMat Mlocal(freeNodes, freeNodes);
        Mlocal.reserve(NnzPerRow);

        //temporary entry for filling
        dtype entry(0);

#pragma omp for nowait,schedule(static)
        for (k = 0; k < geo.Nelem; ++k) { //iteration over the elements(subintervals)
            Jacobian(geo.element[k], geo.vertices, J);
            //abosulte value of the determinant
            jac = fabs(J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1));

            for (uint i(0); i < triVert; ++i) { //iteration over the boundaries
                for (uint j = i; j < triVert; ++j) {
                    //fetch the values of the boundaries of the subinterval
                    if (geo.isInDirichletBoundary[geo.element[k].vert[i]] == true || geo.isInDirichletBoundary[geo.element[k].vert[j]] == true)
                        continue; //skip the boundaries of the domain - functions equal to nil

                    //determine function index
                    FuncIdx1 = geo.fIdx[ geo.element[k].vert[i] ];
                    FuncIdx2 = geo.fIdx[ geo.element[k].vert[j] ];
                    //integrate
                    integral = 0.0;

                    for (uint s(0); s < QuadOrder; ++s)
                        integral += GD[s][2] * Phi(i, GD[s][0], GD[s][1]) * Phi(j, GD[s][0], GD[s][1]);

                    integral *= jac; //transform the measure (constant)

                    //temporary entry - fill the real and imaginary parts of the matrix
                    // with the same values
                    entry = integral;

                    //diagonal and superdiagonal
                    Mlocal.coeffRef(FuncIdx1, FuncIdx2) += entry;
                    if (FuncIdx1 != FuncIdx2) {
                        //subdiagonal
                        Mlocal.coeffRef(FuncIdx2, FuncIdx1) += entry;
                    }
                }
            }
        }
        //finalize and update
        Mlocal.makeCompressed();
#pragma omp critical (Mupdate)
        {
            //std::cout << "thread " << omp_get_thread_num() << " updating matrix" << std::endl;
            M += Mlocal;
        }
#pragma omp barrier
        //end of omp parallel
    }

}

void StiffnessMatrixSetup(CSRMat& A, const geometry2& geo, const dtype GD[][3], const uint QuadOrder,
                          const uint& freeNodes, const Eigen::VectorXi& NnzPerRow) {
    A.setZero();
    if (A.isCompressed() == true)
        A.uncompress();

#pragma omp parallel shared(A)
    {
        dtype jac = 0.0;
        dtype integral = 0.0;
        uint FuncIdx1 = 0;
        uint FuncIdx2 = 0;
        uint k(0);
        Eigen::Vector2d tmp;
        Eigen::Vector2d dPhi1, dPhi2;
        Eigen::Matrix2d J, Jinv;
        A.setZero();
        //create a local copy of the matrix A to be filled
        CSRMat Alocal(freeNodes, freeNodes);
        Alocal.reserve(NnzPerRow);

        //temporary entry for filling
        dtype entry(0);

#pragma omp for nowait,schedule(static)
        for (k = 0; k < geo.Nelem; ++k) { //iteration over the elements(subintervals)
            Jacobian(geo.element[k], geo.vertices, J);
            //abosulte value of the determinant
            jac = fabs(J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1));
            //compute the inverse transposed jacobian matrix
            InvJacobian(geo.element[k], geo.vertices, Jinv);

            for (uint i(0); i < triVert; ++i) { //iteration over the boundaries
                for (uint j = i; j < triVert; ++j) {
                    //fetch the values of the boundaries of the subinterval
                    if (geo.isInDirichletBoundary[geo.element[k].vert[i]] == true || geo.isInDirichletBoundary[geo.element[k].vert[j]] == true)
                        continue; //skip the boundaries of the domain - functions equal to nil

                    //determine function index
                    FuncIdx1 = geo.fIdx[ geo.element[k].vert[i] ];
                    FuncIdx2 = geo.fIdx[ geo.element[k].vert[j] ];
                    //integrate
                    integral = 0.0;


                    for (uint s(0); s < QuadOrder; ++s) {
                        //g(geo.element[k], geo.vertices, q);
                        dPhi(i, tmp, GD[s][0], GD[s][1]);
                        dPhi1 = Jinv*tmp;
                        dPhi(j, tmp, GD[s][0], GD[s][1]);
                        dPhi2 = Jinv*tmp;
                        integral += GD[s][2] * dPhi1.dot(dPhi2);
                    }
                    integral *= jac; //transform the measure (constant)

                    //temporary entry - fill the real and imaginary parts of the matrix
                    // with the same values
                    entry = integral;


                    Alocal.coeffRef(FuncIdx1, FuncIdx2) += entry;
                    if (FuncIdx1 != FuncIdx2) {
                        //subdiagonal
                        Alocal.coeffRef(FuncIdx2, FuncIdx1) += entry;
                    }
                }
            }
        }
        //finalize and update global stiffness matrix
        Alocal.makeCompressed();
#pragma omp critical (Aupdate)
        {
            //std::cout << "thread " << omp_get_thread_num() << " updating matrix" << std::endl;
            A += Alocal;
        }
#pragma omp barrier
        //end of omp parallel
    }

}

void DampingMatrixSetup(CSRMat& C, const geometry2& geo, const dtype GD[][3], const uint QuadOrder,
                        const uint& freeNodes, const Eigen::VectorXi& NnzPerRow) {

    C.setZero();
    if (C.isCompressed() == true)
        C.uncompress();

#pragma omp parallel shared(C)
    {
        dtype jac(0);
        dtype integral(0);
        uint FuncIdx1(0);
        uint FuncIdx2(0);
        uint k(0);
        Eigen::Vector2d tmp;
        Eigen::Vector2d dPhi1, dPhi2;
        Eigen::Matrix2d J, Jinv;
        CSRMat Clocal(freeNodes, freeNodes);
        Clocal.reserve(NnzPerRow);

#pragma omp for nowait,schedule(static)
        for (k = 0; k < geo.Nelem; ++k) { //iteration over the elements(subintervals)
            Jacobian(geo.element[k], geo.vertices, J);
            //abosulte value of the determinant
            jac = fabs(J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1));
            //compute the inverse transposed jacobian matrix
            InvJacobian(geo.element[k], geo.vertices, Jinv);


            for (uint i(0); i < triVert; ++i) { //iteration over the boundaries
                for (uint j = i; j < triVert; ++j) {
                    //fetch the values of the boundaries of the subinterval
                    if (geo.isInDirichletBoundary[geo.element[k].vert[i]] == true || geo.isInDirichletBoundary[geo.element[k].vert[j]] == true)
                        continue; //skip the boundaries of the domain - functions equal to nil

                    //determine function index
                    FuncIdx1 = geo.fIdx[ geo.element[k].vert[i] ];
                    FuncIdx2 = geo.fIdx[ geo.element[k].vert[j] ];
                    //skip the diagonal due to antisymmetry
                    if (FuncIdx1 != FuncIdx2) {
                        //integrate
                        integral = 0.0;

                        //ordering gamess
                        if (FuncIdx1 > FuncIdx2) {
                            for (uint s(0); s < QuadOrder; ++s) {
                                //phi_i*dphi_j
                                dPhi(i, tmp, GD[s][0], GD[s][1]);
                                dPhi1 = Jinv*tmp;
                                //Map the values from reference triangle back to the domain
                                dPhi2(0) = GD[s][0];
                                dPhi2(1) = GD[s][1];
                                g(geo.element[k], geo.vertices, dPhi2);
                                //v = (-y,x) -> xd_y-yd_x
                                tmp(0) = -dPhi2(1);
                                tmp(1) = dPhi2(0);
                                //phi_j*dphi_i will be the same up to a sign!
                                integral += GD[s][2] * Phi(j, GD[s][0], GD[s][1]) * dPhi1.dot(tmp);
                            }
                            integral *= -1;
                        }
                        else {
                            for (uint s(0); s < QuadOrder; ++s) {
                                //phi_i*dphi_j
                                dPhi(j, tmp, GD[s][0], GD[s][1]);
                                dPhi1 = Jinv*tmp;
                                //Map the values from reference triangle back to the domain
                                dPhi2(0) = GD[s][0];
                                dPhi2(1) = GD[s][1];
                                g(geo.element[k], geo.vertices, dPhi2);
                                //v = (-y,x) -> xd_y-yd_x
                                tmp(0) = -dPhi2(1);
                                tmp(1) = dPhi2(0);
                                //phi_j*dphi_i will be the same up to a sign!
                                integral += GD[s][2] * Phi(i, GD[s][0], GD[s][1]) * dPhi1.dot(tmp);
                            }
                        }

                        integral *= jac; //transform the measure (constant)



                        //superdiagonal
                        Clocal.coeffRef(FuncIdx1, FuncIdx2) += integral;
                        //subdiagonal phi_j*dphi_i = -phi_i*dphi_j
                        Clocal.coeffRef(FuncIdx2, FuncIdx1) -= integral;

                    }
                }
            }
        }
        //finalize and update
        Clocal.makeCompressed();
#pragma omp critical (Cupdate)
        {
            //std::cout << "thread " << omp_get_thread_num() << " updating matrix" << std::endl;
            C += Clocal;
        }
#pragma omp barrier
        //end of omp parallel
    }

}

void shift(geometry2 &geo, const dtype& x, const dtype& y) {
    uint k(0);
    for (k = 0; k < geo.Nvert; ++k) {
        geo.vertices[k].x[0] += x;
        geo.vertices[k].x[1] += y;
    }
}

inline dtype radius(const Eigen::Vector2d& x) {
    return x.norm();
}

inline dtype gamma(const dtype& omega, const Eigen::Vector2d& x) {
    return 1.0 / std::sqrt(1.0 - SQU(x.norm() * omega / c));
}

inline dtype gammaSquared(const dtype& omega, const Eigen::Vector2d& x) {
    return 1.0 / (1.0 - SQU(x.norm() * omega / c));
}

void MassMatrixSetup_Full(CSRMat& M, const geometry2& geo, const dtype GD[][3], const uint QuadOrder,
                          const uint& freeNodes, const Eigen::VectorXi& NnzPerRow, const dtype& omega, const dtype& IOR) {

    M.setZero();
    if (M.isCompressed() == true)
        M.uncompress();
#pragma omp parallel shared(M)
    {
        dtype jac(0.0);
        dtype integral(0.0);
        dtype gammaFactor(0);
        uint FuncIdx1(0);
        uint FuncIdx2(0);
        uint k(0);

        Eigen::Vector2d x;
        Eigen::Matrix2d J;
        CSRMat Mlocal(freeNodes, freeNodes);
        Mlocal.reserve(NnzPerRow);

#pragma omp for nowait,schedule(static)
        for (k = 0; k < geo.Nelem; ++k) { //iteration over the elements(subintervals)
            Jacobian(geo.element[k], geo.vertices, J);
            //abosulte value of the determinant
            jac = fabs(J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1));

            for (uint i(0); i < triVert; ++i) { //iteration over the boundaries
                for (uint j = i; j < triVert; ++j) {
                    //fetch the values of the boundaries of the subinterval
                    if (geo.isInDirichletBoundary[geo.element[k].vert[i]] == true || geo.isInDirichletBoundary[geo.element[k].vert[j]] == true)
                        continue; //skip the boundaries of the domain - functions equal to nil

                    //determine function index
                    FuncIdx1 = geo.fIdx[ geo.element[k].vert[i] ];
                    FuncIdx2 = geo.fIdx[ geo.element[k].vert[j] ];
                    //integrate
                    integral = 0.0;

                    for (uint s(0); s < QuadOrder; ++s) {
                        //map the nodes to the triangle
                        x(0) = GD[s][0];
                        x(1) = GD[s][1];
                        g(geo.element[k], geo.vertices, x);
                        //compute the gamma factor
                        gammaFactor = (SQU(IOR) - SQU(radius(x) * omega / c)) * gammaSquared(omega, x);
                        //add to the integral
                        integral += GD[s][2] * gammaFactor * Phi(i, GD[s][0], GD[s][1]) * Phi(j, GD[s][0], GD[s][1]);
                    }

                    integral *= jac; //transform the measure (constant)


                    //diagonal and superdiagonal
                    Mlocal.coeffRef(FuncIdx1, FuncIdx2) += integral;
                    if (FuncIdx1 != FuncIdx2) {
                        //subdiagonal
                        Mlocal.coeffRef(FuncIdx2, FuncIdx1) += integral;
                    }
                }
            }
        }
        //finalize and update
        Mlocal.makeCompressed();
#pragma omp critical (Mupdate)
        {
            //std::cout << "thread " << omp_get_thread_num() << " updating matrix" << std::endl;
            M += Mlocal;
        }
#pragma omp barrier
        //end of omp parallel
    }
}

void StiffnessMatrixSetup_Auxiliary(CSRMat& A, const geometry2& geo, const dtype GD[][3], const uint QuadOrder,
                                    const uint& freeNodes, const Eigen::VectorXi& NnzPerRow) {

    A.setZero();
    if (A.isCompressed() == true)
        A.uncompress();

#pragma omp parallel shared(A)
    {
        dtype jac = 0.0;
        dtype integral = 0.0;
        uint FuncIdx1 = 0;
        uint FuncIdx2 = 0;
        uint k(0);
        Eigen::Vector2d x, alpha;
        Eigen::Vector2d dPhi1, dPhi2;
        Eigen::Matrix2d J, Jinv;


        //create a local copy of the matrix A to be filled
        CSRMat Alocal(freeNodes, freeNodes);
        Alocal.reserve(NnzPerRow);

#pragma omp for nowait,schedule(static)
        for (k = 0; k < geo.Nelem; ++k) { //iteration over the elements(subintervals)
            Jacobian(geo.element[k], geo.vertices, J);
            //abosulte value of the determinant
            jac = fabs(J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1));
            //compute the inverse transposed jacobian matrix
            InvJacobian(geo.element[k], geo.vertices, Jinv);

            for (uint i(0); i < triVert; ++i) { //iteration over the boundaries
                for (uint j = i; j < triVert; ++j) {
                    //fetch the values of the boundaries of the subinterval
                    if (geo.isInDirichletBoundary[geo.element[k].vert[i]] == true || geo.isInDirichletBoundary[geo.element[k].vert[j]] == true)
                        continue; //skip the boundaries of the domain - functions equal to nil

                    //determine function index
                    FuncIdx1 = geo.fIdx[ geo.element[k].vert[i] ];
                    FuncIdx2 = geo.fIdx[ geo.element[k].vert[j] ];
                    //integrate
                    integral = 0.0;


                    for (uint s(0); s < QuadOrder; ++s) {

                        //compute the d_phi flow
                        x(0) = GD[s][0];
                        x(1) = GD[s][1];
                        g(geo.element[k], geo.vertices, x);
                        alpha(0) = -x(1);
                        alpha(1) = x(0);
                        //compute the gradients
                        dPhi(i, x, GD[s][0], GD[s][1]);
                        dPhi1 = Jinv*x;

                        dPhi(j, x, GD[s][0], GD[s][1]);
                        dPhi2 = Jinv*x;
                        integral += GD[s][2] * alpha.dot(dPhi1) * alpha.dot(dPhi2);
                    }
                    integral *= jac; //transform the measure (constant)


                    Alocal.coeffRef(FuncIdx1, FuncIdx2) += integral;
                    if (FuncIdx1 != FuncIdx2) {
                        //subdiagonal
                        Alocal.coeffRef(FuncIdx2, FuncIdx1) += integral;
                    }
                }
            }
        }
        //finalize and update global stiffness matrix
        Alocal.makeCompressed();
#pragma omp critical (Aupdate)
        {
            //std::cout << "thread " << omp_get_thread_num() << " updating matrix" << std::endl;
            A += Alocal;
        }
#pragma omp barrier
        //end of omp parallel
    }

}

dtype L2FuncErr(const cVec& vals, const geometry2& geo, const dtype GD[][3], const uint QuadOrder) {
    dtype NrmSq(0), integral(0), jac(0);
    cdtype arg(0);
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
            //|psi|^2
            integral += GD[s][2] * SQU( std::abs(arg) );
        }
        integral *= jac;
        //add to the total
        NrmSq += integral;
    }
    //return the ||f||_L^2 = sqrt( ||f||^2 )
    return std::sqrt(NrmSq);
}

void saveEV(const std::string Fileprefix, const dtype omega, const uint Neig, const geometry2& geo, cdtype** const EV,
            const bool writeMesh, const uint freeNodes, const dtype GD[][3], const uint QuadOrder) {
    std::ofstream savefile;
    std::string Fname, tmpString;
    

    if(writeMesh == true){
        Fname = "Plot_Vertices_" + Fileprefix + ".dat";
        
        savefile.open(Fname.c_str());
        //write the vertices of the mesh
        if (savefile.is_open()) {
            for (uint i(0); i < geo.Nvert; ++i) {
                savefile << geo.vertices[i].x[0] << '\t' << geo.vertices[i].x[1] << std::endl;
            }
            //close the file
            savefile.close();
        }

        Fname = "Plot_Elements_" + Fileprefix + ".dat";
        savefile.open(Fname.c_str());

        //write the elements
        if (savefile.is_open()) {
            for (uint i(0); i < geo.Nelem; ++i)
                savefile << geo.element[i].vert[0] << '\t' << geo.element[i].vert[1] << '\t' << geo.element[i].vert[2] << std::endl;

            savefile.close();
        }
    }

    
    //helper vector
    cVec eigenmode(freeNodes);
    dtype L2Norm(0);
    //Normalize the eigenvectors w.r.t. the L^2 norm
    for(uint i(0); i < Neig; ++i){
        //copy the entries into the vector
        eigenmode.setZero();
        for(uint j(0); j<freeNodes; ++j)
            eigenmode(j) = EV[j][i];
        
        //normalize - note that the last two parameters are defined globally in GaussQuadParam.h
        L2Norm = L2FuncErr(eigenmode, geo, GD, QuadOrder);
        eigenmode /= L2Norm;
        //copy back
        for(uint j(0); j<freeNodes; ++j)
             EV[j][i] = eigenmode(j);
    }
    
    
    
    std::stringstream converter;
    converter.setf(ios::fixed);
    converter.precision(5);
    converter << omega;
    converter >> tmpString;
    
    
    Fname = "Eigenvalues_" + Fileprefix + "_omega_"+tmpString+".dat";
    //finally store the eigenvalues
    savefile.open(Fname.c_str());

    uint idx(0), j(0);
    if (savefile.is_open()) {
        for (uint i(0); i < geo.Nvert; ++i) {
            //homogeneous dirichlet boundary
            if (geo.isInDirichletBoundary[ i ]) {
                for (j = 0; j < Neig; ++j)
                    savefile << 0 << '\t';
                savefile << std::endl;
                //next vertex
                continue;
            }

            idx = geo.fIdx[i];
            for (j = 0; j < Neig; ++j) 
                savefile << std::setprecision(DBL_DIG) << EV[idx][j] << '\t';
            savefile << std::endl;
            
        }
        
        savefile.close();
    }

}