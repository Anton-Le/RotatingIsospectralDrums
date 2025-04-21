# Supporting Files #    {#supportfiles}

## The Parameter File ##

The file \ref parameters.txt contains a list of key-value pairs of parameters
for the program.
This file is used for convenience
because, generally, only one or two parameters change between simulations
and it would be a pain to supply all parameters for each simulation on the CLI.

Each key-value pair should be specified in a new line and none should occur twice.


## The Mesh Files ##

The mesh is encoded in three files. A vertex file (e.g. \ref ID1_1_Coord.dat )
containing a list of all the vertices of a mesh with the coordinates
of each vertex given in a separate line. The order of the vertices
in this file matters for their row doubles as their index.

The element file (e.g. \ref ID1_1_Elements.dat ) contains a list of index triplets.
Each triplet must be given in a separate line because the line number
doubles as the index of the triangular element.
A triplet consists of the indices of the vertices which define the element.

The boundary file (e.g. \ref ID1_1_Boundary.dat). Therein each line represents
a directed segment of the boundary. The order of the vertex indices defines
the orientation of the segment.
The segments of the Dirichlet and Neumann boundaries must form a closed boundary chain.



 \file parameters.txt

 The file containing a list of key-value pairs defining the 
 parameters of the simulation. This file is used for convenience
 because, generally, only one or two parameters change between simulations
 and it would be a pain to supply all parameters for each simulation on the CLI.

 \verbinclude parameters.txt

