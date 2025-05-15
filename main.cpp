#include "Utils.hpp"
#include <iostream>
#include <sstream>
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace PolyhedralLibrary;

int main()
{
    string path = "/home/appuser/Data/progettoPCS2025/Platonic_solids/Icosaedro";
    
    // Mesh di partenza e quella generata
    PolyhedralMesh icosahedron;
    PolyhedralMesh geodetic;

    // Numero di segmenti da usare per la suddivisione
    int num_segments = 3;

    // Importa la mesh dell'icosaedro
    if (!ImportMesh(path, icosahedron)) {
        cerr << "Error importing the Icosahedron." << endl;
        return 1;
    }

    cout << "Icosahedron successfully imported." << endl;
    cout << "Number of vertices: " << icosahedron.NumCell0Ds << endl;
    cout << "Number of edges: " << icosahedron.NumCell1Ds << endl;
    cout << "Number of faces: " << icosahedron.NumCell2Ds << endl;

    // Genera il solido geodetico
    if (!GenerateGeodeticPolyhedronType1(icosahedron, geodetic, num_segments)) {
        cerr << "Error generating the geodetic solid." << endl;
        return 1;
    }

    cout << "Geodetic solid successfully generated with num_segments = " << num_segments << endl;
    cout << "Number of vertices: " << geodetic.NumCell0Ds << endl;
    cout << "Number of edges: " << geodetic.NumCell1Ds << endl;
    cout << "Number of faces: " << geodetic.NumCell2Ds << endl;

    return 0;
}

