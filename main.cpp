#include "Utils.hpp"
#include <iostream>
#include <sstream>
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace PolyhedralLibrary;

int main()
{
    string path = "/home/appuser/Data/ProgettoPCS2025/Platonic_solids";

    // Mesh di partenza e quella generata
    PolyhedralMesh inputMesh;
    PolyhedralMesh geodetic;

    // Parametri per il solido geodetico
    int num_segments = 3;

    // Parametri inutilizzati
    char p = ' ', q = ' ', b = ' ', c = ' ';

    string solidType;
    cout << "Which solid do you want to generate? (tetrahedron, octahedron, icosahedron): ";
    cin >> solidType;

    if (!PolyhedralChoice(path, inputMesh, p, q, b, c, solidType)) {
        cerr << "Error importing selected solid." << endl;
        return 1;
    }

    if (!GeodeticPolyhedron(inputMesh, geodetic, num_segments)) {
        cerr << "Error generating the geodetic solid." << endl;
        return 1;
    }

    cout << "Geodetic solid successfully generated with num_segments = " << num_segments << endl;
    cout << "Number of vertices: " << geodetic.NumCell0Ds << endl;
    cout << "Number of edges: " << geodetic.NumCell1Ds << endl;
    cout << "Number of faces: "   << geodetic.NumCell2Ds << endl;

    return 0;
}




