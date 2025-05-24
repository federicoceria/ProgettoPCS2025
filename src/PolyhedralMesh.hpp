#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary {
// NOTA IMPORTANTE: MODIFICO I VETTORI DI VETTORI IN MATRICI, PER ELIMINATRE IL SEGFAULT. LASCIO IL VECCHIO CODICE COMMENTATO A FIANCO.
struct PolyhedralMesh
{
	// Celle 0D
    int NumCell0Ds = 0;
	vector<int> Cell0DsId = {};
	MatrixXd Cell0DsCoordinates = {}; // vector<Vector3d> Cell0DsCoordinates = {};
	
	
	// Celle 1D
    int NumCell1Ds = 0;
	vector<int> Cell1DsId = {};
	MatrixXi Cell1DsVertices = {}; //vector<Vector2i> Cell1DsVertices = {};
	
	
	// Celle 2D
    int NumCell2Ds = 0;
	vector<int> Cell2DsId = {};
	vector<int> Cell2DsNumVertices = {};
	vector<int> Cell2DsNumEdges = {};
	vector<vector<int>> Cell2DsVertices = {};  // vertici e spigoli delle celle bidimensionali
	vector<vector<int>> Cell2DsEdges = {}; 
};

}
