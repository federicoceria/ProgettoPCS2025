#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary {

struct PolyhedralMesh
{
	// Celle 0D
    unsigned int NumCell0Ds = 0;
	vector<unsigned int> Cell0DsId = {};
	MatrixXd Cell0DsCoordinates = {}; 
	
	
	// Celle 1D
    unsigned int NumCell1Ds = 0;
	vector<unsigned int> Cell1DsId = {};
	MatrixXi Cell1DsVertices = {};
	
	
	// Celle 2D
    unsigned int NumCell2Ds = 0;
	vector<unsigned int> Cell2DsId = {};
	vector<unsigned int> Cell2DsNumVertices = {};
	vector<unsigned int> Cell2DsNumEdges = {};
	vector<vector<unsigned int>> Cell2DsVertices = {};  // vertici e spigoli delle celle bidimensionali
	vector<vector<unsigned int>> Cell2DsEdges = {};
		
	double epsilon = 1.0e-8; 
};

}
