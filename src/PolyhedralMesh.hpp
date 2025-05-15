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
	vector<unsigned int> Cell0DsMarker = {}; 
	vector<Vector3d> Cell0DsCoordinates = {}; 
	
	
	// Celle 1D
    unsigned int NumCell1Ds = 0;
	vector<unsigned int> Cell1DsId = {};
	vector<unsigned int> Cell1DsMarker = {};
	vector<Vector2i> Cell1DsVertices = {};
	
	
	// Celle 2D
    unsigned int NumCell2Ds = 0;
	vector<unsigned int> Cell2DsId = {};
	vector<vector<unsigned int>> Cell2DsVertices = {};
	vector<vector<unsigned int>> Cell2DsEdges = {};
	
	
	// Celle 3D
	unsigned int NumCell3Ds = 0;
	vector<unsigned int> Cell3DsId = {};
	vector<vector<unsigned int>> Cell3DsVertices = {};  
	vector<vector<unsigned int>> Cell3DsEdges = {};
	vector<vector<unsigned int>> cell3DsFaces = {};
	
    // MatrixXd Points;
	// MatrixXi Segments;
	
	const double epsilon = 1.0e-8; 
};

}