#pragma once

#include "PolyhedralMesh.hpp"
//#include "UCDUtilities.hpp" questo file non è ancora presente nella src, quindi lo ho commentato
#include <limits>

using namespace std;

namespace PolyhedralLibrary
{
bool ParameterSelection(string& path, PolyhedralMesh& mesh);

bool PolyhedralChoice( string& path,
						PolyhedralMesh& mesh,
						const char& p,
						const char& q, 
						const char& b,
						const char& c,
						bool& walk);

bool ImportMesh( string& path, PolyhedralMesh& mesh);

bool ImportCell0Ds( string& path, PolyhedralMesh& mesh);

bool ImportCell1Ds( string& path, PolyhedralMesh& mesh);

bool ImportCell2Ds( string& path, PolyhedralMesh& mesh);

/*bool GeodeticPolyhedron(const PolyhedralMesh& PlatonicPolyhedron, PolyhedralMesh& GeodeticSolid, const int& num_segments);

bool CheckDuplicatesVertex(const vector<Vector3d>& coords, const Vector3d& point, int current_id, int& duplicate_id);

bool CheckDuplicatesEdge(const vector<Vector2i>& edges, int v1, int v2, int& current_edge_id);

bool GenerateGoldbergClassI(int p, int q, int b, int c, PolyhedralMesh& GoldbergSolid);

/*
bool CheckFaces(PolyhedralMesh& mesh);

bool CheckMarker0Ds(PolyhedralMesh& mesh);

bool CheckMarker1Ds(PolyhedralMesh& mesh);

bool ExpPoints(PolyhedralMesh& mesh, const string& FilePath);

bool ExpSegments(PolyhedralMesh& mesh, const string& FilePath);
*/
}