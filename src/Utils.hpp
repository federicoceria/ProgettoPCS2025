#pragma once

#include "PolyhedralMesh.hpp"
#include "UCDUtilities.hpp"
#include <limits>

using namespace std;

namespace PolyhedralLibrary
{
bool ImportVector(const string& path, PolyhedralMesh& mesh);

bool PolyhedralChoice(const string& path,
						PolyhedralMesh& mesh,
						const char& p,
						const char& q, 
						const char& b,
						const char & c,
						bool& walk);

bool ImportMesh(const string& path, PolyhedralMesh& mesh);

bool ImportCell0Ds(const string& path, PolyhedralMesh& mesh);

bool ImportCell1Ds(const string& path, PolyhedralMesh& mesh);

bool ImportCell2Ds(const string& path, PolyhedralMesh& mesh);

bool GeodeticPolyhedron1(const PolyhedralMesh& PlatonicPolyhedron, PolyhedralMesh& GeodeticSolid, const int& num_segments);

bool CheckDuplicatesVertex(const std::vector<Vector3d>& coords, const Vector3d& point, int current_id, int& duplicate_id);

bool CheckDuplicatesEdge(const std::vector<Vector2i>& edges, int v1, int v2, int& current_edge_id);

bool GoldbergClassI(int p, int q, int b, int c, PolyhedralMesh& GoldbergSolid);

/*
bool CheckFaces(PolyhedralMesh& mesh);

bool CheckMarker0Ds(PolyhedralMesh& mesh);

bool CheckMarker1Ds(PolyhedralMesh& mesh);

bool ExpPoints(PolyhedralMesh& mesh, const string& FilePath);

bool ExpSegments(PolyhedralMesh& mesh, const string& FilePath);
*/
}