#pragma once

#include "PolyhedralMesh.hpp"
//#include "UCDUtilities.hpp" 
#include <limits>

using namespace std;

namespace PolyhedralLibrary
{
bool ParameterSelection(string& path, PolyhedralMesh& mesh);

bool PolyhedralChoice(string& path,
						PolyhedralMesh& mesh,
						const char& p,
						const char& q, 
						const char& b,
						const char& c,
						bool& walk);

bool ImportMesh(const string& path, PolyhedralMesh& mesh);

bool ImportCell0Ds(const string& path, PolyhedralMesh& mesh);

bool ImportCell1Ds(const string& path, PolyhedralMesh& mesh);

bool ImportCell2Ds(const string& path, PolyhedralMesh& mesh);

bool GeodeticPolyhedron(const PolyhedralMesh& Platonic, PolyhedralMesh& Geodetic, const int& segments);

void GenerateTriangles(const PolyhedralMesh& Platonic, PolyhedralMesh& Geodetic, map<array<int, 4>, int>& coefficients, int segments, int& edges_id, int& faces_id);

void DualMesh(PolyhedralMesh& InputMesh, PolyhedralMesh& DualMesh);

void Sort_Faces(const vector<int>& UnsortedFaces, vector<int>& SortedFaces, const PolyhedralMesh& Mesh);

bool CheckVertices(const MatrixXd& mat, const Vector3d& vec, int& matSize, int& duplicate_id);
//bool CheckVertices(const vector<Vector3d>& coords, const Vector3d& point, int current_id, int& duplicate_id);

bool CheckEdges(const MatrixXi& mat, int& dimension, const int& v1, const int& v2, int& existing_edge_id);
// bool CheckEdges(const vector<Vector2i>& edges, int v1, int v2, int& current_edge_id, int& existing_edge_id);

void Projection(PolyhedralMesh& mesh);

/*
bool GenerateGoldbergClassI(int p, int q, int b, int c, PolyhedralMesh& Goldberg);

bool CheckFaces(PolyhedralMesh& mesh);

bool CheckMarker0Ds(PolyhedralMesh& mesh);

bool CheckMarker1Ds(PolyhedralMesh& mesh);

bool ExpPoints(PolyhedralMesh& mesh, const string& FilePath);

bool ExpSegments(PolyhedralMesh& mesh, const string& FilePath);
*/
}