#pragma once

#include "PolyhedralMesh.hpp"
#include <limits>

using namespace std;

namespace PolyhedralLibrary
{
bool ImportMesh(const string& path, PolyhedralMesh& mesh);

bool ImportCell0Ds(const string& path, PolyhedralMesh& mesh);

bool ImportCell1Ds(const string& path, PolyhedralMesh& mesh);

bool ImportCell2Ds(const string& path, PolyhedralMesh& mesh);

void GeodeticPolyhedron(const PolyhedralMesh& Platonic, PolyhedralMesh& Geodetic, const int& segments);

void GenerateTriangles(const PolyhedralMesh& Platonic, PolyhedralMesh& Geodetic, map<array<int, 4>, int>& coefficients, int segments, int& edges_id, int& faces_id);

void DualMesh(const PolyhedralMesh& StartPolyhedron, PolyhedralMesh& DualPolyhedron);

void Sort_Faces(const vector<int>& UnsortedFaces, vector<int>& SortedFaces, const PolyhedralMesh& Mesh);

bool CheckVertices(const MatrixXd& mesh, const Vector3d& point, int& dimension, int& duplicate_id);

bool CheckEdges(const MatrixXi& verts, const int& v1, const int& v2, int& dimension, int& existing_edge_id);

void Projection(PolyhedralMesh& mesh);

bool ExportPolyhedralData(const PolyhedralMesh& mesh);

bool ShortestPath(const PolyhedralMesh& mesh, const int& start, const int& end, vector<int>& path, MatrixXd& W);   // prima del vettore c'erano double& length, int& NumPath

bool isInteger(const string& str);

}