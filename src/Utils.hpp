#pragma once

#include "PolyhedralMesh.hpp"
#include "UCDUtilities.hpp"
#include <limits>

using namespace std;

namespace PolyhedralLibrary
{
bool ImportMesh(const string& path, PolyhedralMesh& mesh);

bool ImportCell0Ds(const string& path, PolyhedralMesh& mesh);

bool ImportCell1Ds(const string& path, PolyhedralMesh& mesh);

bool ImportCell2Ds(const string& path, PolyhedralMesh& mesh);

bool ImportCell3Ds(const string& path, PolyhedralMesh& mesh);

bool CheckLength(PolyhedralMesh& mesh);

bool CheckAreas(PolyhedralMesh& mesh);

bool CheckFaces(PolyhedralMesh& mesh);

bool CheckMarker0Ds(PolyhedralMesh& mesh);

bool CheckMarker1Ds(PolyhedralMesh& mesh);

bool ExpPoints(PolyhedralMesh& mesh, const string& FilePath);

bool ExpSegments(PolyhedralMesh& mesh, const string& FilePath);

}