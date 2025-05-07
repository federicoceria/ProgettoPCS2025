#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolygonalLibrary {

struct PolygonalMesh
{
    unsigned int NumCell0Ds;
    unsigned int NumCell1Ds;
    unsigned int NumCell2Ds;

    vector<unsigned int> Cell0DsId;
    vector<unsigned int> Cell1DsId;

    Eigen::MatrixXd Cell0DsCoordinates;
    Eigen::MatrixXi Cell1DsExtrema;

    vector<vector<unsigned int>> Cell2DsVertices;
	vector<vector<unsigned int>> Cell2DsEdges;

	//I use a map to store markers, which associates the marker (of type unsigned int) to the list of points that have that marker
    map<unsigned int, list<unsigned int>> Cell0DMarkers;
	map<unsigned int, list<unsigned int>> Cell1DMarkers;
	map<unsigned int, list<unsigned int>> Cell2DMarkers;
};

}