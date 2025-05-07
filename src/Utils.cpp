#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/Eigen"

using namespace Eigen;
using namespace std;

namespace PolygonalLibrary
{
bool ImportMesh(PolygonalMesh& mesh)
{
	if(!ImportCell0Ds(mesh))
        return false;

    if(!ImportCell1Ds(mesh))
        return false;

    if(!ImportCell2Ds(mesh))
        return false;

	edges_length(mesh);
	check_area(mesh);

	return true;
} 

bool ImportCell0Ds(PolygonalMesh& mesh)
{
	ifstream file("Cell0Ds.csv");

    if(file.fail())
        return false;
	
	//I create a list of strings to put all the lines of the file
    list<string> listLines;

    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // I remove the header of the file from the list
    listLines.pop_front();
	
	// I count how many rows there are (equivalent to the number of 0D points present)
    mesh.NumCell0Ds = listLines.size();

    if (mesh.NumCell0Ds == 0)
    {
        cerr << "There is no cell 0D" << endl;
        return false;
    }
	
	// Allocate memory for the vector for all ID point
    mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
	// I create a matrix of zeros to save the coordinates of the points. The third line of the matrix that corresponds to a possible z coordinate will always be zero
    mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);

    for (string& line : listLines)
    {
		// I replace the ; with just spaces
		replace(line.begin(), line.end(), ';', ' ');
		// I convert the string to a stringstream
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
		
		// The matrix Cell0DsCoordinates has 3 rows and each id column represent a different point
        converter >>  id >> marker >> mesh.Cell0DsCoordinates(0, id) >> mesh.Cell0DsCoordinates(1, id);
        
		mesh.Cell0DsId.push_back(id);

		// If the marker doesn't belong to the mesh boundary, I store it in the map
        if(marker != 0)
		{
			// If I have seen this marker before, I memorize the new point to the list
            // If this is the first time I see this marker, I create a new list that contains this point
			auto it = mesh.Cell0DMarkers.find(marker);
			if(it != mesh.Cell0DMarkers.end())
				mesh.Cell0DMarkers[marker].push_back(id);
			else
				mesh.Cell0DMarkers.insert({marker, {id}});
		}

    }

    return true;
}

// Same procedure as the previous function
bool ImportCell1Ds(PolygonalMesh& mesh)
{
	ifstream file("Cell1Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    listLines.pop_front();

    mesh.NumCell1Ds = listLines.size();

    if (mesh.NumCell1Ds == 0)
    {
        cerr << "There is no cell 1D" << endl;
        return false;
    }

    mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
	// I create a matrix of zero with 2 rows and as many columns as there are edges
    mesh.Cell1DsExtrema = Eigen::MatrixXi::Zero(2, mesh.NumCell1Ds);

    for (string& line : listLines)
    {
		replace(line.begin(), line.end(), ';', ' ');
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
		
		// The origin and the end of the edge are put in the Cell1DsExtrema matrix. The columns indicate the edge ID
        converter >>  id >> marker >>  mesh.Cell1DsExtrema(0, id) >>  mesh.Cell1DsExtrema(1, id);
        mesh.Cell1DsId.push_back(id);

		if(marker != 0)
		{
			auto it = mesh.Cell1DMarkers.find(marker);
			if(it != mesh.Cell1DMarkers.end())
				mesh.Cell1DMarkers[marker].push_back(id);
			else
				mesh.Cell1DMarkers.insert({marker, {id}});
		}
        
    }

    return true;
}

bool ImportCell2Ds(PolygonalMesh& mesh)
{
    ifstream file;
    file.open("./Cell2Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    listLines.pop_front();

    mesh.NumCell2Ds = listLines.size();

    if (mesh.NumCell2Ds == 0)
    {
        cerr << "There is no cell 2D" << endl;
        return false;
    }

    mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds);

    for (string& line : listLines)
    {
		replace(line.begin(), line.end(), ';', ' ');
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        unsigned int numvertices;
        unsigned int numedges;

        converter >> id >> marker >> numvertices;

		if(marker != 0)
		{
			auto it = mesh.Cell2DMarkers.find(marker);
			if(it != mesh.Cell2DMarkers.end())
				mesh.Cell2DMarkers[marker].push_back(id);
			else
				mesh.Cell2DMarkers.insert({marker, {id}});
		}

		//I create a vector of vertices, allocate memory, read them one by one and add them to the vector finally save in the mesh
        vector<unsigned int> vector_vertices;
        vector_vertices.reserve(numvertices);
        for(unsigned int i=0; i < numvertices; i++)
        {   
			unsigned int vertex;
            converter >> vertex;
			vector_vertices.push_back(vertex);
        }

        mesh.Cell2DsVertices.push_back(vector_vertices);
 
		//Same procedure for the edges of each polygon
        vector<unsigned int> vector_edges;
        vector_edges.reserve(numedges);
        for(unsigned int j = 0; j < numedges; j++)
        {
			unsigned int edge;
            converter >> edge;
			vector_edges.push_back(edge);
        }

        mesh.Cell2DsEdges.push_back(vector_edges);

	}

    return true;
}

bool edges_length(PolygonalMesh& mesh)
{

	for(unsigned int i = 0; i < mesh.NumCell1Ds; i++)
	{
		//I take the indices of the 2 points that make up edge number i
		unsigned int origin_vertice = mesh.Cell1DsExtrema(0, i);
		unsigned int end_vertice = mesh.Cell1DsExtrema(1, i);
		
		//I get the coordinates of the point
		double X_origin = mesh.Cell0DsCoordinates(0, origin_vertice);
		double Y_origin = mesh.Cell0DsCoordinates(1, origin_vertice);
		double X_end = mesh.Cell0DsCoordinates(0, end_vertice);
		double Y_end = mesh.Cell0DsCoordinates(1, end_vertice);	
			
		//I use the euclidean distance formula
		double length = sqrt(pow(X_origin - X_end, 2) + pow(Y_origin - Y_end, 2));
		
		//I set a tollerance
		if(length < 1e-16)
		{
			cout << "Test failed: the edge " << i << ", has length equal to zero." << endl;
			return false;
		}
		
	}

	cout << "All edges have non-zero length" << endl;
	return true;
}

bool check_area(PolygonalMesh& mesh)
{
	//I iterate over each polygon
	for(unsigned int i = 0; i < mesh.NumCell2Ds; i++)
	{
		double area = 0.0;
		unsigned int n = mesh.Cell2DsVertices[i].size();
		//I iterate over all vertices of the polygon
		for(unsigned int j = 0; j < n; j++)
		{
			//I take two vertices of the polygon
			unsigned int index_v1 = mesh.Cell2DsVertices[i][j];
			unsigned int index_v2 = mesh.Cell2DsVertices[i][(j+1) % n];
			
			//I get the coordinates of the vertices
			double X_v1 = mesh.Cell0DsCoordinates(0, index_v1);
			double Y_v1 = mesh.Cell0DsCoordinates(1, index_v1);
			double X_v2 = mesh.Cell0DsCoordinates(0, index_v2);
			double Y_v2 = mesh.Cell0DsCoordinates(1, index_v2);
			
			area += (X_v1 * Y_v2) - ( X_v2 * Y_v1);
		}
		
		area = 0.5 * abs(area);
	
		//I set a tollerance
		if(area < 1e-16)
		{
			cout << "Test failed: the polygon" << i << " has equal to zero" << endl;
			return false;
		}

	}

	cout << "All polygons have non-zero area." << endl;
	return true;
}

}
