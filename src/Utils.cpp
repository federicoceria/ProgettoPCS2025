#include "Utils.hpp"
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

/*double EdgeLength(const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2)
{
	return sqrt(pow(x2 - x1 , 2) + pow(y2 - y1 , 2) + pow(z2 - z1, 2));	
}

// ****************************************************************************
*/

namespace PolyhedralLibrary
{
bool PolyhedralChoice(const string& path, 
						PolyhedralMesh& mesh, 
						const char& p, 
						const char& q, 
						const char& b, 
						const char& c,
						bool& walk)
{
	string ans;
	string polihedron;
	string filePath;
	
	if (p=='3')
	{
		switch(q)
		{
			case '3':
				polihedron = "/Tetraedro";
				break;
			case '4':
				polihedron = "/Ottaedro";
				break;
			case '5':
				polihedron = "/Icosaedro";
				break;
			default:
				return false;
		}
		
		filePath = path + polihedron;
		ImportMesh(filePath, mesh);
		
		// Nota: qui andrà inserita la funzione che costruisce il poliedro GEODETICO
		// fornendo la mesh, i parametri b e c (per la triangolazione) e il percorso;
		
		// if (q=='3') Goldberg: duale del poliedro geodetico.
		
		return true;
	}
	
	return true;
}

//*******************************************************************************************+
	
bool ImportMesh(const string& path, PolyhedralMesh& mesh)
{
	if(!ImportCell0Ds(path, mesh))
    {
        cerr << "File Cell0Ds.csv not found" << endl;
        return false;
    }

    if(!ImportCell1Ds(path, mesh))
    {
        cerr << "File Cell1Ds.csv not found" << endl;
        return false;
    }

    if(!ImportCell2Ds(path, mesh))
    {
        cerr << "File Cell2Ds.csv not found" << endl;
        return false;
    }

	return true;
}

//****************************************************************************

bool ImportVector(const string& path, PolyhedralMesh& mesh)
{
	vector<string> v;
	char p;
	char q;
	char b;
	char c;
	char Id1;
	char Id2;
	char tmp;
	bool walk = false;
	
	cout << "Insert each parameter when asked and push enter to confirm." << endl;
	cout << "Insert p: ";
	cin >> p;
	cout << "Insert q: ";
	cin >> q;
	cout << "Insert b: ";
	cin >> b;
	cout << "Insert c: ";
	cin >> c;
	
	cout << "If you don't want to evaluate the shortest path between two vertices, please digit 'n' and enter it in the next two insertions."; 
	cout << "Insert Id1: ";
	cin >> Id1;
	cout << "Insert Id2: ";
	cin >> Id2;
	
	if(Id1 == 'n' && Id2 == 'n')
	{
		PolyhedralChoice(path, mesh, p, q, b, c, walk);
	}
	else if(Id1 != 'n' && Id2 != 'n')
	{
		walk = true;
		PolyhedralChoice(path, mesh, p, q, b, c, walk);
	}
	else
	{
		cerr << "Error: data not valid" << endl;
	}
	return true;
}

// ***************************************************************************

bool ImportCell0Ds(const string& path, PolyhedralMesh& mesh)
{
	string filePath = path + "/Cell0Ds.csv";
	ifstream file0(filePath);

    if (!file0)
    {
        return false;
    }
    
	list<string> lines;
    string line;
    while(getline(file0,line))
    {
        lines.push_back(line);
    }
    lines.pop_front();
    
	mesh.NumCell0Ds = lines.size();
    mesh.Cell0DsCoordinates.reserve(mesh.NumCell0Ds);
    mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
	
    Vector3d Coordinates;
    char tmp;
    unsigned int Id = 0;
	
	for(const auto& l : lines)
	{
		stringstream ss(l);
		ss >> Id >> tmp >> Coordinates(0) >> tmp >> Coordinates(1) >> tmp >> Coordinates(2);
		mesh.Cell0DsId.push_back(Id) ;
		mesh.Cell0DsCoordinates.push_back(Coordinates);
	}
	
	
    return true; 
}

// ***************************************************************************

bool ImportCell1Ds(const string& path, PolyhedralMesh& mesh)
{
	string filePath = path + "/Cell1Ds.csv";
	ifstream file1(filePath);

    if (!file1)
    {
        return false;
    }
    
    list<string> lines;
    string line;
    while(getline(file1,line))
    {
        lines.push_back(line);
    }

    lines.pop_front();
    mesh.NumCell1Ds = lines.size();
    mesh.Cell1DsVertices.reserve(mesh.NumCell1Ds);
    mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
    
    Vector2i Vertices;
    char tmp;
    unsigned int Id;
	
	for(const auto& l : lines)
	{
		stringstream ss(l);
		ss >> Id >> tmp >> Vertices(0) >> tmp >> Vertices(1);
		
		mesh.Cell1DsId.push_back(Id);
		mesh.Cell1DsVertices.push_back(Vertices);
	}
	
    return true;
}

// ***************************************************************************

bool ImportCell2Ds(const string& path, PolyhedralMesh& mesh)
{
	string filePath = path + "/Cell2Ds.csv";
	ifstream file2(filePath);

    if (!file2)
    {
        return false;
    }
    
    list<string> lines;
    string line;
    while(getline(file2,line))
    {
        lines.push_back(line);
    }

    lines.pop_front();
    mesh.NumCell2Ds = lines.size();
    mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsId.reserve(mesh.NumCell2Ds);
    
    char tmp;
    unsigned int Id;
    
	for(const auto& l : lines)
	{
		stringstream ss(l);
		ss >> Id >> tmp;

		vector<unsigned int> Vertices(3);
		for(unsigned int i = 0; i < 3; i++)
		{
		    ss >> Vertices[i] >> tmp;
		}
		
		vector<unsigned int> Edges(3);
		for(unsigned int i = 0; i < 3; i++)
		{
			ss >> Edges[i] >> tmp;
		}
		
		mesh.Cell2DsVertices.push_back(Vertices);
		mesh.Cell2DsEdges.push_back(Edges);
		mesh.Cell2DsId.push_back(Id);
	}
		
	return true;
}
/*
// *******************************************************************************

bool CheckLength(PolyhedralMesh& mesh)
{
	for(size_t i = 0; i < mesh.NumCell1Ds; i++)
	{
		double length = EdgeLength(
		mesh.Cell0DsCoordinates[mesh.Cell1DsVertices[i][0]](0),		 
		mesh.Cell0DsCoordinates[mesh.Cell1DsVertices[i][0]](1),		 
		mesh.Cell0DsCoordinates[mesh.Cell1DsVertices[i][1]](0),		 
		mesh.Cell0DsCoordinates[mesh.Cell1DsVertices[i][1]](1));	
		
		if(length < mesh.epsilon)
		{
			cout << "Il lato: " << i << " ha lunghezza non valida\n";
			return false;
		}
	}
	
	return true;
}

// ***************************************************************************

bool CheckAreas(PolyhedralMesh& mesh)
{
	double x1;
	double x2;
	double y1;
	double y2;
	for (unsigned int i = 0; i < mesh.NumCell2Ds; i++)
	{
		vector<unsigned int> Edges = mesh.Cell2DsEdges[i];
		double Area = 0.0;
		for (unsigned int j = 0; j < Edges.size(); j++)
		{
			const unsigned int Vertex1 = mesh.Cell1DsVertices[Edges[j]][0];
			const unsigned int Vertex2 = mesh.Cell1DsVertices[Edges[j]][1];

			x1 = mesh.Cell0DsCoordinates[Vertex1](0);
			y1 = mesh.Cell0DsCoordinates[Vertex1](1);
			x2 = mesh.Cell0DsCoordinates[Vertex2](0);
			y2 = mesh.Cell0DsCoordinates[Vertex2](1);
			
			Area = Area + (x1*y2 - y1*x2); // Formula per calcolare l'area: 1/2 * abs(sum( x(i)*y(i+1) - x(i+1)*y(i) ))
		}
		Area = 0.5*abs(Area);

		if (Area < mesh.epsilon)
		{ 
			cout << "Il poligono: " << i << " ha area nulla" << endl;
		
			return false;
		}

	}
	return true;
}

// ***************************************************************************

bool CheckFaces(PolyhedralMesh& mesh)
{
	// da verificare faces.edges[e].end == faces.edges[(e+1)%E].origin;
	// da verificare faces.vertices[e] == faces.edges[e].origin;
	
	return true;
}

// ****************************************************************************

bool CheckMarker0Ds(PolyhedralMesh& mesh)
{
	double x;
	double y;
	
	for(const auto& i : mesh.Cell0DsId)
	{
		x = mesh.Cell0DsCoordinates[i](0);
		y = mesh.Cell0DsCoordinates[i](1);
		if(mesh.Cell0DsMarker[i] == 0)
		{
			if( abs(x-1.0) < mesh.epsilon || abs(y-1.0) < mesh.epsilon ||
				x < mesh.epsilon || y < mesh.epsilon)
			{
				cout << "Condizione non soddisfatta per i = " << i << endl;
				return false;
			}
		}
		else
		{
			if( !(abs(x-1.0) < mesh.epsilon || abs(y-1.0) < mesh.epsilon ||
				x < mesh.epsilon || y < mesh.epsilon))
			{
				cout << "Condizione non soddisfatta per i = " << i << endl;
				return false;
			}
		}
	}
	return true;
}

// ***************************************************************************

bool CheckMarker1Ds(PolyhedralMesh& mesh)
{
	unsigned int count = 0;
	
	for(const auto& i : mesh.Cell1DsVertices)
	{
		vector<double> x;
		vector<double> y;
		
		for(const auto& j : i)
		{
			x.push_back(mesh.Cell0DsCoordinates[j](0));
			y.push_back(mesh.Cell0DsCoordinates[j](1));
		}
			
		if(mesh.Cell1DsMarker[count] == 0)
		{
			if( (abs(x[0]-1.0) < mesh.epsilon || abs(y[0]-1.0) < mesh.epsilon || x[0] < mesh.epsilon || y[0] < mesh.epsilon) 
				&& 
				(abs(x[1]-1.0) < mesh.epsilon || abs(y[1]-1.0) < mesh.epsilon || x[1] < mesh.epsilon || y[1] < mesh.epsilon) )
			{
				cout << "Condizione non soddisfatta per i = " << count << endl;
				return false;
			}
		}
		else
		{
			if( !((abs(x[0]-1.0) < mesh.epsilon || abs(y[0]-1.0) < mesh.epsilon || x[0] < mesh.epsilon || y[0] < mesh.epsilon) 
				|| 
				!(abs(x[1]-1.0) < mesh.epsilon || abs(y[1]-1.0) < mesh.epsilon || x[1] < mesh.epsilon || y[1] < mesh.epsilon)))
			{
				cout << "Condizione non soddisfatta per i = " << count << endl;
				return false;
			}
		}
		
		count ++;
	}
	return true;
}

// ***************************************************************************

bool ExpPoints(PolyhedralMesh& mesh, const string& FilePath)
{
	mesh.Points.resize(3, mesh.NumCell0Ds);
	for(size_t i = 0; i < mesh.NumCell0Ds; i++)
	{
		mesh.Points(0,i) = mesh.Cell0DsCoordinates[i][0];
		mesh.Points(1,i) = mesh.Cell0DsCoordinates[i][1];
		mesh.Points(2,i) = 0.0;
	}
	
	VectorXi Materials0Ds(mesh.NumCell0Ds);
	for (size_t i = 0; i < mesh.NumCell0Ds; ++i)
	{
		Materials0Ds[i] = static_cast<int>(mesh.Cell0DsMarker[i]);
	}
	Gedim::UCDUtilities utilities;
	utilities.ExportPoints(FilePath, mesh.Points, {}, Materials0Ds);
	
	return true;
}

// ***************************************************************************

bool ExpSegments(PolyhedralMesh& mesh, const string& FilePath)
{
	mesh.Segments.resize(2, mesh.NumCell1Ds);
	for(size_t i = 0; i < mesh.NumCell1Ds; i++)
	{
		mesh.Segments(0,i) = mesh.Cell1DsVertices[i][0];
		mesh.Segments(1,i) = mesh.Cell1DsVertices[i][1];
	}
	
	VectorXi Materials1Ds(mesh.Cell1DsMarker.size());
	for (size_t i = 0; i < mesh.Cell1DsMarker.size(); ++i)
	{
		Materials1Ds[i] = static_cast<int>(mesh.Cell1DsMarker[i]);
	}
	
	Gedim::UCDUtilities utilities;
	utilities.ExportSegments(FilePath, mesh.Points, mesh.Segments, {}, {}, Materials1Ds);
	
	return true;
} */
}



/* mesh.NumCell0Ds = linesT.size();
    mesh.Cell0DsCoordinates.reserve(mesh.NumCell0Ds);
    mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
	
*/