#include "Utils.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include "PolyhedralMesh.hpp"

using namespace std;

/*double EdgeLength(const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2)
{
	return sqrt(pow(x2 - x1 , 2) + pow(y2 - y1 , 2) + pow(z2 - z1, 2));	
}

// ****************************************************************************
*/

namespace PolyhedralLibrary{


bool PolyhedralChoice( string& path, 
						PolyhedralMesh& mesh, 
						const char& p, 
						const char& q, 
						const char& b, 
						const char& c,
						bool& walk)
{
	// string ans;  serve??
	string polyhedron;
	string filePath;
	
	if (p=='3')
	{
		switch(q)
		{
			case '3':
				polyhedron = "/Tetraedro";
				break;
			case '4':
				polyhedron = "/Ottaedro";
				break;
			case '5':
				polyhedron = "/Icosaedro";
				break;
			default:
				return false;
		}
		
		path = path + polyhedron;
		//ImportMesh(filePath, mesh);
		
		return true;
	}
	
	
	return false; 
}
//*******************************************************************************************
	/*verifica che tutti i file siano aperti*/
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

bool ParameterSelection(string& path, PolyhedralMesh& mesh)
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
		if (!PolyhedralChoice(path, mesh, p, q, b, c, walk))
			return false;
	}
	else if(Id1 != 'n' && Id2 != 'n')
	{
		walk = true;
		if (!PolyhedralChoice(path, mesh, p, q, b, c, walk))
			return false;
	}
	else
	{
		cerr << "Error: data not valid" << endl;
		return false;
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

// ***************************************************************************

bool GeodeticPolyhedron(const PolyhedralMesh& Platonic, PolyhedralMesh& Geodetic, const int& segments)
{
    int points_id = 0;   // inizializziamo a zero tutti gli id del poliedro
    int edges_id = 0;
    int faces_id = 0;
	int duplicate_id = 0;

	// Calcolo dei punti da generare per ogni faccia del poliedro
    int total_points = Platonic.NumCell2Ds * ((segments + 1) * (segments + 2) / 2);/*somma Gaussiana: 
	in poche parole, la formula corrisponde alla somma dei numeri da 1 a segments +1, ove segments=b. */

	// allochiamo spazio in memoria per gli Id e le coordinate delle celle 0D del poliedro geodetico.
    Geodetic.Cell0DsId.reserve(total_points);
    Geodetic.Cell0DsCoordinates = MatrixXd::Zero(3, total_points);

    map<array<int, 4>, int> point_coefficients;/*dichiara una mappa chiamata coefficients dove:
	la chiave è un array di 4 interi (std::array<int, 4>), mentre il valore associato a ogni chiave è un intero*/

	// Numero max di spigoli da generare (si considera il caso con il numero massimo, cioè l'icosaedro)
    int total_edges = 30 * segments * segments;
    Geodetic.Cell1DsId.reserve(total_edges);
    Geodetic.Cell1DsVertices.reserve(total_edges);
	
	// Numero max di facce triangolari da generare (si considera il caso con il numero massimo, cioè l'icosaedro)
    int total_faces = 20 * segments * segments;
    Geodetic.Cell2DsId.reserve(total_faces);
    Geodetic.Cell2DsVertices.resize(total_faces);
    Geodetic.Cell2DsEdges.resize(total_faces);


	// Ciclo su ogni faccia del poliedro
    for (const auto& j : Platonic.Cell2DsId)
    {
		// Estrazione dei 3 vertici dalla faccia corrente (tramite l'estrazione della j-esima colonna) e salvataggio di essi in vettori dinamici
        Vector3d Vert1 = Platonic.Cell0DsCoordinates.col(Platonic.Cell2DsVertices[j][0]);
        Vector3d Vert2 = Platonic.Cell0DsCoordinates.col(Platonic.Cell2DsVertices[j][1]);
        Vector3d Vert3 = Platonic.Cell0DsCoordinates.col(Platonic.Cell2DsVertices[j][2]);


		// Genera i punti interni alla faccia con suddivisione baricentrica
        for (int i = 0; i <= segments; i++)   /* per ale e fede: ho cambiato da ++i a ++j (semanticamente è la stessa cosa, ++i è più efficiente per tipi complessi, 
												ma con vettori e int non cambia niente, quindi teniamo i piedi per terra :) !!!!!!!DA CANCELLARE!!!!!!! */
        {
            for (int j = 0; j <= i; j++)	/* ciclo in i: per ogni i da 0 a b genera i punti interni */
            {
                int a = segments - i;  /* Dato un triangolo con vertici: Vert1, Vert2, Vert3, qualsiasi punto P nel piano del triangolo può essere scritto come: 
										P = a * Vert1 + b * Vert2 + c*Vert3, dove i pesi sommano a 1, e se il punto coincide con uno dei tre vertici, allora il 
										peso di tale vertice è 1 e gli altri due sono 0 (ad esempio, la prima iterazione del ciclo ha a=segments, b=c=0, dunque 
										stiamo generando Vert1. */      
				int b = i - j;
                int c = j;

				/* Calcolo del punto tramite combinazione lineare (a e b diventano double per evitare problemi con la divisione); dividiamo per segments in modo da avere
				la normalizzazione (infatti a, b e c sommano a segments, mentre dovrebbero sommare a 1) */
                Vector3d Point = (double(a) / segments) * Vert1 + (double(b) / segments) * Vert2 + (double(c) / segments) * Vert3;

				// Chiave della mappa: coefficienti + ID della faccia
                array<int, 4> coeffs;
				coeffs[0] = a;
				coeffs[1] = b;
				coeffs[2] = c;
				coeffs[3] = Id;

				/* Il procedimento implementato non garantisce che non ci sia sovrapposizione di punti: per questo, si richiama la funzione CheckDuplicatesVertex */ 
                if (!CheckDuplicatesVertex(Geodetic.Cell0DsCoordinates, Point, points_id, duplicate_id))
                {
                    coefficients[coeffs] = points_id;  // assegno all'array l'id corrente (all'inizio, points_id è inizializzato a 0) 
                    Geodetic.Cell0DsId.push_back(points_id);   // aggiungo l'array in coda agli id del poliedro
                    Geodetic.Cell0DsCoordinates.push_back(Point);  // aggiungo il punto in coda al poliedro geodetico (si può anche pensare di inserire le coordinate una per una?)
                    Geodetic.NumCell0Ds++;   // incremento il contatore di celle 0D del poliedro geodetico e quello dei punti
                    points_id++;
                }
                else
                {
                    coefficients[coeffs] = duplicate_id;  // lo aggiungiamo tra i punti duplicati
                }
            }
        }
    }

	// Costruzione dei triangoli all'interno di ogni faccia originale
    for (const auto& id : Platonic.Cell2DsId)
    {
        for (int i = 0; i < segments; i++)
        {
            for (int j = 0; j < segments - i; j++)
            {
                // Triangolo “a punta in su”
                int v1 = coefficients[{segments - i - j, i, j, id}];       // richiamiamo tre punti adiacenti per fare la triangolazione: nella prima iterazione, 
				// v1 è coefficients[{segments, 0, 0, id}], v2 è coefficients[{segments -1, 0, 1, id}], v3 è coefficients[{segments -1, 1, 0, id}]
                int v2 = coefficients[{segments - i - (j + 1), i, j + 1, id}];
                int v3 = coefficients[{segments - (i + 1) - j, i + 1, j, id}];

				// generazione delle facce: si parte da queste perché così dopo gli spigoli sanno già dove attaccarsi senza il rischio di duplicati
                Geodetic.NumCell2Ds++;    // geodetic è una struct di tipo PolyhedralMesh, dunque ha NumCell2Ds (inizialmente uguale a 0)
				Geodetic.Cell2DsId.push_back(faces_id);
				Geodetic.Cell2DsNumVertices[faces_id] = 3;  // inizializziamo le facce a triangolari, così dopo cicleremo direttamente su quelle
				Geodetic.Cell2DsNumEdges[faces_id] = 3;
				vector<int> Vertices = {Vertex1, Vertex2, Vertex3};
				Geodetic.Cell2DsVertices[faces_id] = Vertices;
				Geodetic.Cell2DsEdges[faces_id].resize(3);  // potevamo già crearle di dimensione 3? secondo ale si
				
				//crezione delle facce della triangolazione
                for (int k = 0; k < 3; k++)
                {
                    int v_start = Geodetic.Cell2DsVertices[faces_id][k];
                    int v_end;
					
					if (k == 2)
						v_end = Geodetic.Cell2DsVertices[faces_id][0];  // impongo che il terzo spigolo si ricolleghi al primo: condizione di chiusura della faccia verificata.
					else
						v_end = Geodetic.Cell2DsVertices[faces_id][k+1];  // altrimenti, lo spigolo legato alla k-esima faccia si ricollega al vertice iniziale della faccia successiva
					
                    if (!CheckDuplicatesEdge(Geodetic.Cell1DsVertices, v_start, v_end, edges_id))   // check se lo spigolo non è già esistente
                    {
                        Geodetic.Cell1DsId.push_back(edges_id);
                        Geodetic.Cell1DsVertices.push_back(Vector2i(v_start, v_end));
						GeodeticSolid.Cell2DsEdges[faces_id][k] = edges_id;  // andiamo a inserire gli id dei vertici
                        Geodetic.NumCell1Ds++;
                        edges_id++;
                    }
                }

                faces_id++;

                // Triangolo “a punta in giù”
                if (i > 0)     // per i=0 non funziona: i triangoli a punta in giù non toccano i vertici della faccia che stanno triangolando, ma solo punti interni dei lati
                {
                    int Vert4 = coefficients[{num_segments - (i - 1) - j, i - 1, j, static_cast<int>(id)}];
					/*i - 1: stai “salendo” di una riga nella griglia triangolare
					j: stessa colonna
					num_segments - (i - 1) - j: mantiene la condizione a + b + c = num_segments
					In pratica, Vert4 è il vertice interno che chiude il triangolo capovolto (v. sotto*/
					Geodetic.NumCell2Ds++;
                    Geodetic.Cell2DsId.push_back(faces_id);
                    Geodetic.Cell2DsNumVertices[faces_id] = 3;
					Geodetic.Cell2DsNumEdges[faces_id] = 3;
					Vertices[2] = Vert4;
					Geodetic.Cell2DsVertices[faces_id] = Vertices;
                    Geodetic.Cell2DsEdges[faces_id].resize(3);

                    for (int k = 0; k < 3; k++)   // stessa struttura di prima
                    {
                        int v_start = Geodetic.Cell2DsVertices[faces_id][k];
                        int v_end = Geodetic.Cell2DsVertices[faces_id][(k + 1) % 3];
						
						if ( k == 2 )
							v_end = Geodetic.Cell2DsVertices[faces_id][0];
						else
							v_end = Geodetic.Cell2DsVertices[faces_id][k+1];
						
                        if (!CheckDuplicatesEdge(Geodetic.Cell1DsVertices, v_start, v_end, edges_id))
                        {
                            Geodetic.Cell1DsId.push_back(edges_id);
                            Geodetic.Cell1DsVertices.push_back(Vector2i(v_start, v_end));
							GeodeticSolid.Cell2DsEdges[faces_id][k] = edges_id;
                            Geodetic.NumCell1Ds++;
                            edges_id++;
                        }
                    }

                    Geodetic.NumCell2Ds++;
                    faces_id++;
                }
            }
        }
    }

    return true;
}

// ***************************************************************************

bool CheckDuplicatesVertex(const std::vector<Vector3d>& coords, const Vector3d& point, int current_id, int& duplicate_id)
{
    for (int i = 0; i < current_id; ++i)
    {
        if ((coords[i] - point).norm() < 1e-8)
        {
            duplicate_id = i;
            return true;
        }
    }
    return false;
}

// ***************************************************************************

bool CheckDuplicatesEdge(const std::vector<Vector2i>& edges, int v1, int v2, int& current_edge_id)
{
    for (int i = 0; i < current_edge_id; ++i)
    {
        int a = edges[i][0];
        int b = edges[i][1];
        if ((a == v1 && b == v2) || (a == v2 && b == v1))
        {
            return true;
        }
    }
    return false;
}

bool GenerateGoldbergClassI(int p, int q, int b, int c, PolyhedralMesh& GoldbergSolid) 
{
	return true;
}

// ***************************************************************************

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
