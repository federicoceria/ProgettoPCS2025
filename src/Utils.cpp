#include "Utils.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <queue>
#include "PolyhedralMesh.hpp"
#include <cctype>
#include <string>
#include "UCDUtilities.hpp"

using namespace std;

namespace PolyhedralLibrary
{


/************************************************************************************************/

// Verifica che tutti i file siano aperti correttamente, salva i dati nella mesh
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


/************************************************************************************************/


//se riesce ad prire il file .csv allora procede con il salvataggio dei dati all'interno della mesh
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
    mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
	mesh.Cell0DsCoordinates = MatrixXd::Zero(3,mesh.NumCell0Ds);
	
    char tmp;
    int Id;
	
	for(const auto& l : lines)
	{
		stringstream ss(l);
		ss >> Id >> tmp >> mesh.Cell0DsCoordinates(0,Id) >> tmp >> mesh.Cell0DsCoordinates(1,Id) >> tmp >> mesh.Cell0DsCoordinates(2,Id);
		mesh.Cell0DsId.push_back(Id);
	}
	
    return true; 
}

/************************************************************************************************/


//se riesce ad prire il file .csv allora procede con il salvataggio dei dati all'interno della mesh
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
    mesh.Cell1DsVertices = MatrixXi::Zero(2, mesh.NumCell1Ds); 
    mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
    
    char tmp;
    int Id;
	
	for(const auto& l : lines)
	{
		stringstream ss(l);
		ss >> Id >> tmp >> mesh.Cell1DsVertices(0,Id) >> tmp >> mesh.Cell1DsVertices(1,Id);
		// ss >> Id >> tmp >> Vertices(0) ...
		mesh.Cell1DsId.push_back(Id);
		// mesh.Cell1DsVertices.push_back(Vertices);
	}
	
    return true;
}

/************************************************************************************************/


//se riesce ad prire il file .csv allora procede con il salvataggio dei dati all'interno della mesh
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
    int Id;
    
	for(const auto& l : lines)
	{
		stringstream ss(l);
		ss >> Id >> tmp;

		vector<int> Vertices(3);
		for(int i = 0; i < 3; i++)
		{
		    ss >> Vertices[i] >> tmp;
		}
		
		vector<int> Edges(3);
		for(int i = 0; i < 3; i++)
		{
			ss >> Edges[i] >> tmp;
		}
		
		mesh.Cell2DsVertices.push_back(Vertices);
		mesh.Cell2DsEdges.push_back(Edges);
		mesh.Cell2DsId.push_back(Id);
	}
		
	return true;
}

/************************************************************************************************/


// esporta i dati contenuti all'interno della mesh
bool ExportPolyhedralData(const PolyhedralMesh& mesh)
{
    {
        ofstream vertexFile("Cell0Ds.txt");
        if (!vertexFile.is_open())
            return false;

        vertexFile << "Id;X;Y;Z" << endl;

        for (int i = 0; i < mesh.NumCell0Ds; i++) {
            vertexFile << mesh.Cell0DsId[i] << ";"
                       << mesh.Cell0DsCoordinates(0, i) << ";"
                       << mesh.Cell0DsCoordinates(1, i) << ";"
                       << mesh.Cell0DsCoordinates(2, i) << endl;
        }
    }

    {
        ofstream edgeFile("Cell1Ds.txt");
        if (!edgeFile.is_open())
            return false;

        edgeFile << "Id;Origin;End" << endl;

        for (int i = 0; i < mesh.NumCell1Ds; i++) {
            edgeFile << mesh.Cell1DsId[i] << ";"
                     << mesh.Cell1DsVertices(0, i) << ";"
                     << mesh.Cell1DsVertices(1, i) << endl;
        }
    }

    {
        ofstream faceFile("Cell2Ds.txt");
        if (!faceFile.is_open())
            return false;

        faceFile << "Id;NumVertices;VerticesId;NumEdges;EdgesId" << endl;

        for (int i = 0; i < mesh.NumCell2Ds; i++) {
            faceFile << mesh.Cell2DsId[i] << ";"
                     << mesh.Cell2DsNumVertices[i];

            for (int j = 0; j < mesh.Cell2DsNumVertices[i]; j++)
                faceFile << ";" << mesh.Cell2DsVertices[i][j];

            faceFile << ";" << mesh.Cell2DsNumEdges[i];

            for (int j = 0; j < mesh.Cell2DsNumEdges[i]; j++)
                faceFile << ";" << mesh.Cell2DsEdges[i][j];

            faceFile << endl;
        }
    }

    {
        ofstream polyhedronFile("Cell3Ds.txt");
        if (!polyhedronFile.is_open())
            return false;

        polyhedronFile << "Id;NumVertices;VerticesId;NumEdges;EdgesId;NumFaces;FacesId" << endl;

        for (int i = 0; i < mesh.NumCell3Ds; i++) {
            polyhedronFile << mesh.Cell3DsId[i] << ";"
                           << mesh.Cell3DsNumVertices[i];

            for (int j = 0; j < mesh.Cell3DsNumVertices[i]; j++)
                polyhedronFile << ";" << mesh.Cell3DsVertices[i][j];

            polyhedronFile << ";" << mesh.Cell3DsNumEdges[i];

            for (int j = 0; j < mesh.Cell3DsNumEdges[i]; j++)
                polyhedronFile << ";" << mesh.Cell3DsEdges[i][j];

            polyhedronFile << ";" << mesh.Cell3DsNumFaces[i];

            for (int j = 0; j < mesh.Cell3DsNumFaces[i]; j++)
                polyhedronFile << ";" << mesh.Cell3DsFaces[i][j];

            polyhedronFile << endl;
        }
    }

    return true;
}

/**************************************************************************************************/


// Questa funzione prende in input una mesh e ne proietta i vertici sulla sfera di raggio 1 
void Projection(PolyhedralMesh& mesh)
{
	for(int i = 0; i < mesh.NumCell0Ds; i++)
	{
		//double Norm = (mesh.Cell0DsCoordinates[i]).norm();
		double Norm = (mesh.Cell0DsCoordinates.col(i)).norm();
		//mesh.Cell0DsCoordinates[i][0] = mesh.Cell0DsCoordinates[i][0]/Norm;
		mesh.Cell0DsCoordinates(0,i) = mesh.Cell0DsCoordinates(0,i)/Norm;
		mesh.Cell0DsCoordinates(1,i) = mesh.Cell0DsCoordinates(1,i)/Norm;
		mesh.Cell0DsCoordinates(2,i) = mesh.Cell0DsCoordinates(2,i)/Norm;
	}
}

/***************************************************************************************************/


// Questa funzione genera il geodetico a partire da un poliedro platonico, scelto dall'utente grazie aii dati presi in input
void GeodeticPolyhedron(const PolyhedralMesh& Platonic, PolyhedralMesh& Geodetic, const int& segments) 
{
	int point_id = 0;			
	int edge_id = 0;		
	int face_id = 0;
	
	
	int total_points = (Platonic.NumCell2Ds)*((segments + 1) * (segments + 2) / 2);	
	Geodetic.Cell0DsId.reserve(total_points);
	Geodetic.Cell0DsCoordinates = MatrixXd::Zero(3,total_points);


	/* allochiamo spazio in memoria per gli Id e i vertici delle celle 1D del poliedro geodetico.
		Numero max di spigoli da generare (si considera il caso con il numero massimo, cioè l'icosaedro)*/
	int total_edges = 30 * segments * segments; //questa formula è scritta nel pdf
	Geodetic.Cell1DsVertices = MatrixXi::Zero(2, total_edges);  
	Geodetic.Cell1DsId.reserve(total_edges);


	/* allochiamo spazio in memoria per gli Id ,i lati e i vertici delle celle 2D del poliedro geodetico.
		Numero max di facce triangolari da generare (si considera il caso con il numero massimo, cioè l'icosaedro)*/
	int total_faces = 20 * segments * segments; //questa formula è scritta nel pdf*/
	Geodetic.Cell2DsId.reserve(total_faces);

	Geodetic.Cell2DsEdges.resize(total_faces);
	Geodetic.Cell2DsNumEdges.resize(total_faces);

	Geodetic.Cell2DsNumVertices.resize(total_faces);
	Geodetic.Cell2DsVertices.resize(total_faces);
	
	

	map<array<int, 4>, int> coefficients;/*dichiara una mappa chiamata coefficients dove:
	la chiave è un array di 4 interi (std::array<int, 4>), mentre il valore associato a ogni chiave è un intero
	questa mappa associa ad ogni punto creato in points_id con un vettore di 4 interi, i cui 4 interi corrispondono ad a,b,c, Idcell2d, ovvero l'Id della faccia*/

	int duplicate_id = 0;

	for (const auto& id : Platonic.Cell2DsId) {
		
		/*per ogni faccia salvo i punti dei vertici, che mi servirannopoi per fare la combinazione convessa, e generare
		gli altri vertici della triangolazione*/
		Vector3d v1 = Platonic.Cell0DsCoordinates.col(Platonic.Cell2DsVertices[id][0]);
		Vector3d v2 = Platonic.Cell0DsCoordinates.col(Platonic.Cell2DsVertices[id][1]);
		Vector3d v3 = Platonic.Cell0DsCoordinates.col(Platonic.Cell2DsVertices[id][2]);
		
		/*con questo for genero appunto i vertici della triangolazione*/
		for (int i = 0; i <= segments; i++) {
			for (int j = 0; j <= i; j++) {
				
				/* Dato un triangolo con vertici: v1, v2, v3, qualsiasi punto P nel piano del triangolo può essere scritto come: 
				P = a * v1 + b * v2 + c*v3, dove i pesi sommano a 1. ora per non utilizzare dei double ma degli int (più semplici da iterare), utilizziamo la
				divisione per segments. a questo punto ossiamo disegnare tutte i punti presenti su di una faccia del poliedro attraverso
				la combinazioni dei valori a,b,c moltiplicate per i vertici */
				int a = segments - i;
				int b = i - j;
				int c = j;
				
				/*procedo appunto con la combinazione lineare cnvessa dei vertici della faccia del poliedro platonico
				per creare i nuovi vertici apartenenti al poliedro geodetico */
				Vector3d Point = double(a)/segments*v1 + double(b)/segments*v2 + double(c)/segments*v3;
				
				array<int, 4> coeffs; // Chiave della mappa: coefficienti + ID della faccia
				coeffs[0] = a;
				coeffs[1] = b;
				coeffs[2] = c;
				coeffs[3] = id;
				
				/* controllo che Point non sia già salvato in Geodetic.Cell0DsCoordinates, nel caso in cui non lo fosse
				lo associo al punto già esistente grazie alla mappa */ 
				if (!CheckVertices(Geodetic.Cell0DsCoordinates, Point, point_id, duplicate_id)){

					coefficients[coeffs] = point_id;

					Geodetic.Cell0DsId.push_back(point_id);
					Geodetic.Cell0DsCoordinates.col(point_id) = Point;

					point_id++;
					Geodetic.NumCell0Ds++;
					
				}
				else
					coefficients[coeffs] = duplicate_id;
			}
		}
	}
	
	Geodetic.Cell0DsCoordinates.conservativeResize(3, Geodetic.NumCell0Ds);
	Projection(Geodetic);
	
	/* costruisco adesso attraverso questi cicli gli spigoli della triangolazione, e le facce delimitate dai vari spigoli
	salvo il tutto poi all'interno di geodetic mesh*/
	for  (const auto& id : Platonic.Cell2DsId){
		for (int i = 0; i < segments; i++){
			for (int j = 0; j < segments - i; j++){
				
				int v1 = coefficients[{i,segments-i-j,j,id}];
				int v2 = coefficients[{i,segments-i-(j+1),j+1,id}];
				int v3 = coefficients[{i+1,segments-(i+1)-j,j,id}];
				vector<int> Vertex_face = {v1, v2, v3};
				
				Geodetic.NumCell2Ds++;

				Geodetic.Cell2DsId.push_back(face_id);

				Geodetic.Cell2DsNumVertices[face_id] = 3;
				Geodetic.Cell2DsVertices[face_id] = Vertex_face;

				Geodetic.Cell2DsNumEdges[face_id] = 3;
				Geodetic.Cell2DsEdges[face_id].resize(3);
				
				// edges
				for (int k = 0; k < 3; k++) {
					int v_start = Geodetic.Cell2DsVertices[face_id][k];
					int v_end;
					if ( k == 2 )
						v_end = Geodetic.Cell2DsVertices[face_id][0];
					else
						v_end = Geodetic.Cell2DsVertices[face_id][k+1];
					
					if(!CheckEdges(Geodetic.Cell1DsVertices, v_start, v_end, edge_id, duplicate_id)){
						Geodetic.NumCell1Ds++;

						Geodetic.Cell1DsId.push_back(edge_id);

						Geodetic.Cell1DsVertices(0, edge_id) = v_start;
						Geodetic.Cell1DsVertices(1, edge_id) = v_end;

						Geodetic.Cell2DsEdges[face_id][k] = edge_id;

						edge_id++;
						
					}
					else
						Geodetic.Cell2DsEdges[face_id][k] = Geodetic.Cell1DsId[duplicate_id];
				}
				face_id++;
				
				/*la combinazione di vertici utilizzata nel codice sopra non è sufficiente per la creazione di tutte 
				le possibili facce, serve infatti considerare anche i triangoli che hanno la punta rivolta verso il 
				basso, e che esistono solo per i>0 (stiamo infatti concettualmente "salendo" di 1 nella griglia triangolare)*/
				if(i > 0){
					
					int v4 = coefficients[{i-1,segments-(i-1)-(j+1),(j+1),id}];
					
					
					Geodetic.NumCell2Ds++;
					Geodetic.Cell2DsId.push_back(face_id);

					
					Vertex_face[2] = v4; /*il vertice va salvato proprio nella posizone 2!*/
					

					Geodetic.Cell2DsNumVertices[face_id] = 3;
					Geodetic.Cell2DsVertices[face_id] = Vertex_face;


					Geodetic.Cell2DsNumEdges[face_id] = 3;
					Geodetic.Cell2DsEdges[face_id].resize(3);
					

					for (int k = 0; k < 3; k++) {
						int v_start = Geodetic.Cell2DsVertices[face_id][k];
						int v_end;
						if ( k == 2 )
							v_end = Geodetic.Cell2DsVertices[face_id][0];
						else
							v_end = Geodetic.Cell2DsVertices[face_id][k+1];
						
						if (!CheckEdges(Geodetic.Cell1DsVertices, v_start, v_end, edge_id, duplicate_id)){
							
							Geodetic.NumCell1Ds++;

							Geodetic.Cell1DsId.push_back(edge_id);

							Geodetic.Cell1DsVertices(0, edge_id) = v_start;
							Geodetic.Cell1DsVertices(1, edge_id) = v_end;

							Geodetic.Cell2DsEdges[face_id][k] = edge_id;
							
							edge_id++;
							
						}
						else
							Geodetic.Cell2DsEdges[face_id][k] = Geodetic.Cell1DsId[duplicate_id];
					}

					face_id++;

				}
				
			}
		}
	}
	/*eseguo qualche resize per sistemare la mesh*/
	Geodetic.Cell1DsVertices.conservativeResize(2, Geodetic.NumCell1Ds);

	Geodetic.Cell2DsNumVertices.resize(Geodetic.NumCell2Ds);
	Geodetic.Cell2DsVertices.resize(Geodetic.NumCell2Ds);


	Geodetic.Cell2DsNumEdges.resize(Geodetic.NumCell2Ds);
	Geodetic.Cell2DsEdges.resize(Geodetic.NumCell2Ds);
	
	/*in conclusione salvo i risultati ottenuti dentro la mesh anche per la Cell3D*/
	Geodetic.NumCell3Ds++;
	Geodetic.Cell3DsId.push_back(0);


	Geodetic.Cell3DsNumVertices.push_back(Geodetic.NumCell0Ds);
	Geodetic.Cell3DsVertices.push_back(Geodetic.Cell0DsId);


	Geodetic.Cell3DsNumEdges.push_back(Geodetic.NumCell1Ds);
	Geodetic.Cell3DsEdges.push_back(Geodetic.Cell1DsId);


	Geodetic.Cell3DsNumFaces.push_back(Geodetic.NumCell2Ds);
	Geodetic.Cell3DsFaces.push_back(Geodetic.Cell2DsId);
}

/************************************************************************************************/

//
void DualMesh(const PolyhedralMesh& StartPolyhedron, PolyhedralMesh& DualPolyhedron) 
{
	int baricenter_id = 0;
	int face_id = 0;
	int edge_id = 0;
	
	//Il Duale ha un numero di facce uguale al numero di vertici del poliedro di partenza
	DualPolyhedron.NumCell0Ds = StartPolyhedron.NumCell2Ds;
	
	//Il Duale ha lo stesso numero di edges del poliedro di partenza, grazie alla formula di Eulero, qualunque sia la varietà su cui si fa la mesh :)
	DualPolyhedron.NumCell1Ds = StartPolyhedron.NumCell1Ds;
	
	//Il Duale ha un numero di facce uguale al numero di vertici del poliedro di partenza
	DualPolyhedron.NumCell2Ds = StartPolyhedron.NumCell0Ds;
	DualPolyhedron.Cell0DsId.reserve(DualPolyhedron.NumCell0Ds);
	DualPolyhedron.Cell0DsCoordinates = MatrixXd::Zero(3,DualPolyhedron.NumCell0Ds);	
	
	DualPolyhedron.Cell1DsVertices = MatrixXi::Zero(2, DualPolyhedron.NumCell1Ds);
	DualPolyhedron.Cell1DsId.reserve(DualPolyhedron.NumCell1Ds);
	
	DualPolyhedron.Cell2DsId.reserve(DualPolyhedron.NumCell2Ds);
	DualPolyhedron.Cell2DsEdges.resize(DualPolyhedron.NumCell2Ds);
	DualPolyhedron.Cell2DsNumEdges.resize(DualPolyhedron.NumCell2Ds);
	DualPolyhedron.Cell2DsVertices.resize(DualPolyhedron.NumCell2Ds);
	DualPolyhedron.Cell2DsNumVertices.resize(DualPolyhedron.NumCell2Ds);
	int duplicate_id = 0;
	//Mappa che associa all'id della faccia l'id del baricentro corrispondente, da usare se cambiassimo gli id dopo, per ora sono uguali
	map<int, int> Faces_bar;
	for (const auto& id : StartPolyhedron.Cell2DsId) {
		Vector3d baricenter;
		
		// salvo in 3 vettori le coordinate dei vertici della faccia corrente del solido platonico
		Vector3d Vertex1 = StartPolyhedron.Cell0DsCoordinates.col(StartPolyhedron.Cell2DsVertices[id][0]);
		Vector3d Vertex2 = StartPolyhedron.Cell0DsCoordinates.col(StartPolyhedron.Cell2DsVertices[id][1]);
		Vector3d Vertex3 = StartPolyhedron.Cell0DsCoordinates.col(StartPolyhedron.Cell2DsVertices[id][2]);
		
		// salvo le coordinate del baricentro
		baricenter = (1.0/3.0)*Vertex1 + (1.0/3.0)*Vertex2 + (1.0/3.0)*Vertex3;
		DualPolyhedron.Cell0DsId.push_back(baricenter_id);
		// salvo le coordinate dentro alle Cell0DsCoordinates del poliedro Duale
		DualPolyhedron.Cell0DsCoordinates(0, id) = baricenter(0);
		DualPolyhedron.Cell0DsCoordinates(1,id) = baricenter(1);
		DualPolyhedron.Cell0DsCoordinates(2, id) = baricenter(2);
		
		//Associo all'id della faccia l'id del baricentro nella mappa, per ora ridondante ma poi sarà meglio
		Faces_bar[id] = baricenter_id;
		baricenter_id ++;
	}

	for(const auto& vertex_id: StartPolyhedron.Cell0DsId){
		
		//Vettore che contiene le facce che hanno il vertice in comune
		vector<int> VertexFaces;
		for(const auto& face_id: StartPolyhedron.Cell2DsId){
			for(int j = 0; j < 3; j++){
				if (StartPolyhedron.Cell2DsVertices[face_id][j] == vertex_id){
					//Se la faccia a cui sono arrivato contiene il vertice, allora la aggiungo al vettore 
					//delle facce avente quel vertice in comune
					VertexFaces.push_back(face_id);
					break;
				}
			}
		}
		
		//PROBLEMA: il vettore VertexFaces contiene tutte le facce adiacenti a un vertice ma NON è ordinato in modo sensato, 
		//per costruire gli edges devo ordinarlo in modo che ogni faccia nel vettore abbia come successiva la faccia adiacente, 
		//ovvero quella che ha l'edge in comune con la faccia corrente, per ordinare questo vettore chiamo la funzione OrderFaces
		//il vettore ordered_Faces è passato per riferimento, in modo che venga aggiornato dalla funzione order_Faces.
		vector<int> ordered_faces;
		Sort_Faces(VertexFaces, ordered_faces, StartPolyhedron);
		
		//la valenza del vertice è pari alla lunghezza del vettore di facce che condividono il vertice dato 
		//ATTENZIONE: Questa parte non è superflua, perché le valenze NON SONO sempre 3 per il generico solido geodetico!!!
		int valence = ordered_faces.size();
		vector<int> New_vertices;
		
		//qui associo a ogni faccia del poliedro di partenza l'id del vertice nel duale corrispondente
		for(const auto& VertexFace: ordered_faces)
			New_vertices.push_back(Faces_bar[VertexFace]);
		
		
		DualPolyhedron.Cell2DsId.push_back(face_id);
		DualPolyhedron.Cell2DsVertices[face_id] = New_vertices;
		DualPolyhedron.Cell2DsEdges[face_id].resize(valence);
		
		DualPolyhedron.Cell2DsNumVertices[face_id] = valence;
		DualPolyhedron.Cell2DsNumEdges[face_id] = valence;
		
		//Questa è la creazione degli edges, praticamente identica a quella per il solido geodetico, 
		//con la differenza che il vettore di vertici della faccia ha tanti elementi quanti la valenza del 
		//vertice, e a parte questo dettaglio è tutto uguale!
		for (int k = 0; k < valence; k++) {
			int originVertex = DualPolyhedron.Cell2DsVertices[face_id][k];
			int endVertex;
			if ( k == valence-1 )
				endVertex = DualPolyhedron.Cell2DsVertices[face_id][0];
			else
				endVertex = DualPolyhedron.Cell2DsVertices[face_id][k+1];
			if(!CheckEdges(DualPolyhedron.Cell1DsVertices, originVertex, endVertex, edge_id, duplicate_id)){
				DualPolyhedron.Cell1DsId.push_back(edge_id);
				DualPolyhedron.Cell1DsVertices(0, edge_id) = originVertex;
				DualPolyhedron.Cell1DsVertices(1, edge_id) = endVertex;
				DualPolyhedron.Cell2DsEdges[face_id][k] = edge_id;
				edge_id++;
			}
			else
				DualPolyhedron.Cell2DsEdges[face_id][k] = DualPolyhedron.Cell1DsId[duplicate_id];
		}
		face_id++;
	}

	Projection(DualPolyhedron);
	
	DualPolyhedron.Cell1DsVertices.conservativeResize(2, DualPolyhedron.NumCell1Ds);
	DualPolyhedron.Cell2DsNumVertices.resize(DualPolyhedron.NumCell2Ds);
	DualPolyhedron.Cell2DsNumEdges.resize(DualPolyhedron.NumCell2Ds);
	DualPolyhedron.Cell2DsVertices.resize(DualPolyhedron.NumCell2Ds);
	DualPolyhedron.Cell2DsEdges.resize(DualPolyhedron.NumCell2Ds);
	
	// GENERAZIONE POLIEDRO
	DualPolyhedron.NumCell3Ds++;
	DualPolyhedron.Cell3DsId.push_back(0);
	DualPolyhedron.Cell3DsNumVertices.push_back(DualPolyhedron.NumCell0Ds);
	DualPolyhedron.Cell3DsNumEdges.push_back(DualPolyhedron.NumCell1Ds);
	DualPolyhedron.Cell3DsNumFaces.push_back(DualPolyhedron.NumCell2Ds);
	DualPolyhedron.Cell3DsVertices.push_back(DualPolyhedron.Cell0DsId);
	DualPolyhedron.Cell3DsEdges.push_back(DualPolyhedron.Cell1DsId);
	DualPolyhedron.Cell3DsFaces.push_back(DualPolyhedron.Cell2DsId);
}
		
/************************************************************************************************/


// Questa funzione controlla che uno spigolo non esista già 
bool CheckEdges(const MatrixXi& verts, const int& v1, const int& v2, int& edge_id, int& existing_edge_id)
{
    for (int i = 0; i < edge_id; i++)
    {
        int a = verts(0,i);
        int b = verts(1,i);
        if ((a == v1 && b == v2) || (a == v2 && b == v1)) 
		// Questo if controlla se la coppia di variabili a e b è uguale alla coppia v1 e v2, in qualunque ordine.
        {
			existing_edge_id = i;
            return true;
        }
    }
    return false;
}

/************************************************************************************************/


// Questa funzione controlla che un vertice non esista già 
bool CheckVertices(const MatrixXd& mesh, const Vector3d& point, int& point_id, int& duplicate_id) //mesh corrisponde a Geodetic.Cell0Dscoordinates
{
	double eps = 1e-8;
    for (int i = 0; i < point_id; i++)
    {
        if ((mesh.col(i) - point).norm() < eps)
        {
            duplicate_id = i;
            return true;
        }
    }
    return false;
}

/************************************************************************************************/

void Sort_Faces(const vector<int>& UnsortedFaces, vector<int>& SortedFaces, const PolyhedralMesh& Mesh) 
{
    // Creiamo una copia delle facce da ordinare per lavorare su di esse
    vector<int> RemainingFaces = UnsortedFaces; // Facce ancora da ordinare.
	int CurrentFace = RemainingFaces[0];
	
    // Selezioniamo la prima faccia arbitrariamente e la inseriamo nel vettore ordinato
    SortedFaces.push_back(CurrentFace); //(RemainingFaces.front());
	// Rimozione della faccia perchè già stata ordinata
    RemainingFaces.erase(RemainingFaces.begin());

    // Iteriamo fino a quando tutte le facce sono state ordinate
    while (!RemainingFaces.empty()) {
        // int CurrentFace = SortedFaces.back();  // Ultima faccia ordinata
        vector<int> CurrentEdges = Mesh.Cell2DsEdges[CurrentFace]; // Recuperiamo i suoi spigoli per confrontarli con le altre facce

        bool Match = false; 
		// int BestMatch = -1; // valore di default
        int BestMatchIndex;

		int RemainingFacesSize = RemainingFaces.size();
        // Scorriamo le facce rimaste per trovare quella più adiacente
        for (int i = 0; i < RemainingFacesSize; i++) {
            vector<int> CandidateEdges = Mesh.Cell2DsEdges[RemainingFaces[i]]; // Ottenimento degli spigoli
            //int SharedEdges = 0;  // Contatore degli spigoli condivisi

            /* Controlla quanti spigoli hanno in comune
			Confronta ogni spigolo della faccia attuale (CurrentEdges) con quelli della faccia candidata (CandidateEdges).
			Se trova un match, incrementa SharedEdges.*/
            for (int Edge1 : CurrentEdges) {
                for (int Edge2 : CandidateEdges) {
                    if (Edge1 == Edge2) {
                        Match = true;
						BestMatchIndex = i;
						//SharedEdges++;
						break;
                    }
                }
				if(Match)
					break;
            }
			if(Match)
				break;
			/*Seleziono la miglior faccia
        	Se questa faccia ha più spigoli in comune, aggiorniamo la migliore scelta*/
            //if (SharedEdges > 0 && (BestMatch == -1 || SharedEdges > Mesh.Cell2DsEdges[BestMatch].size())) {
			/*if(Match)
			{
				BestMatch = RemainingFaces[i];
                BestMatchIndex = i;
            }*/
        }

        if(Match) 
		{
		CurrentFace = RemainingFaces[BestMatchIndex];
        SortedFaces.push_back(CurrentFace);
        RemainingFaces.erase(RemainingFaces.begin() + BestMatchIndex);
		}
    }
}

/************************************************************************************************/

bool ShortestPath(const PolyhedralMesh& mesh, const int& start, const int& end, vector<int>& path, MatrixXd& W)   // prima del vettore c'erano double& length, int& NumPath
{
		
	/*if (start >= mesh.NumCell0Ds || end >= mesh.NumCell0Ds || start < 0 || end < 0)
	{
		cerr << "The inserted Ids are not valid." << endl;
		return false;
	}*/
	
	// generazione della lista di adiacenza, poiché è tutto indicizzato sequenzialmente, 
	// conviene usare un vector di vector anziché un vector di liste
	vector<vector<int>> adjacency_list;
	adjacency_list.reserve(mesh.NumCell0Ds);
	
	// iterando su tutti gli id dei punti, 
	// per ciascun punto itero in tutti gli id dei segmenti (origin,end) e guardo quali hanno per estremo quel punto.
	// L'estremo che è diverso dal punto in questione lo appendo alla "lista" di adiacenza.
	
	for(int i = 0; i < mesh.NumCell0Ds; i++)
	{
		vector<int> adj;
		for(const auto& edge : mesh.Cell1DsId)
		{
			int Origin = mesh.Cell1DsVertices(0,edge);
			int End = mesh.Cell1DsVertices(1,edge);
			
			if (Origin == i)
				adj.push_back(End);
			else if(End == i)
				adj.push_back(Origin);
		}
		adjacency_list.push_back(adj);
	}
	
	//creo la matrice delle distanza
	for(size_t i = 0; i < adjacency_list.size(); i++)
	{
		for(const auto& adj: adjacency_list[i])
		{
			W(i,adj) = (mesh.Cell0DsCoordinates.col(adj)-mesh.Cell0DsCoordinates.col(i)).norm();
		}
	}
	
	// algoritmo di Dijkstra per esplorare il grafo, pred è un vettore
	// ausiliario usato per ricostruire il percorso

	vector<int> pred(mesh.NumCell0Ds, -1);
	vector<double> dist(mesh.NumCell0Ds, 1000.0);
	priority_queue<pair<int, double>, vector<pair<int, double>>, greater<pair<int, double>>> PQ;
	
	pred[start] = start;
	dist[start] = 0;
	
	for(int i = 0; i < mesh.NumCell0Ds; i++)
		PQ.push({i, dist[i]});
	while(!PQ.empty())
	{
		int u = PQ.top().first;
		int p = PQ.top().second;
		PQ.pop();
		if (dist[u] < p)
			continue;
		for(const auto& w : adjacency_list[u])
		{
			if( dist[w] > dist[u] + W(u,w) ) 
			{
				dist[w] = dist[u] + W(u,w);
				pred[w] = u;
				PQ.push({w, dist[w]});
			}
		}
	}
	
	// path contiene gli id dei vertici che compongono il cammino minimo 
	// al contrario, perché sono id provenienti dal vettore pred
	
	
	int v = end;
	/*if (pred[v] == -1) {
		std::cerr << "No path exists from " << start << " to " << end << std::endl;
		return false;
	}*/

	while (v != start) {
		path.push_back(v);
		v = pred[v];
	}
	path.push_back(start);
	//std::reverse(path.begin(), path.end()); // opzionale, se vuoi da start → end
	
	int size = path.size();
	cout << "Shortest path between the vertices " << start << " and " << end << ":" << endl;
	for (int i=0; i < size-1; i++)
	{
		cout << path[i] << ",";
	}
	cout << path[size-1] << "." << endl;
	
	return true;
}

/*******************************************************************************************/

bool isInteger(const string& str)
{
	if (str.empty())
		return false;
	
	int L = str.length();
	for (int i = 0; i < L; i++)
	{
		if (!isdigit(str[i]))
			return false;
	}	
	return true;
}

/********************************************************************************************/

}

// PER ALE: CANCELLA TUTTO; LA GRAFFA DI CHIUSURA DEL NAMESPACE E' GIA' QUELLA DI RIGA 867, QUINDI NON SERVE AGGIUNGERE NULLA. DA QUI IN GIU' E' TUTTO INUTILE

/*bool ExpPath(PolyhedralMesh& mesh, vector<int> path, double& length, int& NumPath, MatrixXd& W)
{
	vector<double> PathPointsProperties(mesh.NumCell0Ds, 0.0);
		for (const auto& point : path)
			PathPointsProperties[point] = 1.0;

			
		Gedim::UCDProperty<double> ShortPathProperty;
		ShortPathProperty.Label = "shortest path";
		ShortPathProperty.UnitLabel = "";
		ShortPathProperty.Size = PathPointsProperties.size();
		ShortPathProperty.NumComponents = 1;
		ShortPathProperty.Data = PathPointsProperties.data();  


		vector<Gedim::UCDProperty<double>> PointsProperties;
		PointsProperties.push_back(ShortPathProperty);
	
		Gedim::UCDUtilities utilities;
		utilities.ExportPoints("./Cell0Ds.inp",
								mesh.Cell0DsCoordinates,
								PointsProperties);
	
		vector<int> pathEdges; 
		vector<double> PathEdgesProperties(mesh.NumCell1Ds, 0.0);

		for (size_t i = 0; i < path.size()-1; i++){
			int v1 = path[i];
			int v2 = path[i+1];
			for(const auto& edge: mesh.Cell1DsId){
				if ((mesh.Cell1DsVertices(0, edge) == v1 && mesh.Cell1DsVertices(1,edge) == v2) || 
					(mesh.Cell1DsVertices(0, edge) == v2 && mesh.Cell1DsVertices(1,edge) == v1)){
						pathEdges.push_back(edge);
						PathEdgesProperties[edge] = 1.0;
					}
			}	
		}
		
		for(size_t i = 0; i < path.size()-1; i++)
			length += W(path[i],path[i+1]);
		NumPath = pathEdges.size();
		
		
		ShortPathProperty.Label = "shortest path";
		ShortPathProperty.UnitLabel = "";
		ShortPathProperty.Size = PathEdgesProperties.size();
		ShortPathProperty.NumComponents = 1;
		ShortPathProperty.Data = PathEdgesProperties.data();  
	
		vector<Gedim::UCDProperty<double>> EdgesProperties;
		EdgesProperties.push_back(ShortPathProperty);
		utilities.ExportSegments("./Cell1Ds.inp",
								mesh.Cell0DsCoordinates,
								mesh.Cell1DsVertices,
								{},
								EdgesProperties);
		return true;
	}
		
	
/*bool GenerateGoldbergClassI(int p, int q, int b, int c, PolyhedralMesh& Goldberg) 
{
	return true;
}

// ***************************************************************************


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