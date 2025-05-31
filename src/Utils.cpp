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


bool PolyhedralChoice(string& path) // c'erano anche PolyhedralMesh& mesh, const char& b, const char& c, bool& walk, const char& p, const char& q
{
	string polyhedron;
	string filePath;
	string p_prov;
	string q_prov;
	
	cout << "Enter each parameter when asked." << endl;
	cout << "Insert p: " << endl;
	cin >> p_prov;
	cout << "Insert q: " << endl;
	cin >> q_prov;
	
	int p = stoi(p_prov);
	int q = stoi(q_prov);
	
	if (p == 3)
	{
		switch(q)
		{
			case 3:
				polyhedron = "/Tetraedro";
				break;
			case 4:
				polyhedron = "/Ottaedro";
				break;
			case 5:
				polyhedron = "/Icosaedro";
				break;
			default:
				return false;
		}
		
		path = path + polyhedron;
		
		return true;
	}
	
	
	return false; 
}

/************************************************************************************************/

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

/************************************************************************************************/

/*bool ParameterSelection(string& path, PolyhedralMesh& mesh)
{
	// vector<string> v;
	char p;
	char q;
	char b;
	char c;
	char Id1;
	char Id2;
	char tmp;
	bool walk = false;
	
	cout << "Insert each parameter when asked and push enter to confirm." << endl;
	cout << "Insert p: " << endl;
	cin >> p;
	cout << "Insert q: " << endl;
	cin >> q;
	cout << "Insert b: " << endl;
	cin >> b;
	cout << "Insert c: " << endl;
	cin >> c;
	
	cout << "If you don't want to evaluate the shortest path between two vertices, please digit 'n' and enter it in the next two insertions." << endl; 
	cout << "Insert Id1: " << endl;
	cin >> Id1;
	cout << "Insert Id2: " << endl;
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
*/
/************************************************************************************************/

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
    //mesh.Cell0DsCoordinates.reserve(mesh.NumCell0Ds);
    mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
	mesh.Cell0DsCoordinates = MatrixXd::Zero(3,mesh.NumCell0Ds);
	
    // Vector3d Coordinates;
    char tmp;
    int Id;
	
	for(const auto& l : lines)
	{
		stringstream ss(l);
		// ss >> Id >> tmp >> Coordinates(0) >> tmp >> Coordinates(1) >> tmp >> Coordinates(2);
		ss >> Id >> tmp >> mesh.Cell0DsCoordinates(0,Id) >> tmp >> mesh.Cell0DsCoordinates(1,Id) >> tmp >> mesh.Cell0DsCoordinates(2,Id);
		mesh.Cell0DsId.push_back(Id);
		//mesh.Cell0DsCoordinates.push_back(Coordinates);
	}
	
    return true; 
}

/************************************************************************************************/

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
    mesh.Cell1DsVertices = MatrixXi::Zero(2, mesh.NumCell1Ds); //.reserve(mesh.NumCell1Ds);
    mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
    
    // Vector2i Vertices;
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



bool GeodeticPolyhedron(const PolyhedralMesh& Platonic, PolyhedralMesh& Geodetic, const int& segments)
{
    int points_id = 0;   // inizializziamo a zero tutti gli id del poliedro
    int edges_id = 0;
    int faces_id = 0;
	int duplicate_id = 0;

	/*NOTA!!!!!! A QUESTO PUNTO CI SARà DA FARE UN IF, IN MODO CHE IN BASE ALLA q FORNITA SI RISERVI PER IL VETTORE LO SPAZIO CORRETTO!!!!*/
	//QUESTO è IL CASO DELL'ICOESAEDRO (PER ORA), SI DOVRà CONTROLLARE LA CORRETTEZZA DEGLI INPUT

	// Calcolo il numero dei punti da generare per il poliedro in generale,
	//moltiplicando il numero di punti su una faccia per il numero di facce
    int total_points = Platonic.NumCell2Ds * ((segments + 1) * (segments + 2) / 2);/*somma Gaussiana: 
	in poche parole, la formula corrisponde alla somma dei numeri da 1 a segments +1, ove segments=b.  */

	// allochiamo spazio in memoria per gli Id e le coordinate delle celle 0D del poliedro geodetico.
    Geodetic.Cell0DsId.reserve(total_points);
    Geodetic.Cell0DsCoordinates = MatrixXd::Zero(3,total_points); // Geodetic.Cell0DsCoordinates.reserve(total_points);


    map<array<int, 4>, int> coefficients;/*dichiara una mappa chiamata coefficients dove:
	la chiave è un array di 4 interi (std::array<int, 4>), mentre il valore associato a ogni chiave è un intero*/
	//questa mappa associa ad ogni punto creato in points_id con un vettore di 4 interi, i cui 4 interi corrispondono ad a,b,c, Idcell2d, ovvero l'Id della faccia


	// allochiamo spazio in memoria per gli Id e i vertici delle celle 1D del poliedro geodetico.
	// Numero max di spigoli da generare (si considera il caso con il numero massimo, cioè l'icosaedro)
    int total_edges = 30 * segments * segments; //questa formula è scritta nel pdf
    Geodetic.Cell1DsId.reserve(total_edges); 
    Geodetic.Cell1DsVertices = MatrixXi::Zero(2, total_edges); // Geodetic.Cell1DsVertices.reserve(total_edges); 


	// allochiamo spazio in memoria per gli Id ,i lati e i vertici delle celle 2D del poliedro geodetico.
	// Numero max di facce triangolari da generare (si considera il caso con il numero massimo, cioè l'icosaedro)
    int total_faces = 20 * segments * segments; //questa formula è scritta nel pdf
    Geodetic.Cell2DsId.reserve(total_faces); 
    Geodetic.Cell2DsVertices.resize(total_faces);
    Geodetic.Cell2DsEdges.resize(total_faces);
	Geodetic.Cell2DsNumEdges.resize(total_faces);
	Geodetic.Cell2DsNumVertices.resize(total_faces);


	//QUA CI ANDREBBE LA FINE DEGLI IF

	// Ciclo su ogni faccia del poliedro (Platonic)
    for (const auto& id : Platonic.Cell2DsId) // QUESTO FOR è DA MODIFICARE, METTERE ALTRO PARAMETRO AL POSTO DI ID, E SCRIVERE POI DENTRO
											//QUALCOSA PER L'INCREMENTO DEL ID DELLA FACCIA
    {
		// Estrazione dei 3 vertici dalla faccia corrente e salvataggio di essi in vettori       
		//Vector3d v1 = Platonic.Cell0DsCoordinates[Platonic.Cell2DsVertices[id][0]]; 
        Vector3d v1 = Platonic.Cell0DsCoordinates.col(Platonic.Cell2DsVertices[id][0]);
		Vector3d v2 = Platonic.Cell0DsCoordinates.col(Platonic.Cell2DsVertices[id][1]);
        Vector3d v3 = Platonic.Cell0DsCoordinates.col(Platonic.Cell2DsVertices[id][2]); // tra le parentesi [] c'è un punto                                             


		// Genera i punti interni alla faccia con suddivisione baricentrica
        for (int i = 0; i <= segments; i++)  //<=segments perchè è il numero (b) di divisioni che vogliamo fare su ogni spigolo
        {
            for (int j = 0; j <= i; j++)	/* ciclo in i: per ogni i da 0 a b genera i punti interni */
            {
                int a = segments - i;  /* Dato un triangolo con vertici: v1, v2, v3, qualsiasi punto P nel piano del triangolo può essere scritto come: 
										P = a * v1 + b * v2 + c*v3, dove i pesi sommano a 1. ora per non utilizzare dei double ma degli int, utilizziamo la
										divisione per segments. a questo punto ossiamo disegnare tutte i punti presenti su di una faccia del poliedro attraverso
										la combinazioni dei valori a,b,c moltiplicate per i vertici */      
				int b = i - j;
                int c = j;
		/*for (const auto& id : Platonic.Cell2DsId) METODO UTILIZZABILEL PER COSTRUIRE I PUNTI DELLA TRIANGOLARIZZAZIONE
    {
        for (int i = 0; i < segments; i++)
        {
            for (int j = 0; j < segments - i; j++)
            {
                int a = segments -i -j;
				int b = i;
				int c = j;
*/

				/* Calcolo del punto tramite combinazione lineare (a e b diventano double per evitare problemi con la divisione); dividiamo per segments in modo da avere
				la normalizzazione (infatti a, b e c sommano a segments, mentre dovrebbero sommare a 1) */
                Vector3d Point = (double(a) / segments) * v1 + (double(b) / segments) * v2 + (double(c) / segments) * v3; /*Ho così creato il vettore
				di coordinate corrispondente al punto creato (funzione di a,b,c,IDfaccia)*/ 

				// Chiave della mappa: coefficienti + ID della faccia
                array<int, 4> coeffs;
				coeffs[0] = a;
				coeffs[1] = b;
				coeffs[2] = c;
				coeffs[3] = id;

				/* Il procedimento implementato non garantisce che non ci sia sovrapposizione di punti: per questo, si richiama la funzione CheckDuplicatesVertex,
				questa controlla che il vertice non sia duplicato. se è duplicato  */ 
                if (!CheckVertices(Geodetic.Cell0DsCoordinates, Point, points_id, duplicate_id))
                {
                    coefficients[coeffs] = points_id;  // assegno all'array l'id corrente (all'inizio, points_id è inizializzato a 0) 
                    Geodetic.Cell0DsId.push_back(points_id);   // aggiungo il nuovo punto (points_id) in coda agli id del poliedro
                    // Geodetic.Cell0DsCoordinates.push_back(Point);  // aggiungo il punto in coda al poliedro geodetico (si può anche pensare di inserire le coordinate una per una?mi sembra decisamente inutile ma volendo...(ale))
                    Geodetic.Cell0DsCoordinates(0,points_id) = Point[0];
					Geodetic.Cell0DsCoordinates(1,points_id) = Point[1];
					Geodetic.Cell0DsCoordinates(2,points_id) = Point[2];
					Geodetic.NumCell0Ds++;   // incremento il contatore di celle 0D del poliedro geodetico e quello dei punti
                    points_id++;
                }
                else
                {
                    coefficients[coeffs] = duplicate_id;  // lo aggiungiamo tra i punti duplicati, l'id corretto è stato preso dentro alla funzione CheckVertices
                } //A QUESTO PUNTO HO CREATO TUTTI I PUNTI PER LA TRIANGOLAZIONE
            }
        }
    }
	
	Geodetic.Cell0DsCoordinates.conservativeResize(3, Geodetic.NumCell0Ds);  

	// Normalizzazione degli elementi del poliedro. ATTENZIONE PERCHè TUTTI I VERTICI DEL GEODETUCO DEVONO ESSERE NORMALIZZATI!! (dovrebbe essere già così)
	Projection(Geodetic);
	// Triangolazione
	GenerateTriangles(Platonic, Geodetic, coefficients, segments, edges_id, faces_id);

	return true;
}

/************************************************************************************************/

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

/************************************************************************************************/

void GenerateTriangles(const PolyhedralMesh& Platonic, PolyhedralMesh& Geodetic, map<array<int, 4>, int>& coefficients, int segments, int& edges_id, int& faces_id)
{
	for (const auto& id : Platonic.Cell2DsId)
    {
        for (int i = 0; i < segments; i++)
        {
            for (int j = 0; j < segments - i; j++)
            {

				//Secondo me c'è un errore, perchè se poniamo un punto di estremo, come un punto che ha c=segments (b), allora poi v2 avrebbe coordinate come
				//un punto che ha c=segments+1, il che è impossibile. no l'errore non c'è perchè il for non implica il <=, ma solo il <.
				//Potrebbe essere risolto imponendo una condizione su a --> if ( (a=)segments -i -j >0)... allora esegui codice...-->
				//--> vogliamo evitare che prenda i vertici che giacciono sullo spigolo destro del triangolo, per quei vertici infatti non ci sono
				// triangolazioni da fare.
	
                // Triangolo “a punta in su”
                int v1 = coefficients[{segments - i - j, i, j, id}];       // richiamiamo tre punti adiacenti per fare la triangolazione: nella prima iterazione, 
				// v1 è coefficients[{segments, 0, 0, id}], v2 è coefficients[{segments -1, 0, 1, id}], v3 è coefficients[{segments -1, 1, 0, id}]
                int v2 = coefficients[{segments - i - (j + 1), i, j + 1, id}];
                int v3 = coefficients[{segments - (i + 1) - j, i + 1, j, id}];

				//int v2 = coefficients[{segments - i - j - 1, i, j + 1, id}]; metodo equivalente (Non hai cambiato quasi niente (ale))
                //int v3 = coefficients[{segments - i - 1 - j, i + 1, j, id}];
				
				/*int v1 = coefficients[{segments - i - j, i, j, id}];    
                int v2 = coefficients[{segments - i - j - 1, i, j + 1, id}];
                int v3 = coefficients[{segments - (i + 1) - j, i + 1, j, id}];*/

				//ogni volta che seleziono questi tre vertici

				int existing_edge_id = 0;

				// generazione delle facce: si parte da queste perché così dopo gli spigoli sanno già dove attaccarsi senza il rischio di duplicati
                Geodetic.NumCell2Ds++;    // geodetic è una struct di tipo PolyhedralMesh, dunque ha NumCell2Ds (inizialmente uguale a 0)
				Geodetic.Cell2DsId.push_back(faces_id);
				Geodetic.Cell2DsNumVertices[faces_id] = 3;  // inizializziamo le facce a triangolari, così dopo cicleremo direttamente su quelle
				Geodetic.Cell2DsNumEdges[faces_id] = 3;		//infatti ogni faccia creata è sempre triangolare, ed ha sempre 3 lati e 3 angoli
				vector<int> Vertices = {v1, v2, v3};
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
					
                    if (!CheckEdges(Geodetic.Cell1DsVertices, v_start, v_end, edges_id, existing_edge_id))   // check se lo spigolo non è già esistente
                    {
                        Geodetic.Cell1DsId.push_back(edges_id);
                        //Geodetic.Cell1DsVertices.push_back(Vector2i(v_start, v_end));
						Geodetic.Cell1DsVertices(0,edges_id) = v_start;
						Geodetic.Cell1DsVertices(1,edges_id) = v_end;
						Geodetic.Cell2DsEdges[faces_id][k] = edges_id;  // andiamo a inserire gli id dei vertici
                        Geodetic.NumCell1Ds++;
                        edges_id++;
                    }
					else
						Geodetic.Cell2DsEdges[faces_id][k] = existing_edge_id; 
                }
				
                faces_id++;

                // Triangolo “a punta in giù”
				/*io avrei messo 
				if(a>0 && c>1){ ----> non serve fare così perchè nel for usa il < e non il <=, come ho già precedentemente fatto notare
				//v1= coefficients[{a,b,c,Id}]
				//v2 = coefficients[{a-1,b+1,c,Id}]
				v4= coefficients[{a,b+1,c-1}]

				Geodetic.NumCell2Ds++;
                    Geodetic.Cell2DsId.push_back(faces_id);
                    Geodetic.Cell2DsNumVertices[faces_id] = 3;
					Geodetic.Cell2DsNumEdges[faces_id] = 3;
					//vector<int> Vertices = {v1, v2, v3};
					Vertices[2] = v4; // va bene solo questo secondo i miei calcoli.
					Geodetic.Cell2DsVertices[faces_id] = Vertices;
				}*/

				//secondo me questa condizione continua ad essere sbagliata, perchè no esclude che dei vertici adiacenti sul lato possano essere presi
                if (i > 0)     // per i=0 non funziona: i triangoli a punta in giù non toccano i vertici della faccia che stanno triangolando, ma solo punti interni dei lati
                {
                    int v4 = coefficients[{segments - (i - 1) - j, i - 1, j, id}];
/* NOTA X ALE (anche3 abbastanza inutile, fatta la sera a mente poco lucida)(per come è fatto il mio disegno io avrei
if (j (=c (?)) >0)
allora v4 = coefficients[{a,b+1,c-1,Id}] 
ove a = segments -i -j;
	b = i 
	c=j, ove la c indica il livello (altezza) sulla faccia triangolarizzata, ma ficnhè non la immmagino come ha fatto lui per capire le differenze
	non la cambio, anche se dovrebbe funzionare ugualmente dovessi riscrivere il tutto coerentemente)*/


					/*i - 1: stai “salendo” di una riga nella griglia triangolare
					j: stessa colonna
					segments - (i - 1) - j: mantiene la condizione a + b + c = segments
					In pratica, v4 è il vertice interno che chiude il triangolo capovolto (v. sotto*/
					Geodetic.NumCell2Ds++;
                    Geodetic.Cell2DsId.push_back(faces_id);
                    Geodetic.Cell2DsNumVertices[faces_id] = 3;
					Geodetic.Cell2DsNumEdges[faces_id] = 3;
					Vertices[2] = v4;
					Geodetic.Cell2DsVertices[faces_id] = Vertices;
                    Geodetic.Cell2DsEdges[faces_id].resize(3);

                    for (int k = 0; k < 3; k++)   // stessa struttura di prima
                    {
                        int v_start = Geodetic.Cell2DsVertices[faces_id][k];
                        //int v_end = Geodetic.Cell2DsVertices[faces_id][(k + 1) % 3]; CHI L'HA FATTA QUESTA COSA CON MODULO 3? NOI? IO NO (ale)
																					//riassume l'if e l'else fatti sotto.
						int v_end;
						
						if (k == 2)
							v_end = Geodetic.Cell2DsVertices[faces_id][0];
						else
							v_end = Geodetic.Cell2DsVertices[faces_id][k+1];
						
                        if (!CheckEdges(Geodetic.Cell1DsVertices, v_start, v_end, edges_id, existing_edge_id))
                        {
                            Geodetic.Cell1DsId.push_back(edges_id);
                            // Geodetic.Cell1DsVertices.push_back(Vector2i(v_start, v_end));
							Geodetic.Cell1DsVertices(0, edges_id) = v_start;
							Geodetic.Cell1DsVertices(1, edges_id) = v_end;
							Geodetic.Cell2DsEdges[faces_id][k] = edges_id;
                            Geodetic.NumCell1Ds++;
                            edges_id++;
                        }
						else
							Geodetic.Cell2DsEdges[faces_id][k] = existing_edge_id; 
                    }

                    // Geodetic.NumCell2Ds++;
                    faces_id++;
                }
            }
        }
    }
	Geodetic.Cell1DsVertices.conservativeResize(2, Geodetic.NumCell1Ds);
	Geodetic.Cell2DsNumVertices.resize(Geodetic.NumCell2Ds);
	Geodetic.Cell2DsNumEdges.resize(Geodetic.NumCell2Ds);
	Geodetic.Cell2DsVertices.resize(Geodetic.NumCell2Ds);
	Geodetic.Cell2DsEdges.resize(Geodetic.NumCell2Ds);
}

/************************************************************************************************/

void DualMesh(PolyhedralMesh& InputMesh, PolyhedralMesh& DualMesh)
{
    // Identificatori iniziali per il duale
    int center_id = 0;  // Identificatore per i baricentri (nuovi vertici)
    int new_face_id = 0;  // Identificatore per le nuove facce
    int new_edge_id = 0;  // Identificatore per gli spigoli
	int existing_edge_id = 0; // Variabile per salvare l'edge esistente

    // Il duale ha tanti vertici quante facce ha il poliedro originale
    DualMesh.NumCell0Ds = InputMesh.NumCell2Ds;
    DualMesh.Cell0DsId.resize(DualMesh.NumCell0Ds);
    //DualMesh.Cell0DsCoordinates.resize(DualMesh.NumCell0Ds);
	DualMesh.Cell0DsCoordinates = MatrixXd::Zero(3,DualMesh.NumCell0Ds);

    // Gli spigoli nel duale rimangono in numero uguale a quelli del poliedro originale
    DualMesh.NumCell1Ds = InputMesh.NumCell1Ds;
	DualMesh.Cell1DsId.resize(DualMesh.NumCell1Ds);
	// DualMesh.Cell1DsVertices.resize(DualMesh.NumCell1Ds);

    // Il duale ha tante facce quante sono i vertici del poliedro originale
    DualMesh.NumCell2Ds = InputMesh.NumCell0Ds;
	DualMesh.Cell1DsVertices = MatrixXi::Zero(2, DualMesh.NumCell1Ds);
    DualMesh.Cell2DsId.resize(DualMesh.NumCell2Ds);
    DualMesh.Cell2DsVertices.resize(DualMesh.NumCell2Ds);
    DualMesh.Cell2DsEdges.resize(DualMesh.NumCell2Ds);

    // Mappa per associare le facce del poliedro ai baricentri nel duale (tramite gli id)
    map<int, int> FaceCenters;
	// Itera su tutte le facce del poliedro
    for (const auto& face_id : InputMesh.Cell2DsId) {
        Vector3d centro;  // Inizializza il baricentro

		Vector3d V1 = InputMesh.Cell0DsCoordinates.col(InputMesh.Cell2DsVertices[face_id][0]);
		Vector3d V2 = InputMesh.Cell0DsCoordinates.col(InputMesh.Cell2DsVertices[face_id][1]);
		Vector3d V3 = InputMesh.Cell0DsCoordinates.col(InputMesh.Cell2DsVertices[face_id][2]);
        
		// Itera su tutti i vertici della faccia corrente per calcolarne il baricentro
        /*for (size_t j = 0; j < InputMesh.Cell2DsVertices[face_id].size(); j++) {
			// Somma le coordinate di ciascun vertice della faccia
            centro_id += InputMesh.Cell0DsCoordinates[InputMesh.Cell2DsVertices[face_id][j]];
        }*/

		// Divide la somma delle coordinate per il numero di vertici per ottenere la media
        //centro_id /= InputMesh.Cell2DsVertices[face_id].size(); 

		centro = (1.0/3.0)*V1 + (1.0/3.0)*V2 + (1.0/3.0)*V3;
        // Salva l'ID del baricentro nel duale
        DualMesh.Cell0DsId.push_back(center_id);
        // DualMesh.Cell0DsCoordinates[center_id] = centro_id;
		DualMesh.Cell0DsCoordinates(0, face_id) = centro(0);
		DualMesh.Cell0DsCoordinates(1, face_id) = centro(1);
		DualMesh.Cell0DsCoordinates(2, face_id) = centro(2);
		
        // Associa l'ID del baricentro alla faccia originale
        FaceCenters[face_id] = center_id;

		// Incrementa il contatore per il prossimo baricentro
        center_id++;
    }

	// Ora creiamo le facce del duale utilizzando i baricentri trovati sopra
	// Itera attraverso tutti gli id dei vertici del poliedro originale, perchè ogni faccia nel poliedro duale corrisponde ad un vertice dell'originale
    for (const auto& vertex_id : InputMesh.Cell0DsId) {
        vector<int> AdjacentFaces;  //vettore che conterrà gli ID delle facce adiacenti al vertice corrente

        // Trova tutte le facce che contengono questo vertice
		// Itera sulle facce del poliedro originale
        for (const auto& face_id : InputMesh.Cell2DsId) {
			// Scorre tutti i vertici della faccia corrente (face_id)
            for (int j = 0; j < 3; j++) {
				/* Confronta vertex_id con ogni vertice della faccia corrente,
				se vero face_id contiene vertex_id, allora aggiungo l'id della faccia (face_id) ad AdjacentFaces*/ 
                if (InputMesh.Cell2DsVertices[face_id][j] == vertex_id) {
                    AdjacentFaces.push_back(face_id);
                    break;
                }
            }
        }

        // Le facce adiacenti non sono necessariamente ordinate! Le ordiniamo con una funzione dedicata
        vector<int> SortedFaces;
        Sort_Faces(AdjacentFaces, SortedFaces, InputMesh);

		int SortedFacesSize = SortedFaces.size();     // dava warning perchè .size() fornisce un tipo size_t, mentre k è int. così è a posto.
        // Creiamo una nuova faccia nel duale, che sarà composta dai baricentri delle facce trovate sopra
        vector<int> NewDualVertices; // conterrà i vertici della nuova faccia nel duale
		//Itera attraverso gli id di tutte le facce ordinate
        for (const auto& sorted_face : SortedFaces)
            NewDualVertices.push_back(FaceCenters[sorted_face]); // usa gli id dei baricentri

		// Assegna un ID alla nuova faccia nel poliedro duale
        DualMesh.Cell2DsId.push_back(new_face_id); //DualMesh.Cell2DsId[new_face_id] = new_face_id;
		// Salva i vertici della nuova faccia nel poliedro duale
        DualMesh.Cell2DsVertices[new_face_id] = NewDualVertices;
		// Alloco spazio per gli spigoli
        DualMesh.Cell2DsEdges[new_face_id].resize(SortedFacesSize);

        // Itera su tutti i vertici della faccia nel duale per determinare gli spigoli
		for (int k = 0; k < SortedFacesSize; k++) {
			int d_start = DualMesh.Cell2DsVertices[new_face_id][k];  // Vertice iniziale dell'edge
			// Se siamo all'ultimo vertice della faccia, il vertice finale deve essere il primo per chiudere la faccia
			int d_end = DualMesh.Cell2DsVertices[new_face_id][k+1];//[(k + 1) % SortedFaces.size()];

			// Verifica che l'edge non sia un duplicato e, se necessario, lo aggiunge
			if (!CheckEdges(DualMesh.Cell1DsVertices, d_start, d_end, new_edge_id, existing_edge_id)) {
				// Se l'edge non esiste, lo aggiungiamo al poliedro duale
				DualMesh.Cell1DsId.push_back(new_edge_id);  // Assegna un nuovo ID allo spigolo
				//DualMesh.Cell1DsVertices.push_back(Vector2i(d_start, d_end));  // Salva gli estremi dello spigolo
				DualMesh.Cell1DsVertices(0, new_edge_id) = d_start;
				DualMesh.Cell1DsVertices(1, new_edge_id) = d_end;
				DualMesh.Cell2DsEdges[new_face_id][k] = new_edge_id;  // Associa lo spigolo alla faccia attuale
				new_edge_id++;  // Incrementa l'ID per il prossimo spigolo
			} else {
				// Se l'edge già esiste, usa l'ID di quello già presente
				DualMesh.Cell2DsEdges[new_face_id][k] = existing_edge_id;
			}
		}

        new_face_id++;
    }

    // Proiezione sulla sfera per mantenere una rappresentazione geometrica corretta
    Projection(DualMesh);
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

//bool CheckVertices(const vector<Vector3d>& coords, const Vector3d& point, int current_id, int& duplicate_id)  // prima il penultimo parametro era int current_id)
bool CheckVertices(const MatrixXd& mesh, const Vector3d& point, int& dimension, int& duplicate_id)
{
	double eps = 1e-8;
    for (int i = 0; i < dimension; i++)
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

// bool CheckEdges(const vector<Vector2i>& edges, int v1, int v2, int& current_edge_id, int& existing_edge_id) //v1 e v2 sono i vertici dello spigolo con il current_edge_id
bool CheckEdges(const MatrixXi& verts, const int& v1, const int& v2, int& dimension, int& existing_edge_id)
{
    for (int i = 0; i < dimension; i++)  // i < current_edge_id
    {
        int a = verts(0,i);//int a = edges[i][0];
        int b = verts(1,i);// int b = edges[i][1];
        if ((a == v1 && b == v2) || (a == v2 && b == v1)) 
		// //Questo if controlla se la coppia di variabili a e b è uguale alla coppia v1 e v2, in qualunque ordine.
        {
			existing_edge_id = i;
            return true;
        }
    }
    return false;
}

/************************************************************************************************/

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
}
