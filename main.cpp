#include "Utils.hpp"
#include <iostream>
#include <sstream>
#include "PolyhedralMesh.hpp"
#include "UCDUtilities.hpp" 

using namespace std;
using namespace PolyhedralLibrary;

int main()
{
    PolyhedralMesh Platonic; //si chiamava mesh;
    PolyhedralMesh Geodetic;
	PolyhedralMesh Goldberg;
    string polyhedron;
    string path = "/home/appuser/Data/ProgettoPCS2025/Platonic_solids";
    string p_prov;
	string q_prov;
	string b_prov;
	string c_prov;
	string id1_prov;
	string id2_prov;

    /*prende in input path e mesh, chiede con un cin i parametri di schlafli, e in base a quelli, grazie 
    alla funzione PolyhedralChoice, modifica il path*/
    /*if(!PolyhedralChoice(path))  // if(!ParameterSelection(path, Platonic)); ho tolto Platonic perché di fatto polyhedralchoice non la usa mai
	{
		cerr << "An error occurred while selecting the parameters: in order to generate a geodetic solid, p has to be equal to 3 and q has to be chosen between the values 3,4,5." << endl;
		return 1;
	}

    /*con questa funzione, che al suo interno chiama anche Cell0/1/2Ds, salvo i dati nella mesh.
    Inoltre questa controlla anche che il salvataggio avvenga correttamente*/

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
				return 1;
                break;
		}

        path = path + polyhedron;
    }

    if(!ImportMesh(path, Platonic))
	{
		cerr << "An error occurred while importing the mesh" << endl;
		return 1;
	}
    
    /*ricordiamoci di controllare gli input mettendo l'errore se b e c !=0*/
	cout << "Insert b: " << endl;
	cin >> b_prov;
	cout << "Insert c: " << endl;
	cin >> c_prov;
	
    int c = stoi(b_prov);
    int b = stoi(c_prov);
	
    if(b > 0 && c == 0)
    {
        GeodeticPolyhedron(Platonic, Geodetic, b);
		//DualMesh(Geodetic, Goldberg);
        cout << "Geodetic solid successfully generated with num_segments = " << b << endl;
    }
    else if(b == 0 && c > 0)
    {
        GeodeticPolyhedron(Platonic, Geodetic, c);
		//DualMesh(Geodetic, Goldberg);
        cout << "Geodetic solid successfully generated with num_segments = " << c << endl;
    }
	else
    {
        cerr << "Error: Invalid parameters for geodetic solid generation." << endl;
        return 1;
    }
    
    // GENERAZIONE DEL POLIEDRO DUALE SE q = 3
    if (q == 3) 
    {
        DualMesh(Geodetic, Goldberg);
        Geodetic = Goldberg;
        cout << "Generated a Goldberg polyhedron {3+, 3} for segments (" << b << ", " << c << ")." << endl;
    } 
    else
    {
        cout << "Generated a geodetic polyhedron {3, " << q << "+} for segments (" << b << ", " << c << ")." << endl;
    }

	cout << "Shortest Path: please, insert the starting id. If you do not want to evaluate it, please enter without inserting anything." << endl;
	cin >> id1_prov;
	cout << "Now, please insert the ending id. As the previous one, if you do not want to evaluate it, please enter without inserting anything." << endl;
	cin >> id2_prov;
	
	int starting_id;
	int ending_id;
	
	if(id1_prov == "" && id2_prov == "")
	{
		cout << "The shortest path will not be evaluated." << endl;
	}
	else if (id1_prov != "" && id2_prov == "")
	{
		starting_id = stoi(id1_prov);
		ending_id = stoi(id2_prov);
	}
	
	double length = 0.0;
	int PathEdges = 0;
	vector<int> visited;
	ShortestPath(Geodetic, starting_id, ending_id, length, PathEdges, visited);
    
	// Esportazione file TXT
    if (!ExportPolyhedralData(Geodetic)) {
        cerr << "Error: Failed to export polyhedral data." << endl;
        return 1;
    }

	// Esportazione per Paraview (solo vertici e spigoli)
    Gedim::UCDUtilities utilities;
    utilities.ExportPoints("./Cell0Ds.inp", Geodetic.Cell0DsCoordinates);
    utilities.ExportSegments("./Cell1Ds.inp", Geodetic.Cell0DsCoordinates, Geodetic.Cell1DsVertices);

    cout << "Export completed successfully!" << endl;
    return 0;
}
	/*Gedim::UCDUtilities utilities;	
    utilities.ExportPoints("./Cell0Ds.inp",
                           Geodetic.Cell0DsCoordinates);

    utilities.ExportSegments("./Cell1Ds.inp",
								Geodetic.Cell0DsCoordinates,
								Geodetic.Cell1DsVertices); 

    utilities.ExportPoints("./Cell0DsDual.inp",
                           Goldberg.Cell0DsCoordinates);

    utilities.ExportSegments("./Cell1DsDual.inp",
								Goldberg.Cell0DsCoordinates,
								Goldberg.Cell1DsVertices); 
    
*/


    
    














    //! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
    /*string solidType;
    cout << "Which solid do you want to generate? (tetrahedron, octahedron, icosahedron): ";
    cin >> solidType;

    if (!PolyhedralChoice(path, inputMesh, p, q, b, c, walk)) {
        cerr << "Error importing selected solid." << endl;
        return 1;
    }

    if (!GeodeticPolyhedron(inputMesh, Geodetic, num_segments)) {
        cerr << "Error generating the Geodetic solid." << endl;
        return 1;
    }    questa parte è commentata in quanto da problemi con l'esecuzione del codice, c'è un problema con
        la firma della funzione, ovvero coi parametri che le vengono passati, solidtype dovrebbe essere
        un valore booleano, ma a quel punto perde di senso. per me solidtype si potrebbe anche eliminare
        e si potrebbero riscrivere queste righe in modo dierso, ma non lo faccio in quanto non lo ho scritto
        io e non vorre che servisse nell'idea di qualcuno a fare qualcosa.
    //! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

    cout << "Geodetic solid successfully generated with num_segments = " << num_segments << endl;
    cout << "Number of vertices: " << Geodetic.NumCell0Ds << endl;
    cout << "Number of edges: " << Geodetic.NumCell1Ds << endl;
    cout << "Number of faces: "   << Geodetic.NumCell2Ds << endl;*/

 /*   cout << "andate tutti al concerto dei POLIFONICI" << endl;
    return 0;
}*/




