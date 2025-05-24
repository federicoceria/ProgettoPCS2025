#include "Utils.hpp"
#include <iostream>
#include <sstream>
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace PolyhedralLibrary;

int main()
{
    PolyhedralMesh Platonic; //si chiamava mesh;
    PolyhedralMesh Geodetic;
	PolyhedralMesh Goldberg;
    string path = "/home/appuser/Data/ProgettoPCS2025/Platonic_solids";
	string b_prov;
	string c_prov;

    /*prende in input path e mesh, chiede con un cin i parametri di schlafli, e in base a quelli, grazie 
    alla funzione PolyhedralChoice, modifica il path*/
    if(!ParameterSelection(path, Platonic))
	{
		cerr << "An error occurred while selecting the parameters" << endl;
		return 1;
	}

    /*con questa funzione, che al suo interno chiama anche Cell0/1/2Ds, salvo i dati nella mesh.
    Inoltre questa controlla anche che il salvataggio avvenga correttamente*/
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
	/*
    if(b > 0 && c == 0)
    {
        GeodeticPolyhedron(Platonic, Geodetic, b);
		DualMesh(Geodetic, Goldberg);
        cout << "Geodetic solid successfully generated with num_segments = " << b << endl;
    }
    else if(b == 0 && c > 0)
    {
        GeodeticPolyhedron(Platonic, Geodetic, c);
		DualMesh(Geodetic, Goldberg);
        cout << "Geodetic solid successfully generated with num_segments = " << c << endl;
    }
    //else
        //seconda classe
	
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

    cout << "andate tutti al concerto dei POLIFONICI" << endl;
    return 0;
}




