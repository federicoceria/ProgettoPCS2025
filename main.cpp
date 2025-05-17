#include "Utils.hpp"
#include <iostream>
#include <sstream>
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace PolyhedralLibrary;

int main()
{
    PolyhedralMesh mesh;
    string path = "/home/appuser/Data/ProgettoPCS2025/Platonic_solids";

    /*poliedral choice modifica il path per referenza in base alla q che ho passato per parametro
    PolyhedralChoice(const string& path, 
						PolyhedralMesh& mesh, 
						const char& p, 
						const char& q, 
						const char& b, 
						const char& c,
						bool& walk);*/ 

    /*prende in input path e mesh, chiede con un cin i paramtri di schlafi, e in base a quelli, grazie 
    alla funzione polyhedralchoice, modifico il path*/
    ParameterSelection(path, mesh);

    /*con questa funzione, che al suo interno chiama anche Cell0/1/2Ds, salvo i dati nella mesh.
    Inoltre questa controlla anche che il salvataggio avvenga correttamente*/
    ImportMesh(path, mesh);


    
    














    //! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
    /*string solidType;
    cout << "Which solid do you want to generate? (tetrahedron, octahedron, icosahedron): ";
    cin >> solidType;

    if (!PolyhedralChoice(path, inputMesh, p, q, b, c, walk)) {
        cerr << "Error importing selected solid." << endl;
        return 1;
    }

    if (!GeodeticPolyhedron(inputMesh, geodetic, num_segments)) {
        cerr << "Error generating the geodetic solid." << endl;
        return 1;
    }    questa parte è commentata in quanto da problemi con l'esecuzione del codice, c'è un problema con
        la firma della funzione, ovvero coi parametri che le vengono passati, solidtype dovrebbe essere
        un valore booleano, ma a quel punto perde di senso. per me solidtype si potrebbe anche eliminare
        e si potrebbero riscrivere queste righe in modo dierso, ma non lo faccio in quanto non lo ho scritto
        io e non vorre che servisse nell'idea di qualcuno a fare qualcosa.
    //! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

    cout << "Geodetic solid successfully generated with num_segments = " << num_segments << endl;
    cout << "Number of vertices: " << geodetic.NumCell0Ds << endl;
    cout << "Number of edges: " << geodetic.NumCell1Ds << endl;
    cout << "Number of faces: "   << geodetic.NumCell2Ds << endl;*/

    cout << "francesca barra" << endl;
    return 0;
}




