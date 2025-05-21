#include "Utils.hpp"
#include <iostream>
#include <sstream>
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace PolyhedralLibrary;

int main()
{
    PolyhedralMesh mesh;
    PolyhedralMesh Geodetic;
    string path = "/home/appuser/Data/ProgettoPCS2025/Platonic_solids";

    /*prende in input path e mesh, chiede con un cin i parametri di schlafi, e in base a quelli, grazie 
    alla funzione polyhedralchoice, modifico il path*/
    ParameterSelection(path, mesh);

    /*con questa funzione, che al suo interno chiama anche Cell0/1/2Ds, salvo i dati nella mesh.
    Inoltre questa controlla anche che il salvataggio avvenga correttamente*/
    ImportMesh(path, mesh);
    
    /*ricordiamoci di controllare gli imput mettendo l'errore se b e c !=0* /
    if(b==0 && c!=0)
    {
        GeodeticPolyhedron(mesh, Geodetic, c);
        cout << "Geodetic solid successfully generated with num_segments = " << c << endl;
        cout << "Number of vertices: " << Geodetic.NumCell0Ds << endl;
        cout << "Number of edges: " << Geodetic.NumCell1Ds << endl;
        cout << "Number of faces: "   << Geodetic.NumCell2Ds << endl;
    }
    else    //sarà da mettere così (c==0 && b!=0)
    {
        GeodeticPolyhedron(mesh, Geodetic, b);
        cout << "Geodetic solid successfully generated with num_segments = " << b << endl;
        cout << "Number of vertices: " << Geodetic.NumCell0Ds << endl;
        cout << "Number of edges: " << Geodetic.NumCell1Ds << endl;
        cout << "Number of faces: "   << Geodetic.NumCell2Ds << endl;
    }*/
    //else
        //seconda classe


    



    
    














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




