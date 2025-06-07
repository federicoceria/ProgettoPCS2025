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

	int p;
	int q;
	
    cout << "Insert p: " << endl;
	cin >> p_prov;
	if(!isInteger(p_prov))
	{
		cerr << "The inserted value for p is not valid." << endl;
		return 1;
	}
	else
		p = stoi(p_prov);
	
	cout << "Insert q: " << endl;
	cin >> q_prov;
	if(!isInteger(q_prov))
	{
		cerr << "The inserted value for q is not valid." << endl;
		return 1;
	}
	else
		q = stoi(q_prov);


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
	else
	{
		cerr << "The inserted value for p should be equal to 3. It is not possible to generate a geodetic polyhedron with p = " << p << "." << endl;
		return 1;
	}
	
    if(!ImportMesh(path, Platonic))
	{
		cerr << "An error occurred while importing the mesh" << endl;
		return 1;
	}
    
    /*ricordiamoci di controllare gli input mettendo l'errore se b e c !=0*/
	
	int b;
	int c;
	
	cout << "Insert b: " << endl;
	cin >> b_prov;
	if(!isInteger(b_prov))
	{
		cerr << "The inserted value for b is not valid." << endl;
		return 1;
	}
	else
		b = stoi(b_prov);
	
	cout << "Insert c: " << endl;
	cin >> c_prov;
	if(!isInteger(c_prov))
	{
		cerr << "The inserted value for c is not valid." << endl;
		return 1;
	}
	else
		c = stoi(c_prov);
    
	
    if(b > 0 && c == 0)
    {
        GeodeticPolyhedron(Platonic, Geodetic, b);
        cout << "Geodetic solid of the first class successfully generated with num_segments = " << b << endl;
    }
    else if(b == 0 && c > 0)
    {
        GeodeticPolyhedron(Platonic, Geodetic, c);
        cout << "Geodetic solid of the first class successfully generated with num_segments = " << c << endl;
    }
	else
    {
        cerr << "Error: Invalid parameters for type 1 geodetic solid generation." << endl;
        return 1;
    }
    
    // GENERAZIONE DEL POLIEDRO DUALE SE q = 3
    if (q == 3) 
    {
        DualMesh(Geodetic, Goldberg);
        //fatto per evitare di fare degli if dopo
        Geodetic = Goldberg;
        cout << "Goldberg polyhedron {3+, 3} for segments (" << b << ", " << c << ") generated succesfully." << endl;
    } 
    else
    {
        cout << "Geodetic polyhedron {3, " << q << "+} for segments (" << b << ", " << c << ") generated succesfully." << endl;
    }
	
	
	int starting_id;
	int ending_id;
	
	cout << "Shortest Path: please, insert the starting id. If you do not want to evaluate it, please enter the letter n." << endl;
	cin >> id1_prov;
	cout << "Now, please insert the ending id. As the previous one, if you do not want to evaluate it, please enter the letter n." << endl;
	cin >> id2_prov;

	if(id1_prov == "n" && id2_prov == "n")
	{
		cout << "The shortest path will not be evaluated." << endl;
	}
	else if (id1_prov != "n" && id2_prov != "n")
	{
		if(!isInteger(id1_prov))    // check if the two ids are integers.
		{
			cerr << "The inserted value for the starting Id is not valid." << endl;
			return 1;
		}
		else
			starting_id = stoi(id1_prov);
		
		if(!isInteger(id2_prov))
		{
			cerr << "The inserted value for the ending Id is not valid." << endl;
			return 1;
		}
		else
			ending_id = stoi(id2_prov);
		
		// check if the ids are out of range for the geodetic polyhedron
		if (starting_id >= Geodetic.NumCell0Ds || ending_id >= Geodetic.NumCell0Ds || starting_id < 0 || ending_id < 0)
		{
		cerr << "The inserted Ids are not valid." << endl;
		return 1;
		}
		
        double length = 0.0;
        int NumPath = 0;
        vector<int> visited;
		MatrixXd W = MatrixXd::Zero(Geodetic.NumCell0Ds, Geodetic.NumCell0Ds);
        ShortestPath(Geodetic, starting_id, ending_id, visited, W);    //  length, NumPath,

        cout << "ho eseguito l'if" << endl;


        vector<double> PathPointsProperties(Geodetic.NumCell0Ds, 0.0);
		for (const auto& point : visited)
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
								Geodetic.Cell0DsCoordinates,
								PointsProperties);
	
		vector<int> pathEdges; 
		vector<double> PathEdgesProperties(Geodetic.NumCell1Ds, 0.0);

		for (size_t i = 0; i < visited.size()-1; i++){
			int v1 = visited[i];
			int v2 = visited[i+1];
			for(const auto& edge: Geodetic.Cell1DsId){
				if ((Geodetic.Cell1DsVertices(0, edge) == v1 && Geodetic.Cell1DsVertices(1,edge) == v2) || 
					(Geodetic.Cell1DsVertices(0, edge) == v2 && Geodetic.Cell1DsVertices(1,edge) == v1))
					{
						pathEdges.push_back(edge);
						PathEdgesProperties[edge] = 1.0;
					}
			}	
		}
		
		for(size_t i = 0; i < visited.size()-1; i++)
			length += W(visited[i],visited[i+1]);
		NumPath = pathEdges.size();
		
		
		ShortPathProperty.Label = "shortest path";
		ShortPathProperty.UnitLabel = "";
		ShortPathProperty.Size = PathEdgesProperties.size();
		ShortPathProperty.NumComponents = 1;
		ShortPathProperty.Data = PathEdgesProperties.data();  
	
		vector<Gedim::UCDProperty<double>> EdgesProperties;
		EdgesProperties.push_back(ShortPathProperty);
		utilities.ExportSegments("./Cell1Ds.inp",
								Geodetic.Cell0DsCoordinates,
								Geodetic.Cell1DsVertices,
								{},
								EdgesProperties);
    }   
    else
	{
        cerr << "Error: inserted IDs are not valid." << endl;
    }

    
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

 /*   cout << "andate tutti al concerto dei POLIFONICI" << endl;*/

/*	return 0;
}  */




