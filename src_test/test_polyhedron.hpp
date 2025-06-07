#pragma once

#include <iostream>
#include <gtest/gtest.h>
#include <algorithm>
#include "Utils.hpp"

using namespace PolyhedralLibrary;

const int segments = 5;

// questa calcola il numero di vertici che hanno una certa valenza data in input
int VertexDegree(int& ExpectedDegree, const std::vector<int>& Vertices, const std::vector<std::vector<int>>& Faces) {
	int NumVertexOfDegree = 0;
	for (const auto& Vertex : Vertices) {
		int VertexDegree = 0;
		for (const auto& listOfVertex : Faces ) //listOfVertex=vettore dei vertici che formano una faccia.
			if (std::find(listOfVertex.begin(), listOfVertex.end(), Vertex) != listOfVertex.end())
				VertexDegree++;
		if (VertexDegree == ExpectedDegree)
			NumVertexOfDegree++;
	}
	
	return NumVertexOfDegree;
}		

//SE NON SERVE CANCELLARE!!!
/*array<int, 3> SolidProperties(const int& V, const int& E, const int& F, const int& TriangulationParameter) {
	int Exp_V = 2*segments*segments +2;
	int Exp_E = ;
	int Exp_F = F * ( 3 * segments * segments + 3 * segments);
	array<int, 3> solid_properties = {Exp_V, Exp_E, Exp_F};
	return solid_properties;*/
 
//testa se la funzione geodetica associata al tetraedro funziona
TEST(TestGeodeticPolyhedron, TestTetrahedronTriangulation)
{
	PolyhedralMesh Platonic;
	if (!ImportMesh("../Platonic_solids/Tetraedro/", Platonic))
		FAIL() << "Something went wrong during the importation of the mesh: Platonic";

	PolyhedralMesh Geodetic;
	GeodeticPolyhedron(Platonic, Geodetic, segments);
	
	int T = segments*segments;
	int ExpectedVertices = 2*T + 2;
	int ExpectedEdges = 6*T;
	int ExpectedFaces = 4*T;
	
	EXPECT_EQ(Geodetic.NumCell0Ds, ExpectedVertices);
	EXPECT_EQ(Geodetic.NumCell1Ds, ExpectedEdges);
	EXPECT_EQ(Geodetic.NumCell2Ds, ExpectedFaces);
	
	int Degree1 = 3;
	int Degree2 = 6;
	int ExpectedNumVertexDegree1 = 4;
	int ExpectedNumVertexDegree2 = 2*(T-1);
	
	int NumVertexOfDegree1 = VertexDegree(Degree1, Geodetic.Cell0DsId, Geodetic.Cell2DsVertices);
	int NumVertexOfDegree2 = VertexDegree(Degree2, Geodetic.Cell0DsId, Geodetic.Cell2DsVertices);
	
	EXPECT_EQ(ExpectedNumVertexDegree1, NumVertexOfDegree1);
	EXPECT_EQ(ExpectedNumVertexDegree2, NumVertexOfDegree2);
}

//testa se la funzione geodetica associata al' ottaedro funziona
TEST(TestGeodeticPolyhedron, TestOctahedronTriangulation)
{
	PolyhedralMesh Platonic;
	if (!ImportMesh("../Platonic_solids/Ottaedro/", Platonic))
		FAIL() << "Something went wrong during the importation of the mesh: Platonic";

	PolyhedralMesh Geodetic;
	GeodeticPolyhedron(Platonic, Geodetic, segments);
	
	int T = segments*segments;
	int ExpectedVertices = 4*T + 2;
	int ExpectedEdges = 12*T;
	int ExpectedFaces = 8*T;
	
	EXPECT_EQ(Geodetic.NumCell0Ds, ExpectedVertices);
	EXPECT_EQ(Geodetic.NumCell1Ds, ExpectedEdges);
	EXPECT_EQ(Geodetic.NumCell2Ds, ExpectedFaces);
	
	int Degree1 = 4;
	int Degree2 = 6;
	int ExpectedNumVertexDegree1 = 6;
	int ExpectedNumVertexDegree2 = 4*(T-1);
	
	int NumVertexOfDegree1 = VertexDegree(Degree1, Geodetic.Cell0DsId, Geodetic.Cell2DsVertices);
	int NumVertexOfDegree2 = VertexDegree(Degree2, Geodetic.Cell0DsId, Geodetic.Cell2DsVertices);
	
	EXPECT_EQ(ExpectedNumVertexDegree1, NumVertexOfDegree1);
	EXPECT_EQ(ExpectedNumVertexDegree2, NumVertexOfDegree2);
}
//testa se la funzione geodetica associata all' icoesaedro funziona
TEST(TestGeodeticPolyhedron, TestIcosahedronTriangulation)
{
	PolyhedralMesh Platonic;
	if (!::ImportMesh("../Platonic_solids/Icosaedro/", Platonic))
	FAIL() << "Something went wrong during the importation of the mesh: Platonic";
	PolyhedralMesh Geodetic;
	GeodeticPolyhedron(Platonic, Geodetic, segments);
	
	int T = segments*segments;
	int ExpectedVertices = 10*T + 2;
	int ExpectedEdges = 30*T;
	int ExpectedFaces = 20*T;
	
	EXPECT_EQ(Geodetic.NumCell0Ds, ExpectedVertices);
	EXPECT_EQ(Geodetic.NumCell1Ds, ExpectedEdges);
	EXPECT_EQ(Geodetic.NumCell2Ds, ExpectedFaces);
	
	int Degree1 = 5;
	int Degree2 = 6;
	int ExpectedNumVertexDegree1 = 12;
	int ExpectedNumVertexDegree2 = 10*(T-1);
	
	int NumVertexOfDegree1 = VertexDegree(Degree1, Geodetic.Cell0DsId, Geodetic.Cell2DsVertices);
	int NumVertexOfDegree2 = VertexDegree(Degree2, Geodetic.Cell0DsId, Geodetic.Cell2DsVertices);
	
	EXPECT_EQ(ExpectedNumVertexDegree1, NumVertexOfDegree1);
	EXPECT_EQ(ExpectedNumVertexDegree2, NumVertexOfDegree2);
}

//FUNZIONE TEST FEDE
TEST(TestDualPolyhedron, TestType1)
{
	PolyhedralMesh Platonic;
	if (!::ImportMesh("../Platonic_solids/Tetraedro/", Platonic))
		FAIL() << "Something went wrong during the creation of the platonic polyhedron mesh";

	PolyhedralMesh Geodetic;
	GeodeticPolyhedron(Platonic, Geodetic, segments);
	
	PolyhedralMesh DualPolyhedron;
	DualMesh(Geodetic, DualPolyhedron);
	
	int T = segments*segments;
	int ExpectedVertices = 4*T;
	int ExpectedFaces = 2*T+2;
	
	EXPECT_EQ(DualPolyhedron.NumCell0Ds, ExpectedVertices);
	EXPECT_EQ(DualPolyhedron.NumCell2Ds, ExpectedFaces);
}

//PER CAPIRE QUESTE DEVI VEDERE IL FUNZIONAMENTO DI ORDERED_FACES-FEDE
TEST(TestSortFaces, Test_unordered)
{
	PolyhedralMesh Platonic;
	if (!::ImportMesh("../Platonic_solids/Ottaedro/", Platonic))
		FAIL() << "Something went wrong during the creation of the platonic polyhedron mesh";
	
	vector<int> unsorted_faces = {4, 5, 3, 2};
	vector<int> sorted_faces;
	
	Sort_Faces(unsorted_faces, sorted_faces, Platonic);
	vector<int>expected_sorted_faces  = {4, 5, 2, 3};
	
	EXPECT_EQ(sorted_faces, expected_sorted_faces);
}

//PER CAPIRE QUESTE DEVI VEDERE IL FUNZIONAMENTO DI ORDERED_FACES-FEDE
TEST(TestSortFaces, Test_sorted)
{
	PolyhedralMesh Platonic;
	if (!::ImportMesh("../Platonic_solids/Icosaedro/", Platonic))
		FAIL() << "Something went wrong during the creation of the platonic polyhedron mesh";
	
	vector<int> unsorted_faces = {10, 11, 12, 17, 16};
	vector<int> sorted_faces;
	
	Sort_Faces(unsorted_faces, sorted_faces, Platonic);
	vector<int>expected_sorted_faces  = {10, 11, 12, 17, 16};
	
	EXPECT_EQ(sorted_faces, expected_sorted_faces);
}

//FUNZIONE GABRI
TEST(TestShortestPath, ShortestPath)
{
	
	PolyhedralMesh Platonic;
	if (!::ImportMesh("../Platonic_solids/Ottaedro/", Platonic))
		FAIL() << "Something went wrong during the creation of the platonic polyhedron mesh";
	
	PolyhedralMesh Geodetic;
	GeodeticPolyhedron(Platonic, Geodetic, 2);
	
	double path_length;
	int edges_in_path;
	vector<int> path_vertices;	
	MatrixXd W = MatrixXd::Zero(Geodetic.NumCell0Ds, Geodetic.NumCell0Ds);
	
	if(!ShortestPath(Geodetic, 3, 7, path_vertices, W))
		FAIL() << "Something went wrong during the execution of ShortestPath function";
	
	vector<int> expected_path = {7,3};
	
	EXPECT_EQ(path_vertices, expected_path);
	/*EXPECT_EQ(edges_in_path, 1);
	EXPECT_NEAR(0.765367, path_length, 1e-6);*/
}
//ANCORA GABRI
TEST(TestShortestPath, ShortestPathDual)
{
	PolyhedralMesh Platonic;
	PolyhedralMesh DualPolyhedron;
	
	if (!::ImportMesh("../Platonic_solids/Tetraedro/", Platonic))
		FAIL() << "Something went wrong during the creation of the platonic polyhedron mesh";
	
	PolyhedralMesh Geodetic;
	
	GeodeticPolyhedron(Platonic, Geodetic, 3);
	DualMesh(Geodetic, DualPolyhedron);
	
	double path_length;
	int number_edges_in_path;
	vector<int> path_vertices;	
	MatrixXd W = MatrixXd::Zero(DualPolyhedron.NumCell0Ds, DualPolyhedron.NumCell0Ds);
	
	if(!ShortestPath(DualPolyhedron, 4, 8, path_vertices, W))
		FAIL() << "Something went wrong during the execution of ShortestPath function";
	
	vector<int>expected_path = {8, 3, 4};
	
	EXPECT_EQ(path_vertices, expected_path);
	/*EXPECT_EQ(number_edges_in_path, 3);
	EXPECT_NEAR(1.316185, path_length,1e-6);*/
}
