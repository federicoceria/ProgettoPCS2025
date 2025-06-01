/* #pragma once

#include <iostream>
#include <gtest/gtest.h>
#include <algorithm>
#include "Utils.hpp"



const int segments = 5;

// questa calcola il numero di vertici che hanno una certa valenza
int VertexDegree(int& ExpectedDegree, const std::vector<int>& Vertices, const std::vector<std::vector<int>>& Faces) {
	int NumVertexOfExpectedDegree = 0;
	for (const auto& Vertex : Vertices) {
		int VertexDegree = 0;
		for (const auto& listOfVertex : Faces ) //listOfVertex=vettore dei vertici che formano una faccia.
			if (std::find(listOfVertex.begin(), listOfVertex.end(), Vertex) != listOfVertex.end())
				VertexDegree++;
		if (VertexDegree == ExpectedDegree)
			NumVertexOfExpectedDegree++;
	}
	
	return NumVertexOfExpectedDegree;
}		
/*NON SO A COSA SERVA STA ROBA
array<int, 3> SolidProperties(const int& V, const int& E, const int& F, const int& TriangulationParameter) {
	int Exp_V = V + E*(2*segments-1) + F * ( (3.0*segments*segments)/2.0 - (3.0*segments/2.0) + 1);
	int Exp_E = E * 2 * segments + F * ( (9.0*segments*segments)/2.0 + (3.0*segments/2.0) );
	int Exp_F = F * ( 3 * segments * segments + 3 * segments);
	array<int, 3> solid_properties = {Exp_V, Exp_E, Exp_F};
	return solid_properties;
} * /

TEST(TestGeodeticPolyhedron, TestTetrahedronTriangulation)
{
    //TIPO QUA SOTTO CI ANDREBBE POLYHEDRALMESH E NON POLIEDRON
	PolyhedronMesh PlatonicPolyhedron;
	if (!FileManagement::ImportPolyhedronMesh(PlatonicPolyhedron, "../Platonic_solids/Tetraedro/"))
		FAIL() << "Something went wrong during the creation of the platonic polyhedron mesh";

	PolyhedronMesh GeodeticPolyhedron;
	Generation::GeodeticSolidType1(PlatonicPolyhedron, GeodeticPolyhedron, segments);
	
	int T = segments*segments;
	int ExpectedVertices = 2*T + 2;
	int ExpectedEdges = 6*T;
	int ExpectedFaces = 4*T;
	
	EXPECT_EQ(GeodeticPolyhedron.NumCell0Ds, ExpectedVertices);
	EXPECT_EQ(GeodeticPolyhedron.NumCell1Ds, ExpectedEdges);
	EXPECT_EQ(GeodeticPolyhedron.NumCell2Ds, ExpectedFaces);
	
	int Degree1 = 3;
	int Degree2 = 6;
	int ExpectedNumDegree1 = 4;
	int ExpectedNumDegree2 = 2*(T-1);
	
	int NumDegree1 = VertexDegree(Degree1, GeodeticPolyhedron.Cell0DsId, GeodeticPolyhedron.Cell2DsVertices);
	int NumDegree2 = VertexDegree(Degree2, GeodeticPolyhedron.Cell0DsId, GeodeticPolyhedron.Cell2DsVertices);
	
	EXPECT_EQ(ExpectedNumDegree1, NumDegree1);
	EXPECT_EQ(ExpectedNumDegree2, NumDegree2);
}

TEST(TestGeodeticPolyhedron, TestOctahedronTriangulation)
{
	PolyhedronMesh PlatonicPolyhedron;
	if (!FileManagement::ImportPolyhedronMesh(PlatonicPolyhedron, "../Platonic_solids/Ottaedro/"))
		FAIL() << "Something went wrong during the creation of the platonic polyhedron mesh";

	PolyhedronMesh GeodeticPolyhedron;
	Generation::GeodeticSolidType1(PlatonicPolyhedron, GeodeticPolyhedron, segments);
	
	int T = segments*segments;
	int ExpectedVertices = 4*T + 2;
	int ExpectedEdges = 12*T;
	int ExpectedFaces = 8*T;
	
	EXPECT_EQ(GeodeticPolyhedron.NumCell0Ds, ExpectedVertices);
	EXPECT_EQ(GeodeticPolyhedron.NumCell1Ds, ExpectedEdges);
	EXPECT_EQ(GeodeticPolyhedron.NumCell2Ds, ExpectedFaces);
	
	int Degree1 = 4;
	int Degree2 = 6;
	int ExpectedNumDegree1 = 6;
	int ExpectedNumDegree2 = 4*(T-1);
	
	int NumDegree1 = VertexDegree(Degree1, GeodeticPolyhedron.Cell0DsId, GeodeticPolyhedron.Cell2DsVertices);
	int NumDegree2 = VertexDegree(Degree2, GeodeticPolyhedron.Cell0DsId, GeodeticPolyhedron.Cell2DsVertices);
	
	EXPECT_EQ(ExpectedNumDegree1, NumDegree1);
	EXPECT_EQ(ExpectedNumDegree2, NumDegree2);
}

TEST(TestGeodeticPolyhedron, TestIcosahedronTriangulation)
{
	PolyhedronMesh PlatonicPolyhedron;
	if (!FileManagement::ImportPolyhedronMesh(PlatonicPolyhedron, "../Platonic_solids/Icosaedro/"))
		FAIL() << "Something went wrong during the creation of the platonic polyhedron mesh";

	PolyhedronMesh GeodeticPolyhedron;
	Generation::GeodeticSolidType1(PlatonicPolyhedron, GeodeticPolyhedron, segments);
	
	int T = segments*segments;
	int ExpectedVertices = 10*T + 2;
	int ExpectedEdges = 30*T;
	int ExpectedFaces = 20*T;
	
	EXPECT_EQ(GeodeticPolyhedron.NumCell0Ds, ExpectedVertices);
	EXPECT_EQ(GeodeticPolyhedron.NumCell1Ds, ExpectedEdges);
	EXPECT_EQ(GeodeticPolyhedron.NumCell2Ds, ExpectedFaces);
	
	int Degree1 = 5;
	int Degree2 = 6;
	int ExpectedNumDegree1 = 12;
	int ExpectedNumDegree2 = 10*(T-1);
	
	int NumDegree1 = VertexDegree(Degree1, GeodeticPolyhedron.Cell0DsId, GeodeticPolyhedron.Cell2DsVertices);
	int NumDegree2 = VertexDegree(Degree2, GeodeticPolyhedron.Cell0DsId, GeodeticPolyhedron.Cell2DsVertices);
	
	EXPECT_EQ(ExpectedNumDegree1, NumDegree1);
	EXPECT_EQ(ExpectedNumDegree2, NumDegree2);
}

TEST(TestDualPolyhedron, TestType1)
{
	PolyhedronMesh PlatonicPolyhedron;
	if (!FileManagement::ImportPolyhedronMesh(PlatonicPolyhedron, "../Platonic_solids/Tetraedro/"))
		FAIL() << "Something went wrong during the creation of the platonic polyhedron mesh";

	PolyhedronMesh GeodeticPolyhedron;
	Generation::GeodeticSolidType1(PlatonicPolyhedron, GeodeticPolyhedron, segments);
	
	PolyhedronMesh DualPolyhedron;
	Generation::Dual(GeodeticPolyhedron, DualPolyhedron);
	
	int T = segments*segments;
	int ExpectedVertices = 4*T;
	int ExpectedFaces = 2*T+2;
	
	EXPECT_EQ(DualPolyhedron.NumCell0Ds, ExpectedVertices);
	EXPECT_EQ(DualPolyhedron.NumCell2Ds, ExpectedFaces);
}

//PER CAPIRE QUESTE DEVI VEDERE IL FUNZIONAMENTO DI ORDERED_FACES
TEST(TestOrderFaces, Test_unordered)
{
	PolyhedronMesh PlatonicPolyhedron;
	if (!FileManagement::ImportPolyhedronMesh(PlatonicPolyhedron, "../Platonic_solids/Ottaedro/"))
		FAIL() << "Something went wrong during the creation of the platonic polyhedron mesh";
	
	vector<int> unordered_faces = {4, 5, 3, 2};
	vector<int> ordered_faces;
	
	OrderFaces(unordered_faces, ordered_faces, PlatonicPolyhedron);
	vector<int>expected_ordered_faces  = {4, 5, 2, 3};
	
	EXPECT_EQ(ordered_faces, expected_ordered_faces);
}

//PER CAPIRE QUESTE DEVI VEDERE IL FUNZIONAMENTO DI ORDERED_FACES
TEST(TestOrderFaces, Test_ordered)
{
	PolyhedronMesh PlatonicPolyhedron;
	if (!FileManagement::ImportPolyhedronMesh(PlatonicPolyhedron, "../Platonic_solids/Icosaedro/"))
		FAIL() << "Something went wrong during the creation of the platonic polyhedron mesh";
	
	vector<int> unordered_faces = {10, 11, 12, 17, 16};
	vector<int> ordered_faces;
	
	OrderFaces(unordered_faces, ordered_faces, PlatonicPolyhedron);
	vector<int>expected_ordered_faces  = {10, 11, 12, 17, 16};
	
	EXPECT_EQ(ordered_faces, expected_ordered_faces);
}


/*
TEST(TestShortestPath, ShortestPathOnType1)
{
	
	PolyhedronMesh PlatonicPolyhedron;
	if (!FileManagement::ImportPolyhedronMesh(PlatonicPolyhedron, "../Platonic_solids/Ottaedro/"))
		FAIL() << "Something went wrong during the creation of the platonic polyhedron mesh";
	
	PolyhedronMesh GeodeticPolyhedron;
	Generation::GeodeticSolidType1(PlatonicPolyhedron, GeodeticPolyhedron, 2);
	
	double path_length;
	int number_edges_in_path;
	vector<int> path_vertices;	
	
	if(!Generation::ShortestPath(GeodeticPolyhedron, 3, 7, path_length, number_edges_in_path, path_vertices))
		FAIL() << "Something went wrong during the execution of ShortestPath function";
	
	vector<int>expected_path = {7,3};
	
	EXPECT_EQ(path_vertices, expected_path);
	EXPECT_EQ(number_edges_in_path, 1);
	EXPECT_NEAR(0.765367, path_length, 1e-6);
}

TEST(TestShortestPath, ShortestPathOnDual)
{
	
	PolyhedronMesh PlatonicPolyhedron;
	PolyhedronMesh DualPolyhedron;
	
	if (!FileManagement::ImportPolyhedronMesh(PlatonicPolyhedron, "../Platonic_solids/Tetraedro/"))
		FAIL() << "Something went wrong during the creation of the platonic polyhedron mesh";
	
	PolyhedronMesh GeodeticPolyhedron;
	
	Generation::GeodeticSolidType1(PlatonicPolyhedron, GeodeticPolyhedron, 3);
	Generation::Dual(GeodeticPolyhedron, DualPolyhedron);
	
	double path_length;
	int number_edges_in_path;
	vector<int> path_vertices;	
	
	if(!Generation::ShortestPath(DualPolyhedron, 21, 27, path_length, number_edges_in_path, path_vertices))
		FAIL() << "Something went wrong during the execution of ShortestPath function";
	
	vector<int>expected_path = {27, 25, 26, 21};
	
	EXPECT_EQ(path_vertices, expected_path);
	EXPECT_EQ(number_edges_in_path, 3);
	EXPECT_NEAR(1.316185, path_length,1e-6);
}
*/