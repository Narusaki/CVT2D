#include "..\CVT2D\CVT2D.h"

#include <iostream>
#include <fstream>
#include <random>

using namespace std;

template< typename T >
Line ConstructHalfPlane(T p0, T p1)
{
	auto A = p0.y() - p1.y();
	auto B = p1.x() - p0.x();
	auto C = p0.x()*p1.y() - p1.x()*p0.y();
	return Line(CGAL::to_double(A), CGAL::to_double(B), CGAL::to_double(C));
}

bool LoadBoundary(const char *fileName, Nef_polyhedron &boundaryNef)
{
	int processRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
	ifstream input(fileName);
	if (!input)
	{
		if (processRank == 0)
			cout << "Cannot open boundary file!" << endl;
		return false;
	}

	// loading boundary information
	TPolygon_2 outLine;
	list< TPolygon_2 > partitionedOuter, partitionedInner;
	list< TPolygon_2 > boundaryOuter, boundaryInner;
	string curLine;
	bool isOuter = true;
	while (getline(input, curLine))
	{
		if (curLine[0] == '#')
		{
			if (!outLine.is_empty())
				if (isOuter) boundaryOuter.push_back(outLine);
				else boundaryInner.push_back(outLine);
				outLine.clear();
				isOuter = curLine == "#Outer";
		}
		else
		{
			double x, y;
			stringstream sin;
			sin << curLine; sin >> x >> y;
			outLine.push_back(TPoint_2(x, y));
		}
	}
	if (!outLine.is_empty())
	{
		if (isOuter) boundaryOuter.push_back(outLine);
		else boundaryInner.push_back(outLine);
	}

	// partitioning...
	if (processRank == 0)
		cout << "Partitioning..." << endl;
	for (auto &Outer : boundaryOuter)
	{
		Traits partition_traits;
		Validity_traits validity_traits;
		list< TPolygon_2 > curPartition;
		CGAL::optimal_convex_partition_2(Outer.vertices_begin(), Outer.vertices_end(), back_inserter(curPartition), partition_traits);
		assert(CGAL::partition_is_valid_2(Outer.vertices_begin(), Outer.vertices_end(), curPartition.begin(), curPartition.end(), validity_traits));
		partitionedOuter.insert(partitionedOuter.end(), curPartition.begin(), curPartition.end());
	}
	for (auto &Inner : boundaryInner)
	{
		Traits partition_traits;
		Validity_traits validity_traits;
		list< TPolygon_2 > curPartition;
		CGAL::optimal_convex_partition_2(Inner.vertices_begin(), Inner.vertices_end(), back_inserter(curPartition), partition_traits);
		assert(CGAL::partition_is_valid_2(Inner.vertices_begin(), Inner.vertices_end(), curPartition.begin(), curPartition.end(), validity_traits));
		partitionedInner.insert(partitionedInner.end(), curPartition.begin(), curPartition.end());
	}
	if (processRank == 0)
	{
		cout << "Partition result: " << endl;
		cout << "Outer: " << partitionedOuter.size() << ", Inner: " << partitionedInner.size() << endl;
	}

	boundaryNef = Nef_polyhedron(Nef_polyhedron::EMPTY);
	for (auto &Outer : partitionedOuter)
	{
		auto curPolygon = Nef_polyhedron(Nef_polyhedron::COMPLETE);
		for (auto v = Outer.vertices_begin(); v != Outer.vertices_end(); ++v)
		{
			auto vNext = v; ++vNext;
			if (vNext == Outer.vertices_end()) vNext = Outer.vertices_begin();
			curPolygon *= Nef_polyhedron(ConstructHalfPlane(*v, *vNext));
		}
		boundaryNef += curPolygon;
	}
	for (auto &Inner : partitionedInner)
	{
		auto curPolygon = Nef_polyhedron(Nef_polyhedron::COMPLETE);
		for (auto v = Inner.vertices_begin(); v != Inner.vertices_end(); ++v)
		{
			auto vNext = v; ++vNext;
			if (vNext == Inner.vertices_end()) vNext = Inner.vertices_begin();
			curPolygon *= Nef_polyhedron(ConstructHalfPlane(*v, *vNext));
		}
		boundaryNef -= curPolygon;
	}

	return true;
}

void GenInitGenerators(const Nef_polyhedron& boundaryNef, vector< Point_2 > &generators, int Ng)
{
	generators.clear();

	// construct bounding box of boundary
	typedef Nef_polyhedron::Explorer Explorer;
	Explorer E = boundaryNef.explorer();
	auto faceIter = E.faces_begin(); ++faceIter; ++faceIter;
	double xmin = 1e30, xmax = -1e30, ymin = 1e30, ymax = -1e30;

	for (; faceIter != E.faces_end(); ++faceIter)
	{
		Explorer::Halfedge_around_vertex_const_circulator heIter(E.halfedge(faceIter));
		auto heStartIter = heIter;
		do
		{
			if (E.is_standard(E.target(heIter)) && E.is_standard(E.source(heIter)))
			{
				auto p = E.point(E.source(heIter));
				xmin = min(CGAL::to_double(xmin), CGAL::to_double(p.x()));
				xmax = max(CGAL::to_double(xmax), CGAL::to_double(p.x()));
				ymin = min(CGAL::to_double(ymin), CGAL::to_double(p.y()));
				ymax = max(CGAL::to_double(ymax), CGAL::to_double(p.y()));
			}
			heIter = heIter->next();
		} while (heIter != heStartIter);
	}

	default_random_engine gen;
	uniform_real_distribution<double> distX(xmin, xmax);
	uniform_real_distribution<double> distY(ymin, ymax);
	for (int i = 0; i < Ng; ++i)
	{
		auto p = Nef_polyhedron::Point(distX(gen), distY(gen));
		if (boundaryNef.contains(boundaryNef.locate(p)))
			generators.push_back(Point_2(p));
		else --i;
	}
}

void LoadInitGenerator(const char *fileName, vector< Point_2 > &generators)
{
	ifstream input(fileName);
	if (!input)
	{
		cout << "Cannot open initial generators file!" << endl;
		return;
	}
	generators.clear();
	double x, y;
	while (input >> x >> y)
		generators.push_back(Point_2(x, y));
	input.close();
}

int main(int argc, char **argv)
{
	if (argc < 3)
	{
		cout << "USAGE: [.exe] [.boundary] [#generator] [initGenerators]" << endl;
		return -1;
	}

	MPI_Init(&argc, &argv);

	string directory;
	if (argc == 3)
		directory = ".";
	else
	{
		directory = argv[3];
		if (directory.find("\\") == string::npos) directory = ".";
		else directory = directory.substr(0, directory.rfind("\\"));
	}
		
	int processRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
	Nef_polyhedron boundaryNef;
	if (!LoadBoundary(argv[1], boundaryNef))
	{
		if (processRank == 0)
			cout << "Cannot load boundary file!" << endl;
		return -2;
	}
	int Ng = stoi(argv[2]);
	vector< Point_2 > generators;
	if (argc == 3)
		GenInitGenerators(boundaryNef, generators, Ng);
	else
		LoadInitGenerator(argv[3], generators);

	CCVT2D cvt;
	cvt.isSilent = false;
	cvt.AssignBoundary(boundaryNef);
	cvt.AssignGeneratorNum(Ng);
	cvt.AssignInitGenerators(generators);

	cvt.SetMaxIteration(100);
	cvt.SetMinEnergyChange(0.00000);
	cvt.SetOutputDirectory(directory);

	cvt.Execute();

	if (processRank == 0)
	{
		ofstream output(directory + "\\finalState.txt");
		cvt.PrintGenerators(output);
		output.close();
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
	return 0;
}