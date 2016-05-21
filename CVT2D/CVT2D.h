// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the CVT2D_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// CVT2D_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef CVT2D_EXPORTS
#define CVT2D_API __declspec(dllexport)
#else
#define CVT2D_API __declspec(dllimport)
#endif

#include <iostream>
#include <fstream>
#include <cassert>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel			K;
typedef CGAL::Delaunay_triangulation_2<K>							DT;
typedef K::Point_2					Point_2;
typedef K::Segment_2				Segment_2;
typedef K::Ray_2					Ray_2;
typedef K::Line_2					Line_2;

// This class is exported from the CVT2D.dll
class CVT2D_API CCVT2D {
public:
	CCVT2D();
	~CCVT2D();

	void AssignBoundary();

	void AssignGeneratorNum(int genNum_);
	void AssignInitGenerators();

	void SetMaxIteration(int maxIter_);
	void SetMinMove(double minMove_);


	void Execute();
	
private:
	DT dt;
};
