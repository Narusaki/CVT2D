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
#include <vector>
#include <list>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Iso_rectangle_2.h>

#include <CGAL/Boolean_set_operations_2.h>

#include <CGAL/Exact_rational.h>
#include <CGAL/Extended_cartesian.h>
#include <CGAL/Nef_polyhedron_2.h>

// includes for convex partition
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Partition_is_valid_traits_2.h>
#include <CGAL/polygon_function_objects.h>
#include <CGAL/partition_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <cassert>

// typedefs for Voronoi diagram
typedef CGAL::Exact_predicates_exact_constructions_kernel			K;
typedef CGAL::Delaunay_triangulation_2<K>							DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT, AT, AP>                                  VD;

// typedefs for generic geometry
typedef K::Point_2													Point_2;
typedef K::Segment_2												Segment_2;
typedef K::Ray_2													Ray_2;
typedef K::Line_2													Line_2;
typedef VD::Site_2													Site_2;
typedef VD::Face_handle												Face_handle;
typedef VD::Halfedge_handle											Halfedge;
typedef VD::Vertex_handle											Vertex_handle;
typedef CGAL::Polygon_2<K>											Polygon_2;
typedef CGAL::Polygon_with_holes_2<K>								Polygon_with_holes_2;

// typedefs for boolean operations
typedef CGAL::Exact_rational RT;
typedef CGAL::Extended_cartesian<RT> Extended_kernel;
typedef CGAL::Nef_polyhedron_2<Extended_kernel> Nef_polyhedron;
typedef Nef_polyhedron::Line  Line;

// typedefs for convex partition
typedef CGAL::Partition_traits_2<K>                         Traits;
typedef CGAL::Is_convex_2<Traits>                           Is_convex_2;
typedef Traits::Polygon_2                                   TPolygon_2;
typedef Traits::Point_2                                     TPoint_2;
typedef CGAL::Partition_is_valid_traits_2 < Traits, Is_convex_2 > Validity_traits;
typedef CGAL::Creator_uniform_2<int, TPoint_2>               Creator;

// This class is exported from the CVT2D.dll
class CVT2D_API CCVT2D {
public:
	CCVT2D();
	~CCVT2D();

	void AssignBoundary(const Nef_polyhedron &boundary_);

	void AssignGeneratorNum(int genNum_);
	void AssignInitGenerators(const std::vector< Point_2 > &generators_);

	void SetMaxIteration(int maxIter_);
	void SetMinMove(double minMove_);

	void Execute();

	void PrintGenerators(std::ostream &output);

private:
	template< typename T >
	Line ConstructHalfPlane(T p0, T p1)
	{
		auto A = p0.y() - p1.y();
		auto B = p1.x() - p0.x();
		auto C = p0.x()*p1.y() - p1.x()*p0.y();
		return Line(CGAL::to_double(A), CGAL::to_double(B), CGAL::to_double(C));
	}

	Point_2 CalcCentroidOfPolygon(const Polygon_2& polygon);

private:
	VD vd;
	std::vector< Point_2 > generators;
	Nef_polyhedron boundaryNef;

	int maxIter;
};
