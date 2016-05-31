// CVT2D.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "CVT2D.h"

using namespace std;

// This is the constructor of a class that has been exported.
// see CVT2D.h for the class definition
CCVT2D::CCVT2D()
{
	return;
}

CCVT2D::~CCVT2D()
{

}

void CCVT2D::AssignBoundary(const Nef_polyhedron &boundary_)
{
	boundaryNef = boundary_;
}

void CCVT2D::AssignGeneratorNum(int genNum_)
{
	generators.clear();
	generators.resize(genNum_);
}

void CCVT2D::AssignInitGenerators(const vector< Point_2 > &generators_)
{
	generators = generators_;
}

void CCVT2D::SetMaxIteration(int maxIter_)
{
	maxIter = maxIter_;
}

void CCVT2D::SetMinMove(double minMove_)
{

}


void CCVT2D::Execute()
{
	for (int i = 0; i < maxIter; ++i)
	{
		cout << "Iteration " << i << endl;
		// construct delaunay triangulation
		vd.clear();
		for (auto generator : generators)
			vd.insert(Site_2(generator.x(), generator.y()));
		
		string fileName = "iter_";
		fileName += to_string(i) + ".txt";
		ofstream output(fileName);
		PrintGenerators(output);
		output.close();

		generators.clear();

		double maxMoveDist = 0.0;
		for (VD::Site_iterator siteIter = vd.sites_begin(); siteIter != vd.sites_end(); ++siteIter)
		{
			// construct each Voronoi cell
			VD::Locate_result locRes = vd.locate(*siteIter);
			Face_handle *cell = get<Face_handle>(&locRes);
			VD::Ccb_halfedge_circulator halfEdgeCirculator = (*cell)->ccb();
			/*while (!halfEdgeCirculator->has_source() && halfEdgeCirculator->has_target()) ++halfEdgeCirculator;*/
			auto startHalftEdgeCirculator = halfEdgeCirculator;

			Nef_polyhedron vorCellNef(Nef_polyhedron::COMPLETE);
			do {
				if (halfEdgeCirculator->has_source() && halfEdgeCirculator->has_target())
				{
					vorCellNef *= Nef_polyhedron(ConstructHalfPlane(halfEdgeCirculator->source()->point(), halfEdgeCirculator->target()->point()), Nef_polyhedron::EXCLUDED);
				}
				else if (halfEdgeCirculator->has_source() && !halfEdgeCirculator->has_target())
				{
					// calculate direction
					auto delaunayEdge = halfEdgeCirculator->dual();
					auto v0 = delaunayEdge.first->vertex(delaunayEdge.first->ccw(delaunayEdge.second))->point();
					auto v1 = delaunayEdge.first->vertex(delaunayEdge.first->cw(delaunayEdge.second))->point();
					if (v0 != *siteIter) swap(v0, v1);
					K::Vector_2 direction(-(v1 - v0).y(), (v1 - v0).x());
					Point_2 endPoint = halfEdgeCirculator->source()->point() + direction;
					vorCellNef *= Nef_polyhedron(ConstructHalfPlane(Point_2(halfEdgeCirculator->source()->point()), endPoint), Nef_polyhedron::EXCLUDED);
					// calculate an endpoint that is far enough
				}
				else if (!halfEdgeCirculator->has_source() && halfEdgeCirculator->has_target())
				{
					// calculate direction
					auto delaunayEdge = halfEdgeCirculator->dual();
					auto v0 = delaunayEdge.first->vertex(delaunayEdge.first->ccw(delaunayEdge.second))->point();
					auto v1 = delaunayEdge.first->vertex(delaunayEdge.first->cw(delaunayEdge.second))->point();
					if (v0 != *siteIter) swap(v0, v1);
					K::Vector_2 direction((v1 - v0).y(), -(v1 - v0).x());
					Point_2 startPoint = halfEdgeCirculator->target()->point() + direction;
					vorCellNef *= Nef_polyhedron(ConstructHalfPlane(startPoint, Point_2(halfEdgeCirculator->target()->point())), Nef_polyhedron::EXCLUDED);
					// calculate an endpoint that is far enough
				}
				else
				{
					auto delaunayEdge = halfEdgeCirculator->dual();
					auto v0 = delaunayEdge.first->vertex(delaunayEdge.first->ccw(delaunayEdge.second))->point();
					auto v1 = delaunayEdge.first->vertex(delaunayEdge.first->cw(delaunayEdge.second))->point();
					if (v0 != *siteIter) swap(v0, v1);
					K::Vector_2 direction(-(v1 - v0).y(), (v1 - v0).x());
					Point_2 startPoint((v0.x() + v1.x()) / 2.0, (v0.y() + v1.y()) / 2.0);
					Point_2 endPoint = startPoint + direction;
					vorCellNef *= Nef_polyhedron(ConstructHalfPlane(startPoint, endPoint));
				}
				++halfEdgeCirculator;
			} while (halfEdgeCirculator != startHalftEdgeCirculator);

			vorCellNef *= boundaryNef;

			// now the Voronoi cell is in the format of Nef_polygon
			// to calculate its centroid we can simply record the outline-boundaries and holes, 
			// and do a weighted addition&substraction
			typedef Nef_polyhedron::Explorer Explorer;
			Explorer E = vorCellNef.explorer();
			auto faceIter = E.faces_begin(); ++faceIter; ++faceIter;
			Explorer::Halfedge_around_vertex_const_circulator heIter(E.halfedge(faceIter));
			auto heStartIter = heIter;

			Point_2 curCentroid(0.0, 0.0);
			K::FT totalArea = 0.0;

			for (; faceIter != E.faces_end(); ++faceIter)
			{
				Polygon_2 vorBoundary;
				do
				{
					if (E.is_standard(E.target(heIter)) && E.is_standard(E.source(heIter)))
					{
						auto p = E.point(E.source(heIter));
						vorBoundary.push_back(Point_2(CGAL::to_double(p.x()), CGAL::to_double(p.y())));
					}
					heIter = heIter->next();
				} while (heIter != heStartIter);

				Point_2 c = CalcCentroidOfPolygon(vorBoundary);
				K::FT curArea = (E.mark(faceIter) ? vorBoundary.area() : -vorBoundary.area());
				c = Point_2(curArea * c.x(), curArea * c.y());
				
				curCentroid = Point_2(curCentroid.x() + c.x(), curCentroid.y() + c.y());
				totalArea += curArea;
			}

			// move the generator to centroid
			if (totalArea == 0.0) continue;
			generators.push_back(Point_2(curCentroid.x() / totalArea, curCentroid.y() / totalArea));
// 			cout << totalArea << endl;
// 			cout << "site: " << siteIter->x() << " " << siteIter->y();
// 			cout << " -- >";
// 			cout << generators.back().x() << " " << generators.back().y() << endl;

			double curMoveDist = sqrt(CGAL::to_double((*siteIter - generators.back()).squared_length()));
			maxMoveDist = max(maxMoveDist, curMoveDist);
		}
		cout << "Max move dist: " << maxMoveDist << ", generator left: " << generators.size() << endl;
	}
}

void CCVT2D::PrintGenerators(ostream &output)
{
	for (auto p : generators)
	{
		output << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << endl;
	}
}

Point_2 CCVT2D::CalcCentroidOfPolygon(const Polygon_2& polygon)
{
	Polygon_2::FT x = 0.0, y = 0.0;
	for (auto vert = polygon.vertices_begin(); vert != polygon.vertices_end(); ++vert)
	{
		auto vertNext = vert; ++vertNext; 
		if (vertNext == polygon.vertices_end()) vertNext = polygon.vertices_begin();
		Polygon_2::FT x0 = vert->x(), y0 = vert->y();
		Polygon_2::FT x1 = vertNext->x(), y1 = vertNext->y();
		x += (x0 + x1)*(x0*y1 - x1*y0);
		y += (y0 + y1)*(x0*y1 - x1*y0);
	}
	Polygon_2::FT A = polygon.area();
	return Point_2(x / 6.0 / A, y / 6.0 / A);
}