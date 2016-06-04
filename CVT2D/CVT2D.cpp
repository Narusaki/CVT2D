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

void CCVT2D::SetMinEnergyChange(double minEnergyChange_)
{
	minEnergyChange = minEnergyChange_;
}


void CCVT2D::Execute()
{
	K::FT prevEnergy = 0.0;
	for (int i = 0; i < maxIter; ++i)
	{
		if (!isSilent)
			cout << "Iteration " << i << endl;
		// construct delaunay triangulation
		vd.clear();
		for (auto generator : generators)
			vd.insert(Site_2(generator.x(), generator.y()));
		
		if (!isSilent)
		{
			string fileName = directory + "\\iter_";
			fileName += to_string(i) + ".txt";
			ofstream output(fileName);
			PrintGenerators(output);
			output.close();
		}

		generators.clear();

		double maxMoveDist = 0.0;
		K::FT energy = 0.0;
		int cnt = 0;
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

			Point_2 curCentroid(0.0, 0.0);
			K::FT totalArea = 0.0;

			for (; faceIter != E.faces_end(); ++faceIter)
			{
				Explorer::Halfedge_around_vertex_const_circulator heIter(E.halfedge(faceIter));
				auto heStartIter = heIter;
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

				auto curEnergy = CalcCellEnergy(*siteIter, vorBoundary);
				if (E.mark(faceIter)) energy += curEnergy; else energy -= curEnergy;
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
		if (!isSilent)
		{
			cout << "Max move dist: " << maxMoveDist << ", generator left: " << generators.size()
				<< ", energy: " << CGAL::to_double(energy) << endl;
		}
		if (fabs(CGAL::to_double((energy - prevEnergy) / max(energy, prevEnergy))) < minEnergyChange)
			break;
		prevEnergy = energy;
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

K::FT CCVT2D::CalcCellEnergy(const Point_2 &center, const Polygon_2 &poly)
{
	Point_2 p0, p1, q0, q1;
	K::FT integral = 0.0;

	// partition the cell into triangles
	CDT cdt;
	cdt.insert_constraint(poly.vertices_begin(), poly.vertices_end(), true);
	mark_domains(cdt);

	// rotate each triangle to align the longest edge with x-axis
	for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
	{
		if (!fit->info().in_domain()) continue;
	
		Point_2 ps[3];
		for (int i = 0; i < 3; ++i)
			ps[i] = fit->vertex(i)->point();

		K::Vector_2 dir = ps[1] - ps[0];
		
		K::FT cosA = CGAL::to_double(dir.x()) / sqrt(CGAL::to_double(dir.squared_length()));
		if (cosA > 1.0) cosA = 1.0; else if (cosA < -1.0) cosA = -1.0;
		K::FT sinA = sqrt(1.0 - CGAL::to_double(cosA*cosA));
		if (dir.y() < 0.0) sinA *= -1.0;

		Point_2 ps_[3];
		for (int i = 0; i < 3; ++i)
			ps_[i] = Point_2(cosA * ps[i].x() + sinA * ps[i].y(), -sinA * ps[i].x() + cosA * ps[i].y());
		Point_2 center_(cosA * center.x() + sinA * center.y(), -sinA * center.x() + cosA * center.y());

		// integrate
		K::FT curIntegral = CalcEquation(center_, ps_[2], ps_[0], ps_[2], ps_[1], ps_[0].y(), ps_[2].y());
		integral += curIntegral > 0.0 ? curIntegral : -curIntegral;
	}

	return integral;
}

K::FT CCVT2D::CalcEquation(const Point_2 &center,
	const Point_2 &p0, const Point_2 &p1,
	const Point_2 &q0, const Point_2 &q1,
	const K::FT &y0, const K::FT &y1)
{
	K::FT A = p1.y() - p0.y(), B = p0.x() - p1.x(), C = p1.x()*p0.y() - p0.x()*p1.y();
	K::FT D = q1.y() - q0.y(), E = q0.x() - q1.x(), F = q1.x()*q0.y() - q0.x()*q1.y();

	K::FT a0 = -E / D, b0 = -F / D, a1 = -B / A, b1 = -C / A;
	return CalcSubEquation(center, a0, b0, a1, b1, y0, y1);
}

double CCVT2D::CalcEquation2(const Point_2 &center,
	const Point_2 &p0, const Point_2 &p1,
	const Point_2 &q0, const Point_2 &q1,
	const K::FT &y0, const K::FT &y1)
{
	K::FT A = p1.y() - p0.y(), B = p0.x() - p1.x(), C = p1.x()*p0.y() - p0.x()*p1.y();
	K::FT D = q1.y() - q0.y(), E = q0.x() - q1.x(), F = q1.x()*q0.y() - q0.x()*q1.y();

	K::FT a0 = -E / D, b0 = -F / D, a1 = -B / A, b1 = -C / A;

	Symbolic x("x"), y("y");
	Symbolic f = integrate(
		integrate(
		(x - CGAL::to_double(center.x()))*(x - CGAL::to_double(center.x())) + (y - CGAL::to_double(center.y()))*(y - CGAL::to_double(center.y())), 
		x, 
		CGAL::to_double(a0)*y + CGAL::to_double(b0), CGAL::to_double(a1)*y + CGAL::to_double(b1)), 
		y);     // => 1/2*x^(2)+x
	/*Symbolic f = integrate(integrate(1, x, -sqrt(1-y*y), sqrt(1-y*y)), y);     // => 1/2*x^(2)+x*/
	return fabs(f[y == CGAL::to_double(y1)] - f[y == CGAL::to_double(y0)]);
}

K::FT CCVT2D::CalcSubEquation(const Point_2 &center,
	const K::FT &a0, const K::FT &b0, const K::FT &a1, const K::FT &b1,
	const K::FT &y0, const K::FT &y1)
{
	return CalcSubEquation(center, a1, b1, y0, y1) - CalcSubEquation(center, a0, b0, y0, y1);
}

K::FT CCVT2D::CalcSubEquation(const Point_2 &center, const K::FT &a, const K::FT &b, const K::FT &y0, const K::FT &y1)
{
	return CalcSubEquation(center, a, b, y1) - CalcSubEquation(center, a, b, y0);
}

K::FT CCVT2D::CalcSubEquation(const Point_2 &center, const K::FT &a, const K::FT &b, const K::FT &y)
{
	return (a*a*a / 3.0 + a)*y*y*y*y / 4.0 +
		(a*a*(b - center.x()) + (b - 2.0*center.y()*a))*y*y*y / 3.0 +
		(a*(b - center.x())*(b - center.x()) + a*center.y()*center.y() - 2.0*center.y()*b)*y*y / 2.0 +
		((b - center.x())*(b - center.x())*(b - center.x())/3.0 + b*center.y()*center.y()) * y;
}

void CCVT2D::mark_domains(CDT& ct, CDT::Face_handle start, int index, std::list<CDT::Edge>& border)
{
	if (start->info().nesting_level != -1){
		return;
	}
	std::list<CDT::Face_handle> queue;
	queue.push_back(start);
	while (!queue.empty()){
		CDT::Face_handle fh = queue.front();
		queue.pop_front();
		if (fh->info().nesting_level == -1){
			fh->info().nesting_level = index;
			for (int i = 0; i < 3; i++){
				CDT::Edge e(fh, i);
				CDT::Face_handle n = fh->neighbor(i);
				if (n->info().nesting_level == -1){
					if (ct.is_constrained(e)) border.push_back(e);
					else queue.push_back(n);
				}
			}
		}
	}
}

void CCVT2D::mark_domains(CDT& cdt)
{
	for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
		it->info().nesting_level = -1;
	}
	std::list<CDT::Edge> border;
	mark_domains(cdt, cdt.infinite_face(), 0, border);
	while (!border.empty()){
		CDT::Edge e = border.front();
		border.pop_front();
		CDT::Face_handle n = e.first->neighbor(e.second);
		if (n->info().nesting_level == -1){
			mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
		}
	}
}