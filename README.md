# CVT2D
A Lloyd algorithm implementation constructing Centroidal Voronoi Diagram in 2D. It uses the Voronoi\_diagram\_2 and Nef\_Polyhedron\_2 package in CGAL and works on arbitrary boundary contraints (i.e. whether it is convex or concave, connected or multi-connected, genus-1 or multi-genus)

## Example

<div align="center">
<img src="example/circle_400.000000_8.png" width="200" align="center"/>
<img src="example/circle_400.000000_100.PNG" width="200" align="center"/>
<img src="example/ring_400.000000_200.000000_8.png" width="200" align="center"/>
<img src="example/ring_400.000000_200.000000_100.png" width="200" align="center"/>
<br>
<caption align="bottom">Figure 1. Results.</caption>
</div>

<div align="center">
<img src="example/shrimp_10.PNG" width="200" align="center"/>
<img src="example/shrimp_20.PNG" width="200" align="center"/>
<img src="example/shrimp_50.PNG" width="200" align="center"/>
<img src="example/shrimp_100.PNG" width="200" align="center"/>
<br>
<caption align="bottom">Figure 2. Results on a shrimp-shaped constraint with 10, 20, 50 and 100 generators.</caption>
</div>
