#include "delaunay_mesh.h"

template <typename T>
double DelaunayMesh::DistanceSquared(const T &p0, const T &p1) {
  double v[3] = {p0.x() - p1.x(), p0.y() - p1.y(), p0.z() - p1.z()};
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

// Specializing template functions for point with different CGAL kernels.
template double DelaunayMesh::DistanceSquared(const DelaunayMesh::Delaunay::Point &p0, const DelaunayMesh::Delaunay::Point &p1);
template double DelaunayMesh::DistanceSquared(const DelaunayMesh::TrianglePrimitive::Point &p0, const DelaunayMesh::TrianglePrimitive::Point &p1);


void DelaunayMesh::AddPoint(const DelaunayMesh::DPoint &p, const DelaunayMesh::VertexInfo &vi) {
  auto v = d.nearest_vertex(p);
  // Found a nearest vertex within merge_threshold distance.
  if (d.is_vertex(v) && DistanceSquared(p, v->point()) < merge_threshold) {
    auto vp = v->point();
    int n = v->info().n;
    d.move(v, Delaunay::Point(
          (p.x() + vp.x() * n) / (n + 1),
          (p.y() + vp.y() * n) / (n + 1),
          (p.z() + vp.z() * n) / (n + 1)));
    v->info().n++;
    // Merges camera list;
    v->info().cams.insert(vi.cams.begin(), vi.cams.end());
  } else {
    auto v = d.insert(p);
    v->info() = vi;
  }

}

DelaunayMesh::FacetAndNormal::FacetAndNormal(Facet f) {
  this->f = f;
  auto &ch = f.first;
  auto &v0 = ch->vertex((f.second + 1) % 4)->point();
  auto &v1 = ch->vertex((f.second + 2) % 4)->point();
  auto &v2 = ch->vertex((f.second + 3) % 4)->point();
  auto v10 = v1 - v0;
  auto v21 = v2 - v1;
  n = CGAL::cross_product(v10, v21);
  n = n / sqrt(n.squared_length());
  if (f.second % 2 == 0)
    n = -n;
}

// Utility function to convert a custom point type to CGAL point type.
DelaunayMesh::Point DelaunayMesh::TrianglePrimitive::convert(const DelaunayMesh::DPoint &p) const {
  return DelaunayMesh::Point(p[0], p[1], p[2]);
}
    
// On-the-fly conversion from the internal data to the CGAL types.
DelaunayMesh::Datum DelaunayMesh::TrianglePrimitive::datum() const {
  const auto &v0 = m_pt->f.first->vertex((m_pt->f.second + 1) % 4);
  const auto &v1 = m_pt->f.first->vertex((m_pt->f.second + 2) % 4);
  const auto &v2 = m_pt->f.first->vertex((m_pt->f.second + 3) % 4);
  return Datum(convert(v0->point()), convert(v1->point()), convert(v2->point()));
}

// Returns a reference point which must be on the primitive.
DelaunayMesh::Point DelaunayMesh::TrianglePrimitive::reference_point() const { 
  const auto &p = m_pt->f.first->vertex((m_pt->f.second + 1) % 4)->point();
  return DelaunayMesh::Point(p[0], p[1], p[2]);
}

// Builds AABB tree for fast ray intersection. Return the number of triangles.
int DelaunayMesh::BuildAABBTree() {
  triangles.clear();
  tree.clear();
  for (auto it = d.finite_facets_begin(); it != d.finite_facets_end(); it++) {
    triangles.push_back(FacetAndNormal(*it));
  }
  tree.insert(triangles.begin(),triangles.end());
  return triangles.size();
}

void DelaunayMesh::RayIntersect(const DelaunayMesh::KSC::Ray_3 &ray, std::list<Tree::Intersection_and_primitive_id<KSC::Ray_3>::Type> &intersections) {
  tree.all_intersections(ray, std::back_inserter(intersections));

}

// Assigns IDs to all tetrahedrons after Delaunay is complete. Returns number of cells.
int DelaunayMesh::AssignTetrahedronIds() {
  int idCount = 0;
  for (auto it = d.all_cells_begin(); it != d.all_cells_end(); it++) {
    it->info().id = idCount++;
  }
  return idCount;
}
