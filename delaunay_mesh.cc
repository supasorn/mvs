#include "delaunay_mesh.h"


void DelaunayMesh::AddPoint(const DelaunayMesh::DPoint &p, const DelaunayMesh::VertexInfo &vi) {
  auto v = d.nearest_vertex(p);
  // Found a nearest vertex within merge_threshold distance.
  if (d.is_vertex(v) && (p - v->point()).squared_length() < merge_threshold) {
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

inline int DelaunayMesh::IsSameSide(DelaunayMesh::KSC::Vector_3 &a, const DelaunayMesh::K::Vector_3 &b) {
    return a.x() * b.x() + a.y() * b.y() + a.z() * b.z() > 0;
}

int DelaunayMesh::EqualZero(double &a) {
  return std::abs(a) < 1e-10;
}

void DelaunayMesh::AssignCost(Point &point, Point &camera, double cost) {
  typedef std::pair<const DelaunayMesh::FacetAndNormal*, double> sortedFacets; 
  struct Comparator {
    static bool compare(sortedFacets &a, sortedFacets &b) { return a.second < b.second; }
  };

  double cameraToPointDist = (camera - point).squared_length();
  KSC::Ray_3 ray_query(camera, point);
  auto toCamera = camera - point;

  std::list<IntersectReturn> inters;
  RayIntersect(ray_query, inters);

  // Store intersected facets and its distance from the depth point.
  std::vector<sortedFacets> dists;
  for (auto &inter: inters) {

    // Only retrieve point intersection here because segment intersection
    // always has triangle's normal perpedicular to the ray.
    if (const Point *p = boost::get<Point>(&(inter.first))) {
      const DelaunayMesh::Facet &f = inter.second->f;
      double dist = (*p - camera).squared_length() - cameraToPointDist;
      dists.push_back(std::make_pair(inter.second, dist));
    } 
  }
  sort(dists.begin(), dists.end(), Comparator::compare);

  for (int i = 0; i < dists.size(); i++) {
    if (EqualZero(dists[i].second)) continue;

    int id[2];
    getIncidentTetrahedrons(dists[i].first, id[0], id[1]);
    // True iff normal and ray from point to camera are on the same side.
    int ss = IsSameSide(toCamera, dists[i].first->n);

    if (i == 0) {
      g->add_tweights(nodes[id[!ss]], cost, 0);
    }

    if (dists[i].second < 0) {
      g->add_edge(nodes[id[!ss]], nodes[id[ss]], cost, 0);
    } else { 
      g->add_tweights(nodes[id[!ss]], 0, cost);
      break;
    }
  }
}

void DelaunayMesh::getIncidentTetrahedrons(const DelaunayMesh::FacetAndNormal *f, int &id0, int &id1) {
  id0 = f->f.first->info().id; 
  id1 = d.mirror_facet(f->f).first->info().id;
}

void DelaunayMesh::SaveObj(std::string name) {
  std::vector<DPoint> _points;
  std::vector<Eigen::Vector3i> _faces;
  for (auto &triangle : triangles) {
    int id[2];
    getIncidentTetrahedrons(&triangle, id[0], id[1]);
    if (IsInsideSurface(id[0]) != IsInsideSurface(id[1])) {
      auto &f = triangle.f;
      int base = _points.size() + 1;
      for (int i = 1; i <= 3; i++) {
        auto v = f.first->vertex((f.second + i) % 4)->point();
        _points.push_back(v);
      }
      if ((f.second % 2) ^ IsInsideSurface(id[0]))
        _faces.push_back(Eigen::Vector3i(base, base+1, base+2));
      else
        _faces.push_back(Eigen::Vector3i(base, base+2, base+1));

    }
  }
  FILE *fo = fopen(name.c_str(), "w");
  for (int i = 0; i < _points.size(); i++) {
    fprintf(fo, "v %lf %lf %lf\n", _points[i][0], _points[i][1], _points[i][2]);
  }
  for (int i = 0; i < _faces.size(); i++) { 
    fprintf(fo, "f %d %d %d\n", _faces[i][0], _faces[i][1], _faces[i][2]);
  }
  fclose(fo);
}

double CosineOfCircumsphere(DelaunayMesh::Facet f, DelaunayMesh::K::Vector_3 &normal) {
  auto &ch = f.first;
  auto center = CGAL::circumcenter(ch->vertex(0)->point(), ch->vertex(1)->point(), ch->vertex(2)->point(), ch->vertex(3)->point());
  auto &p0 = ch->vertex((f.second + 1) % 4)->point();
  auto radius = center - p0;
  // This is dot product.
  double height = normal * radius;
  double length = sqrt(radius.squared_length()); 
  return height / length;
}

int DelaunayMesh::IsInsideSurface(int id) {
  return g->what_segment(nodes[id]) == Graph::SINK;
}

double timestamp() {
  timeval start; 
  gettimeofday(&start, NULL);
  return ((start.tv_sec) + start.tv_usec/1000000.0);
}

void DelaunayMesh::ExtractSurface(std::vector<Point> &camera) {
  BuildAABBTree();
  AssignTetrahedronIds();

  nodes.resize(d.number_of_cells());

	g = new Graph();
  for (int i = 0; i < nodes.size(); i++) {
    nodes[i] = g->add_node();
  }

  for (auto it = d.finite_vertices_begin(); it != d.finite_vertices_end(); it++) {
    auto p = it->point();
    std::set<int> &cams = it->info().cams;
    //printf("%f %f %f %d\n", p.x(), p.y(), p.z(), cams.size());

    Point point = Point(p.x(), p.y(), p.z());
    double cost = cams.size();
    for (int camId : cams) {
      AssignCost(point, camera[camId], cost);
    }
  }


  for (auto &triangle : triangles) {
    // Surface quality. Higher -> smoother.
    double lambda = 20;
    double surfaceQuality = lambda * (1 - std::min(
        CosineOfCircumsphere(triangle.f, triangle.n), 
        -CosineOfCircumsphere(d.mirror_facet(triangle.f), triangle.n)));

    int id[2];
    getIncidentTetrahedrons(&triangle, id[0], id[1]);
    g->add_edge(nodes[id[0]], nodes[id[1]], surfaceQuality, surfaceQuality);

  }

  printf("Solving min-cut\n");
  double start = timestamp();
  Graph::flowtype flow = g->maxflow();
  printf("flow = %lf\n", flow);
  printf("done\n");
  printf("Time %f\n", timestamp() - start);

}
