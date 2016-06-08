#include "delaunay_mesh.h"

// *********    TO BE REMOVED     ***********
// ********* Start Profiling code ***********
std::map<std::string, double> ticmap;
std::map<std::string, double> ticamap; // accumulate
double timestamp() {
  timeval start; 
  gettimeofday(&start, NULL);
  return ((start.tv_sec) + start.tv_usec/1000000.0);
}

void tic(std::string name) {
  ticmap[name] = timestamp();
}

double toc(std::string name) {
  if (ticmap.find(name) != ticmap.end()) {
    double ret = timestamp() - ticmap[name];
    printf("Time %s: %f\n", name.c_str(), ret);
    return ret;
  }
  return -1;
}

void toca(std::string name) {
  if (ticmap.find(name) != ticmap.end()) {
    if (ticamap.find(name) == ticamap.end())
      ticamap[name] = timestamp() - ticmap[name];
    else
      ticamap[name] += timestamp() - ticmap[name];
  }
}

void tocaList() {
  printf("Profiling\n");
  for (std::map<std::string, double>::iterator i = ticamap.begin(); i != ticamap.end(); i++) {
    printf("  %s : %f\n", i->first.c_str(), i->second);
  }
}
// ********* End Profiling code ***********


// Adds depth point to the Delaunay.
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

// Assigns IDs to all vertices and tetrahedrons after Delaunay is complete. Returns number of cells.
void DelaunayMesh::AssignIds() {
  int idCount = 0;
  for (auto it = d.all_cells_begin(); it != d.all_cells_end(); it++) {
    it->info().id = idCount++;
  }

  idCount = 0;
  for (auto it = d.all_vertices_begin(); it != d.all_vertices_end(); it++) {
    it->info().id = idCount++;
  }
}

// Check if dot product between a and b > 0. It's a bit annoying because a and
// b use different CGAL kernel, and a - b isn't overloaded in CGAL.
inline int IsSameSide(DelaunayMesh::KSC::Vector_3 &a, const DelaunayMesh::K::Vector_3 &b) {
    return a.x() * b.x() + a.y() * b.y() + a.z() * b.z() > 0;
}

// Adds an edge for node id0 to source and sink.
void DelaunayMesh::AddCostUnary(int id0, double src, double snk) {
  omp_set_lock(&lock[id0]);
  cost_unary[id0].first += src;
  cost_unary[id0].second += snk;
  omp_unset_lock(&lock[id0]);
}

// Adds an edge from and to node id0 to id1 with cost0 and cost1 respectively.
void DelaunayMesh::AddCostBinary(int id0, int id1, double cost0, double cost1) {
  omp_set_lock(&lock[id0]);
  if (cost_binary[id0].find(id1) == cost_binary[id0].end())
    cost_binary[id0][id1] = cost0;
  else
    cost_binary[id0][id1] += cost0;
  omp_unset_lock(&lock[id0]);

  omp_set_lock(&lock[id1]);
  if (cost_binary[id1].find(id0) == cost_binary[id1].end())
    cost_binary[id1][id0] = cost1;
  else
    cost_binary[id1][id0] += cost1;
  omp_unset_lock(&lock[id1]);
}

// Adds a set of edges for the network flow for point and camera with cost.
void DelaunayMesh::AssignCost(Point &point, Point &camera, double cost) {
  typedef std::pair<const DelaunayMesh::FacetAndNormal*, double> sortedFacets; 

  double cameraToPointDist = (camera - point).squared_length();
  KSC::Ray_3 ray_query(camera, point);
  auto toCamera = camera - point;

  tic("AABB");
  std::list<IntersectReturn> inters;
  tree.all_intersections(ray_query, std::back_inserter(inters));
  toca("AABB");

  const FacetAndNormal* minFacet[2] = {0, 0};
  double minVal[2] = {1e10, 1e10};
  int id[2];
  tic("Add");
  for (auto &inter: inters) {
    // Only retrieve point intersection here because segment intersection
    // always has triangle's normal perpedicular to the ray.
    if (const Point *p = boost::get<Point>(&(inter.first))) {
      const DelaunayMesh::Facet &f = inter.second->f;
      double dist = (*p - camera).squared_length() - cameraToPointDist;
      if (std::abs(dist) < 1e-10) continue;

      if (dist > 0) {
        if (!minFacet[0] || dist < minVal[0]) {
          minFacet[0] = inter.second;
          minVal[0] = dist;
        }
      } else {
        if (!minFacet[1] || dist < minVal[1]) {
          minFacet[1] = inter.second;
          minVal[1] = dist;
        }
        GetIncidentTetrahedrons(inter.second, id[0], id[1]);
        // True iff normal and ray from point to camera are on the same side.
        int ss = IsSameSide(toCamera, inter.second->n);
        AddCostBinary(id[!ss], id[ss], cost, 0);
      }
    } 
  }
  toca("Add");

  for (int i = 0; i < 2; i++) {
    if (minFacet[i]) {
      GetIncidentTetrahedrons(minFacet[i], id[0], id[1]);
      AddCostUnary(id[!IsSameSide(toCamera, minFacet[i]->n)], i * cost, (1 - i) * cost);
    }
  }
}

// Returns id0, and id1 of tetrahedra incident to facet f.
void DelaunayMesh::GetIncidentTetrahedrons(const DelaunayMesh::FacetAndNormal *f, int &id0, int &id1) {
  id0 = f->f.first->info().id; 
  id1 = d.mirror_facet(f->f).first->info().id;
}

// Saves mesh to obj file.
void DelaunayMesh::SaveObj(std::string name) {
  std::map<int, int> used;
  std::vector<Vertex_handle> sorted;
  std::vector<Eigen::Vector3i> faces;

  for (auto &triangle : triangles) {
    int id[2];
    GetIncidentTetrahedrons(&triangle, id[0], id[1]);
    if (IsInsideSurface(id[0]) != IsInsideSurface(id[1])) {
      auto &f = triangle.f;
      int sign = (f.second % 2) ^ IsInsideSurface(id[0]) ? 1 : -1;
      faces.push_back(Eigen::Vector3i());
      for (int i = 1; i <= 3; i++) {
        Vertex_handle vh = f.first->vertex((4 + f.second + i * sign) % 4);
        int vid = vh->info().id;
        if (used.find(vid) == used.end()) {
          sorted.push_back(vh);
          used[vid] = sorted.size();
        }
        faces.back()[i-1] = used[vid];
      }
    }
  }

  FILE *fo = fopen(name.c_str(), "w");
  for (Vertex_handle &vh : sorted) {
    auto &p = vh->point();
    fprintf(fo, "v %lf %lf %lf\n", p.x(), p.y(), p.z());
  }

  for (auto &face : faces) {
    fprintf(fo, "f %d %d %d\n", face[0], face[1], face[2]);
  }
  fclose(fo);
}

// Returns the cosine angle between the facet with inward-facing normal and 
// circimscribed sphere in the tetrahedron.
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

// Returns true if tetrahedron with this id is inside the mesh.
int DelaunayMesh::IsInsideSurface(int id) {
  return g->what_segment(nodes[id]) == Graph::SINK;
}

// Runs the main algorithm. Extract mesh surface given cameras' centers.
void DelaunayMesh::ExtractSurface(std::vector<Point> &camera) {
  BuildAABBTree();
  AssignIds();

  nodes.resize(d.number_of_cells());
  cost_unary.resize(d.number_of_cells());
  cost_binary.resize(d.number_of_cells());
  lock.resize(d.number_of_cells());

	g = new Graph();
  for (int i = 0; i < nodes.size(); i++) {
    nodes[i] = g->add_node();
    omp_init_lock(&lock[i]);
  }

  double t = timestamp();
  std::vector<Vertex_handle> vs(d.number_of_vertices());
  auto vit = vs.begin();
  for (auto it = d.finite_vertices_begin(); it != d.finite_vertices_end(); it++, vit++) 
    *vit = it;

  printf("Constructing Graph\n");
  printf("  Surface Visibility ... "); fflush(stdout);

#pragma omp parallel for
  for (int i = 0; i < vs.size(); i++) {
    auto p = vs[i]->point();
    std::set<int> &cams = vs[i]->info().cams;

    Point point = Point(p.x(), p.y(), p.z());
    double cost = cams.size();
    for (int camId : cams) {
      AssignCost(point, camera[camId], cost);
    }
  }
  printf("Done %f\n", timestamp() - t);


  printf("  Surface Quality ... "); fflush(stdout);
  t = timestamp();
  for (auto &triangle : triangles) {
    // Surface quality. Higher -> smoother.
    //double lambda = 0.5;
    double sq = surface_quality_lambda * (1 - std::min(
        CosineOfCircumsphere(triangle.f, triangle.n), 
        -CosineOfCircumsphere(d.mirror_facet(triangle.f), triangle.n)));

    int id[2];
    GetIncidentTetrahedrons(&triangle, id[0], id[1]);
    AddCostBinary(id[0], id[1], sq, sq);
  }
  printf("Done %f\n", timestamp() - t);
  tocaList();

  for (int i = 0; i < nodes.size(); i++) {
    for (auto const &k : cost_binary[i]) {
      if (k.first > i) {
        g->add_edge(nodes[i], nodes[k.first], cost_binary[i][k.first], cost_binary[k.first][i]);
      }
    }
  }
  for (int i = 0; i < nodes.size(); i++) {
    g->add_tweights(nodes[i], cost_unary[i].first, cost_unary[i].second);
  }

  printf("Solving min-cut ... "); fflush(stdout);
  t = timestamp();
  Graph::flowtype flow = g->maxflow();
  printf("Done %f\n", timestamp() - t);
  printf("flow = %lf\n", flow);

  for (int i = 0; i < nodes.size(); i++) {
    omp_destroy_lock(&lock[i]);
  }
}
