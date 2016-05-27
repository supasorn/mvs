#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <iostream>
#include <set>
#include "delaunay_mesh.h"
#include "viz.h"
#include "maxflow/v2_adjacency_list/graph.h"

DelaunayMesh dm;
using namespace std;

typedef DelaunayMesh::DPoint DPoint;
typedef DelaunayMesh::Cell_handle Cell_handle;
typedef DelaunayMesh::Tree::Primitive_id Primitive_id;
typedef DelaunayMesh::VertexInfo VertexInfo;

typedef DelaunayMesh::KSC::Ray_3 Ray;
typedef DelaunayMesh::K::Vector_3 Vector_3;
typedef DelaunayMesh::KSC::Line_3 Line;
typedef DelaunayMesh::KSC::Point_3 KPoint;
typedef DelaunayMesh::KSC::Segment_3 Segment;
typedef DelaunayMesh::KSC::Triangle_3 Triangle;
typedef DelaunayMesh::KSC::Intersect_3 Intersect;
typedef DelaunayMesh::Tree::Intersection_and_primitive_id<DelaunayMesh::KSC::Ray_3>::Type IntersectReturn;


vector<Mat> camera;

vector<Graph::node_id> nodes;
Graph *g;

float rand1(float c = 1) {
  return c * (2.0 * rand() / RAND_MAX - 1);
}
void printPoint(DelaunayMesh::DPoint &p) {
  printf("[%f %f %f]\n", p.x(), p.y(), p.z());
}


void keyboard(unsigned char key, int x, int y) {
  DelaunayMesh::Delaunay &d = dm.d;
  float c = 0.5;
  auto p = d.finite_vertices_begin()->point();
  if (key == 'a') {
    d.move(d.finite_vertices_begin(), DPoint(p.x()+0.1, p.y(), p.z()));
  } else if (key == 's') {
    d.move(d.finite_vertices_begin(), DPoint(p.x()-0.1, p.y(), p.z()));
  }
}

//void printFacet(Delaunay::Finite_facets_iterator &f) {
void printFacet(DelaunayMesh::Delaunay::Facet f) {
  auto v0 = f.first->vertex((f.second + 1) % 4);
  auto v1 = f.first->vertex((f.second + 2) % 4);
  auto v2 = f.first->vertex((f.second + 3) % 4);
  printPoint(v0->point());
  printPoint(v1->point());
  printPoint(v2->point());
}


KPoint p0(0, 0, 0);
KPoint p1(2, 3, 4);
//KPoint p1(1, 0, 0);
vector<pair<Cell_handle, int> > intersected;

void display() {
  DelaunayMesh::Delaunay &d = dm.d;

  glBegin(GL_LINES);

  glColor4ub(255, 255, 255, 80);
  for (auto it = d.finite_edges_begin(); it != d.finite_edges_end(); it++) {
    auto p0 = it->first->vertex(it->second)->point();
    auto p1 = it->first->vertex(it->third)->point();
    glVertex3f(p0.x(), p0.y(), p0.z());
    glVertex3f(p1.x(), p1.y(), p1.z());
  }
  glEnd();

  for (auto it = d.finite_vertices_begin(); it != d.finite_vertices_end(); it++) {
    auto p = it->point();
    glPointSize(1 + 0.25 * it->info().n);
    glBegin(GL_POINTS);
    glColor3ub(0, 255, 0);
    glVertex3f(p.x(), p.y(), p.z());
    glEnd();
  }

  for (auto it = d.finite_facets_begin(); it != d.finite_facets_end(); it++) {
    const DelaunayMesh::Facet &f = *it;

    int id0 = f.first->info().id;
    int id1 = d.mirror_facet(f).first->info().id;
    if (g->what_segment(nodes[id0]) != g->what_segment(nodes[id1])) {
      glColor4ub(28, 159, 213, 80);
      glBegin(GL_TRIANGLES);
      auto v0 = f.first->vertex((f.second + 1) % 4)->point();
      auto v1 = f.first->vertex((f.second + 2) % 4)->point();
      auto v2 = f.first->vertex((f.second + 3) % 4)->point();
      glVertex3f(v0.x(), v0.y(), v0.z());
      glVertex3f(v1.x(), v1.y(), v1.z());
      glVertex3f(v2.x(), v2.y(), v2.z());
      glEnd();
    }
  }

  /*
  glPointSize(10);
  glBegin(GL_POINTS);
  glColor3ub(0, 255, 0);
  glVertex3f(p0[0], p0[1], p0[2]);
  glEnd();


  glBegin(GL_LINES);
  glColor3ub(255, 255, 255);
  glVertex3f(p0[0], p0[1], p0[2]);
  glVertex3f(p1[0], p1[1], p1[2]);
  glEnd();
  for (int i = 0; i < intersected.size(); i++) {
    glBegin(GL_TRIANGLES);
    glColor4ub(28, 159, 213, 80);
    //glColor4ub(237, 155, 13, 80);
    auto p = intersected[i].first->vertex((intersected[i].second + 1) % 4)->point();
    glVertex3f(p[0], p[1], p[2]);
    p = intersected[i].first->vertex((intersected[i].second + 2) % 4)->point();
    glVertex3f(p[0], p[1], p[2]);
    p = intersected[i].first->vertex((intersected[i].second + 3) % 4)->point();
    glVertex3f(p[0], p[1], p[2]);
    glEnd();
  }*/

}
void init() {
  glDepthFunc(GL_ALWAYS);
  glEnable(GL_POINT_SMOOTH);

  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}


void RebuildTree() {
  dm.BuildAABBTree();

  Ray ray_query(p0 ,p1);


  std::list<IntersectReturn> inters;
  dm.RayIntersect(ray_query, inters);

  intersected.clear();
  for (auto &inter : inters) {
    if (boost::get<KPoint>(&(inter.first))) {
      printf("point\n");
      const DelaunayMesh::Facet *f = &inter.second->f;

      cout << inter.first << endl;
      cout << f->first->info().id << endl;
      intersected.push_back(make_pair(f->first, f->second));
    } else if (boost::get<Segment>(&(inter.first))) {
      printf("segment\n");
      cout << inter.first;
      cout << inter.second->f.first->info().id << endl;
    }
  }

  /*
  for (auto i = primitives.begin(); i != primitives.end(); i++) {
    auto triangle = *i;
    intersected.push_back(make_pair(triangle->first, triangle->second));
  }*/

}
void keyboard2(unsigned char key, int x, int y) {
  DelaunayMesh::Delaunay &d = dm.d;
  if (key == 'a') {
    float c = 0.5;
    dm.AddPoint(DPoint(rand1(c), rand1(c), rand1(c)), VertexInfo());
    RebuildTree();
  }
}
void test2() {
  float c = 0.5;
  DelaunayMesh::Delaunay &d = dm.d;
  dm.AddPoint(DPoint(0, 0, 0), VertexInfo());
  dm.AddPoint(DPoint(c, 0, 0), VertexInfo());
  dm.AddPoint(DPoint(0, c, 0), VertexInfo());
  dm.AddPoint(DPoint(0, 0, c), VertexInfo());
  dm.AddPoint(DPoint(c, c, c), VertexInfo());


  int id = 0;
  for (auto it = d.all_cells_begin(); it != d.all_cells_end(); ++it) {
    it->info().id = id++;
  }

  for (auto it = d.finite_facets_begin(); it != d.finite_facets_end(); it++) {
    //printf("%d\n", d.triangle(face));
    printFacet(*it);
    printFacet(d.mirror_facet(*it));
    printf("%d - %d\n", it->first->info().id, d.mirror_facet(*it).first->info().id);

    printf("\n");
  }


  RebuildTree();

  /*
  KPoint p0(0.1, 0.1, 0.1);
  KPoint p1(0.2, 1, 0.1);
  Ray ray_query(p0 ,p1);

  std::list<Primitive_id> primitives;
  tree.all_intersected_primitives(ray_query, std::back_inserter(primitives));
  printf("%d\n", primitives.size());
  for (auto i = primitives.begin(); i != primitives.end(); i++) {
    auto triangle = *i;
    //auto &p0 = triangle->m_pa;
    //auto &p1 = triangle->m_pb;
    //auto &p2 = triangle->m_pc;
    //printf("%f %f %f\n", p0->m_x, p0->m_y, p0->m_z);
    //printf("%f %f %f\n", p1->m_x, p1->m_y, p1->m_z);
    //printf("%f %f %f\n", p2->m_x, p2->m_y, p2->m_z);
  }*/

  Viz::setDisplayCallback(display);
  Viz::setKeyboardCallback(keyboard2);
  Viz::setInitCallback(init);
  Viz::startWindow(800, 800);
}

//int SameSide(DelaunayMesh::KSC::Vector_3 ray, const Vector_3 &normal) {
//int SameSide(DelaunayMesh::KSC::Vector_3 ray, const DelaunayMesh::K::Vector_3 &normal) {
  //return DelaunayMesh::KSC::Vector_3(ray.x(), ray.y(), ray.z()) * DelaunayMesh::KSC::Vector_3(normal.x(), normal.y(), normal.z()) > 0;
//}

void iterateAllRays() {
  DelaunayMesh::Delaunay &d = dm.d;

  typedef pair<const DelaunayMesh::FacetAndNormal*, double> sortedFacets; 
  struct Comparator {
    static bool compare(sortedFacets &a, sortedFacets &b) { return a.second < b.second; }
  };

  nodes.resize(d.number_of_cells());
	g = new Graph();

  for (int i = 0; i < nodes.size(); i++) {
    nodes[i] = g->add_node();
  }

  for (auto it = d.finite_vertices_begin(); it != d.finite_vertices_end(); it++) {
    auto p = it->point();
    set<int> &cams = it->info().cams;
    //printf("%f %f %f %d\n", p.x(), p.y(), p.z(), cams.size());

    KPoint p1 = KPoint(p.x(), p.y(), p.z());
    for (int camId : cams) {
      KPoint p0 = KPoint(camera[camId].at<double>(0, 0),
          camera[camId].at<double>(1, 0),
          camera[camId].at<double>(2, 0));

      double cameraToPointDist = DelaunayMesh::DistanceSquared(p0, p1);
      Ray ray_query(p0, p1);
      auto toCamera = p0 - p1;

      std::list<IntersectReturn> inters;
      dm.RayIntersect(ray_query, inters);

      // Store intersected facets and its distance from the depth point.
      vector<sortedFacets> dists;
      for (auto &inter: inters) {

        // Only retrieve point intersection here because segment intersection
        // always has triangle's normal perpedicular to the ray.
        if (const KPoint *p = boost::get<KPoint>(&(inter.first))) {
          const DelaunayMesh::Facet &f = inter.second->f;
          double dist = DelaunayMesh::DistanceSquared(*p, p0) - cameraToPointDist;
          //printf(" tetra1: %d\n", f.first->info().id);
          //printf(" tetra2: %d\n", d.mirror_facet(f).first->info().id);
          //printf("%f %f %f\n", p->x(), p->y(), p->z());
          //printf("%f\n", dist);
          //cout << inter.first << endl;
          //cout << f->first->info().id << endl;

          //intersected.push_back(make_pair(f.first, f.second));
          
          //cout << " - " << f.second << ":" << b << endl;
          //printf("[%f %f %f]\n", tt.x(), tt.y(), tt.z());
          //printf("infinite %d\n", inf);
          //printFacet(*f);
          dists.push_back(make_pair(inter.second, dist));
        } 
      }
      //printf("%d\n", camId);
      sort(dists.begin(), dists.end(), Comparator::compare);

      // ------ 0
      auto &normal = dists[0].first->n;
      // True iff normal and ray from point to camera are on the same side.
      int ss = toCamera.x() * normal.x() + toCamera.y() * normal.y() + toCamera.z() * toCamera.z() > 0;
      int odd = dists[0].first->f.second % 2;

      if (ss ^ odd) {
        g->add_tweights(nodes[d.mirror_facet(dists[0].first->f).first->info().id], 100, 0);
      } else {
        g->add_tweights(nodes[dists[0].first->f.first->info().id], 100, 0);
      }
      // ------ 0
      //
      int i;
      for (i = 0; i < dists.size(); i++) {
        if (dists[i].second == 0) {
          while (i < dists.size() && dists[i++].second == 0);
          break;
        }
        //printf("  %e %d\n", dists[i].second, dists[i].second == 0);
        auto &f = dists[i].first->f; // Facet
        auto &ch = f.first; // Cell handle 
        const int &id = f.second; // Vertex id.

        //auto diff = ch->vertex(id)->point() - ch->vertex((id + 1) % 4)->point();
        //int b = SameSide(diff, dists[i].first->n);
        int inf = d.is_infinite(ch->vertex(id));


        auto &normal = dists[i].first->n;
        int ss = toCamera.x() * normal.x() + toCamera.y() * normal.y() + toCamera.z() * toCamera.z() > 0;
        int odd = dists[i].first->f.second % 2;

        int id0 = dists[i].first->f.first->info().id;
        int id1 = d.mirror_facet(dists[i].first->f).first->info().id;

        if (ss ^ odd) {
          g->add_edge(nodes[id0], nodes[id1], 100, 0);
        } else {
          g->add_edge(nodes[id0], nodes[id1], 0, 100);
        }
        //printf("  %e %d\n", dists[i].second, b);
      }
      if (i < dists.size()) {
        auto &normal = dists[i].first->n;
        // True iff normal and ray from point to camera are on the same side.
        int ss = toCamera.x() * normal.x() + toCamera.y() * normal.y() + toCamera.z() * toCamera.z() > 0;
        int odd = dists[i].first->f.second % 2;

        if (ss ^ odd) {
          g->add_tweights(nodes[d.mirror_facet(dists[i].first->f).first->info().id], 0, 100);
        } else {
          g->add_tweights(nodes[dists[i].first->f.first->info().id], 0, 100);
        }
      }
    }
  }
  printf("Solving min-cut\n");
	Graph::flowtype flow = g -> maxflow();
  int count = 0;
  for (int i = 0; i < nodes.size(); i++) {
    printf("%d %d\n", i, g->what_segment(nodes[i]) == Graph::SOURCE);
    count += g->what_segment(nodes[i]) == Graph::SOURCE;
  }
  printf("count = %d\n", count);
}

void test3() {
  DelaunayMesh::Delaunay &d = dm.d;

  FILE *fi = fopen("point.bin", "rb");
  int n;
  // Write the number of cameras, and cameras' centers.
  fread(&n, sizeof(int), 1, fi);
  printf("#cameras = %d\n", n);
  camera.resize(n);
  for (int i = 0; i < n; i++) {
    camera[i] = Mat(3, 1, CV_64F);
    fread(&camera[i].at<double>(0, 0), sizeof(double), 3, fi);
  }

  fread(&n, sizeof(int), 1, fi);
  printf("#points = %d\n", n);
  float p[3], c[3];
  int vc[2];
  for (int i = 0; i < n; i++) {
    fread(p, sizeof(float), 3, fi);
    fread(c, sizeof(float), 3, fi);
    fread(vc, sizeof(int), 2, fi);

    if (i % 100 < 90) continue;

    if (p[0] == 0 && p[1] == 0 && p[2] == 0)
      exit(0);

    dm.AddPoint(DPoint(p[0], p[1], p[2]), VertexInfo(vc[0], vc[1]));

    if (rand() % 60 == 0) {
      p0 = KPoint(p[0], p[1], p[2]);
      p1 = KPoint(
          camera[vc[0]].at<double>(0, 0),
          camera[vc[0]].at<double>(1, 0),
          camera[vc[0]].at<double>(2, 0));
    }

  }

  //RebuildTree();
  dm.BuildAABBTree();
  dm.AssignTetrahedronIds();

  iterateAllRays();

  Viz::setDisplayCallback(display);
  Viz::setKeyboardCallback(keyboard2);
  Viz::setInitCallback(init);
  Viz::startWindow(800, 800);

  fclose(fi);
  printf("Done reading\n");
}

void testGraph() {
  Graph::node_id nodes[2];
	Graph *g = new Graph();
	nodes[0] = g -> add_node();
	nodes[1] = g -> add_node();
	g -> add_tweights(nodes[0], 100, 1);
	g -> add_tweights(nodes[1], 1, 100);
	//g -> add_edge(nodes[0], nodes[1], 3, 4);

	Graph::flowtype flow = g -> maxflow();

	printf("Flow = %d\n", flow);
	printf("Minimum cut:\n");
	if (g->what_segment(nodes[0]) == Graph::SOURCE)
		printf("node0 is in the SOURCE set\n");
	else
		printf("node0 is in the SINK set\n");
	if (g->what_segment(nodes[1]) == Graph::SOURCE)
		printf("node1 is in the SOURCE set\n");
	else
		printf("node1 is in the SINK set\n");
 
  exit(0);
	delete g;
}
int main() {
  //testGraph();
  //testRay();
  //test2();
  test3();

  return 0;
}
