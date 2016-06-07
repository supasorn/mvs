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
//#include "gflags/gflags.h"
//#include "glog/logging.h"

DelaunayMesh dm;
using namespace std;

//#include <gflags/gflags.h>
//using namespace gflags;

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


#include "debug.h"

vector<KPoint> camera;

vector<int> nodesOut;
DEFINE_string(guy, "George_W_Bush", "a");

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

  glPointSize(4);
  glBegin(GL_POINTS);
  for (int i = 0; i < camera.size(); i++) {
    glColor3ub(255, 0, 0);
    //glVertex3f(camera[i].at<double>(0, 0), camera[i].at<double>(1, 0), camera[i].at<double>(2, 0));
    glVertex3f(camera[i].x(), camera[i].y(), camera[i].z());
  }
  glEnd();

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
    //if (g->what_segment(nodes[id0]) != g->what_segment(nodes[id1])) {
    if (nodesOut[id0] != nodesOut[id1]) {
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
  glDisable(GL_CULL_FACE);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void keyboard2(unsigned char key, int x, int y) {
  DelaunayMesh::Delaunay &d = dm.d;
  if (key == 'a') {
    float c = 0.5;
    dm.AddPoint(DPoint(rand1(c), rand1(c), rand1(c)), VertexInfo());
    //RebuildTree();
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


  //RebuildTree();

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



void ExtractSurface() {
  for (auto &triangle : dm.triangles) {
    int id[2];
    dm.GetIncidentTetrahedrons(&triangle, id[0], id[1]);
    if (nodesOut[id[0]] != nodesOut[id[1]]) {
      addFacet(triangle.f, nodesOut[id[0]]);
    }
  }
  outObj();
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
    //camera[i] = Mat(3, 1, CV_64F);
    double t[3];
    fread(t, sizeof(double), 3, fi);
    camera[i] = KPoint(t[0], t[1], t[2]);
  }

  fread(&n, sizeof(int), 1, fi);
  printf("#points = %d\n", n);
  float p[3], c[3];
  int vc[2];
  for (int i = 0; i < n; i++) {
    fread(p, sizeof(float), 3, fi);
    fread(c, sizeof(float), 3, fi);
    fread(vc, sizeof(int), 2, fi);

    if (i % 100 < 70) continue;

    dm.AddPoint(DPoint(p[0], p[1], p[2]), VertexInfo(vc[0], vc[1]));

  }
  fclose(fi);
  printf("Done reading\n");

  //RebuildTree();
  //dm.BuildAABBTree();
  //dm.AssignTetrahedronIds();
  //iterateAllRays();
  dm.ExtractSurface(camera);
  nodesOut.resize(dm.d.number_of_cells());
  int count = 0;
  for (int i = 0; i < nodesOut.size(); i++) {
    nodesOut[i] = dm.IsInsideSurface(i);
    count += nodesOut[i]; 
    //printf("%d\n", nodesOut[i]);
  }
  printf("count = %d, %d\n", count, nodesOut.size() - count);
  dm.SaveObj("debug2.obj");


  //ExtractSurface();

  Viz::setDisplayCallback(display);
  Viz::setKeyboardCallback(keyboard2);
  Viz::setInitCallback(init);
  Viz::startWindow(800, 800);

}

void testGraph() {
  Graph::node_id nodes[2];
	Graph *g = new Graph();
	nodes[0] = g -> add_node();
	nodes[1] = g -> add_node();
	g -> add_tweights(nodes[0], 100, 1);
	g -> add_tweights(nodes[1], 1, 100);
  g -> add_edge(nodes[0], nodes[1], 3, 8);
  g -> add_edge(nodes[0], nodes[1], 3, 0);

	Graph::flowtype flow = g -> maxflow();

	printf("Flow = %f\n", flow);
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
int main(int argc, char** argv) {
  VLOG (1) << "eys";
  google::ParseCommandLineFlags(&argc, &argv, false);
  google::InitGoogleLogging(argv[0]);
  printf("%s\n", FLAGS_guy.c_str());
  //testGraph();
  //test2();
  test3();
  //test4();

  return 0;
}
