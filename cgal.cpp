#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/IO/Color.h>
#include <iostream>
#include <set>

#include "viz.h"




class VertexInfo {
  public:
  unsigned char col[3];
  std::set<int> cams;
  int n;

  VertexInfo() {
    n = 1;
  }

  VertexInfo(unsigned char col0, unsigned char col1, unsigned char col2) {
    n = 1;
    col[0] = col0;
    col[1] = col1;
    col[2] = col2;
  }
};

class TetraInfo {
  public:
  int id;
};

class DelaunayMesh {
public:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K> Vb;
  typedef CGAL::Triangulation_cell_base_with_info_3<TetraInfo, K> Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
  typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;

  Delaunay d;

  double DistanceSquared(const Delaunay::Point &p0, const Delaunay::Point &p1) {
    double v[3] = {p0.x() - p1.x(), p0.y() - p1.y(), p0.z() - p1.z()};
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  }



  void AddPoint(const Delaunay::Point &p, const VertexInfo &vi) {
    auto v = d.nearest_vertex(p);
    // Found a nearest vertex within merge_threshold distance.
    if (d.is_vertex(v) && DistanceSquared(p, v->point()) < merge_threshold) {
      // TODO(Aek): merge camera list
      auto vp = v->point();
      int n = v->info().n;
      d.move(v, Delaunay::Point(
            (p.x() + vp.x() * n) / (n + 1),
            (p.y() + vp.y() * n) / (n + 1),
            (p.z() + vp.z() * n) / (n + 1)));
      v->info().n++;
    } else {
      auto v = d.insert(p);
      v->info() = vi;
    }

  }

private:
  double merge_threshold = 0.005;
};

typedef DelaunayMesh::Delaunay::Point   DPoint;
DelaunayMesh dm;


using namespace std;


void display() {
  DelaunayMesh::Delaunay &d = dm.d;
  glBegin(GL_LINES);
    glColor3ub(100, 100, 100);
  for (auto it = d.finite_edges_begin(); it != d.finite_edges_end(); it++) {
    //printf("%d\n", d.triangle(face));
    auto p0 = it->first->vertex(it->second)->point();
    auto p1 = it->first->vertex(it->third)->point();
    glVertex3f(p0.x(), p0.y(), p0.z());
    glVertex3f(p1.x(), p1.y(), p1.z());
  }
  glEnd();

  for (auto it = d.finite_vertices_begin(); it != d.finite_vertices_end(); it++) {
    auto p = it->point();
    glPointSize(4 * it->info().n);
    glBegin(GL_POINTS);
    glColor3ub(it->info().col[0], it->info().col[1], it->info().col[2]);
    glVertex3f(p.x(), p.y(), p.z());
    glEnd();
  }

}

float rand1(float c = 1) {
  return c * (2.0 * rand() / RAND_MAX - 1);
}
void printPoint(DPoint &p) {
  printf("%f %f %f\n", p.x(), p.y(), p.z());
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

void testDelaunay() {
  DelaunayMesh::Delaunay &d = dm.d;
  srand(0);
  float c = 0.5;
  for (int i = 0; i < 5; i++) {
    d.insert(DPoint(rand1(c), rand1(c), rand1(c)));
  }
  Viz::setDisplayCallback(display);
  Viz::setKeyboardCallback(keyboard);
  Viz::startWindow(800, 800);
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
void keyboard2(unsigned char key, int x, int y) {
  DelaunayMesh::Delaunay &d = dm.d;
  if (key == 'a') {
    float c = 0.5;
    dm.AddPoint(DPoint(rand1(c), rand1(c), rand1(c)), VertexInfo(255, 0, 0));
  }
}

void test2() {
  float c = 0.5;
  DelaunayMesh::Delaunay &d = dm.d;
  dm.AddPoint(DPoint(0, 0, 0), VertexInfo(255, 0, 0));
  dm.AddPoint(DPoint(c, 0, 0), VertexInfo(255, 0, 0));
  dm.AddPoint(DPoint(0, c, 0), VertexInfo(255, 0, 0));
  dm.AddPoint(DPoint(0, 0, c), VertexInfo(255, 0, 0));
  dm.AddPoint(DPoint(c, c, c), VertexInfo(255, 0, 0));
  //d.insert(DPoint(0, 0, 0));
  //d.insert(DPoint(c, 0, 0));
  //d.insert(DPoint(0, c, 0));
  //d.insert(DPoint(0, 0, c));
  //d.insert(DPoint(c, c, c));


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

  Viz::setDisplayCallback(display);
  Viz::setKeyboardCallback(keyboard2);
  Viz::startWindow(800, 800);
}
int main() {
  DelaunayMesh::Delaunay &d = dm.d;

  test2();
  testDelaunay();
  d.insert(DPoint(0,0,0));
  d.insert(DPoint(1,0,0));
  d.insert(DPoint(0,1,0));
  d.insert(DPoint(0,0,1));
  
  //auto v = d.nearest_vertex(DPoint(0.9, 0, 0));
  //DPoint p = v->point();
  //printf("%f %f %f\n", p.x(), p.y(), p.z());

  for (auto it = d.finite_facets_begin(); it != d.finite_facets_end(); it++) {
    //printf("%d\n", d.triangle(face));
    auto v0 = it->first->vertex((it->second + 1) % 4);
    auto v1 = it->first->vertex((it->second + 2) % 4);
    auto v2 = it->first->vertex((it->second + 3) % 4);
    printPoint(v0->point());
    printPoint(v1->point());
    printPoint(v2->point());
    printf("\n");
  }

  return 0;
}

