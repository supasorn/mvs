#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/IO/Color.h>
#include <iostream>

#include "viz.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K>                      Delaunay;
typedef Delaunay::Point   DPoint;


using namespace std;

Delaunay d;

void display() {
  glBegin(GL_POINTS);
  for (auto it = d.finite_vertices_begin(); it != d.finite_vertices_end(); it++) {
    auto p = it->point();
    glVertex3f(p.x(), p.y(), p.z());
  }
  glEnd();

  glBegin(GL_LINES);
  for (auto it = d.finite_edges_begin(); it != d.finite_edges_end(); it++) {
    //printf("%d\n", d.triangle(face));
    auto p0 = it->first->vertex(it->second)->point();
    auto p1 = it->first->vertex(it->third)->point();
    glVertex3f(p0.x(), p0.y(), p0.z());
    glVertex3f(p1.x(), p1.y(), p1.z());
  }
  glEnd();

}

float rand1(float c = 1) {
  return c * (2.0 * rand() / RAND_MAX - 1);
}
void printPoint(DPoint &p) {
  printf("%f %f %f\n", p.x(), p.y(), p.z());
}


void keyboard(unsigned char key, int x, int y) {
  float c = 0.5;
  auto p = d.finite_vertices_begin()->point();
  if (key == 'a') {
    d.move(d.finite_vertices_begin(), DPoint(p.x()+0.1, p.y(), p.z()));
  } else if (key == 's') {
    d.move(d.finite_vertices_begin(), DPoint(p.x()-0.1, p.y(), p.z()));
  }
}

void testDelaunay() {
  srand(0);
  float c = 0.5;
  for (int i = 0; i < 10; i++) {
    d.insert(DPoint(rand1(c), rand1(c), rand1(c)));
  }

  Viz::setDisplayCallback(display);
  Viz::setKeyboardCallback(keyboard);
  Viz::startWindow(800, 800);
}
int main() {

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

