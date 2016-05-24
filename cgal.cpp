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
typedef DelaunayMesh::Tds::Cell_handle Cell_handle;
DelaunayMesh dm;


using namespace std;



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

//void printFacet(Delaunay::Finite_facets_iterator &f) {
void printFacet(DelaunayMesh::Delaunay::Facet f) {
  auto v0 = f.first->vertex((f.second + 1) % 4);
  auto v1 = f.first->vertex((f.second + 2) % 4);
  auto v2 = f.first->vertex((f.second + 3) % 4);
  printPoint(v0->point());
  printPoint(v1->point());
  printPoint(v2->point());
}

//typedef CGAL::Epick K;
typedef CGAL::Simple_cartesian<double> K; // Could use Epick, but Simple_cartesian is supposed to be fastest.
// custom point type
struct My_point {
    double m_x;
    double m_y;
    double m_z;
    My_point(const double x,
        const double y,
        const double z)
        : m_x(x), m_y(y), m_z(z) {}
};
// custom triangle type with
// three pointers to points
struct My_triangle {
    My_point *m_pa;
    My_point *m_pb;
    My_point *m_pc;
    My_triangle(My_point *pa,
        My_point *pb,
        My_point *pc)
        : m_pa(pa), m_pb(pb), m_pc(pc) {}
};

struct TetraTriangle {
  Cell_handle ch;
  int vi;
  TetraTriangle(Cell_handle c, int vi) : ch(c), vi(vi) {}

};
// the custom triangles are stored into a vector
//typedef std::vector<My_triangle>::const_iterator Iterator;
typedef std::vector<TetraTriangle>::const_iterator Iterator;
// The following primitive provides the conversion facilities between
// the custom triangle and point types and the CGAL ones
struct My_triangle_primitive {
public:
    // this is the type of data that the queries returns. For this example
    // we imagine that, for some reasons, we do not want to store the iterators
    // of the vector, but raw pointers. This is to show that the Id type
    // does not have to be the same as the one of the input parameter of the 
    // constructor.
    //typedef const My_triangle* Id;
    typedef const TetraTriangle* Id;
    // CGAL types returned
    typedef K::Point_3    Point; // CGAL 3D point type
    typedef K::Triangle_3 Datum; // CGAL 3D triangle type
private:
    Id m_pt; // this is what the AABB tree stores internally
public:
    My_triangle_primitive() {} // default constructor needed
    // the following constructor is the one that receives the iterators from the 
    // iterator range given as input to the AABB_tree
    My_triangle_primitive(Iterator it)
        : m_pt(&(*it)) {}
    const Id& id() const { return m_pt; }
    // utility function to convert a custom 
    // point type to CGAL point type.
    Point convert(const My_point *p) const {
        return Point(p->m_x,p->m_y,p->m_z);
    }

    Point convert(const DPoint &p) const {
        return Point(p[0], p[1], p[2]);
    }
    /*
    // on the fly conversion from the internal data to the CGAL types
    Datum datum() const {
        return Datum(convert(m_pt->m_pa),
            convert(m_pt->m_pb),
            convert(m_pt->m_pc));
    }
    // returns a reference point which must be on the primitive
    Point reference_point() const
    { return convert(m_pt->m_pa); }*/
    Datum datum() const {
      auto v0 = m_pt->ch->vertex((m_pt->vi + 1) % 4);
      auto v1 = m_pt->ch->vertex((m_pt->vi + 2) % 4);
      auto v2 = m_pt->ch->vertex((m_pt->vi + 3) % 4);
      return Datum(convert(v0->point()), convert(v1->point()), convert(v2->point()));
    }
    Point reference_point() const { 
      auto v0 = m_pt->ch->vertex((m_pt->vi + 1) % 4);
      auto p = v0->point();
      return Point(p[0], p[1], p[2]);
    }
};


typedef K::FT FT;
typedef K::Ray_3 Ray;
typedef K::Line_3 Line;
typedef K::Point_3 KPoint;
typedef K::Triangle_3 Triangle;

//typedef std::list<Triangle>::iterator Iterator;
//typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
//typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
//typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

typedef CGAL::AABB_traits<K, My_triangle_primitive> My_AABB_traits;
typedef CGAL::AABB_tree<My_AABB_traits> Tree;

typedef Tree::Primitive_id Primitive_id;

/*
void testRay() {
  My_point a(1, 0, 0);
  My_point b(0, 1, 0);
  My_point c(0, 0, 1);
  My_point d(0, 0, 0);
  std::vector<My_triangle> triangles;
  triangles.push_back(My_triangle(&a, &b, &c));
  triangles.push_back(My_triangle(&a, &b, &d));
  triangles.push_back(My_triangle(&a, &d, &c));

  //KPoint a(1.0, 0.0, 0.0);
  //KPoint b(0.0, 1.0, 0.0);
  //KPoint c(0.0, 0.0, 1.0);
  //KPoint d(0.0, 0.0, 0.0);
  //std::list<Triangle> triangles;
  //triangles.push_back(Triangle(a,b,c));
  //triangles.push_back(Triangle(a,b,d));
  //triangles.push_back(Triangle(a,d,c));
  //
  Tree tree(triangles.begin(),triangles.end());

  KPoint p0(0.4, 0.01, 0.4);
  KPoint p1(0.5, 1, 0.5);
  Ray ray_query(p0 ,p1);

  std::list<Primitive_id> primitives;
  tree.all_intersected_primitives(ray_query, std::back_inserter(primitives));
  printf("%d\n", primitives.size());
  for (auto i = primitives.begin(); i != primitives.end(); i++) {
    auto triangle = *i;
    //auto &p0 = triangle->vertex(0);
    //auto &p1 = triangle->vertex(1);
    //auto &p2 = triangle->vertex(2);
    //printf("%f %f %f\n", p0[0], p0[1], p0[2]);
    //printf("%f %f %f\n", p1[0], p1[1], p1[2]);
    //printf("%f %f %f\n", p2[0], p2[1], p2[2]);
    auto &p0 = triangle->m_pa;
    auto &p1 = triangle->m_pb;
    auto &p2 = triangle->m_pc;
    printf("%f %f %f\n", p0->m_x, p0->m_y, p0->m_z);
    printf("%f %f %f\n", p1->m_x, p1->m_y, p1->m_z);
    printf("%f %f %f\n", p2->m_x, p2->m_y, p2->m_z);

  }


  std::cout << tree.number_of_intersected_primitives(ray_query)
    << " intersections(s) with ray query" << std::endl;
  // compute closest point and squared distance
  exit(0);
}*/

Tree tree;

KPoint p0(0.1, 0.1, 0.1);
KPoint p1(2, 3, 4);
vector<pair<Cell_handle, int> > intersected;

void display() {
  DelaunayMesh::Delaunay &d = dm.d;
  
  glBegin(GL_LINES);
  glColor4ub(255, 255, 255, 80);
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
  }

}
void init() {
  glDepthFunc(GL_ALWAYS);
  glEnable(GL_POINT_SMOOTH); 

  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void RebuildTree() {
  DelaunayMesh::Delaunay &d = dm.d;
  std::vector<TetraTriangle> triangles;
  for (auto it = d.finite_facets_begin(); it != d.finite_facets_end(); it++) {
    triangles.push_back(TetraTriangle(it->first, it->second));
  }
  tree.clear();
  tree.insert(triangles.begin(),triangles.end());

  Ray ray_query(p0 ,p1);

  std::list<Primitive_id> primitives;
  tree.all_intersected_primitives(ray_query, std::back_inserter(primitives));
  
  intersected.clear();
  for (auto i = primitives.begin(); i != primitives.end(); i++) {
    auto triangle = *i;
    intersected.push_back(make_pair(triangle->ch, triangle->vi));
  }
  
}
void keyboard2(unsigned char key, int x, int y) {
  DelaunayMesh::Delaunay &d = dm.d;
  if (key == 'a') {
    float c = 0.5;
    dm.AddPoint(DPoint(rand1(c), rand1(c), rand1(c)), VertexInfo(255, 0, 0));
    RebuildTree();
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
int main() {
  //testRay();
  test2();

  return 0;
}

