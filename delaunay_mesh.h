#ifndef MVSA_MESH_EXTRACTION_DELAUNAY_MESH_
#define MVSA_MESH_EXTRACTION_DELAUNAY_MESH_

#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <Eigen/Core>
#include <glog/logging.h>

#include <set>
#include <sys/time.h>
#include <time.h>
#include "omp.h"

#include "maxflow/v2_adjacency_list/graph.h"
// This class is used to perform mesh extraction. The input is a set of point cloud.
// Each point has a world-coordinate and a list of camera visible from that point.
// The output is a visibility-consistent triangle mesh. 
class DelaunayMesh {
public:
  // Vertex info class for used with Delaunay. This stores number of merged points
  // and a list of camaera visible from this vertex.
  struct VertexInfo {
    std::set<int> cams;
    int n;
    int id;

    VertexInfo() : n(1) {}
    VertexInfo(int cam0, int cam1) : n(1) { 
      cams.insert(cam0); 
      cams.insert(cam1);
    }
  };

  // Tetrahedron info class for used with Delaunay. It stores the id of this 
  // tetrahedron.
  struct TetraInfo {
    int id;
  };

  // Could use Exact_predicates_inexact_constructions_kernel for AABB, but 
  // Simple_cartesian is supposed to be fastest.
  typedef CGAL::Simple_cartesian<double> KSC; 
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

  typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K> Vb;
  typedef CGAL::Triangulation_cell_base_with_info_3<TetraInfo, K> Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
  typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
  typedef DelaunayMesh::Tds::Cell_handle Cell_handle;
  typedef DelaunayMesh::Tds::Vertex_handle Vertex_handle;
  typedef DelaunayMesh::Delaunay::Point   DPoint;
  typedef Delaunay::Facet Facet;

  struct FacetAndNormal {
    Facet f;
    // Unnormalized normal used for figuring out whether a camera ray intersects
    // front or back facet.
    K::Vector_3 n;

    // Assigns f and compute n.
    FacetAndNormal(Facet f);
  };


  typedef std::vector<FacetAndNormal>::const_iterator Iterator;

  // The following primitive provides the conversion facilities between
  // the custom triangle and point types and the CGAL ones.
  struct TrianglePrimitive {
  typedef KSC::Point_3    Point; // CGAL 3D point type.
  typedef KSC::Triangle_3 Datum; // CGAL 3D triangle type.
  public:
    // This is the type of data that the queries returns. 
    typedef const FacetAndNormal* Id;
    // CGAL types returned.
  private:
    Id m_pt; // This is what the AABB tree stores internally.
  public:
    TrianglePrimitive() {} 

    // The following constructor is the one that receives the iterators from the 
    // iterator range given as input to the AABB_tree.
    TrianglePrimitive(Iterator it) : m_pt(&(*it)) {}
    const Id& id() const { return m_pt; }

    // Utility function to convert a custom point type to CGAL point type.
    Point convert(const DPoint &p) const;
    
    // On-the-fly conversion from the internal data to the CGAL types.
    Datum datum() const; 

    // Returns a reference point which must be on the primitive.
    Point reference_point() const;
  };

  typedef CGAL::AABB_traits<KSC, TrianglePrimitive> My_AABB_traits;
  typedef CGAL::AABB_tree<My_AABB_traits> Tree;
  typedef Tree::Primitive_id Primitive_id;
  typedef TrianglePrimitive::Point Point;
  typedef TrianglePrimitive::Datum Datum;
  typedef Tree::Intersection_and_primitive_id<KSC::Ray_3>::Type IntersectReturn;

  // Adds depth point to the Delaunay.
  void AddPoint(const Delaunay::Point &p, const VertexInfo &vi);

  // Runs the main algorithm. Extract mesh surface given cameras' centers.
  void ExtractSurface(std::vector<Point> &camera);

  // Returns true if tetrahedron with this id is inside the mesh.
  int IsInsideSurface(int id);

  // Saves mesh to obj file.
  void SaveObj(std::string name);

  Delaunay d;
  Tree tree;

  // Triangle in AABB. This vector MUST persist throughout AABB usage. 
  std::vector<FacetAndNormal> triangles;

private:
  // If a new vertex is closer than merge_threshold to an existing vertex, 
  // merge the two vertices. Otherwise, add the new vertex into Delaunay.
  // TODO(supasorn): make a setter or option for this.
  double merge_threshold = 0.00001;
  //double merge_threshold = 0.00005;
  //double merge_threshold = 0.0000001;
  //
  double surface_quality_lambda = 10;

  // Graph-cut nodes.
  std::vector<Graph::node_id> nodes;
  Graph *g;

  // cost_binary[i][j] is the edge cost from node i to node j.
  std::vector<std::map<int, double> > cost_binary;

  // cost_unary[i] = (src, sink) contains the edge costs from node i
  // to source and sink.
  std::vector<std::pair<double, double> > cost_unary;

  std::vector<omp_lock_t> lock;

  // Adds an edge for node id0 to source and sink.
  void AddCostUnary(int id0, double src, double snk);

  // Adds an edge from and to node id0 and id1 with cost0 and cost1 respectively.
  void AddCostBinary(int id0, int id1, double cost0, double cost1);

  // Adds a set of edges for the network flow for point and camera with cost.
  void AssignCost(Point &point, Point &camera, double cost); 

  // Builds AABB tree for fast ray intersection. Return the number of triangles.
  int BuildAABBTree();

  // Assigns IDs to all vertices and tetrahedrons after Delaunay is complete. Returns number of cells.
  void AssignIds();

  // Returns id0, and id1 of tetrahedra incident to facet f.
  void GetIncidentTetrahedrons(const FacetAndNormal *f, int &id0, int &id1);
};

#endif
