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
#include <set>

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

  // Could use Epick kernel, but Simple_cartesian is supposed to be fastest.
  typedef CGAL::Simple_cartesian<double> KSC; 
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K> Vb;
  typedef CGAL::Triangulation_cell_base_with_info_3<TetraInfo, K> Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
  typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
  typedef DelaunayMesh::Tds::Cell_handle Cell_handle;
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

  // Adds depth point to the Delaunay.
  void AddPoint(const Delaunay::Point &p, const VertexInfo &vi);

  // Builds AABB tree for fast ray intersection. Return the number of triangles.
  int BuildAABBTree();
  void RayIntersect(const DelaunayMesh::KSC::Ray_3 &ray, std::list<Tree::Intersection_and_primitive_id<KSC::Ray_3>::Type> &intersections);
  int AssignTetrahedronIds();

  template <typename T>
  static double DistanceSquared(const T &p0, const T &p1);

  Delaunay d;
  Tree tree;
private:
  // If a new vertex is closer than merge_threshold to an existing vertex, 
  // merge the two vertices. Otherwise, add the new vertex into Delaunay.
  // TODO(supasorn): make a setter or option for this.
  double merge_threshold = 0.00001;

  // Triangle in AABB. This vector MUST persist throughout AABB usage. 
  std::vector<FacetAndNormal> triangles;

};

#endif
