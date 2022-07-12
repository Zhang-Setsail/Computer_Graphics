#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Task 1.
    int point_number = points.size();
    std::vector<Vector2D> final_points;
    if (point_number == 1) {
      return points;
    }
    Vector2D mid_point;
    for (int i = 0; i < point_number-1; i++) {
      mid_point = (1 - t) * points[i] + t * points[i + 1];
      final_points.push_back(mid_point);
    }
    // final point number = init point number - 1
    return final_points;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Task 2.
    int point_number = points.size();
    std::vector<Vector3D> final_points;
    if (point_number == 1) {
      return points;
    }
    Vector3D mid_point;
    for (int i = 0; i < point_number-1; i++) {
      mid_point = (1 - t) * points[i] + t * points[i + 1];
      final_points.push_back(mid_point);
    }
    // final point number = init point number - 1
    return final_points;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Task 2.
    if (points.size() == 1) {
      return points[0];
    }
    std::vector<Vector3D> points_less;
    points_less = evaluateStep(points, t);
    return evaluate1D(points_less, t);
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Task 2.
    std::vector<Vector3D> points_on_v;
    Vector3D point_on_u;
    Vector3D point_on_uv;
    int n = controlPoints.size();
    for (int i = 0; i < n ; i++) {
      point_on_u = evaluate1D(controlPoints[i], u);
      points_on_v.push_back(point_on_u);      
    }
    point_on_uv = evaluate1D(points_on_v, v);
    return point_on_uv;
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Task 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.

    Vector3D sum_weighted_normal(0.0);
    HalfedgeCIter h = halfedge();
    do {
      if (!h->face()->isBoundary()) 
      {
        Vector3D vertex0 = position;
        Vector3D vertex1 = h->next()->vertex()->position;
        Vector3D vertex2 = h->next()->next()->vertex()->position;
        Vector3D v0_v1, v0_v2;
        double area;
        Vector3D area_weighted_normal;

        v0_v1 = vertex1 - vertex0;
        v0_v2 = vertex2 - vertex0;
        area = (cross(v0_v1,v0_v2).norm())/2;
        area_weighted_normal = area * h->face()->normal();
        sum_weighted_normal = sum_weighted_normal + area_weighted_normal;
      }
      h = h->twin()->next();    
    }while (h != halfedge());


    // //Test
    // cout<<sum_weighted_normal.unit()<<endl;
    // Vector3D test3d;
    // test3d = sum_weighted_normal;
    // test3d.normalize();
    // cout<<test3d<<endl;
    return sum_weighted_normal.unit();

  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Task 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
    if (!e0->isBoundary())
    {

      HalfedgeIter f1_h0 = e0->halfedge();
      HalfedgeIter f1_h1 = f1_h0->next();
      HalfedgeIter f1_h2 = f1_h1->next();
      HalfedgeIter f2_h0 = f1_h0->twin();
      HalfedgeIter f2_h1 = f2_h0->next();
      HalfedgeIter f2_h2 = f2_h1->next();
      HalfedgeIter f1_h1_twin = f1_h1->twin();
      HalfedgeIter f1_h2_twin = f1_h2->twin();
      HalfedgeIter f2_h1_twin = f2_h1->twin();
      HalfedgeIter f2_h2_twin = f2_h2->twin();
      VertexIter v0 = f1_h0->vertex();
      VertexIter v1 = f2_h0->vertex();
      VertexIter v2 = f1_h2->vertex();
      VertexIter v3 = f2_h2->vertex();
      EdgeIter f1_e1 = f1_h1->edge();
      EdgeIter f1_e2 = f1_h2->edge();
      EdgeIter f2_e1 = f2_h1->edge();
      EdgeIter f2_e2 = f2_h2->edge();
      FaceIter f1 = f1_h0->face();
      FaceIter f2 = f2_h0->face();

      f1_h0->setNeighbors(f1_h1,f2_h0,v2,e0,f1);
      f1_h1->setNeighbors(f1_h2,f2_h2_twin,v3,f2_e2,f1);
      f1_h2->setNeighbors(f1_h0,f1_h1_twin,v1,f1_e1,f1);
      f2_h0->setNeighbors(f2_h1,f1_h0,v3,e0,f2);
      f2_h1->setNeighbors(f2_h2,f1_h2_twin,v2,f1_e2,f2);
      f2_h2->setNeighbors(f2_h0,f2_h1_twin,v0,f2_e1,f2);
      f1_h1_twin->setNeighbors(f1_h1_twin->next(),f1_h2,v2,f1_e1,f1_h1_twin->face());
      f1_h2_twin->setNeighbors(f1_h2_twin->next(),f2_h1,v0,f1_e2,f1_h2_twin->face());
      f2_h1_twin->setNeighbors(f2_h1_twin->next(),f2_h2,v3,f2_e1,f2_h1_twin->face());
      f2_h2_twin->setNeighbors(f2_h2_twin->next(),f1_h1,v1,f2_e2,f2_h2_twin->face());
      v0->halfedge() = f2_h2;
      v1->halfedge() = f1_h2;
      v2->halfedge() = f2_h1;
      v3->halfedge() = f1_h1;
      e0->halfedge() = f1_h0;
      f1_e1->halfedge() = f1_h2;
      f1_e2->halfedge() = f2_h1;
      f2_e1->halfedge() = f2_h2;
      f2_e2->halfedge() = f1_h1;
      f1->halfedge() = f1_h0;
      f2->halfedge() = f2_h0;
    }
    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Task 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
    if (!e0->isBoundary())
    {
      HalfedgeIter h0 = e0->halfedge();
      HalfedgeIter h1 = h0->next();
      HalfedgeIter h2 = h1->next();
      HalfedgeIter h3 = h0->twin();
      HalfedgeIter h4 = h3->next();
      HalfedgeIter h5 = h4->next();

      HalfedgeIter h6 = h1->twin();
      HalfedgeIter h7 = h2->twin();
      HalfedgeIter h8 = h4->twin();
      HalfedgeIter h9 = h5->twin();

      VertexIter v0 = h0->vertex();
      VertexIter v1 = h3->vertex();
      VertexIter v2 = h2->vertex();
      VertexIter v3 = h5->vertex();

      EdgeIter e1 = h1->edge();
      EdgeIter e2 = h2->edge();
      EdgeIter e3 = h4->edge();
      EdgeIter e4 = h5->edge();

      FaceIter f1 = h0->face();
      FaceIter f2 = h3->face();

      //New things
      HalfedgeIter h10 = newHalfedge();
      HalfedgeIter h11 = newHalfedge();
      HalfedgeIter h12 = newHalfedge();
      HalfedgeIter h13 = newHalfedge();
      HalfedgeIter h14 = newHalfedge();
      HalfedgeIter h15 = newHalfedge();

      VertexIter v = newVertex();
      
      EdgeIter e5 = newEdge();
      EdgeIter e6 = newEdge();
      EdgeIter e7 = newEdge();

      FaceIter f3 = newFace();
      FaceIter f4 = newFace();

      h0->setNeighbors(h1,h3,v,e0,f1);
      h1->setNeighbors(h2,h6,v1,e1,f1);
      h2->setNeighbors(h0,h11,v2,e5,f1);
      h3->setNeighbors(h4,h0,v1,e0,f2);
      h4->setNeighbors(h5,h15,v,e7,f2);
      h5->setNeighbors(h3,h9,v3,e4,f2);
      h6->setNeighbors(h6->next(),h1,v2,e1,h6->face());
      h7->setNeighbors(h7->next(),h12,v0,e2,h7->face());
      h8->setNeighbors(h8->next(),h14,v3,e3,h8->face());
      h9->setNeighbors(h9->next(),h5,v1,e4,h9->face());
      h10->setNeighbors(h11,h13,v0,e6,f3);
      h11->setNeighbors(h12,h2,v,e5,f3);
      h12->setNeighbors(h10,h7,v2,e2,f3);
      h13->setNeighbors(h14,h10,v,e6,f4);
      h14->setNeighbors(h15,h8,v0,e3,f4);
      h15->setNeighbors(h13,h4,v3,e7,f4);

      v->position = 0.5 * (v0->position + v1->position);
      v->isNew = 1;
      v0->halfedge() = h10;
      v1->halfedge() = h1;
      v2->halfedge() = h12;
      v3->halfedge() = h5;
      v->halfedge() = h0;
      e0->halfedge() = h0;
      e1->halfedge() = h1;
      e2->halfedge() = h12;
      e3->halfedge() = h14;
      e4->halfedge() = h5;
      e5->halfedge() = h2;
      e6->halfedge() = h10;
      e7->halfedge() = h4;
      e0->isNew = 0;
      e6->isNew = 0;
      e5->isNew = 1;
      e7->isNew = 1;
      f1->halfedge() = h0;
      f2->halfedge() = h3;
      f3->halfedge() = h10;
      f4->halfedge() = h13;
      return v;
    }
    else if (e0->isBoundary())
    {
      HalfedgeIter h0;
      if (e0->halfedge()->isBoundary() && (!e0->halfedge()->twin()->isBoundary()))
      {
        h0 = e0->halfedge();
      }
      else if ((!e0->halfedge()->isBoundary()) && (e0->halfedge()->twin()->isBoundary()))
      {
        h0 = e0->halfedge()->twin();
      }
      else 
      {
        return VertexIter();
      }
      
      HalfedgeIter h1 = h0->next();
      HalfedgeIter h2 = h1->next();
      HalfedgeIter h3 = h0->twin();
      HalfedgeIter h4 = h3->next();
      HalfedgeIter h5 = h4->next();

      HalfedgeIter h6 = h1->twin();
      HalfedgeIter h7 = h2->twin();
      HalfedgeIter h8 = h4->twin();
      HalfedgeIter h9 = h5->twin();

      VertexIter v0 = h0->vertex();
      VertexIter v1 = h3->vertex();
      VertexIter v2 = h2->vertex();
      VertexIter v3 = h5->vertex();

      EdgeIter e1 = h1->edge();
      EdgeIter e2 = h2->edge();
      EdgeIter e3 = h4->edge();
      EdgeIter e4 = h5->edge();

      FaceIter f1 = h0->face();
      FaceIter f2 = h3->face();
      //New things
      HalfedgeIter h10 = newHalfedge();
      HalfedgeIter h11 = newHalfedge();
      HalfedgeIter h12 = newHalfedge();
      HalfedgeIter h13 = newHalfedge();

      VertexIter v = newVertex();
      
      EdgeIter e5 = newEdge();
      EdgeIter e6 = newEdge();

      FaceIter f3 = newFace();

      h0->setNeighbors(h1,h3,v,e0,f1);
      h1->setNeighbors(h2,h6,v1,e1,f1);
      h2->setNeighbors(h10,h7,v0,e2,f1);
      h3->setNeighbors(h4,h0,v1,e0,f2);
      h4->setNeighbors(h5,h13,v,e6,f2);
      h5->setNeighbors(h3,h9,v3,e4,f2);
      h6->setNeighbors(h6->next(),h1,v2,e1,h6->face());
      h7->setNeighbors(h7->next(),h2,v0,e2,h7->face());
      h8->setNeighbors(h8->next(),h12,v3,e3,h8->face());
      h9->setNeighbors(h9->next(),h5,v1,e4,h9->face());
      h10->setNeighbors(h0,h11,v0,e5,f1);
      h11->setNeighbors(h12,h10,v,e5,f3);
      h12->setNeighbors(h13,h8,v0,e3,f3);
      h13->setNeighbors(h11,h4,v3,e6,f3);

      v->position = 0.5 * (v0->position + v1->position);
      v->isNew = 1;
      v0->halfedge() = h10;
      v1->halfedge() = h1;
      v2->halfedge() = h2;
      v3->halfedge() = h5;
      v->halfedge() = h0;

      e0->halfedge() = h0;
      e1->halfedge() = h1;
      e2->halfedge() = h2;
      e3->halfedge() = h12;
      e4->halfedge() = h5;
      e5->halfedge() = h10;
      e6->halfedge() = h4;
      
      e0->isNew = 0;
      e5->isNew = 0;
      e6->isNew = 1;

      f1->halfedge() = h0;
      f2->halfedge() = h3;
      f3->halfedge() = h11;

      return v;
    }
    

    return VertexIter();
  }
  



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Task 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    
    // 4. Flip any new edge that connects an old and new vertex.

    // 5. Copy the new vertex positions into final Vertex::position.

    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd();e++) {
      e->isNew = 0;
      HalfedgeIter AB = e->halfedge();
      HalfedgeIter BC = AB->next();
      HalfedgeIter CA = BC->next();
      HalfedgeIter BA = AB->twin();
      HalfedgeIter AD = BA->next();
      HalfedgeIter DB = AD->next();

      VertexIter A = AB->vertex();
      VertexIter B = BA->vertex();
      VertexIter C = CA->vertex();
      VertexIter D = DB->vertex();

      e->newPosition = 0.375 * (A->position + B->position) + 0.125 * (C->position + D->position);
      // don't use 3/8 and 1/8
    }

    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) { 
      v->isNew = 0;
      HalfedgeIter h = v->halfedge();
      Vector3D neighbor_sum(0.0);
      int n = v->degree();
      float u;
      for (int i = 0; i < n; i++)
      {
        neighbor_sum = neighbor_sum + h->twin()->vertex()->position;
        h = h->twin()->next();
      }
      if (n == 3)
      {
        u = 3.0 / 16.0;
      }
      else {
        u = 3.0/ (8.0 * (float)n);
      }
      v->newPosition = (1.0 - n * u) * v->position + u * neighbor_sum;
    }
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
      VertexIter v1 = e->halfedge()->vertex();
      VertexIter v2 = e->halfedge()->twin()->vertex();
      bool need_split = !(v1->isNew||v2->isNew);
      if(need_split) {
        VertexIter v = mesh.splitEdge(e);
        v->newPosition = e->newPosition;
      }
    }
    int time = 0;
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
      VertexIter v1 = e->halfedge()->vertex();
      VertexIter v2 = e->halfedge()->twin()->vertex();
      bool one_new_vertex;
      one_new_vertex = (v1->isNew && !(v2->isNew))||(v2->isNew && !(v1->isNew));
      bool need_flip = e->isNew && one_new_vertex;
      if(need_flip) {
        mesh.flipEdge(e);
        time = time + 1;
      }
    }
    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
      v->position = v->newPosition;
      v->isNew = 0;
    }
  }
}
