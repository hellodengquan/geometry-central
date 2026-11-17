// Tests for the surface heat-method distance solver.
//
// These tests pin down the *meaning* of the returned scalar field rather than
// exact geodesic distances. The heat method is a discretization: it approximates
// geodesic distance up to an error that shrinks under mesh refinement. Numerical
// criteria below are therefore expressed relative to the mesh diameter, and are
// chosen to separate ordinary discretization error from genuine mapping bugs
// (wrong source, wrong datum, NaN, information leaking between components).
//
// Some tests document boundaries where the current implementation offers no
// guarantee (empty source sets, zero-area faces, the robust Laplacian on meshes
// with boundary). Those assert the documented behaviour explicitly rather than
// silently accepting a numerical anomaly as a distance.
//
// All meshes are built procedurally; no external mesh assets are used.

#include "geometrycentral/surface/heat_method_distance.h"

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace {

using MeshPtr = std::unique_ptr<ManifoldSurfaceMesh>;
using GeomPtr = std::unique_ptr<VertexPositionGeometry>;
using MeshGeom = std::tuple<MeshPtr, GeomPtr>;

// A self-contained connected piece used to assemble multi-component meshes.
struct ComponentSpec {
  std::vector<Vector3> verts;
  std::vector<std::vector<size_t>> faces;
};

// Concatenate several disconnected components into one mesh, re-indexing faces.
// Components are appended in order, so tests can rely on the resulting vertex
// index ranges.
MeshGeom makeUnion(const std::vector<ComponentSpec>& comps) {
  std::vector<Vector3> verts;
  std::vector<std::vector<size_t>> faces;
  for (const ComponentSpec& c : comps) {
    size_t offset = verts.size();
    for (const Vector3& p : c.verts) verts.push_back(p);
    for (const std::vector<size_t>& f : c.faces) {
      std::vector<size_t> g;
      for (size_t i : f) g.push_back(offset + i);
      faces.push_back(g);
    }
  }
  return makeManifoldSurfaceMeshAndGeometry(faces, verts);
}

ComponentSpec unitRightTriangle(double x0 = 0., double y0 = 0.) {
  ComponentSpec c;
  c.verts = {{x0, y0, 0.}, {x0 + 1., y0, 0.}, {x0, y0 + 1., 0.}};
  c.faces = {{0, 1, 2}};
  return c;
}

// A flat w x h rectangle (two triangles) placed at (x0, y0).
ComponentSpec rectangle(double x0, double y0, double w, double h) {
  ComponentSpec c;
  c.verts = {{x0, y0, 0.}, {x0 + w, y0, 0.}, {x0 + w, y0 + h, 0.}, {x0, y0 + h, 0.}};
  c.faces = {{0, 1, 2}, {0, 2, 3}};
  return c;
}

// Flat (n+1) x (n+1) vertex grid over [0,n]x[0,n]; vertex (i,j) has index
// j*(n+1)+i. Each unit cell is split into two triangles.
MeshGeom makeUnitGrid(int n) {
  std::vector<Vector3> verts;
  for (int j = 0; j <= n; j++) {
    for (int i = 0; i <= n; i++) verts.push_back(Vector3{(double)i, (double)j, 0.});
  }
  auto idx = [&](int i, int j) { return (size_t)(j * (n + 1) + i); };
  std::vector<std::vector<size_t>> faces;
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      faces.push_back({idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)});
      faces.push_back({idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)});
    }
  }
  return makeManifoldSurfaceMeshAndGeometry(faces, verts);
}

size_t gridIndex(int n, int i, int j) { return (size_t)(j * (n + 1) + i); }

// Two triangular "rulers" glued along a unit-width strip. The same connectivity
// is embedded either flat (isometric to the plane) or as an accordion that folds
// in 3D every segment. Both embeddings are intrinsically identical, so any
// intrinsic distance solver must return the same field on both.
MeshGeom makeAccordionStrip(int n, bool folded) {
  std::vector<Vector3> spine;
  if (!folded) {
    for (int i = 0; i <= n; i++) spine.push_back(Vector3{(double)i, 0., 0.});
  } else {
    // Unit steps alternating +x, +z, -x, -z: intrinsic length n, tiny ambient extent.
    std::vector<Vector3> dirs = {{1., 0., 0.}, {0., 0., 1.}, {-1., 0., 0.}, {0., 0., -1.}};
    spine.push_back(Vector3{0., 0., 0.});
    for (int i = 1; i <= n; i++) spine.push_back(spine[i - 1] + dirs[(i - 1) % 4]);
  }
  std::vector<Vector3> verts;
  for (int j = 0; j <= 1; j++) {
    for (int i = 0; i <= n; i++) verts.push_back(spine[i] + Vector3{0., (double)j, 0.});
  }
  auto idx = [&](int i, int j) { return (size_t)(j * (n + 1) + i); };
  std::vector<std::vector<size_t>> faces;
  for (int i = 0; i < n; i++) {
    faces.push_back({idx(i, 0), idx(i + 1, 0), idx(i + 1, 1)});
    faces.push_back({idx(i, 0), idx(i + 1, 1), idx(i, 1)});
  }
  return makeManifoldSurfaceMeshAndGeometry(faces, verts);
}

// Closed octahedron; no boundary, so the robust-Laplacian path is usable.
MeshGeom makeOctahedron() {
  std::vector<Vector3> verts = {{1., 0., 0.}, {-1., 0., 0.}, {0., 1., 0.},
                                {0., -1., 0.}, {0., 0., 1.}, {0., 0., -1.}};
  std::vector<std::vector<size_t>> faces = {{0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4},
                                            {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}};
  return makeManifoldSurfaceMeshAndGeometry(faces, verts);
}

// A single zero-area (collinear) triangle.
MeshGeom makeCollinearTriangle() {
  ComponentSpec c;
  c.verts = {{0., 0., 0.}, {1., 0., 0.}, {2., 0., 0.}};
  c.faces = {{0, 1, 2}};
  return makeUnion({c});
}

// A quad (two triangles) whose apex is almost on the opposite edge. With
// wellShaped=false the area is ~1e-6; the connectivity is identical in both
// variants so the two meshes are directly comparable.
MeshGeom makeApexQuad(bool wellShaped) {
  ComponentSpec c;
  double apexY = wellShaped ? 1.0 : 1e-6;
  c.verts = {{0., 0., 0.}, {2., 0., 0.}, {1., apexY, 0.}, {1., -1., 0.}};
  c.faces = {{0, 1, 2}, {1, 0, 3}};
  return makeUnion({c});
}

// Recompute the in-face distance between the source and each vertex of its face
// purely from edge lengths (law of cosines -> 2D embedding). This independently
// reconstructs the "expected distance to the face vertices" that the
// implementation uses to anchor the zero of the distance field.
double referenceAnchor(IntrinsicGeometryInterface& geom, const SurfacePoint& p) {
  geom.requireEdgeLengths();
  SurfacePoint fp = p.inSomeFace();
  Halfedge h0 = fp.face.halfedge();
  double L0 = geom.edgeLengths[h0.edge()];               // |v0 v1|
  double L1 = geom.edgeLengths[h0.next().edge()];        // |v1 v2|
  double L2 = geom.edgeLengths[h0.next().next().edge()]; // |v2 v0|

  Vector2 P0{0., 0.};
  Vector2 P1{L0, 0.};
  double x = (L2 * L2 + L0 * L0 - L1 * L1) / (2.0 * L0);
  double y = std::sqrt(std::max(0.0, L2 * L2 - x * x));
  Vector2 P2{x, y};

  Vector2 ps = fp.faceCoords.x * P0 + fp.faceCoords.y * P1 + fp.faceCoords.z * P2;

  double anchor = 0.;
  anchor += fp.faceCoords.x * (P0 - ps).norm();
  anchor += fp.faceCoords.y * (P1 - ps).norm();
  anchor += fp.faceCoords.z * (P2 - ps).norm();
  return anchor;
}

// Build the initial condition the solver uses for point sources: a barycentric
// "delta" spread over the vertices of one incident face, summed over all sources.
VertexData<double> sourceRHS(SurfaceMesh& mesh, const std::vector<SurfacePoint>& sources) {
  VertexData<double> rhs(mesh, 0.);
  for (const SurfacePoint& p : sources) {
    SurfacePoint fp = p.inSomeFace();
    Halfedge he = fp.face.halfedge();
    rhs[he.vertex()] += fp.faceCoords.x;
    rhs[he.next().vertex()] += fp.faceCoords.y;
    rhs[he.next().next().vertex()] += fp.faceCoords.z;
  }
  return rhs;
}

double maxAbsDiff(const VertexData<double>& a, const VertexData<double>& b) {
  double worst = 0.;
  for (Vertex v : a.getMesh()->vertices()) worst = std::max(worst, std::abs(a[v] - b[v]));
  return worst;
}

// Compare two fields that live on two *different* mesh objects with matching
// vertex indices (as produced by rebuilding the same fixture). Indexing one
// field with the other's Vertex handles would be invalid: element handles are
// tied to their mesh.
double maxAbsDiffByIndex(const VertexData<double>& a, const VertexData<double>& b) {
  SurfaceMesh& meshA = *a.getMesh();
  SurfaceMesh& meshB = *b.getMesh();
  double worst = 0.;
  for (Vertex v : meshA.vertices()) {
    Vertex w = meshB.vertex(v.getIndex());
    worst = std::max(worst, std::abs(a[v] - b[w]));
  }
  return worst;
}

double fieldRange(const VertexData<double>& d, const std::vector<Vertex>& verts) {
  double lo = std::numeric_limits<double>::max();
  double hi = -std::numeric_limits<double>::max();
  for (Vertex v : verts) {
    lo = std::min(lo, d[v]);
    hi = std::max(hi, d[v]);
  }
  return hi - lo;
}

double meshDiameter(SurfaceMesh& mesh, VertexPositionGeometry& geom) {
  double D = 0.;
  for (Vertex a : mesh.vertices()) {
    for (Vertex b : mesh.vertices()) {
      D = std::max(D, (geom.vertexPositions[a] - geom.vertexPositions[b]).norm());
    }
  }
  return D;
}

std::string fieldToString(const VertexData<double>& d) {
  std::string s;
  for (Vertex v : d.getMesh()->vertices()) s += std::to_string(d[v]) + " ";
  return s;
}

// Pick an interior edge incident on v, so that both adjacent faces exist and
// SurfacePoint::inFace works for either side.
Edge interiorEdgeAt(SurfaceMesh& mesh, Vertex v) {
  for (Edge e : v.adjacentEdges()) {
    if (!e.isBoundary()) return e;
  }
  return Edge();
}

} // namespace

// ============================================================
// =============== Source representation and the zero point
// ============================================================

// A vertex source must give zero exactly at that vertex, a non-negative field,
// and finite values everywhere: the zero set of the field is the source.
TEST(HeatMethodDistanceSuite, VertexSourceIsZeroAtSourceAndNonNegative) {
  {
    MeshGeom mg = makeUnitGrid(10);
    ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
    HeatMethodDistanceSolver solver(*std::get<1>(mg));

    for (size_t k = 0; k < 12; k++) {
      Vertex src = mesh.vertex((k * 13) % mesh.nVertices());
      VertexData<double> d = solver.computeDistance(src);

      EXPECT_TRUE(d.toVector().allFinite()) << "non-finite value for source " << k;
      EXPECT_NEAR(d[src], 0., 1e-12) << "field is not zero at vertex source " << k;
      EXPECT_GE(d.toVector().minCoeff(), -1e-12) << "negative distance for source " << k;
    }
  }

  {
    // Closed surface: no boundary vertices to distract the zero set.
    MeshGeom mg = makeOctahedron();
    ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
    HeatMethodDistanceSolver solver(*std::get<1>(mg));

    VertexData<double> d = solver.computeDistance(mesh.vertex(0));
    EXPECT_TRUE(d.toVector().allFinite());
    EXPECT_NEAR(d[mesh.vertex(0)], 0., 1e-12);
    EXPECT_GE(d.toVector().minCoeff(), -1e-12);
  }
}

// The three entry points (vertex, vertex list, surface point) describe the same
// source and must return the same field.
TEST(HeatMethodDistanceSuite, SourceOverloadsAgree) {
  MeshGeom mg = makeUnitGrid(6);
  ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
  HeatMethodDistanceSolver solver(*std::get<1>(mg));

  Vertex src = mesh.vertex(14);
  VertexData<double> dv = solver.computeDistance(src);
  VertexData<double> dl = solver.computeDistance(std::vector<Vertex>{src});
  VertexData<double> dp = solver.computeDistance(SurfacePoint(src));

  EXPECT_LT(maxAbsDiff(dv, dl), 1e-12);
  EXPECT_LT(maxAbsDiff(dv, dp), 1e-12);
}

// A point in the relative interior of an edge is not a vertex, so no vertex can
// sit at the zero of the field. The solver instead anchors the field so that the
// value interpolated *at the source point* equals the in-face distance to its
// face's vertices; for an edge point that equals the exact along-edge distance
// 2 t (1 - t) L. This test pins that anchoring convention.
TEST(HeatMethodDistanceSuite, EdgeSourceAnchorMatchesInFaceDistance) {
  MeshGeom mg = makeUnitGrid(8);
  ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
  VertexPositionGeometry& geom = *std::get<1>(mg);
  geom.requireEdgeLengths();
  HeatMethodDistanceSolver solver(geom);

  Edge e = interiorEdgeAt(mesh, mesh.vertex(gridIndex(8, 4, 4)));
  ASSERT_TRUE(e != Edge());
  double L = geom.edgeLengths[e];

  for (double t : {0.0, 0.15, 0.3, 0.5, 0.72, 1.0}) {
    SurfacePoint p(e, t);
    VertexData<double> d = solver.computeDistance(p);

    double valueAtSource = p.interpolate(d);
    EXPECT_NEAR(valueAtSource, referenceAnchor(geom, p), 1e-10)
        << "anchor mismatch at t=" << t << " (field: " << fieldToString(d) << ")";
    EXPECT_NEAR(valueAtSource, 2.0 * t * (1.0 - t) * L, 1e-10)
        << "edge source is not anchored at the along-edge distance, t=" << t;

    if (t > 1e-6 && t < 1.0 - 1e-6) {
      EXPECT_GT(d.toVector().minCoeff(), 1e-6)
          << "edge interior source collapsed onto a vertex, t=" << t;
    }
  }
}

// Same anchoring rule for a source in the interior of a face: a positive offset
// of order the edge length, and no vertex at the zero of the field.
TEST(HeatMethodDistanceSuite, FaceSourceAnchorMatchesInFaceDistance) {
  MeshGeom mg = makeUnitGrid(8);
  ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
  VertexPositionGeometry& geom = *std::get<1>(mg);
  geom.requireEdgeLengths();
  HeatMethodDistanceSolver solver(geom);

  for (Face f : {mesh.face(0), mesh.face(70), mesh.face(100)}) {
    for (Vector3 b : {Vector3{1. / 3., 1. / 3., 1. / 3.}, Vector3{0.55, 0.3, 0.15},
                      Vector3{0.1, 0.1, 0.8}}) {
      SurfacePoint p(f, b);
      VertexData<double> d = solver.computeDistance(p);

      EXPECT_NEAR(p.interpolate(d), referenceAnchor(geom, p), 1e-10)
          << "face-source anchor mismatch, face " << f.getIndex() << ", bary " << b;
      EXPECT_GT(p.interpolate(d), 1e-6)
          << "interior face source should not be a zero of the field";
      EXPECT_GE(d.toVector().minCoeff(), -1e-12);
    }
  }
}

// The same geometric point can be written as a vertex, as an edge point (in
// either orientation), or as a face point in any adjacent face. Every legal
// representation must produce the same field; otherwise the "source location"
// is being interpreted inconsistently.
TEST(HeatMethodDistanceSuite, EquivalentRepresentationsProduceSameField) {
  MeshGeom mg = makeUnitGrid(8);
  ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
  VertexPositionGeometry& geom = *std::get<1>(mg);
  geom.requireEdgeLengths();
  HeatMethodDistanceSolver solver(geom);

  // --- vertex, directly and as a face point in every incident face.
  Vertex v = mesh.vertex(gridIndex(8, 4, 3));
  VertexData<double> dv = solver.computeDistance(v);

  std::vector<SurfacePoint> vertexReps = {SurfacePoint(v).inSomeFace()};
  for (Face f : v.adjacentFaces()) vertexReps.push_back(SurfacePoint(v).inFace(f));
  for (const SurfacePoint& rep : vertexReps) {
    EXPECT_LT(maxAbsDiff(dv, solver.computeDistance(rep)), 1e-12)
        << "vertex representation mismatch: " << std::to_string(rep);
  }

  // --- edge point: canonical edge form, halfedge form, reversed halfedge form,
  //     and the equivalent face point in both adjacent faces.
  Edge e = interiorEdgeAt(mesh, v);
  ASSERT_TRUE(e != Edge());
  Halfedge he = e.halfedge();
  for (double t : {0.0, 0.25, 0.5, 0.75, 1.0}) {
    SurfacePoint base(e, t);
    VertexData<double> dBase = solver.computeDistance(base);

    std::vector<SurfacePoint> reps;
    reps.push_back(SurfacePoint(he, t));
    reps.push_back(SurfacePoint(he.twin(), 1.0 - t));
    reps.push_back(base.inFace(he.face()));
    reps.push_back(base.inFace(he.twin().face()));
    for (const SurfacePoint& rep : reps) {
      EXPECT_LT(maxAbsDiff(dBase, solver.computeDistance(rep)), 1e-12)
          << "edge representation mismatch at t=" << t << ": " << std::to_string(rep);
    }
  }

  // --- an edge endpoint is a vertex: edge(t=0) is the tail-vertex source.
  {
    SurfacePoint tailPoint(e, 0.0);
    Vertex tail = he.tailVertex();
    EXPECT_LT(maxAbsDiff(solver.computeDistance(tailPoint), solver.computeDistance(tail)), 1e-12)
        << "edge endpoint t=0 is not consistent with the vertex source";
  }
}

// ============================================================
// =============== Multiple sources
// ============================================================

// With several sources the solver applies one global change of origin. The
// convention is that the mean of (field value at each source minus its anchor)
// is zero; this reduces to the single-source convention for one source. It must
// hold for vertex, edge, and face sources mixed together.
TEST(HeatMethodDistanceSuite, MultiSourceDatumIsMeanOfAnchors) {
  MeshGeom mg = makeUnitGrid(10);
  ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
  VertexPositionGeometry& geom = *std::get<1>(mg);
  geom.requireEdgeLengths();
  HeatMethodDistanceSolver solver(geom);

  // Vertex sources only: each anchor is zero, so the mean of the field at the
  // sources is zero.
  std::vector<Vertex> vSrcs = {mesh.vertex(gridIndex(10, 1, 1)), mesh.vertex(gridIndex(10, 9, 2)),
                               mesh.vertex(gridIndex(10, 4, 9))};
  VertexData<double> dv = solver.computeDistance(vSrcs);

  double vertexSum = 0.;
  for (Vertex v : vSrcs) vertexSum += dv[v];
  EXPECT_NEAR(vertexSum / (double)vSrcs.size(), 0., 1e-10)
      << "mean of the field at vertex sources is not zero";

  // Mixed source types: the same mean-of-(value - anchor) rule must hold.
  Edge e = interiorEdgeAt(mesh, mesh.vertex(gridIndex(10, 5, 5)));
  ASSERT_TRUE(e != Edge());
  std::vector<SurfacePoint> mixed = {SurfacePoint(mesh.vertex(gridIndex(10, 1, 1))), SurfacePoint(e, 0.35),
                                     SurfacePoint(mesh.face(40), Vector3{0.5, 0.3, 0.2})};
  VertexData<double> dm = solver.computeDistance(mixed);

  double residualSum = 0.;
  for (const SurfacePoint& p : mixed) residualSum += p.interpolate(dm) - referenceAnchor(geom, p);
  EXPECT_NEAR(residualSum / (double)mixed.size(), 0., 1e-10)
      << "mixed-source datum is not the mean of the source anchors (values: " << fieldToString(dm) << ")";
}

// The datum does not depend on source ordering, but it is a *weighted* mean: the
// same source listed twice changes the weights. Duplicating every source equally
// is a no-op.
TEST(HeatMethodDistanceSuite, MultiSourceOrderAndDuplicateWeighting) {
  MeshGeom mg = makeUnitGrid(8);
  ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
  HeatMethodDistanceSolver solver(*std::get<1>(mg));

  Vertex A = mesh.vertex(gridIndex(8, 2, 2));
  Vertex B = mesh.vertex(gridIndex(8, 6, 5));

  VertexData<double> dAB = solver.computeDistance(std::vector<SurfacePoint>{SurfacePoint(A), SurfacePoint(B)});
  VertexData<double> dBA = solver.computeDistance(std::vector<SurfacePoint>{SurfacePoint(B), SurfacePoint(A)});
  EXPECT_LT(maxAbsDiff(dAB, dBA), 1e-12) << "source ordering changed the field";

  // Duplicating one source re-weights the datum: (2A+B)/3 != (A+B)/2.
  VertexData<double> dAAB =
      solver.computeDistance(std::vector<SurfacePoint>{SurfacePoint(A), SurfacePoint(A), SurfacePoint(B)});
  EXPECT_GT(maxAbsDiff(dAB, dAAB), 1e-6) << "duplicate source was silently ignored";
  double mean = (dAAB[A] + dAAB[A] + dAAB[B]) / 3.0;
  EXPECT_NEAR(mean, 0., 1e-10) << "multiset datum is not the mean over listed sources";

  // Duplicating *all* sources equally leaves the field unchanged (equal weights).
  VertexData<double> dAABB = solver.computeDistance(
      std::vector<SurfacePoint>{SurfacePoint(A), SurfacePoint(B), SurfacePoint(A), SurfacePoint(B)});
  EXPECT_LT(maxAbsDiff(dAB, dAABB), 1e-12);

  // Repeating a single source any number of times is still the single-source solve.
  VertexData<double> dA = solver.computeDistance(A);
  VertexData<double> dAAA =
      solver.computeDistance(std::vector<SurfacePoint>{SurfacePoint(A), SurfacePoint(A), SurfacePoint(A)});
  EXPECT_LT(maxAbsDiff(dA, dAAA), 1e-12);
}

// The multi-source field should look like the lower envelope of the individual
// single-source fields, up to one additive constant and discretization error.
// This checks that all sources share a single distance datum rather than each
// contributing an independently shifted field.
TEST(HeatMethodDistanceSuite, MultiSourceApproximatesLowerEnvelopeOfSingles) {
  MeshGeom mg = makeUnitGrid(14);
  ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
  VertexPositionGeometry& geom = *std::get<1>(mg);
  HeatMethodDistanceSolver solver(geom);

  double D = meshDiameter(mesh, geom);

  std::vector<Vertex> srcs = {mesh.vertex(gridIndex(14, 2, 3)), mesh.vertex(gridIndex(14, 11, 4)),
                              mesh.vertex(gridIndex(14, 6, 12))};
  std::vector<SurfacePoint> points;
  for (Vertex v : srcs) points.push_back(SurfacePoint(v));
  VertexData<double> multi = solver.computeDistance(points);

  // The field is minimized on the source set.
  Vertex argmin = *mesh.vertices().begin();
  for (Vertex v : mesh.vertices()) {
    if (multi[v] < multi[argmin]) argmin = v;
  }
  bool argminOnSource = false;
  for (Vertex v : srcs) argminOnSource = argminOnSource || (v == argmin);
  EXPECT_TRUE(argminOnSource) << "multi-source field is not minimized at any source";

  // Align each single-source field to zero at its own source, take the envelope,
  // and compare to the multi-source field after removing one best-fit constant.
  std::vector<VertexData<double>> singles;
  for (Vertex v : srcs) singles.push_back(solver.computeDistance(v));

  auto envelopeAt = [&](Vertex v) {
    double envelope = std::numeric_limits<double>::max();
    for (size_t i = 0; i < srcs.size(); i++) {
      envelope = std::min(envelope, singles[i][v] - singles[i][srcs[i]]);
    }
    return envelope;
  };

  double meanOffset = 0.;
  for (Vertex v : mesh.vertices()) meanOffset += multi[v] - envelopeAt(v);
  meanOffset /= (double)mesh.nVertices();

  double worstDeviation = 0.;
  for (Vertex v : mesh.vertices()) {
    worstDeviation = std::max(worstDeviation, std::abs((multi[v] - meanOffset) - envelopeAt(v)));
  }

  // The tolerance covers discretization error, not a per-source datum mismatch
  // (which would be O(distance between sources) ~ D).
  EXPECT_LT(worstDeviation / D, 0.08)
      << "multi-source field deviates from the single-source envelope by " << worstDeviation / D
      << " of the diameter";
}

// ============================================================
// =============== Disconnected components
// ============================================================

// On a mesh with several connected components, a component without a source
// receives no information: the raw (unshifted) solution is exactly zero there,
// and after the global change of origin its reported field is a single constant
// equal to that shift. A constant is not a distance; this test records that
// clearly instead of treating the value as a finite geodesic distance.
TEST(HeatMethodDistanceSuite, SourceFreeComponentIsAUniformDatumOffset) {
  MeshGeom mg = makeUnion({unitRightTriangle(), rectangle(50., 0., 4., 3.)});
  ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
  HeatMethodDistanceSolver solver(*std::get<1>(mg));

  ASSERT_EQ(mesh.nConnectedComponents(), 2);

  std::vector<Vertex> compA = {mesh.vertex(0), mesh.vertex(1), mesh.vertex(2)};
  std::vector<Vertex> compB = {mesh.vertex(3), mesh.vertex(4), mesh.vertex(5), mesh.vertex(6)};

  std::vector<SurfacePoint> sources = {SurfacePoint(compA[0])};
  VertexData<double> d = solver.computeDistance(sources);

  // The raw solution never reaches the source-free component.
  VertexData<double> raw(mesh, solver.computeDistanceRHS(sourceRHS(mesh, sources).toVector()));
  for (Vertex v : compB) {
    EXPECT_NEAR(raw[v], 0., 1e-12) << "information leaked into a source-free component at vertex " << v;
  }

  // The reported field is constant on the source-free component...
  EXPECT_LT(fieldRange(d, compB), 1e-12) << "source-free component did not receive a uniform value";
  // ...and that constant is precisely the global datum shift applied to the
  // sourced component: an arbitrary origin, not a distance.
  double shift = d[compA[0]] - raw[compA[0]];
  EXPECT_NEAR(d[compB[0]], shift, 1e-10) << "source-free constant is not the datum shift";
  EXPECT_GT(std::abs(shift), 1e-3) << "expected a non-zero datum offset in this fixture";

  // Choosing a different source in the sourced component moves the free constant,
  // confirming that it tracks the datum rather than a surface distance.
  VertexData<double> dAlt = solver.computeDistance(std::vector<SurfacePoint>{SurfacePoint(compA[1])});
  EXPECT_GT(std::abs(dAlt[compB[0]] - d[compB[0]]), 1e-3)
      << "free-component value did not follow the datum, which would suggest it is distance-like";

  // The sourced component still looks like a distance field.
  EXPECT_NEAR(d[compA[0]], 0., 1e-12);
  EXPECT_GE(d.toVector().minCoeff(), -1e-12);
}

// A source-free component's offset must not depend on where that component sits
// in space. Rigid motion changes no intrinsic quantity, so the whole field must
// be unchanged; and the field inside the free component must be exactly constant
// regardless of its shape.
TEST(HeatMethodDistanceSuite, SourceFreeComponentOffsetIgnoresItsOwnPlacement) {
  MeshGeom mgOrig = makeUnion({unitRightTriangle(), rectangle(50., 0., 4., 3.)});
  MeshGeom mgTrans = makeUnion({unitRightTriangle(), rectangle(-800., 250., 4., 3.)});
  HeatMethodDistanceSolver solverOrig(*std::get<1>(mgOrig));
  HeatMethodDistanceSolver solverTrans(*std::get<1>(mgTrans));

  VertexData<double> dOrig = solverOrig.computeDistance(std::get<0>(mgOrig)->vertex(0));
  VertexData<double> dTrans = solverTrans.computeDistance(std::get<0>(mgTrans)->vertex(0));

  EXPECT_LT(maxAbsDiffByIndex(dOrig, dTrans), 1e-10)
      << "rigidly moving the source-free component changed the distance field";

  // Reshaping the free component (same connectivity, different edge lengths) must
  // not create any variation inside it: it stays one constant.
  MeshGeom mgReshaped = makeUnion({unitRightTriangle(), rectangle(50., 0., 20., 1.)});
  HeatMethodDistanceSolver solverReshaped(*std::get<1>(mgReshaped));
  VertexData<double> dReshaped = solverReshaped.computeDistance(std::get<0>(mgReshaped)->vertex(0));

  std::vector<Vertex> freeVerts = {std::get<0>(mgReshaped)->vertex(3), std::get<0>(mgReshaped)->vertex(4),
                                   std::get<0>(mgReshaped)->vertex(5), std::get<0>(mgReshaped)->vertex(6)};
  EXPECT_LT(fieldRange(dReshaped, freeVerts), 1e-10)
      << "reshaping the source-free component produced a non-uniform field inside it";
}

// Topology is not the only coupling between components. The heat kernel's time
// step is derived from the *global* mean edge length, so changing the geometry of
// a source-free component perturbs the time step and hence the sourced
// component's field by a small amount, even though heat cannot flow between the
// components. This is a discretization artifact, not a geodesic effect: it is
// small, does not depend on the spatial separation, and does not turn the free
// component's constant into a distance. Record the magnitude so that a future
// change that makes it O(distance) is caught.
TEST(HeatMethodDistanceSuite, FreeComponentShapePerturbsGlobalTimeStepOnly) {
  const int n = 8;
  const size_t nGridVerts = (size_t)(n + 1) * (n + 1);

  auto buildWithFreeComponent = [&](double freeW, double freeH, double sep) {
    ComponentSpec freeComp = rectangle(sep, 0., freeW, freeH);
    // Assemble the sourced grid component by hand so its vertex numbering is
    // known and matches gridIndex().
    ComponentSpec grid;
    for (int j = 0; j <= n; j++) {
      for (int i = 0; i <= n; i++) grid.verts.push_back(Vector3{(double)i, (double)j, 0.});
    }
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < n; i++) {
        grid.faces.push_back({gridIndex(n, i, j), gridIndex(n, i + 1, j), gridIndex(n, i + 1, j + 1)});
        grid.faces.push_back({gridIndex(n, i, j), gridIndex(n, i + 1, j + 1), gridIndex(n, i, j + 1)});
      }
    }
    return makeUnion({grid, freeComp});
  };

  // Solve and copy values out by vertex index, since the mesh is owned by the
  // local MeshGeom and dies at the end of each call.
  auto fieldFor = [&](double freeW, double freeH, double sep) {
    MeshGeom mg = buildWithFreeComponent(freeW, freeH, sep);
    ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
    HeatMethodDistanceSolver solver(*std::get<1>(mg));
    VertexData<double> d = solver.computeDistance(mesh.vertex(0));
    std::vector<double> values;
    for (Vertex v : mesh.vertices()) values.push_back(d[v]);
    return values;
  };

  const double D = std::sqrt(2.) * n;
  std::vector<double> baseline = fieldFor(4., 3., 50.);

  // Moving the free component further away changes nothing at all (no coupling).
  {
    std::vector<double> far = fieldFor(4., 3., 5000.);
    ASSERT_EQ(far.size(), baseline.size());
    for (size_t i = 0; i < nGridVerts; i++) {
      EXPECT_LT(std::abs(far[i] - baseline[i]), 1e-12)
          << "spatial separation of the free component altered the sourced field at vertex " << i;
    }
  }

  // Changing the free component's shape perturbs the sourced region only slightly,
  // through the global time step, and leaves the free component's field constant.
  double worstRel = 0.;
  for (auto shape : std::vector<std::pair<double, double>>{{20., 1.}, {0.5, 0.5}, {100., 50.}}) {
    std::vector<double> d = fieldFor(shape.first, shape.second, 50.);
    ASSERT_EQ(d.size(), baseline.size());

    for (size_t i = 0; i < nGridVerts; i++) {
      worstRel = std::max(worstRel, std::abs(d[i] - baseline[i]) / D);
    }

    double lo = std::numeric_limits<double>::max();
    double hi = -std::numeric_limits<double>::max();
    for (size_t i = nGridVerts; i < d.size(); i++) {
      lo = std::min(lo, d[i]);
      hi = std::max(hi, d[i]);
    }
    EXPECT_LT(hi - lo, 1e-10) << "source-free component became non-constant after reshaping";
  }
  // Well below any real distance scale: this is a time-step artifact, not a
  // spurious path between components. A genuine leak would be O(D).
  EXPECT_LT(worstRel, 0.08)
      << "free-component shape changed the sourced field by " << worstRel << " of the diameter";
}


// When every component carries a source, the same global datum rule applies. For
// two congruent components with corresponding sources the two fields are exact
// mirror images; for asymmetric components the values at the sources still have
// zero mean.
TEST(HeatMethodDistanceSuite, EveryComponentSourcedSharesTheGlobalDatum) {
  {
    MeshGeom mg = makeUnion({unitRightTriangle(), unitRightTriangle(50., 0.)});
    ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
    HeatMethodDistanceSolver solver(*std::get<1>(mg));

    std::vector<SurfacePoint> sources = {SurfacePoint(mesh.vertex(0)), SurfacePoint(mesh.vertex(3))};
    VertexData<double> d = solver.computeDistance(sources);

    EXPECT_NEAR(d[mesh.vertex(0)], 0., 1e-12);
    EXPECT_NEAR(d[mesh.vertex(3)], 0., 1e-12);
    EXPECT_NEAR(d[mesh.vertex(1)], d[mesh.vertex(4)], 1e-12);
    EXPECT_NEAR(d[mesh.vertex(2)], d[mesh.vertex(5)], 1e-12);
  }

  {
    // Asymmetric components: the global datum is still the mean over sources.
    MeshGeom mg = makeUnion({unitRightTriangle(), rectangle(20., 0., 4., 3.)});
    ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
    HeatMethodDistanceSolver solver(*std::get<1>(mg));

    std::vector<SurfacePoint> sources = {SurfacePoint(mesh.vertex(0)), SurfacePoint(mesh.vertex(3))};
    VertexData<double> d = solver.computeDistance(sources);
    EXPECT_NEAR(d[mesh.vertex(0)] + d[mesh.vertex(3)], 0., 1e-10)
        << "global datum is not the mean over sources in different components";
  }
}

// An empty source set has no zero point to anchor to. The current implementation
// produces NaN everywhere (0/0 in the datum normalization). Pin that boundary so
// callers cannot mistake the output for a valid distance field.
TEST(HeatMethodDistanceSuite, EmptySourceSetIsNotAValidDistanceField) {
  MeshGeom mg = makeUnitGrid(4);
  ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
  HeatMethodDistanceSolver solver(*std::get<1>(mg));

  VertexData<double> d = solver.computeDistance(std::vector<SurfacePoint>{});
  EXPECT_FALSE(d.toVector().allFinite())
      << "empty source set returned a finite field; if this becomes supported, update this boundary test";
}

// ============================================================
// =============== Intrinsic invariants and refinement
// ============================================================

// The distance is intrinsic: an isometric unfolding of the same strip must give
// the identical field even though the ambient embedding differs wildly. It also
// confirms the answer follows surface length, not straight-line distance through
// space.
TEST(HeatMethodDistanceSuite, DistanceIsIntrinsicNotAmbient) {
  const int n = 14;
  MeshGeom mgFlat = makeAccordionStrip(n, /*folded=*/false);
  MeshGeom mgFold = makeAccordionStrip(n, /*folded=*/true);
  ManifoldSurfaceMesh& meshFlat = *std::get<0>(mgFlat);
  ManifoldSurfaceMesh& meshFold = *std::get<0>(mgFold);
  VertexPositionGeometry& geomFlat = *std::get<1>(mgFlat);
  VertexPositionGeometry& geomFold = *std::get<1>(mgFold);

  // Sanity: the two embeddings really are intrinsically identical.
  geomFlat.requireEdgeLengths();
  geomFold.requireEdgeLengths();
  ASSERT_EQ(meshFlat.nEdges(), meshFold.nEdges());
  for (size_t i = 0; i < meshFlat.nEdges(); i++) {
    ASSERT_NEAR(geomFlat.edgeLengths[meshFlat.edge(i)], geomFold.edgeLengths[meshFold.edge(i)], 1e-12)
        << "fixture is not an isometric pair";
  }

  HeatMethodDistanceSolver solverFlat(geomFlat);
  HeatMethodDistanceSolver solverFold(geomFold);

  VertexData<double> dFlat = solverFlat.computeDistance(meshFlat.vertex(0));
  VertexData<double> dFold = solverFold.computeDistance(meshFold.vertex(0));

  VertexData<double> dFoldReindexed(meshFlat);
  for (Vertex v : meshFlat.vertices()) dFoldReindexed[v] = dFold[meshFold.vertex(v.getIndex())];

  EXPECT_LT(maxAbsDiff(dFlat, dFoldReindexed), 1e-12)
      << "isometric embeddings produced different distance fields";

  // Walking the strip is much longer than the ambient gap, and the field follows
  // the surface: the far end is ~n, not ~sqrt(2).
  size_t farIndex = gridIndex(n, n, 0);
  double farValue = dFold[meshFold.vertex(farIndex)];
  double ambientGap =
      (geomFold.vertexPositions[meshFold.vertex(0)] - geomFold.vertexPositions[meshFold.vertex(farIndex)]).norm();
  EXPECT_GT(farValue, 0.9 * (double)n) << "distance did not follow the surface";
  EXPECT_GT(farValue, 4.0 * ambientGap) << "distance looks ambient rather than intrinsic";
}

// Scaling the surface by a constant scales distances by the same constant; the
// time step is chosen from the mean edge length, so this checks the
// scale-consistency of that choice.
TEST(HeatMethodDistanceSuite, DistancesScaleEquivariantly) {
  const int n = 8;
  MeshGeom mgUnit = makeUnitGrid(n);

  std::vector<Vector3> verts;
  for (int j = 0; j <= n; j++) {
    for (int i = 0; i <= n; i++) verts.push_back(Vector3{10. * i, 10. * j, 0.});
  }
  std::vector<std::vector<size_t>> faces;
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      faces.push_back({gridIndex(n, i, j), gridIndex(n, i + 1, j), gridIndex(n, i + 1, j + 1)});
      faces.push_back({gridIndex(n, i, j), gridIndex(n, i + 1, j + 1), gridIndex(n, i, j + 1)});
    }
  }
  MeshGeom mgScaled = makeManifoldSurfaceMeshAndGeometry(faces, verts);

  HeatMethodDistanceSolver solverUnit(*std::get<1>(mgUnit));
  HeatMethodDistanceSolver solverScaled(*std::get<1>(mgScaled));
  VertexData<double> dUnit = solverUnit.computeDistance(std::get<0>(mgUnit)->vertex(0));
  VertexData<double> dScaled = solverScaled.computeDistance(std::get<0>(mgScaled)->vertex(0));

  double worstRel = 0.;
  for (Vertex v : std::get<0>(mgUnit)->vertices()) {
    double scaled = 10.0 * dUnit[v];
    double got = dScaled[std::get<0>(mgScaled)->vertex(v.getIndex())];
    worstRel = std::max(worstRel, std::abs(got - scaled) / std::max(1e-12, std::abs(scaled)));
  }
  EXPECT_LT(worstRel, 1e-8) << "distance field is not scale-equivariant";
}

// On a flat, uniformly refined grid the heat method converges to the planar
// Euclidean distance. This checks the *approximation*, with tolerances reflecting
// discretization error that improve under refinement. It intentionally does not
// demand pointwise equality with exact geodesics.
TEST(HeatMethodDistanceSuite, FlatGridApproximationImprovesWithRefinement) {
  double firstRelError = -1.;
  double lastRelError = -1.;

  for (int n : {12, 24}) {
    MeshGeom mg = makeUnitGrid(n);
    ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
    VertexPositionGeometry& geom = *std::get<1>(mg);
    HeatMethodDistanceSolver solver(geom);

    Vertex src = mesh.vertex(gridIndex(n, n / 4, n / 3));
    VertexData<double> d = solver.computeDistance(src);
    Vector3 srcPos = geom.vertexPositions[src];

    double D = meshDiameter(mesh, geom);
    double maxError = 0.;
    for (Vertex v : mesh.vertices()) {
      double exact = (geom.vertexPositions[v] - srcPos).norm();
      maxError = std::max(maxError, std::abs(d[v] - exact));
    }
    double relError = maxError / D;

    EXPECT_LT(relError, 0.08) << "discretization error too large on a refined grid, n=" << n;
    if (firstRelError < 0.) firstRelError = relError;
    lastRelError = relError;
  }

  // Refinement should visibly reduce the error (here by roughly a factor of two).
  EXPECT_LT(lastRelError, 0.7 * firstRelError)
      << "error did not decrease under refinement: " << firstRelError << " -> " << lastRelError;
}

// ============================================================
// =============== Raw right-hand-side contract
// ============================================================

// `computeDistanceRHS` skips the change of origin: its result differs from
// `computeDistance` by a single additive constant. It is also invariant under
// positive scaling of the right-hand side, because the per-face gradient is
// normalized before integration. (It is not additive in the right-hand side.)
TEST(HeatMethodDistanceSuite, RawRHSSolutionIsDeterminedUpToAConstant) {
  MeshGeom mg = makeUnitGrid(8);
  ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
  HeatMethodDistanceSolver solver(*std::get<1>(mg));

  std::vector<SurfacePoint> sources = {SurfacePoint(mesh.vertex(3))};
  VertexData<double> rhs = sourceRHS(mesh, sources);

  VertexData<double> raw(mesh, solver.computeDistanceRHS(rhs.toVector()));
  VertexData<double> shifted = solver.computeDistance(sources);

  double shift = shifted[mesh.vertex(3)] - raw[mesh.vertex(3)];
  double worstDeviation = 0.;
  for (Vertex v : mesh.vertices()) {
    worstDeviation = std::max(worstDeviation, std::abs((shifted[v] - shift) - raw[v]));
  }
  EXPECT_LT(worstDeviation, 1e-10) << "raw solution is not the shifted solution minus one constant";

  // Positive scaling of the right-hand side leaves the raw solution unchanged.
  VertexData<double> rawScaled(mesh, solver.computeDistanceRHS(5.0 * rhs.toVector()));
  double scale = 0.;
  for (Vertex v : mesh.vertices()) scale = std::max(scale, std::abs(raw[v]));
  EXPECT_LT(maxAbsDiff(raw, rawScaled), 1e-8 * std::max(1.0, scale))
      << "raw solution depends on the magnitude of the right-hand side";
}

// ============================================================
// =============== Robust Laplacian and degenerate geometry
// ============================================================

// The robust (tufted/mollified) Laplacian path flips edges to a Delaunay
// triangulation and refuses to flip across a boundary, so the constructor throws
// on any mesh with boundary. This is a current limitation; the test pins it
// rather than assuming robustness applies everywhere.
TEST(HeatMethodDistanceSuite, RobustLaplacianRejectsBoundaryMeshes) {
  MeshGeom mg = makeUnitGrid(4);
  ASSERT_GT(std::get<0>(mg)->nBoundaryLoops(), 0u);
  EXPECT_THROW(HeatMethodDistanceSolver solver(*std::get<1>(mg), 1.0, /*useRobustLaplacian=*/true),
               std::runtime_error);
}

// On a closed mesh the robust and default Laplacians agree: with no boundary and
// no non-manifold vertices there is nothing to repair.
TEST(HeatMethodDistanceSuite, RobustLaplacianAgreesOnClosedManifold) {
  MeshGeom mg = makeOctahedron();
  ManifoldSurfaceMesh& mesh = *std::get<0>(mg);
  VertexPositionGeometry& geom = *std::get<1>(mg);

  HeatMethodDistanceSolver plain(geom);
  HeatMethodDistanceSolver robust(geom, 1.0, /*useRobustLaplacian=*/true);

  VertexData<double> dPlain = plain.computeDistance(mesh.vertex(0));
  VertexData<double> dRobust = robust.computeDistance(mesh.vertex(0));

  EXPECT_LT(maxAbsDiff(dPlain, dRobust), 1e-10)
      << "robust and default Laplacians disagree on a closed manifold: " << fieldToString(dRobust);
}

// A zero-area (degenerate) triangle makes the cotangent weights non-finite, and
// the matrix factorization rejects them. Whether by throwing or by returning a
// non-finite field, the solver must never hand back a plausible finite distance
// for such input. Either outcome is acceptable; the dangerous one is excluded.
TEST(HeatMethodDistanceSuite, ZeroAreaFaceIsNeverSilentlyAccepted) {
  MeshGeom mg = makeCollinearTriangle();
  ManifoldSurfaceMesh& mesh = *std::get<0>(mg);

  bool threw = false;
  bool producedFiniteField = false;
  try {
    HeatMethodDistanceSolver solver(*std::get<1>(mg));
    VertexData<double> d = solver.computeDistance(mesh.vertex(0));
    producedFiniteField = d.toVector().allFinite();
  } catch (const std::exception&) {
    threw = true;
  }

  EXPECT_TRUE(threw || !producedFiniteField)
      << "zero-area triangle produced a finite distance field; this must not be interpreted as a valid distance";
}

// A *nearly* degenerate (sliver) triangle is accepted: the geometry is finite, so
// no exception is raised, but the intrinsic metric is badly conditioned and the
// resulting distances are not trustworthy. The only guarantee is that the output
// is finite. This test documents the magnitude of the failure by contrasting it
// with a well-shaped mesh on identical connectivity.
TEST(HeatMethodDistanceSuite, SliverFacesYieldFiniteButUnreliableValues) {
  MeshGeom mgSliver = makeApexQuad(/*wellShaped=*/false);
  HeatMethodDistanceSolver sliverSolver(*std::get<1>(mgSliver));
  VertexData<double> dSliver = sliverSolver.computeDistance(std::get<0>(mgSliver)->vertex(2));

  EXPECT_TRUE(dSliver.toVector().allFinite()) << "sliver mesh produced non-finite values";
  EXPECT_NEAR(dSliver[std::get<0>(mgSliver)->vertex(2)], 0., 1e-12)
      << "vertex source should still be the zero of the field on a sliver mesh";

  MeshGeom mgWell = makeApexQuad(/*wellShaped=*/true);
  HeatMethodDistanceSolver wellSolver(*std::get<1>(mgWell));
  VertexData<double> dWell = wellSolver.computeDistance(std::get<0>(mgWell)->vertex(2));

  double sliverValue = dSliver[std::get<0>(mgSliver)->vertex(0)];
  double wellValue = dWell[std::get<0>(mgWell)->vertex(0)];

  // The well-shaped control is within ordinary coarse-mesh error of the exact
  // distance from the apex to vertex 0 (here sqrt(2)); it is clearly non-zero.
  EXPECT_NEAR(wellValue, std::sqrt(2.), 0.6) << "well-shaped control mesh is unexpectedly inaccurate";
  EXPECT_GT(wellValue, 0.5);

  // On the sliver the same quantity collapses toward zero: finite, but a
  // numerically broken value rather than a distance. Do not interpret it as one.
  EXPECT_GE(sliverValue, 0.);
  EXPECT_LT(sliverValue, 0.05 * wellValue)
      << "sliver no longer corrupts the field; if the algorithm improves, relax this boundary test";
}
