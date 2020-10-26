#undef NDEBUG

#include "laplacian.h"
#include "mesh.h"
#include <Eigen/SparseCholesky>

using namespace Eigen;
using namespace pmp;
using namespace std;

typedef SparseMatrix<float> SpMat;
typedef PermutationMatrix<Dynamic> Permutation;
typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Mat;

double cotan_weight(const SurfaceMesh &mesh, pmp::Halfedge he) {
  auto points = mesh.get_vertex_property<Point>("v:point");

  // TODO

  return 1;
}

/// Computes the Laplacian matrix in matrix \a L using cotangent weights or the
/// graph Laplacian if useCotWeights==false.
void create_laplacian_matrix(const SurfaceMesh &mesh, SpMat &L,
                             bool useCotWeights) {
  // number of vertices in mesh
  int n = (int)mesh.n_vertices();

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  //tripletList.reserve(n * n);
  int i, j, v_ij;
  for (auto vj : mesh.vertices()) {
    i = vj.idx();
    v_ij = 0;
    for (auto neighbor : mesh.vertices(vj)) {
      j = neighbor.idx();
      tripletList.push_back(T(i, j, 1));
      v_ij++;
    }
    tripletList.push_back(T(i, i, v_ij));
    //cout << i << endl;
  }
  L.setFromTriplets(tripletList.begin(), tripletList.end());
}

/// Computes the permutation putting selected vertices (mask==1) first, and the
/// others at the end. It returns the number of selected vertices.
int create_permutation(const SurfaceMesh &mesh, Permutation &perm) {
  auto masks = mesh.get_vertex_property<int>("v:mask");

  // number of vertices in mesh
  int n = (int)mesh.n_vertices();

  // TODO
  perm.resize(n);
  int i, j;
  int nb_selected_vertices = 0;
  for (auto vj : mesh.vertices()) {
    i = vj.idx();
    if (masks[vj] == 1) {
      perm.indices()[nb_selected_vertices] = i;
      perm.indices()[i] = nb_selected_vertices;
      nb_selected_vertices++;
      //cout << "i: " << i << " --- nb selected vertices: " << nb_selected_vertices << endl;
    } else perm.indices()[i] = i;
  }

  return nb_selected_vertices;
}

/// Performs the poly-harmonic interpolation (order k) over the selected
/// vertices (mask==1) of the vertex attributes u. For each vertex V of index i,
///     if  mask[V]!=1 then u.col(i) is used as input constraints,
///     otherwise, mask[V]==1, and u.col(i) is replaced by the poly-harmonic
///     interpolation of the fixed values.
void poly_harmonic_interpolation(const SurfaceMesh &mesh,
                                 Ref<Mat> u,
                                 int k) {
  // Number of vertices in the mesh
  int n = (int)mesh.n_vertices();

  // 1 - Create the sparse Laplacian matrix
  SpMat L(n,n);
  //cout << "n:" << n << endl;
  //cout << "L size : " << L.size() << endl;
  create_laplacian_matrix(mesh, L, false);
  //cout << "nb non-zeros in L : " << L.nonZeros() << endl;

  // 2 - Create the permutation matrix putting the fixed values at the end,
  //     and the true unknown at the beginning
  Permutation perm;
  int nb_unknowns = create_permutation(mesh, perm);
  cout << "nb unknowns : " << nb_unknowns << endl;
  //for (int i = 0 ; i < n ; i++) cout << "perm[i]: " << perm.indices()[i] << endl;

  // 3 - Apply the permutation to both rows (equations) and columns (unknowns),
  //     i.e., L = P * L * P^-1
  L = L.twistedBy(perm);
  // Perform k-harmonic interpolation
  //for (int i = 1 ; i < k ; i++) L *= L;

  // 4 - solve L * [x^T u^T]^T = 0, i.e., L00 x = - L01 * u
  // solve Ax = b
  u = perm * u;
  SimplicialLDLT<SpMat> solver;
  SpMat A = L.topLeftCorner(nb_unknowns, nb_unknowns);
  solver.compute(A);
  if (solver.info() != Success) {
    cerr << "Cholesky factorisation failed.\n";
    return;
  }
  Mat b = -1 * L.topRightCorner(nb_unknowns, n - nb_unknowns) * u.bottomRows(n - nb_unknowns);
  Mat x = solver.solve(b);
  if (solver.info() != Success) {
    cerr << "Cholesky factorisation failed.\n";
    return;
  }

  // 5 - Copy back the results to u
  for (int i = 0 ; i < nb_unknowns ; i++) {
    for (int j = 0 ; j < 3 ; j++) {
      u(i, j) = x(i, j);
    }
  }
   u = perm.inverse() * u;
}
