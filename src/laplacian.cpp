
#include "laplacian.h"
#include "mesh.h"
#include <Eigen/SparseCholesky>

using namespace Eigen;
using namespace pmp;
using namespace std;

typedef SparseMatrix<float> SpMat;
typedef PermutationMatrix<Dynamic> Permutation;

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
  tripletList.reserve(n * n);

  int i, j, v_ij;
  for (auto vj : mesh.vertices()) {
    i = vj.idx();
    v_ij = 0;
    for (auto neighbor : mesh.vertices(vj)) {
      j = vj.idx();
      tripletList.push_back(T(i, j, -1));
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

  return 0;
}

/// Performs the poly-harmonic interpolation (order k) over the selected
/// vertices (mask==1) of the vertex attributes u. For each vertex V of index i,
///     if  mask[V]!=1 then u.col(i) is used as input constraints,
///     otherwise, mask[V]==1, and u.col(i) is replaced by the poly-harmonic
///     interpolation of the fixed values.
void poly_harmonic_interpolation(const SurfaceMesh &mesh,
                                 Ref<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > u,
                                 int k) {
  // Number of vertices in the mesh
  int n = (int)mesh.n_vertices();

  // 1 - Create the sparse Laplacian matrix
  SpMat L(n,n);
  create_laplacian_matrix(mesh, L, false);
  cout << L.nonZeros() << endl;

  // 2 - Create the permutation matrix putting the fixed values at the end,
  //     and the true unknown at the beginning
  Permutation perm;
  int nb_unknowns = create_permutation(mesh, perm);

  /*
  if (masks[vj]) {

  }
  perm.indices()[???] = ???;
  */

  // 3 - Apply the permutation to both rows (equations) and columns (unknowns),
  //     i.e., L = P * L * P^-1

  // TODO

  // 4 - solve L * [x^T u^T]^T = 0, i.e., L00 x = - L01 * u

  // TODO

  // 5 - Copy back the results to u

  // TODO
}
