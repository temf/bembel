// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_H2MATRIX_H2MULTIPOLE_H_
#define BEMBEL_H2MATRIX_H2MULTIPOLE_H_

namespace Bembel {
namespace H2Multipole {
/**
 *  \ingroup H2Matrix
 *  \brief  Computes the number_of_points Chebychev points.
 */
struct ChebychevRoots {
  ChebychevRoots() {}
  ChebychevRoots(int number_of_points) {
    init_ChebyshevRoots(number_of_points);
  }
  void init_ChebyshevRoots(int n_pts) {
    auto grid_pts = Eigen::ArrayXd::LinSpaced(n_pts, 1, n_pts).reverse();
    double alpha = BEMBEL_PI / (2. * double(n_pts));
    points_ = 0.5 * (alpha * (2 * grid_pts - 1)).cos() + 0.5;
    return;
  }
  Eigen::VectorXd points_;
};
/**
 *  \ingroup H2Matrix
 *  \brief computes the Lagrange polynomials wrt the interpolation
 *         points given by the InterpolationPoints struct wrt.
 *         Newton basis
 **/
template <typename InterpolationPoints>
Eigen::MatrixXd computeLagrangePolynomials(int number_of_points) {
  Eigen::MatrixXd retval =
      Eigen::MatrixXd::Identity(number_of_points, number_of_points);
  Eigen::VectorXd x = InterpolationPoints(number_of_points).points_;
  for (auto i = 0; i < number_of_points; ++i)
    for (auto j = 1; j < number_of_points; ++j)
      for (auto k = number_of_points - 1; k >= j; --k)
        retval(k, i) = (retval(k, i) - retval(k - 1, i)) / (x(k) - x(k - j));
  return retval;
};
/**
 *  \ingroup H2Matrix
 *  \brief evaluates a given polynomial in the Newton basis wrt the
 *         interpolation points at a given location
 **/
template <typename InterpolationPoints>
double evaluatePolynomial(const Eigen::VectorXd &L, double xi) {
  int number_of_points = L.size();
  Eigen::VectorXd x = InterpolationPoints(number_of_points).points_;
  double retval = 0;

  retval = L(number_of_points - 1);

  for (auto i = number_of_points - 2; i >= 0; --i)
    retval = retval * (xi - x(i)) + L(i);

  return retval;
}
/**
 * \ingroup H2Matrix
 * \brief Computes transfer matrices in required order to apply them
 *        all-at-once in a matrix-vector-product, i.e. the order
 *        in the output is [T0 T2 T3 T1].
 *
 *
 */
template <typename InterpolationPoints>
Eigen::MatrixXd computeTransferMatrices(int number_of_points) {
  int np2 = number_of_points * number_of_points;
  Eigen::MatrixXd E(number_of_points, 2 * number_of_points);
  Eigen::MatrixXd T(np2, 4 * np2);
  /// initialize Lagrange polynomials and interpolation nodes
  Eigen::VectorXd x = InterpolationPoints(number_of_points).points_;
  Eigen::MatrixXd L =
      computeLagrangePolynomials<InterpolationPoints>(number_of_points);
  /**
   * initialize values of Lagrange functions on the interpolation points
   * as follows:
   * E(:,0:npts-1) containes the values of the Lagrange polynomials on the
   * interpolation points scaled to [0, 0.5].
   * E(:,npts:2*npts-1) containes the values of the Lagrange polynomials on
   * the interpolation points scaled to [0.5, 1].
   * E(i,j) containes the value of the j.th polynomial on the ith point
   **/
  for (auto j = 0; j < number_of_points; ++j)
    for (auto i = 0; i < number_of_points; ++i) {
      E(i, j) = evaluatePolynomial<InterpolationPoints>(L.col(j), 0.5 * x(i));
      E(i, j + number_of_points) =
          evaluatePolynomial<InterpolationPoints>(L.col(j), 0.5 * x(i) + 0.5);
    }
  // This construction, regard permuation vector, results in the
  // order T0 T2 T3 T1 to apply to the moments
  Eigen::Vector4i permutation;
  permutation << 0, 3, 1, 2;
  for (auto k = 0; k < 4; ++k)
    for (auto i = 0; i < number_of_points; ++i)
      for (auto ii = 0; ii < number_of_points; ++ii)
        for (auto j = 0; j < number_of_points; ++j)
          for (auto jj = 0; jj < number_of_points; ++jj)
            T(j * number_of_points + jj,
              i * number_of_points + ii + np2 * permutation(k)) =
                E(i, j + (k / 2) * number_of_points) *
                E(ii, jj + (k % 2) * number_of_points);

  return T;
}
/**
 *  \ingroup H2Matrix
 *  \brief Computes 1D moment for FMM. All calculations ar performed on [0,1]
 */
template <typename InterpolationPoints, typename LinOp>
struct Moment1D {
  static Eigen::MatrixXd computeMoment1D(const SuperSpace<LinOp> &super_space,
                                         const int cluster_level,
                                         const int cluster_refinements,
                                         const int number_of_points) {
    int n = 1 << cluster_refinements;
    double h = 1. / n;
    int N = 1 << cluster_level;
    double H = 1. / N;
    int polynomial_degree = super_space.get_polynomial_degree();
    int polynomial_degree_plus_one = polynomial_degree + 1;
    int polynomial_degree_plus_one_squared =
        polynomial_degree_plus_one * polynomial_degree_plus_one;
    GaussLegendre<Constants::maximum_quadrature_degree> GL;
    auto Q =
        GL[(int)std::ceil(0.5 * (number_of_points + polynomial_degree - 2))];

    Eigen::VectorXd x = InterpolationPoints(number_of_points).points_;
    Eigen::MatrixXd L =
        computeLagrangePolynomials<InterpolationPoints>(number_of_points);

    Eigen::MatrixXd moment(number_of_points, n * polynomial_degree_plus_one);

    for (auto i = 0; i < number_of_points; ++i) {
      Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar,
                    Eigen::Dynamic, 1>
          intval(polynomial_degree_plus_one, 1);
      for (auto j = 0; j < n; ++j) {
        intval.setZero();
        for (auto k = 0; k < Q.w_.size(); ++k) {
          super_space.addScaledBasis1D(
              &intval,
              Q.w_(k) * std::sqrt(h * H) *
                  evaluatePolynomial<InterpolationPoints>(L.col(i),
                                                          h * (j + Q.xi_(k))),
              Q.xi_(k));
        }
        // we are integrating real stuff, so take the real part
        moment.block(i, j * polynomial_degree_plus_one, 1,
                     polynomial_degree_plus_one) = intval.real().transpose();
      }
    }

    return moment;
  }
};
/**
 *  \ingroup H2Matrix
 *  \brief Computes 1D moment for FMM using derivatives of the basis functions.
 * All calculations ar performed on [0,1].
 */
template <typename InterpolationPoints, typename LinOp>
struct Moment1DDerivative {
  static Eigen::MatrixXd computeMoment1D(const SuperSpace<LinOp> &super_space,
                                         const int cluster_level,
                                         const int cluster_refinements,
                                         const int number_of_points) {
    int n = 1 << cluster_refinements;
    double h = 1. / n;
    int N = 1 << cluster_level;
    double H = 1. / N;
    int polynomial_degree = super_space.get_polynomial_degree();
    int polynomial_degree_plus_one = polynomial_degree + 1;
    int polynomial_degree_plus_one_squared =
        polynomial_degree_plus_one * polynomial_degree_plus_one;
    GaussLegendre<Constants::maximum_quadrature_degree> GL;
    auto Q =
        GL[(int)std::ceil(0.5 * (number_of_points + polynomial_degree - 2))];

    Eigen::VectorXd x = InterpolationPoints(number_of_points).points_;
    Eigen::MatrixXd L =
        computeLagrangePolynomials<InterpolationPoints>(number_of_points);

    Eigen::MatrixXd moment(number_of_points, n * polynomial_degree_plus_one);

    for (auto i = 0; i < number_of_points; ++i) {
      Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar,
                    Eigen::Dynamic, 1>
          intval(polynomial_degree_plus_one, 1);
      for (auto j = 0; j < n; ++j) {
        intval.setZero();
        for (auto k = 0; k < Q.w_.size(); ++k) {
          super_space.addScaledBasis1DDx(
              &intval,
              Q.w_(k) / std::sqrt(h * H) *
                  evaluatePolynomial<InterpolationPoints>(L.col(i),
                                                          h * (j + Q.xi_(k))),
              Q.xi_(k));
        }
        // we are integrating real stuff, so take the real part
        moment.block(i, j * polynomial_degree_plus_one, 1,
                     polynomial_degree_plus_one) = intval.real().transpose();
      }
    }

    return moment;
  }
};
/**
 * \ingroup H2Matrix
 * \brief Computes a single 2D moment for the FMM by tensorisation of the 1D
 * moments
 */
template <typename Mom1D_1, typename Mom1D_2, typename LinOp>
Eigen::MatrixXd moment2DComputer(const SuperSpace<LinOp> &super_space,
                                 const int cluster_level,
                                 const int cluster_refinements,
                                 const int number_of_points) {
  int n = 1 << cluster_refinements;
  auto n2 = n * n;
  int polynomial_degree = super_space.get_polynomial_degree();
  int polynomial_degree_plus_one = polynomial_degree + 1;
  int polynomial_degree_plus_one_squared =
      polynomial_degree_plus_one * polynomial_degree_plus_one;

  // compute 1D moments
  Eigen::MatrixXd moment1D_1 = Mom1D_1::computeMoment1D(
      super_space, cluster_level, cluster_refinements, number_of_points);
  Eigen::MatrixXd moment1D_2 = Mom1D_2::computeMoment1D(
      super_space, cluster_level, cluster_refinements, number_of_points);
  /**
   *  Throughout this code we face the problem of memory serialisation for
   *  the traversel of elements in the element tree. the canonical orders would
   *  either be row major or column major traversal of a given patch
   *  resulting in e.g. element(i,j) = *(elementRoot + i * n + j)
   *  for n = 1 << j. However, to achieve a localisation of matrix blocks if
   *  their corresponding shape functions are geometrically close, we typically
   *  use a hierarchical ordering of the elements. The mapping from hierarchical
   *  to row major ordering is easily facilitated by computing the elements
   *  position in its patch using its llc_. For our special situation of a
   *  balanced quadtree, this mapping is the same for every patch.
   *  In particular, the moments in our construction of the FMM only depend on
   *  the current level of uniform refinement and not on a particular patch.
   *  Next, we explicitly set up this mapping, where index_s determines the
   *  location of a given element along the first coordinate in the reference
   *  domain and index_t its second coordinate in the reference domain.
   *  Note that this piece of code only has to be used once in the entire setup
   *  of the FMM, such that the overhead is quite small.
   **/
  Eigen::VectorXi index_s(n2);
  Eigen::VectorXi index_t(n2);
  index_s.setZero();
  index_t.setZero();
  for (auto j = 0; j < cluster_refinements; ++j) {
    auto inc = 1 << 2 * j;
    for (auto i = 0; i < inc; ++i) {
      index_s(4 * i) = 2 * index_s(i);
      index_s(4 * i + 1) = 2 * index_s(i) + 1;
      index_s(4 * i + 2) = 2 * index_s(i) + 1;
      index_s(4 * i + 3) = 2 * index_s(i);
      index_t(4 * i) = 2 * index_t(i);
      index_t(4 * i + 1) = 2 * index_t(i);
      index_t(4 * i + 2) = 2 * index_t(i) + 1;
      index_t(4 * i + 3) = 2 * index_t(i) + 1;
    }
  }

  // assemble 2D tensor-product moments
  Eigen::MatrixXd moment2D(number_of_points * number_of_points,
                           moment1D_1.cols() * moment1D_2.cols());
  for (auto i = 0; i < number_of_points; ++i)
    for (auto j = 0; j < number_of_points; ++j)
      for (auto k = 0; k < n2; ++k)
        for (auto m1 = 0; m1 < polynomial_degree_plus_one; ++m1)
          for (auto m2 = 0; m2 < polynomial_degree_plus_one; ++m2)
            moment2D(i * number_of_points + j,
                     polynomial_degree_plus_one_squared * k +
                         m1 * polynomial_degree_plus_one + m2) =
                moment1D_1(i, index_s(k) * polynomial_degree_plus_one + m2) *
                moment1D_2(j, index_t(k) * polynomial_degree_plus_one + m1);

  return moment2D;
}
/**
 * \ingroup H2Matrix
 * \brief Computes all 2D moment for the FMM by tensorisation of the 1D
 * moments. Specialice this for your linear operator if you need derivatives on
 * your local shape functions. See e.g. MaxwellSingleLayerOperator for an
 * example.
 */
template <typename InterpolationPoints, typename LinOp>
struct Moment2D {
  static std::vector<Eigen::MatrixXd> compute2DMoment(
      const SuperSpace<LinOp> &super_space, const int cluster_level,
      const int cluster_refinements, const int number_of_points) {
    std::vector<Eigen::MatrixXd> vector_of_moments;
    for (int i = 0;
         i <
         getFunctionSpaceVectorDimension<LinearOperatorTraits<LinOp>::Form>();
         ++i)
      vector_of_moments.push_back(
          moment2DComputer<Moment1D<InterpolationPoints, LinOp>,
                           Moment1D<InterpolationPoints, LinOp>>(
              super_space, cluster_level, cluster_refinements,
              number_of_points));
    return vector_of_moments;
  }
};
/**
 * \ingroup H2Matrix
 * \brief Compute tensor interpolation points on unit square from 1D
 * interpolation points.
 */
Eigen::MatrixXd interpolationPoints2D(const Eigen::VectorXd &x) {
  Eigen::MatrixXd x2(x.size() * x.size(), 2);
  for (auto i = 0; i < x.size(); ++i)
    for (auto j = 0; j < x.size(); ++j) {
      x2.row(i * x.size() + j) = Eigen::Vector2d(x(i), x(j));
    }
  return x2;
}
/**
 * \ingroup H2Matrix
 * \brief Interpolate kernel function on reference domain for FMM.
 */
template <typename LinOp>
Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar, Eigen::Dynamic,
              Eigen::Dynamic>
interpolateKernel(const LinOp &linOp, const SuperSpace<LinOp> &super_space,
                  const Eigen::MatrixXd &x,
                  const Bembel::ElementTreeNode &cluster1,
                  const Bembel::ElementTreeNode &cluster2) {
  SurfacePoint qp1, qp2;
  Eigen::Matrix<
      typename LinearOperatorTraits<LinOp>::Scalar,
      getFunctionSpaceVectorDimension<LinearOperatorTraits<LinOp>::Form>() *
          LinearOperatorTraits<LinOp>::NumberOfFMMComponents,
      getFunctionSpaceVectorDimension<LinearOperatorTraits<LinOp>::Form>() *
          LinearOperatorTraits<LinOp>::NumberOfFMMComponents>
      interpval;
  const int vector_dimension = Bembel::getFunctionSpaceVectorDimension<
      Bembel::LinearOperatorTraits<LinOp>::Form>();
  const int number_of_FMM_components =
      Bembel::LinearOperatorTraits<LinOp>::NumberOfFMMComponents;
  Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar, Eigen::Dynamic,
                Eigen::Dynamic>
      F(vector_dimension * number_of_FMM_components * x.rows(),
        vector_dimension * number_of_FMM_components * x.rows());
  for (int i = 0; i < x.rows(); ++i) {
    super_space.map2surface(cluster1, x.row(i), 1., &qp1);
    for (int j = 0; j < x.rows(); ++j) {
      super_space.map2surface(cluster2, x.row(j), 1., &qp2);
      auto FMM_output = linOp.evaluateFMMInterpolation(qp1, qp2);
      for (int ii = 0; ii < vector_dimension; ++ii)
        for (int jj = 0; jj < vector_dimension; ++jj)
          for (int iii = 0; iii < number_of_FMM_components; ++iii)
            for (int jjj = 0; jjj < number_of_FMM_components; ++jjj)
              F(i + (ii * number_of_FMM_components + iii) * x.rows(),
                j + (jj * number_of_FMM_components + jjj) * x.rows()) =
                  FMM_output(iii + ii * number_of_FMM_components,
                             jjj + jj * number_of_FMM_components);
    }
  }
  return F;
}
/**
 * \ingroup H2Matrix
 * \brief Forward transformation for FMM.
 */
template <typename Scalar>
std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>
forwardTransformation(const Eigen::MatrixXd &moment_matrices,
                      const Eigen::MatrixXd &transfer_matrices, const int steps,
                      const Eigen::Matrix<Scalar, Eigen::Dynamic,
                                          Eigen::Dynamic> &long_rhs_matrix) {
  // get numbers
  int number_of_points = transfer_matrices.rows();
  int number_of_FMM_components = moment_matrices.rows() / number_of_points;
  // apply moment matrices
  std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>
      long_rhs_forward;
  long_rhs_forward.push_back(moment_matrices * long_rhs_matrix);
  // apply transfer matrices
  for (int i = 0; i < steps; ++i) {
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> reshaped =
        Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(
            long_rhs_forward.back().data(), 4 * long_rhs_forward.back().rows(),
            long_rhs_forward.back().cols() / 4);
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> transferred(
        reshaped.rows() / 4, reshaped.cols());

    // iterate over FMM components
    for (int j = 0; j < number_of_FMM_components; ++j) {
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> to_transfer(
          4 * number_of_points, reshaped.cols());
      to_transfer << reshaped.block(
          (j + 0 * number_of_FMM_components) * number_of_points, 0,
          number_of_points, reshaped.cols()),
          reshaped.block((j + 1 * number_of_FMM_components) * number_of_points,
                         0, number_of_points, reshaped.cols()),
          reshaped.block((j + 2 * number_of_FMM_components) * number_of_points,
                         0, number_of_points, reshaped.cols()),
          reshaped.block((j + 3 * number_of_FMM_components) * number_of_points,
                         0, number_of_points, reshaped.cols());
      transferred.block(j * number_of_points, 0, number_of_points,
                        reshaped.cols()) = transfer_matrices * to_transfer;
    }

    long_rhs_forward.push_back(transferred);
  }
  return long_rhs_forward;
}
/**
 * \ingroup H2Matrix
 * \brief Backward transformation for FMM. The content of long_dst_backward is
 * destroyed.
 */
template <typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> backwardTransformation(
    const Eigen::MatrixXd &moment_matrices,
    const Eigen::MatrixXd &transfer_matrices, const int steps,
    std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>
        &long_dst_backward) {
  // get numbers
  int number_of_points = transfer_matrices.rows();
  int number_of_FMM_components = moment_matrices.rows() / number_of_points;
  // apply transfer matrices
  for (int i = steps; i > 0; --i) {
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> transferred(
        4 * long_dst_backward[i].rows(), long_dst_backward[i].cols());
    for (int j = 0; j < number_of_FMM_components; ++j) {
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> prod =
          transfer_matrices.transpose() *
          long_dst_backward[i].block(j * number_of_points, 0, number_of_points,
                                     long_dst_backward[i].cols());
      transferred.block((j + 0 * number_of_FMM_components) * number_of_points,
                        0, number_of_points, prod.cols()) =
          prod.block(0 * number_of_points, 0, number_of_points, prod.cols());
      transferred.block((j + 1 * number_of_FMM_components) * number_of_points,
                        0, number_of_points, prod.cols()) =
          prod.block(1 * number_of_points, 0, number_of_points, prod.cols());
      transferred.block((j + 2 * number_of_FMM_components) * number_of_points,
                        0, number_of_points, prod.cols()) =
          prod.block(2 * number_of_points, 0, number_of_points, prod.cols());
      transferred.block((j + 3 * number_of_FMM_components) * number_of_points,
                        0, number_of_points, prod.cols()) =
          prod.block(3 * number_of_points, 0, number_of_points, prod.cols());
    }
    long_dst_backward[i - 1] +=
        Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(
            transferred.data(), long_dst_backward[i - 1].rows(),
            long_dst_backward[i - 1].cols());
  }
  // apply moment matrices
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> long_dst_matrix =
      moment_matrices.transpose() * long_dst_backward[0];
  return Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(
      long_dst_matrix.data(), long_dst_matrix.size());
}

}  // namespace H2Multipole
}  // namespace Bembel
#endif
