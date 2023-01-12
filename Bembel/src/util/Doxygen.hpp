/**
 * \defgroup AnsatzSpace AnsatzSpace
 * \brief The AnsatzSpace module contains the routines managing the discrete
 * space on the surface of the geometry.
 *
 * This is realised through the four
 * classes SuperSpace , Projector , Glue , and AnsatzSpace . Therein, SuperSpace
 * manages local polynomial bases on every element. Through a transformation
 * matrix generated by the template class Projector which depends on the
 * specializa- tion of LinearOperatorBase and its defined traits, the SuperSpace
 * can be related to B-Spline bases on every patch. To build conforming spaces
 * (in the case of DifferentialForm::DivergenceConforming through continuity of
 * the normal component across patch interfaces, in the case of
 * DifferentialForm::Continuous through global C0-continuity), the template
 * class Glue assembles another transformation matrix to identify degrees of
 * freedom across edges. Then, a coefficient vector in the SuperSpace can be
 * related to one of the smooth B-Spline basis
 **/

/**
 *  \defgroup ClusterTree ClusterTree
 * \brief The ClusterTree module introduces a structure onto the parametric
 * surfaces.
 *
 * Therein, the ClusterTree class takes a Geometry and introduces an
 * element struc- ture on it, in the form of an ElementTree . Each
 * ElementTreeNode of the ElementTree corresponds to an entire parametric
 * mapping (at the trees first level) or to a sub element induced by a recursive
 * refinement strategy. The leaves of the tree correspond to the elements on
 * which the SuperSpace introduces shape functions.
 **/

/**
 *  \defgroup DuffyTrick DuffyTrick
 *  \brief The DuffyTrick module provides quadrature routines for (nearly)
 * singular integrals.
 *
 * Therein, the implementation is ssentially an adaptation of the so-called
 * Sauter-Schwab quadrature rules
 **/

/**
 *  \defgroup Geometry Geometry
 *  \brief This module handles all geometry related concerns.
 *
 *  The Geometry class is the interface between geometry description and the
 * remainder of the code. We provide an implementation of NURBS discretized
 * patches via the Patch class, but the code is easily extensible to other
 * parametric descriptions mapping from [0,1]2 to parts of the geometry, as long
 * as the corresponding methods for point evaluation and evaluation of the
 * pointwise Jacobian are implemented
 **/

/**
 *  \defgroup H2Matrix H2Matrix
 *  \brief The H2Matrix module provides functionality for an efficient
 * compression of the system matrix and reduction of the computational
 *complexity of the matrix-vector multiplication
 **/

/**
 *  \defgroup Helmholtz Helmholtz
 *  \brief The Helmholtz module provides some specializations to solve Helmholtz
 *problems.
 **/

/**
 *  \defgroup HomogenisedLaplace HomogenisedLaplace
 *  \brief The HomogenisedLaplace module provides som specialisations to solve
 *the homogenised Laplace problem.
 */

/**
 *  \defgroup IO IO
 *  \brief The IO module provides input-output functionality.
 **/

/**
 *  \defgroup Laplace Laplace
 *  \brief The Laplace module provides some specializations to solve Laplace
 *problems.
 **/

/**
 *  \defgroup LinearForm LinearForm
 *  \brief The LinearForm template class must be specialized for the assembly of
 *the right hand side.
 **/

/**
 *  \defgroup LinearOperator LinearOperator
 *  \brief Specializations of the LinearOperatorBase class govern the behaviour
 *of most other modules.
 *
 * To provide a valid specialization, methods for kernel
 *evaluation and evaluation of the integrand, i.e., including the test
 *functions, must be provided. Moreover, the corresponding specialization of
 *LinearOperatorTraits must be provided, allowing other classes to determine
 *crucial properties such as the numerical type of the problem (in general
 *double or std::complex<double> ) and the type of discretization, i.e., either
 *DifferentialForm::Continuous , corresponding to a discrete subspace of H1/2 ,
 *DifferentialForm::DivConforming , corresponding to a discrete subspace of
 *H-1/2x(div), or DifferentialForm::Discontinuous , corresponding to a discrete
 *subspace of H-1/2).
 **/

/**
 *  \defgroup Maxwell Maxwell
 *  \brief The Maxwell module provides specializations to solve Maxwell problems
 **/

/**
 *  \defgroup Potential Potential
 *  \brief The Potential module introduces means to evaluate the solution to a
 *PDE via a suitable integral operator taking the unknown of the linear system
 *as input. It relies on a suitable specialization corresponding of the
 *LinearOperatorBase class.
 **/

/**
 *  \defgroup Quadrature Quadrature
 *  \brief The Quadrature module provides quadrature routines for the unit interval/square/cube/.... by means of template recursion.
 **/

/**
 *  \defgroup Spline Spline
 *  \brief The Spline module provides basic routines related to spline function and local polynomials.
 **/
