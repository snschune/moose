//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RayTracingAngularQuadrature.h"

#include "MooseUtils.h"

#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

namespace RayTracingAngularQuadrature
{

void
gaussLegendre(const unsigned int order, std::vector<Real> & x, std::vector<Real> & w)
{
  if (order % 2 == 0)
    mooseError("Order must be odd in gaussLegendre()");

  const std::size_t n = order + 1;
  x.resize(n);
  w.resize(n);
  DenseMatrix<Real> mat(n, n);
  DenseVector<Real> lambda(n);
  DenseVector<Real> lambdai(n);
  DenseMatrix<Real> vec(n, n);

  for (unsigned int i = 1; i < n; ++i)
  {
    Real ri = i;
    mat(i, i - 1) = ri / std::sqrt(((2. * ri - 1.) * (2. * ri + 1.)));
    mat(i - 1, i) = mat(i, i - 1);
  }
  mat.evd_right(lambda, lambdai, vec);

  for (unsigned int i = 0; i < n; ++i)
  {
    x[i] = 0.5 * (lambda(i) + 1.0);
    w[i] = vec(0, i) * vec(0, i);
  }

  // Sort based on the points
  std::vector<std::size_t> sorted_indices(x.size());
  std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
  std::stable_sort(sorted_indices.begin(), sorted_indices.end(), [&x](size_t i1, size_t i2) {
    return x[i1] < x[i2];
  });
  const auto x_copy = x;
  const auto w_copy = w;
  for (std::size_t i = 0; i < x.size(); ++i)
  {
    x[i] = x_copy[sorted_indices[i]];
    w[i] = w_copy[sorted_indices[i]];
  }
}

DenseMatrix<Real>
getRotationMatrix(const Point & direction, const unsigned int dim)
{
  DenseMatrix<Real> rot(3, 3);

  if (dim == 2)
  {
    rot(0, 0) = direction(0);
    rot(0, 1) = -direction(1);
    rot(1, 0) = direction(1);
    rot(1, 1) = direction(0);
    rot(2, 2) = 1;
  }
  else if (dim == 3)
  {
    // Create a local coordinate system around direction
    const Point tx = RayTracingAngularQuadrature::getOrthonormalVector(direction, dim);
    const Point ty = direction.cross(tx);

    // Create rotation matrix and rotate vector omega
    for (unsigned int j = 0; j < 3; ++j)
    {
      rot(j, 0) = tx(j);
      rot(j, 1) = ty(j);
      rot(j, 2) = direction(j);
    }
  }
  else
    ::mooseError("Dimension ", dim, " not supported in aqRoationMatrix()");

  return rot;
}

void
getHalfRange3D(const unsigned int chebyshev_order,
               const unsigned int gauss_legendre_order,
               std::vector<std::pair<Real, Real>> & x,
               std::vector<Real> & w,
               const Real half_angle)
{
  if (half_angle <= 0 || half_angle > M_PI / 2)
    ::mooseError("Half angle out of range in getHalfRange3D(). Must be > 0 and < pi / 2");

  std::vector<Real> chebyshev_x(chebyshev_order);
  std::vector<Real> chebyshev_w(chebyshev_order);

  // Chebysheb quadrature on [0, 2 \pi]
  for (std::size_t j = 0; j < chebyshev_order; ++j)
  {
    chebyshev_x[j] = 2 * (Real)j * M_PI / (Real)chebyshev_order;
    chebyshev_w[j] = 2 * M_PI / (Real)chebyshev_order;
  }

  // Gauss-Legendre quadrature on [0, 1]
  std::vector<Real> gauss_legendre_x;
  std::vector<Real> gauss_legendre_w;
  gaussLegendre(gauss_legendre_order, gauss_legendre_x, gauss_legendre_w);

  // Product quadrature
  const std::size_t product_size = chebyshev_x.size() * gauss_legendre_x.size();
  x.resize(product_size);
  w.resize(product_size);

  std::size_t l = 0;
  for (std::size_t i = 0; i < chebyshev_x.size(); ++i)
    for (std::size_t j = 0; j < gauss_legendre_x.size(); ++j)
    {
      x[l].first = chebyshev_x[i];
      x[l].second = gauss_legendre_x[j] * (1 - std::cos(half_angle)) + std::cos(half_angle);
      w[l] = gauss_legendre_w[j] * chebyshev_w[i];
      ++l;
    }
}

void
getHalfRange2D(const unsigned int gauss_legendre_order,
               std::vector<std::pair<Real, Real>> & x,
               std::vector<Real> & w,
               const Real half_angle)
{
  if (half_angle <= 0 || half_angle > 0.5 * M_PI)
    ::mooseError("Half angle out of range in getHalfRange2D(). Must be > 0 and < pi / 2");

  // Gauss-Legendre quadrature on [0, 1]
  std::vector<Real> gauss_legendre_x;
  std::vector<Real> gauss_legendre_w;
  gaussLegendre(gauss_legendre_order, gauss_legendre_x, gauss_legendre_w);
  x.resize(2 * gauss_legendre_w.size());
  w.resize(2 * gauss_legendre_w.size());

  std::size_t l = 0;
  for (std::size_t j = 0; j < gauss_legendre_w.size(); ++j)
  {
    const auto weight = gauss_legendre_w[j] * half_angle;

    // Positive quad
    x[l].first = 1;
    x[l].second = gauss_legendre_x[j] * half_angle;
    w[l] = weight;
    ++l;

    // Negative quad
    x[l].first = -1;
    x[l].second = -gauss_legendre_x[j] * half_angle;
    w[l] = weight;
    ++l;
  }
}

Point
getOrthonormalVector(const Point & v, const unsigned int dim)
{
  if (MooseUtils::absoluteFuzzyLessEqual(v.norm(), 0))
    ::mooseError("Vector v has norm close to 0 in getOrthonormalVector()");
  if (dim != 2 && dim != 3)
    ::mooseError("Dimension must be 2 or 3 but is provided as ", dim, " in getOrthonormalVector()");

  if (v(0) == 0)
    return Point(1, 0, 0);
  else if (v(1) == 0)
    return Point(0, 1, 0);
  else if (v(2) == 0 && dim == 3)
    return Point(0, 0, 1);

  Point t(-v(1), v(0), 0);
  return t.unit();
}

Point
getDirection(const std::size_t l,
             const DenseMatrix<Real> & rotation_matrix,
             const std::vector<std::pair<Real, Real>> & angles,
             const unsigned int dim)
{
  if (angles.size() <= l)
    ::mooseError("l out of range for the given angles in getDirection()");

  DenseVector<Real> omega(3);
  if (dim == 2)
  {
    const Real theta = angles[l].second;
    if (angles[l].first != 1 && angles[l].first != -1)
      ::mooseError("Invalid angles in getDirection(). You likely used the wrong half range getter "
                   "based on the provided dimension.");

    omega(0) = std::cos(theta);
    omega(1) = std::sin(theta);
  }
  else if (dim == 3)
  {
    const Real phi = angles[l].first;
    const Real mu = angles[l].second;

    omega(0) = sqrt(1 - mu * mu) * cos(phi);
    omega(1) = sqrt(1 - mu * mu) * sin(phi);
    omega(2) = mu;
  }
  else
    ::mooseError("Dimension must be 2 or 3 but is provided as ", dim, " in getDirection()");

  // Rotate
  DenseVector<Real> omega_p(3);
  rotation_matrix.vector_mult(omega_p, omega);

  return Point(omega_p(0), omega_p(1), omega_p(2));
}

}
