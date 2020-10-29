//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MooseTypes.h"

namespace RayTracingAngularQuadrature
{

/**
 * Builds the Gauss-Legendre quadrature over polar angle theta for theta \in
 * [half_range, -half_range] (in rad).
 *
 * In 2D it is more advantageous to integrate over angle as opposed to mu = cos(theta).
 *
 * The abscissae are to be interpreted as angles as follows: theta_j =
 * x[j].second, while x[j].first indicates if it theta is positive or negative (mostly unused).
 *
 * @params gauss_legendre_order The quadrature order (must be odd).
 * @params x To be filled with the abscissae.
 * @params w To be filled with the weights
 * @params half_angle The upper bound of the quadrature range in rad.
 */
void getHalfRange2D(const unsigned int gauss_legendre_order,
                    std::vector<std::pair<Real, Real>> & x,
                    std::vector<Real> & w,
                    const Real half_angle = 0.5 * M_PI);

/**
 * Builds the half range (or less, depending on \p half_angle) Gauss-Legendre-Chebyshev
 * product quadrature.
 *
 * @params chebyhev_order The order for the Chebyshev quadrature (azimuthal). Total
 * number of abcissae in phi.
 * @params gauss_legendre_order The order for the Guass-Legendre quadrature (polar);
 * must be odd.
 * @params x To be filled with the abscissae (phi and mu).
 * @params w To be filled with the product weights.
 * @params half_angle The half-angle upper bound in rad. Defaults to the full half-range.
 */
void getHalfRange3D(const unsigned int chebyshev_order,
                    const unsigned int gauss_legendre_order,
                    std::vector<std::pair<Real, Real>> & x,
                    std::vector<Real> & w,
                    const Real half_angle = 0.5 * M_PI);

/**
 * @return A vector orthogonal to \p v with dimension \p dim.
 */
Point getOrthonormalVector(const Point & v, const unsigned int dim);

/**
 * The angular quadrature can be used to create angular directions w.r.t. to
 * a reference vector.
 *
 * This obtains a rotation matrix for the purposes of rotating the quadrature
 * to be in the coordinate system of said reference vector to be used in
 * getDirection().
 */
DenseMatrix<Real> getRotationMatrix(const Point & normal, const unsigned int dim);

/**
 * Gets the rotated angular quadrature direction associated with index \p l.
 *
 * @param l The direction index (must be < the size of \p angles).
 * @param rotation_matrix The rotation matrix used to rotate the quadrature.
 * @param angles The angles created by getHalfRange2D() or getHalfRange3D().
 * @param dim The dimension.
 */
Point getDirection(const std::size_t l,
                   const DenseMatrix<Real> & rotation_matrix,
                   const std::vector<std::pair<Real, Real>> & angles,
                   const unsigned int dim);

/**
 * Builds the Gauss-Legendre quadrature for order \p order (must be odd).
 * The points are on [0, 1] and the weights sum to 1.
 *
 * Fills the points into \p x and the weights into \p w.
 */
void gaussLegendre(const unsigned int order, std::vector<Real> & x, std::vector<Real> & w);

}
