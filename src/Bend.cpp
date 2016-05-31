/*---------------------------------------------------------------------------*\
 *
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of MiMMO.
 *
 *  MiMMO is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  MiMMO is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with MiMMO. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#include "Bend.hpp"

using namespace mimmo;

/*!Default constructor of Bend
 */
Bend::Bend(){
	m_name = "MiMMO.Bend";
	m_degree.fill({{0,0,0}});
};

/*!Default destructor of Bend
 */
Bend::~Bend(){};

/*!Copy constructor of Bend.
 */
Bend::Bend(const Bend & other):BaseManipulation(other){
	m_coords = other.m_coords;
	m_degree = other.m_degree;
	m_coeffs = other.m_coeffs;
};

/*!Assignement operator of Bend.
 */
Bend & Bend::operator=(const Bend & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_coords = other.m_coords;
	m_degree = other.m_degree;
	m_coeffs = other.m_coeffs;
	return	*this;
};

/*!It gets the coordinates of the degrees of freedom.
 * \return Coordinates of the degrees of freedom.
 */
dvecarr3E*
Bend::getCoords(){
	return(&m_coords);
};

/*!It gets the degrees of polynomial law for each component of displacements of degrees of freedom.
 * \return Degrees of polynomial laws (degree[i][j] = degree of displacement function si = f(xj)).
 */
umatrix33E
Bend::getDegree(){
	return(m_degree);
};

/*!It gets the coefficients of the polynomial laws.
 * \return Coefficients of the polynomial laws. (coeffs[i][j][k] = coefficients aijk of term si = aij * xj^k).
 */
dmat33Evec
Bend::getCoeffs(){
	return(m_coeffs);
};


/*!It sets the coordinates of the degrees of freedom.
 * \param[in] coords Coordinates of the degrees of freedom.
 */
void
Bend::setCoords(dvecarr3E coords){
	m_coords = coords;
};

/*!It sets the degrees of polynomial law for each component of displacements of degrees of freedom.
 * \param[in] degree Degrees of polynomial laws (degree[i][j] = degree of displacement function si = f(xj)).
 */
void
Bend::setDegree(umatrix33E degree){
	m_degree = degree;
};

/*!It sets the degrees of a term of a polynomial law for a component of displacements of degrees of freedom.
 * \param[in] i Components of displacement.
 * \param[in] j Coordinate of the function related to input degree.
 * \param[in] degree Degrees of polynomial laws (degree[i][j] = degree of displacement function si = f(xj)).
 */
void
Bend::setDegree(int i, int j, uint32_t degree){
	m_degree[i][j] = degree;
	m_coeffs[i][j].resize(degree);
};

/*!It sets the coefficients of the polynomial laws.
 * \param[in] coeffs Coefficients of the polynomial laws. (coeffs[i][j][k] = coefficients aijk of term si = aij * xj^k).
 */
void
Bend::setCoeffs(dmat33Evec coeffs){
	m_coeffs = coeffs;
};

/*!It sets the coefficients of the polynomial laws.
 * \param[in] i Components of displacement.
 * \param[in] j Coordinate of the function related to input degree.
 * \param[in] coeffs Coefficients of the polynomial laws. (coeffs[i][j][k] = coefficients aijk of term si = aij * xj^k).
 */
void
Bend::setCoeffs(int i, int j, dvector1D coeffs){
	m_coeffs[i][j] = coeffs;
};

/*!Execution command. It modifies the displacements given by the input manipulation object
 * with the polynomial law.
 * The input has to be set with a dvecarr3E variable (bend it casts the template method
 * getInput to this type) and the result will be of the same type.
 * After exec() the modified displacements are stored in result of base class.
 */
void
Bend::execute(){
	int	ndispl = m_displ.size();
	ndispl = std::min(ndispl, int(m_coords.size()));
	for (int j=0; j<3; j++){
		for (int i=0; i<ndispl; i++){
			for (int z=0; z<3; z++){
				if (m_degree[j][z] > 0){
					for (int k=0; k<m_degree[j][z]+1; k++){
						m_displ[i][j] += pow(m_coords[i][z],(double)k)*m_coeffs[j][z][k];
					}
				}
			}
		}
	}
	return;
};
