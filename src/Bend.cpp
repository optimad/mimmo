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

///*!Default constructor of Bend
// */
Bend::Bend(){};
//
///*!Default destructor of Bend
// */
Bend::~Bend(){};

///*!Copy constructor of Bend.
// */
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
};

/*!It sets the coordinates of the degrees of freedom.
 * \param[in] coords Coordinates of the degrees of freedom.
 */
void
Bend::setCoords(dvecarr3E & coords){
	m_coords = coords;
};

/*!It sets the degrees of polynomial law for each component of displacements of degrees of freedom.
 * \param[in] degrees Degrees of polynomial laws (degree[i][j] = degree of displacement function si = f(xj)).
 */
void
Bend::setDegree(dvecarr3E & degree){
	m_degree = degree;
};

/*!It sets the coefficients of the polynomial laws.
 * \param[in] coeffs Coefficients of the polynomial laws. (coeffs[i][j][k] = coefficients aijk of term si = aij * xj^k).
 */
void
Bend::setCoeffs(dvector3D & coeffs){
	m_coeffs = coeffs;
};

void
Bend::useInfo(){
	if (m_ndeg !=  m_info->m_coords.size() || m_info->m_naxes != 3){
		std::cout << "#degrees : " << m_ndeg << " - #coords : " << m_info->m_coords.size() << std::endl;
		std::cout << "Incoherent Size ---> end of process " << std::endl;
		exit(1001);
	}
	m_coords.resize(m_info->m_coords.size());
	for (int i=0; i<m_ndeg; i++){
		for (int j=0; j<3; j++){
			m_coords[i][j] = m_info->m_coords[i][j];
		}
	}
};

/*!Execution command. It modifies the displacements given by the child manipulation object
 * with the polynomial law. After exec() the original displacements will be permanently modified.
 */
void
Bend::execute(){
	for (int j=0; j<3; j++){
		for (int i=0; i<getNDeg(); i++){
			for (int z=0; z<3; z++){
				if (m_degree[j][z] > 0){
					for (int k=0; k<m_degree[j][z]+1; k++){
						m_displ[i][j] += pow(m_coords[i][z],(double)k)*m_coeffs[j][z][k];
					}
				}
			}
		}
	}
	for (int j=0; j<getNChild(); j++){
		setDisplacementsOut(j, m_displ);
	}
	return;
};
