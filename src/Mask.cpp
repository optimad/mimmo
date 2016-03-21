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
#include "Mask.hpp"

///*!Default constructor of Mask
// */
Mask::Mask(){};
//
///*!Default destructor of Mask
// */
Mask::~Mask(){};

///*!Copy constructor of Mask.
// */
Mask::Mask(const Mask & other):BaseManipulation(other){};

/*!Assignement operator of Mask.
 */
Mask & Mask::operator=(const Mask & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	return(*this);
};

/*!It sets the coordinates of the degrees of freedom.
 * \param[in] coords Coordinates of the degrees of freedom.
 */
void
Mask::setCoords(dvecarr3E & coords){
	m_coords = coords;
};

/*!It sets the limit of coordinates to apply the masking.
 * \param[in] thres Limit of coordinates to apply the masking.
 */
void
Mask::setThresholds(dmatrix32E & thres){
	m_thres = thres;
};

/*!It sets the limit of one coordinate to apply the masking.
 * \param[in] thres Limit of coordinate to apply the masking.
 * \param[in] dir Index of component.
 */
void
Mask::setThresholds(darray2E & thres, int dir){
	m_thres[dir] = thres;
};

/*!It sets the limit of x-coordinate to apply the masking.
 * \param[in] thres Limit of coordinate to apply the masking.
 */
void
Mask::setThresholdx(darray2E & thres){
	m_thres[0] = thres;
};

/*!It sets the limit of y-coordinate to apply the masking.
 * \param[in] thres Limit of coordinate to apply the masking.
 */
void
Mask::setThresholdy(darray2E & thres){
	m_thres[1] = thres;
};

/*!It sets the limit of z-coordinate to apply the masking.
 * \param[in] thres Limit of coordinate to apply the masking.
 */
void
Mask::setThresholdz(darray2E & thres){
	m_thres[2] = thres;
};

/*!It sets the condition to apply the mask
 * (true/false to set to zero the displacements inside/outside the thresholds).
 * \param[in] forward Condition to apply the mask for all the components.
 */
void
Mask::setInside(bool inside){
	for (int i=0; i<3; i++){
		m_inside[i] = inside;
	}
};

/*!It sets the condition to apply the mask
 * (true/false to set to zero the displacements inside/outside the thresholds).
 * \param[in] i Index of component.
 * \param[in] forward Condition to apply the mask for i-th component.
 */
void
Mask::setInside(int i, bool inside){
	if (i >= 0 && i < 3) m_inside[i] = inside;
};

/*!Execution command. It modifies the displacements given by the input manipulation object
 * with the masking conditions.
 * The input has to be set with a dvecarr3E variable (mask it casts the template method getInput
 * to this type) and the result will be of the same type.
 * After exec() the original displacements will be permanently modified.
 */
void
Mask::execute(){
	dvecarr3E displ = *(getInput<dvecarr3E>());
	int	ndispl;
	for (int i=0; i<ndispl; i++){
		if (m_coords[i][0]>m_thres[0][0] && m_coords[i][0]<m_thres[0][1] &&
				m_coords[i][1]>m_thres[1][0] && m_coords[i][1]<m_thres[1][1] &&
				m_coords[i][2]>m_thres[2][0] && m_coords[i][2]<m_thres[2][1]){
			for (int j=0; j<3; j++){
				displ[i][j] = (1-m_inside[j])*displ[i][j];
			}
		}
		else{
			for (int j=0; j<3; j++){
				displ[i][j] = (m_inside[j])*displ[i][j];
			}
		}
	}
	setResult(displ);
	return;
};
