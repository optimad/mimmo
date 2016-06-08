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

using namespace mimmo;

/*!Default constructor of Mask
*/
Mask::Mask(){
	m_name = "MiMMO.Mask";
	m_thres.fill({{0,0}});
	m_inside = {{true, true, true}};
};

/*!Default destructor of Mask
 */
Mask::~Mask(){};

/*!Copy constructor of Mask.
 */
Mask::Mask(const Mask & other):BaseManipulation(other){};

/*!Assignement operator of Mask.
 */
Mask & Mask::operator=(const Mask & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	return(*this);
};

/*! It builds the input/output ports of the object
 */
void
Mask::buildPorts(){
	bool built = true;
	built = (built && createPortIn<dvecarr3E, Mask>(&m_coords, COORDS, 0, {DISPLS, GDISPLS, GLOBAL, LOCAL}));
	built = (built && createPortIn<dvecarr3E, Mask>(&m_displ, DISPLS, 10, {GDISPLS}));
	built = (built && createPortIn<dmatrix32E, Mask>(&m_thres, RANGE, 41));
	built = (built && createPortIn<std::array<bool,3>, Mask>(&m_inside, BOOLS3, 42));
	built = (built && createPortOut<dvecarr3E, Mask>(this, &mimmo::Mask::getCoords, COORDS, 0));
	built = (built && createPortOut<dvecarr3E, Mask>(this, &mimmo::Mask::getDisplacements, DISPLS, 10));
	m_arePortsBuilt = built;
};

/*!It gets the coordinates of points stored in the object.
 * \return Coordinates of points stored in the object.
 */
dvecarr3E
Mask::getCoords(){
	return(m_coords);
};

/*!It gets the displacements of points stored in the object.
 * \return Displacements of points stored in the object.
 */
dvecarr3E
Mask::getDisplacements(){
	return(m_displ);
};

/*!It sets the coordinates of points used by the masking.
 * \param[in] coords Coordinates of points used by the masking.
 */
void
Mask::setCoords(dvecarr3E coords){
	m_coords = coords;
};

/*!It sets the limits of coordinates to apply the masking.
 * \param[in] thres Limits of coordinates to apply the masking.
 */
void
Mask::setThresholds(dmatrix32E thres){
	m_thres = thres;
};

/*!It sets the lower limits of the three coordinates to apply the masking.
 * \param[in] thres Minimum limits of coordinates to apply the masking.
 */
void
Mask::setMinThresholds(darray3E thres){
	for (int dir = 0; dir < 3; ++dir) m_thres[dir][0] = thres[dir];
};

/*!It sets the greater limits of the three coordinates to apply the masking.
 * \param[in] thres Maximum limits of coordinates to apply the masking.
 */
void
Mask::setMaxThresholds(darray3E thres){
	for (int dir = 0; dir < 3; ++dir) m_thres[dir][1] = thres[dir];
};

/*!It sets the limits of one coordinate to apply the masking.
 * \param[in] thres Limits of coordinate to apply the masking.
 * \param[in] dir Index of component.
 */
void
Mask::setThresholds(darray2E thres, int dir){
	m_thres[dir] = thres;
};

/*!It sets the limits of x-coordinate to apply the masking.
 * \param[in] thres Limits of x-coordinate to apply the masking.
 */
void
Mask::setThresholdx(darray2E thres){
	m_thres[0] = thres;
};

/*!It sets the limits of y-coordinate to apply the masking.
 * \param[in] thres Limits of y-coordinate to apply the masking.
 */
void
Mask::setThresholdy(darray2E thres){
	m_thres[1] = thres;
};

/*!It sets the limits of z-coordinate to apply the masking.
 * \param[in] thres Limits of z-coordinate to apply the masking.
 */
void
Mask::setThresholdz(darray2E thres){
	m_thres[2] = thres;
};

/*!It sets the condition to apply the masking
 * (true/false to set to zero the displacements inside/outside the thresholds).
 * \param[in] inside Condition to apply the mask for all the components.
 */
void
Mask::setInside(bool inside){
	for (int i=0; i<3; i++){
		m_inside[i] = inside;
	}
};

/*!It sets the condition to apply the masking
 * (true/false to set to zero the displacements inside/outside the thresholds).
 * \param[in] inside Condition to apply the mask for the three components.
 */
void
Mask::setInside(std::array<bool,3> inside){
	for (int i=0; i<3; i++){
		m_inside[i] = inside[i];
	}
};

/*!It sets the condition to apply the mask
 * (true/false to set to zero the displacements inside/outside the thresholds).
 * \param[in] i Index of component.
 * \param[in] inside Condition to apply the mask for i-th component.
 */
void
Mask::setInside(int i, bool inside){
	if (i >= 0 && i < 3) m_inside[i] = inside;
};

/*!Execution command. It modifies the initial values of m_displ with the masking conditions.
 * After exec() the modified values are stored in the m_displ member.
 */
void
Mask::execute(){
	int	ndispl = m_displ.size();
	ndispl = std::min(ndispl, int(m_coords.size()));
	for (int i=0; i<ndispl; i++){
		if (m_coords[i][0]>m_thres[0][0] && m_coords[i][0]<m_thres[0][1] &&
				m_coords[i][1]>m_thres[1][0] && m_coords[i][1]<m_thres[1][1] &&
				m_coords[i][2]>m_thres[2][0] && m_coords[i][2]<m_thres[2][1]){
			for (int j=0; j<3; j++){
				m_displ[i][j] = (1-m_inside[j])*m_displ[i][j];
			}
		}
		else{
			for (int j=0; j<3; j++){
				m_displ[i][j] = (m_inside[j])*m_displ[i][j];
			}
		}
	}
	return;
};
