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
#include "ClipGeometry.hpp"

using namespace mimmo;

/*!Default constructor of ClipGeometry
*/
ClipGeometry::ClipGeometry(){
	m_name = "MiMMO.ClipGeometry";
	m_plane.fill({{0,0}});
	m_insideOut = false;
	m_patch.reset(nullptr);
};

/*!Default destructor of ClipGeometry
 */
ClipGeometry::~ClipGeometry(){};

/*!Copy constructor of ClipGeometry.
 */
ClipGeometry::ClipGeometry(const ClipGeometry & other){
	*this = other;
};

/*!Assignement operator of ClipGeometry.
 */
ClipGeometry & ClipGeometry::operator=(const ClipGeometry & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_plane = other.m_plane;
	m_insideout = other.m_insideout;
	return(*this);
};

/*! It builds the input/output ports of the object
 */
void
ClipGeometry::buildPorts(){
	bool built = true;
	built = (built && createPortIn<darray4E, ClipGeometry>(this, &mimmo::ClipGeometry::setClipPlane, PortType::M_PLANE, mimmo::pin::containerTAG::ARRAY4, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<bool, ClipGeometry>(this, &mimmo::ClipGeometry::setInsideOut, PortType::M_VALUEB, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::BOOL));
	built = (built && createPortIn<MimmoObject*, ClipGeometry>(this, &mimmo::ClipGeometry::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	built = (built && createPortOut<MimmoObject*, ClipGeometry>(this, &mimmo::ClipGeometry::getClippedPatch, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	m_arePortsBuilt = built;
};

/*!It gets the coordinates of points stored in the object.
 * \return Coordinates of points stored in the object.
 */
dvecarr3E
ClipGeometry::getCoords(){
	return(m_coords);
};

/*!It gets the displacements of points stored in the object.
 * \return Displacements of points stored in the object.
 */
dvecarr3E
ClipGeometry::getDisplacements(){
	return(m_displ);
};

/*!It sets the coordinates of points used by the masking.
 * \param[in] coords Coordinates of points used by the masking.
 */
void
ClipGeometry::setCoords(dvecarr3E coords){
	m_coords = coords;
};

/*!It sets the limits of coordinates to apply the masking.
 * \param[in] thres Limits of coordinates to apply the masking.
 */
void
ClipGeometry::setThresholds(dmatrix32E thres){
	m_thres = thres;
};

/*!It sets the lower limits of the three coordinates to apply the masking.
 * \param[in] thres Minimum limits of coordinates to apply the masking.
 */
void
ClipGeometry::setMinThresholds(darray3E thres){
	for (int dir = 0; dir < 3; ++dir) m_thres[dir][0] = thres[dir];
};

/*!It sets the greater limits of the three coordinates to apply the masking.
 * \param[in] thres Maximum limits of coordinates to apply the masking.
 */
void
ClipGeometry::setMaxThresholds(darray3E thres){
	for (int dir = 0; dir < 3; ++dir) m_thres[dir][1] = thres[dir];
};

/*!It sets the limits of one coordinate to apply the masking.
 * \param[in] thres Limits of coordinate to apply the masking.
 * \param[in] dir Index of component.
 */
void
ClipGeometry::setThresholds(darray2E thres, int dir){
	m_thres[dir] = thres;
};

/*!It sets the limits of x-coordinate to apply the masking.
 * \param[in] thres Limits of x-coordinate to apply the masking.
 */
void
ClipGeometry::setThresholdx(darray2E thres){
	m_thres[0] = thres;
};

/*!It sets the limits of y-coordinate to apply the masking.
 * \param[in] thres Limits of y-coordinate to apply the masking.
 */
void
ClipGeometry::setThresholdy(darray2E thres){
	m_thres[1] = thres;
};

/*!It sets the limits of z-coordinate to apply the masking.
 * \param[in] thres Limits of z-coordinate to apply the masking.
 */
void
ClipGeometry::setThresholdz(darray2E thres){
	m_thres[2] = thres;
};

/*!It sets the condition to apply the masking
 * (true/false to set to zero the displacements inside/outside the thresholds).
 * \param[in] inside Condition to apply the mask for all the components.
 */
void
ClipGeometry::setInside(bool inside){
	for (int i=0; i<3; i++){
		m_inside[i] = inside;
	}
};

/*!It sets the condition to apply the masking
 * (true/false to set to zero the displacements inside/outside the thresholds).
 * \param[in] inside Condition to apply the mask for the three components.
 */
void
ClipGeometry::setInside(std::array<bool,3> inside){
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
ClipGeometry::setInside(int i, bool inside){
	if (i >= 0 && i < 3) m_inside[i] = inside;
};

/*!Execution command. It modifies the initial values of m_displ with the masking conditions.
 * After exec() the modified values are stored in the m_displ member.
 */
void
ClipGeometry::execute(){
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
