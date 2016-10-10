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
#include "RotationGeometry.hpp"

using namespace mimmo;

/*!Default constructor of RotationGeometry
 */
RotationGeometry::RotationGeometry(darray3E origin, darray3E direction){
	m_origin = origin;
	m_direction = direction;
	m_name = "MiMMO.RotationGeometry";
};

/*!Default destructor of RotationGeometry
 */
RotationGeometry::~RotationGeometry(){};

/*!Copy constructor of RotationGeometry.
 */
RotationGeometry::RotationGeometry(const RotationGeometry & other):BaseManipulation(other){
	m_origin = other.m_origin;
	m_direction = other.m_direction;
};

/*!Assignement operator of RotationGeometry.
 */
RotationGeometry & RotationGeometry::operator=(const RotationGeometry & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_origin = other.m_origin;
	m_direction = other.m_direction;
	return(*this);
};


/*! It builds the input/output ports of the object
 */
void RotationGeometry::buildPorts(){
	bool built = true;
	built = (built && createPortIn<darray3E, RotationGeometry>(&m_origin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<darray3E, RotationGeometry>(&m_direction, PortType::M_AXIS, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<double, RotationGeometry>(&m_alpha, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<dvector1D, RotationGeometry>(this, &mimmo::RotationGeometry::setFilter, PortType::M_FILTER, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<MimmoObject*, RotationGeometry>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortOut<dvecarr3E, RotationGeometry>(this, &mimmo::RotationGeometry::getDisplacements, PortType::M_GDISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	m_arePortsBuilt = built;
};

/*!It sets the origin and direction of the rotation axis.
 * \param[in] origin Origin of rotation axis.
 * \param[in] direction Direction of rotation axis.
 */
void
RotationGeometry::setAxis(darray3E origin, darray3E direction){
	m_origin = origin;
	m_direction = direction;
}

/*!It sets the origin of the rotation axis.
 * \param[in] origin Origin of rotation axis.
 */
void
RotationGeometry::setOrigin(darray3E origin){
	m_origin = origin;
}

/*!It sets the direction of the rotation axis.
 * \param[in] direction Direction of rotation axis.
 */
void
RotationGeometry::setDirection(darray3E direction){
	m_direction = direction;
	double L = norm2(m_direction);
	for (int i=0; i<3; i++)
		m_direction[i] /= L;
}

/*!It sets the value of the rotation.
 * \param[in] alpha Value of rotation axis.
 */
void
RotationGeometry::setRotation(double alpha){
	m_alpha = alpha;
}

/*!It sets the filter field to modulate the displacements of the vertices
 * of the target geometry.
 * \param[in] filter Filter field defined on geometry vertices.
 */
void
RotationGeometry::setFilter(dvector1D filter){
	m_filter = filter;
}

/*!
 * Return actual computed displacements field (if any) for the geometry linked.
 * @return  The computed deformation field on the vertices of the linked geometry
 */
dvecarr3E
RotationGeometry::getDisplacements(){
    return m_displ;
};

/*!Execution command. It saves in "rot"-terms the modified axes and origin, by the
 * rotation conditions. This terms can be recovered and passed by a pin to a child object
 * by the related get-methods.
 */
void
RotationGeometry::execute(){

    m_displ.resize(m_geometry->getNVertex());

    //compute coefficients and constant vectors of rodriguez formula
    double a = cos(m_alpha);
    darray3E b =  (1 - cos(m_alpha)) * m_direction;
    double c = sin(m_alpha);

    darray3E point, rotated;
    long ID;
    int idx;
    liimap mapID = m_geometry->getMapDataInv();

    for (auto vertex : m_geometry->getVertices()){
        point = vertex.getCoords();
        ID = vertex.getId();
        idx = mapID[ID];

        point -= m_origin;
        //rodrigues formula
        rotated = a * point +
                  b * dotProduct(m_direction, point) +
                  c * crossProduct(m_direction, point);

        rotated += m_origin;
        point += m_origin;

        m_displ[idx] = (rotated-point)*m_filter[idx];

    }

    return;
};

