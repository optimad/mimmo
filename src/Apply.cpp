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
#include "Apply.hpp"

/*!Default constructor of Apply
 */
Apply::Apply():BaseManipulation(){};

/*! Custom constructor of Apply.
 * It sets the linked geometry and the linked parent of the apply object.
 * \param[in] geometry MiMMO object with the geometry object of the deformation.
 * \param[in] parent Manipulation object linked by the apply object.
 */
Apply::Apply(MimmoObject* geometry, BaseManipulation* child):BaseManipulation(geometry, child){};

/*! Custom constructor of Apply.
 * It sets the linked parent of the apply object.
 * \param[in] parent Manipulation object linked by the apply object.
 */
Apply::Apply(BaseManipulation* parent):BaseManipulation(parent){};

/*!Default destructor of Apply
 */
Apply::~Apply(){};

/*!Copy constructor of Apply.
 */
Apply::Apply(const Apply & other){
	*this = other;
};

/*!Assignement operator of Apply.
 */
Apply & Apply::operator=(const Apply & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
};

/*!Execution command. It applies the deformation given by the parent manipulation
 * to the linked geometry. After exec() the original geometry will be permanently modified.
 */
void
Apply::execute(){
	if (m_geometry == NULL) return;
	dvecarr3E vertex = m_geometry->getVertex();
	long nv = m_geometry->getNVertex();
	for (long i=0; i<nv; i++){
		vertex[i] += m_displ[i];
		m_geometry->modifyVertex(vertex[i], i);
	}
	return;
};
