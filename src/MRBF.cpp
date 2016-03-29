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
 \ *---------------------------------------------------------------------------*/

#include "MRBF.hpp"

using namespace std;


/*! Default Constructor.*/
MRBF::MRBF(){
	m_tol = 0.00001;
	m_name = "MiMMO.MRBF";
};


/*! Default Destructor */
MRBF::~MRBF(){};

/*! Copy Constructor
 *\param[in] other MRBF where copy from
 */
MRBF::MRBF(const MRBF & other){
	*this = other;
};

/*! Copy Operator
 * \param[in] other MRBF where copy from
 */
MRBF & MRBF::operator=(const MRBF & other){
	m_tol = other.m_tol;
	return(*this);
};

/*!It sets the tolerance used in greedy algorithm of bitpit::RBF class to choose
 * a sub-set of the control points.
 * \param[in] tol Tolerance used in greedy algorithm.
 */
void MRBF::setTol(double tol){
	m_tol = tol;
};

/*!It adds a set of points to the control points used in RBF execution.
 * \param[in] nodes Coordinates of control points.
 */
void MRBF::addNodes(dvecarr3E nodes){
	for (auto &node : nodes){
		RBF::addNode(node);
	}
};

/*!It adds a set of points to the control points used in RBF execution by extracting
 * the vertices stored in a MimmoObject container.
 * \param[in] geometry Pointer to MimmoObject that contains the geometry.
 */
void MRBF::addNodes(MimmoObject* geometry){
	int nv = geometry->getNVertex();
	dvecarr3E vertex = geometry->getVertex();
	for (auto &node : vertex){
		RBF::addNode(node);
	}
};

/*!It adds a field of values on the control points to be interpolated by RBF techniques (in this library it is
 * supposed to be a displacements field).
 * \param[in] field Field of displacements of control points.
 */
void MRBF::addField(dvecarr3E field){
	int np = field.size();
	dvector1D f(np);
	for (int i=0; i<3; i++){
		for (int j=0; j<np; j++){
			f[j] = field[j][i];
		}
		RBF::addField(f);
	}
	double maxdispl = 0.0;
	for (int j=0; j<np; j++){
		maxdispl = max(maxdispl, norm2(field[j]));
	}
	setSupportRadius(3*maxdispl);

};

/*!Execution of RBF object. It evaluates the displacements (values) over the point of the
 * linked geometry, given as result of RBF technique implemented in bitpit::RBF base class.
 * The result is stored in the result member of BaseManipulation base class.
 *
 */
void MRBF::execute(){

	greedy(m_tol);

	MimmoObject * container = getGeometry();
	if(container == NULL ) return;

	int nv = container->getNVertex();
	dvecarr3E vertex = container->getVertex();

	dvecarr3E result(nv, darray3E{0,0,0});
	dvector1D displ;
	for(int i=0; i<nv; ++i){
		displ = RBF::evalRBF(vertex[i]);
		for (int j=0; j<3; j++) result[i][j] = displ[j];
	}

	setResult(result);

};


