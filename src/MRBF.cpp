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

// IMPLEMENTATION OF MRBF ***********************************************//
/*
 *	\date			26/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *
 *	\brief Radial Basis Function from point clouds.
 *
 */

/*! Basic Constructor. Doing nothing.*/
MRBF::MRBF(){
	m_tol = 0.00001;
	m_name = "MiMMO.MRBF";
};


/*! Destructor */
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

void MRBF::setTol(double tol){
	m_tol = tol;
};
void MRBF::addNodes(dvecarr3E nodes){
	for (auto &node : nodes){
		RBF::addNode(node);
	}
};

void MRBF::addNodes(MimmoObject* geometry){
	int nv = geometry->getNVertex();
	dvecarr3E vertex = geometry->getVertex();
	for (auto &node : vertex){
		RBF::addNode(node);
	}
};


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

//execute deformation methods
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


