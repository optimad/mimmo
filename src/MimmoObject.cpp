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
#include "MimmoObject.hpp"
#include "Operators.hpp"

using namespace std;

MimmoObject::MimmoObject(int type){
	m_type = type;
	m_geometry = NULL;
	m_internalPatch = false;
}

MimmoObject::MimmoObject(int type, dvecarr3E & vertex, ivector2D * connectivity){
	m_type = type;
	m_internalPatch = true;
//	if (m_type){
//		m_geometry = new Voltri();
//	}else{
//		m_geometry = new Surftri();
//	}
	setVertex(vertex);
//	if (connectivity != NULL)
//		setConnectivity((*connectivity));

}

MimmoObject::MimmoObject(int type, Patch* geometry){
	m_type = type;
	m_geometry = geometry;
	m_internalPatch = false;
}

MimmoObject::~MimmoObject(){
	clear();
};

void
MimmoObject::clear(){
	if (m_internalPatch){
		delete m_geometry;
	}
	m_geometry = NULL;
};

bool
MimmoObject::isEmpty(){
	return (m_geometry == NULL);
};

int
MimmoObject::getType(){
	return m_type;
};


long
MimmoObject::getNVertex(){
	return m_geometry->getVertexCount();
};

long
MimmoObject::getNCells(){
	return m_geometry->getCellCount();
};

dvecarr3E
MimmoObject::getVertex(){
	long nv = getNVertex();
	dvecarr3E result(nv);
	for (int i=0; i<nv; i++){
		result[i] = m_geometry->getVertexCoords(i);
	}
	return result;
};

//ivector2D*
//MimmoObject::getConnectivity(){
//};

Patch*
MimmoObject::getGeometry(){
	return m_geometry;
};

bool
MimmoObject::setVertex(dvecarr3E & vertex){
	if (m_geometry == NULL) return false;
	long nv = vertex.size();
	for (int i=0; i<nv; i++){
		m_geometry->addVertex(index);
		m_geometry->getvertex(index).setCoords(vertex);
	}
	return true;

};

bool
MimmoObject::setVertex(int index, darray3E & vertex){
	if (m_geometry == NULL) return false;
	m_geometry->addVertex(index);
	m_geometry->getvertex(index).setCoords(vertex);
	return true;
};

void
MimmoObject::setConnectivity(ivector2D & connectivity){
	if (m_geometry == NULL) return false;
};

bool
MimmoObject::setGeometry(int type, Patch* geometry){
	if (geometry == NULL) return false;
	m_geometry = geometry;
	m_type = type;
};






