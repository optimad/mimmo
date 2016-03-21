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
#include <iterator>
#include <utility>
#include "BaseManipulation.hpp"

using namespace std;

/*!Default constructor of BaseManipulation.
 * It sets to zero/null each member/pointer.
 */
BaseManipulation::BaseManipulation(){
	m_ndeg 			= 0;
	m_geometry 		= NULL;
	m_relInfo 		= false;
	m_parent.clear();
	m_child.clear();
	m_displ.clear();
	m_result		= NULL;
};

/*!Default destructor of BaseManipulation.
 */
BaseManipulation::~BaseManipulation(){
	clear();
	m_ndeg = 0;
};

/*!Copy constructor of BaseManipulation.
 */
BaseManipulation::BaseManipulation(const BaseManipulation & other){
	m_ndeg 			= other.m_ndeg;
	m_displ 		= other.m_displ;
	m_child 		= other.m_child;
	m_parent 		= other.m_parent;
	m_geometry 		= other.m_geometry;
	m_relInfo 		= other.m_relInfo;
};

/*!Assignement operator of BaseManipulation.
 */
BaseManipulation & BaseManipulation::operator=(const BaseManipulation & other){
	m_ndeg 			= other.m_ndeg;
	m_displ 		= other.m_displ;
	m_parent 		= other.m_parent;
	m_child 		= other.m_child;
	m_geometry 		= other.m_geometry;
	m_relInfo 		= other.m_relInfo;
};

/*!It gets the number of degrees of freedom of the manipulator object.
 * \return #Degrees of freedom.
 */
int
BaseManipulation::getNDeg(){
	return m_ndeg;
};

///*!It gets the number of degrees of freedom of a child of the manipulator object.
// * \param[in] i Index of target child.
// * \return #Degrees of freedom of ìtarget child.
// */
//int
//BaseManipulation::getNDegOut(int i){
//	if (i>m_child.size()-1) return 0;
//	return m_child[i]->getNDeg();
//};

/*!It gets the displacement of the degree of freedom currently stored in the object.
 * \return Displacements of the degrees of freedom.
 */
dvecarr3E&
BaseManipulation::getDisplacements(){
	return m_displ;
};

///*!It gets the displacement of the degree of freedom currently stored in a child of the object.
// * \param[in] i Index of target child.
// * \return Displacements of the degrees of freedom of the target child.
// */
//dvecarr3E*
//BaseManipulation::getDisplacementsOut(int i){
//	if (i>m_child.size()-1) return NULL;
//	return &m_child[i]->getDisplacements();
//};

/*!It gets the number of parent linked to the manipulator object.
 * \return #Parent.
 */
int
BaseManipulation::getNParent(){
	return m_parent.size();
};

/*!It gets the manipulator object linked by this object.
 * \param[in] i Index of target parent.
 * \return Pointer to i-th parent manipulator object.
 */
BaseManipulation*
BaseManipulation::getParent(int i){
	if (i>m_parent.size()-1) return NULL;
	return next(m_parent.begin(), i)->first;
};

/*! Return true if the target is contained in the parent list.
 * \param[in] target BaseManipulation object
 * \param[out] index return actual position in the list(may vary if the parent list is modified). If not in the list return -1.
 * \param[out] result boolean flag.
 */
bool
BaseManipulation::isParent(BaseManipulation * target, int index){
	unordered_map<BaseManipulation *, int>::iterator it;
	it = m_parent.find(target);
	index = -1;
	if(it == m_parent.end()) return false;
	
	index = distance(m_parent.begin(), it);
	return true;
};


/*!It gets the number of children linked to the manipulator object.
 * \return #Children.
 */
int
BaseManipulation::getNChild(){
	return m_child.size();
};

/*!It gets one child object linked by this object.
 * \param[in] i Index of target child.
 * \return Pointer to i-th child manipulator object.
 */
BaseManipulation*
BaseManipulation::getChild(int i){
	if (i>m_child.size()-1) return NULL;
	return next(m_child.begin(), i)->first;
};


/*! Return true if the target is contained in the child list.
 * \param[in] target BaseManipulation object
 * \param[out] index return actual position in the list(may vary if the parent list is modified). If not in the list return -1.
 * \param[out] result boolean flag.
 */
bool
BaseManipulation::isChild(BaseManipulation * target, int index){
	unordered_map<BaseManipulation *, int>::iterator it;
	it = m_child.find(target);
	index = -1;
	if(it == m_child.end()) return false;
	
	index = distance(m_child.begin(), it);
	return true;
};

/*!It gets the geometry linked by the manipulator object.
 * \return Pointer to geometry to be deformed by the manipulator object.
 */
MimmoObject*
BaseManipulation::getGeometry(){
	return m_geometry;
};

/*!It gets if the object is a "release Info" object.
 * \return Is the object a "release Info" object?
 */
bool
BaseManipulation::getReleaseInfo(){
	return m_relInfo;
};

/*!It gets the pointer to the info release by a "release Info" object of the chain.
 * \return Pointer to Info structure.
 */
Info*
BaseManipulation::getInfo(){
	return m_info;
};

/*!It sets the number of degrees of freedom of the manipulator object.
 * \param[in] ndeg #Degrees of freedom.
 */
void
BaseManipulation::setNDeg(int ndeg){
	m_ndeg = ndeg;
	m_displ.resize(m_ndeg);
};

///*!It sets the number of degrees of freedom of a child of the manipulator object.
// * \param[in] i Index of target child.
// * \param[in] ndeg #Degrees of freedom to set in the target child.
// */
//void
//BaseManipulation::setNDegOut(int i, int ndeg){
//	if (i>m_child.size()-1) return;
//	m_child[i]->setNDeg(ndeg);
//};

/*!It sets the displacement of the degree of freedom currently stored in the object.
 * \param[in] displacements Displacements of the degrees of freedom.
 */
void
BaseManipulation::setDisplacements(dvecarr3E displacements){
	m_displ = displacements;
	m_ndeg = m_displ.size();
};

///*!It sets the displacement of the degree of freedom of a child of the object.
// * \param[in] i Index of target child.
// * \param[in] displacements Displacements of the degrees of freedom to set in the target child.
// */
//void
//BaseManipulation::setDisplacementsOut(int i, dvecarr3E & displacements){
//	if (i>m_child.size()-1) return;
//	m_child[i]->setDisplacements(displacements);
//};

/*!It sets the geometry linked by the manipulator object.
 * \param[in] geometry Pointer to geometry to be deformed by the manipulator object.
 */
void
BaseManipulation::setGeometry(MimmoObject* geometry){
	m_geometry = geometry;
};

/*!It sets if the object is a "release Info" object.
 * \param[in] Is the object a "release Info" object?
 */
void
BaseManipulation::setReleaseInfo(bool flag){
	m_relInfo = flag;
};

/*!It clears the displacement of the degree of freedom currently stored in the object.
 */
void
BaseManipulation::clearDisplacements(){
	m_displ.clear();
};

///*!It clears the displacements of the degrees of freedom of the children of the object.
// */
//void
//BaseManipulation::clearDisplacementsOut(){
//	for (int i=0; m_child.size(); i++)	m_child[i]->clearDisplacements();
//};
//
///*!It clears the displacements of the degrees of freedom of a child of the object.
// * \param[in] i Index of target child.
// */
//void
//BaseManipulation::clearDisplacementsOut(int i){
//	if (i>m_child.size()-1) return;
//	m_child[i]->clearDisplacements();
//};

/*!It clears the pointer to the geometry linked by the object.
 */
void
BaseManipulation::unsetGeometry(){
	m_geometry = NULL;
};

void
BaseManipulation::clearResult(){
	m_result.release();
}

void
BaseManipulation::clearInput(){
	m_input.release();
}

/*!It clears the object, by setting to zero/NULL each member/pointer in the object.
 */
void
BaseManipulation::clear(){
	unsetParent();
	unsetChild();
	unsetGeometry();
	clearDisplacements();
	m_relInfo = false;
};

/*!It adds a manipulator object linked by this object.
 * \param[in] parent Pointer to parent manipulator object.
 */
void
BaseManipulation::addParent(BaseManipulation* parent){
	
	if(!m_parent.count(parent)){
		m_parent.insert(pair<BaseManipulation*,int>(parent,1)); //add new parent with counter 1;
	}else{
		m_parent[parent]++; //just incrementing pre-existent parent counter;
	}	
};

/*!It adds a child manipulator object to the children linked by this object.
 * \param[in] child Pointer to child manipulator object.
 */
void
BaseManipulation::addChild(BaseManipulation* child){
	if(!m_child.count(parent)){
		m_child.insert(pair<BaseManipulation*,int>(child,1)); //add new child with counter 1;
	}else{
		m_child[child]++; //just incrementing pre-existent child counter;
	}	
};

/*!It clears the current father object linked by this object.
 * It object is not linked doing nothing.
 * \param[in] parent pointer to BaseManipulation object 
 */
void
BaseManipulation::removeParent(BaseManipulation * parent){
	
	unordered_map<BaseManipulation*, int>::iterator got = m_parent.find(parent);
	if(got != m_parent.end())
		m_parent.erase();
	//need to clear mutual connection
	
};

/*!It clears the current child object linked by this object.
 * It object is not linked doing nothing.
 * \param[in] child pointer to BaseManipulation object  
 */
void
BaseManipulation::removeChild(BaseManipulation * child){
	unordered_map<BaseManipulation*, int>::iterator got = m_child.find(child);
	if(got != m_child.end())
		m_child.erase();
	
	//need to clear mutual connection;
};

/*!It clears all father objects linked by this object. */
void
BaseManipulation::removeAllParent(){
	m_parent.clear();
	//need to clear all mutual connection involving this object
};

/*!It clears all child objects linked by this object.*/
void
BaseManipulation::removeAllChild(){
	m_child.clear();
	//need to clear all mutual connection involving this object
};


/*! Decrement target parent multiplicity, contained in member m_parent.
 * If multiplicity is zero, erase target from list. The method is meant to be used in conjuction 
 * to manual cut off of object pins. 
 * \param[in] parent pointer to BaseManipulation object 
 */
void
BaseManipulation::unsetParent(BaseManipulation * parent){
	
	unordered_map<BaseManipulation*, int>::iterator got = m_parent.find(parent);
	if(got != m_parent.end()){
		m_parent[parent]--;
	}
	if(m_parent[parent] <1)	m_parent.erase();
};

/*! Decrement target child multiplicity, contained in member m_child.
 * If multiplicity is zero, erase target from list. The method is meant to be used in conjuction 
 * to manual cut off of object pins.
 * \param[in] child pointer to BaseManipulation object  
 */
void
BaseManipulation::unsetChild(BaseManipulation * child){
	unordered_map<BaseManipulation*, int>::iterator got = m_child.find(child);
	if(got != m_child.end()){
		m_child[child]--;
	}
	if(m_child[child] <1)	m_child.erase();
};


void
BaseManipulation::removePins(){
	removePinsIn();
	removePinsOut();
}

void
BaseManipulation::removePinsIn(){
	for (int i=0; i<m_pinIn.size(); i++){
		delete m_pinIn[i];
		m_pinIn[i] = NULL;
	}
	m_pinIn.clear();
}

void
BaseManipulation::removePinsOut(){
	for (int i=0; i<m_pinOut.size(); i++){
		delete m_pinOut[i];
		m_pinOut[i] = NULL;
	}
	m_pinOut.clear();
}

void
BaseManipulation::removePinIn(int i){
	if (i<m_pinIn.size()){
		delete m_pinIn[i];
		m_pinIn[i] = NULL;
		m_pinIn.erase(m_pinIn.begin()+i);
	}
}

void
BaseManipulation::removePinOut(int i){
	if (i<m_pinOut.size()){
		delete m_pinOut[i];
		m_pinOut[i] = NULL;
		m_pinOut.erase(m_pinOut.begin()+i);
	}
}

/*!It releases the info, i.e. it creates an Info structure and
 * it sets it by setInfo() method.
 */
void
BaseManipulation::releaseInfo(){
	m_info = new Info();
	setInfo();
}

/*!It recover the info from the Info structure released by a "release Info" object
 * in the chain. It travels through the chain until it finds a not NULL pointer to info
 * or a "release Info" object.
 */
Info*
BaseManipulation::recoverInfo(){
	bool found = false;
	int ich = 0;
	Info* pointInfo = NULL;
	BaseManipulation* pointer;
	if (getReleaseInfo()){
		releaseInfo();
		pointInfo = getInfo();
	}
	while(pointInfo == NULL && ich<m_child.size()){
		pointer = m_child[ich];
		pointInfo = pointer->recoverInfo();
		ich++;
	}
	return pointInfo;
}

/*!It sets the info structure. Virtual method to be eventually
 * implemented in derived manipulation objects.
 * Here it is a void method, i.e. it does nothing.
 */
void
BaseManipulation::setInfo(){
}

/*!It uses the info structure. Virtual method to be eventually
 * implemented in derived manipulation objects.
 * Here it is a void method, i.e. it does nothing.
 */
void
BaseManipulation::useInfo(){
}

/*!Execution command. It runs recoverDisplacements applyFilters and execute.
 * execute is pure virtual and it has to be implemented in a derived class.
 * Temporary Inputs are cleared from class
 */
void
BaseManipulation::exec(){
	m_info = recoverInfo();
	if (m_info != NULL) useInfo();
	execute();
	for (int i=0; i<m_pinOut.size(); i++){
		if (m_pinOut[i]->getLink() != NULL){
			m_pinOut[i]->exec();
		}
	}
	clearInput();
}


//==============================//
//EXTERNAL METHODS				//
//==============================//



