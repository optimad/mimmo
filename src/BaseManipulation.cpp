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
	m_geometry 		= NULL;
	m_pinType		= PinsType::BOTH;
	m_name			= "MiMMO";
};

/*!Default destructor of BaseManipulation.
 */
BaseManipulation::~BaseManipulation(){
	clear();
};

/*!Copy constructor of BaseManipulation.
 */
BaseManipulation::BaseManipulation(const BaseManipulation & other){
	*this = other;
};

/*!Assignement operator of BaseManipulation.
 */
BaseManipulation & BaseManipulation::operator=(const BaseManipulation & other){
	m_geometry 		= other.m_geometry;
	m_parent 		= other.m_parent;
	m_child 		= other.m_child;
	m_pinType		= other.m_pinType;
	m_pinIn			= other.m_pinIn;
	m_pinOut		= other.m_pinOut;
	m_name 			= other.m_name;
	//input and result are not copied (unique pointer of template memebers)
	return (*this);
};

/*!It gets the name of the manipulator object.
 * \return Name of the manipulator object.
 */
string
BaseManipulation::getName(){
	return m_name;
};

/*!It gets the geometry linked by the manipulator object.
 * \return Pointer to geometry to be deformed by the manipulator object.
 */
MimmoObject*
BaseManipulation::getGeometry(){
	return m_geometry;
};

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
 * \return false if not found.
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
 * \return false if not found.
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

PinsType
BaseManipulation::getPinType(){
	return (m_pinType);
}

int
BaseManipulation::getNPinsIn(){
	return (m_pinIn.size());
}

int
BaseManipulation::getNPinsOut(){
	return (m_pinOut.size());
}

/*!It sets the name of the manipulator object.
 * \param[in] name Name of the manipulator object.
 */
void
BaseManipulation::setName(string name){
	m_name = name;
};

/*!It sets the geometry linked by the manipulator object.
 * \param[in] geometry Pointer to geometry to be deformed by the manipulator object.
 */
void
BaseManipulation::setGeometry(MimmoObject* geometry){
	m_geometry = geometry;
};

/*!It clears the pointer to the geometry linked by the object.
 */
void
BaseManipulation::unsetGeometry(){
	m_geometry = NULL;
};

void
BaseManipulation::removePins(){
	removePinsIn();
	removePinsOut();
}

void
BaseManipulation::removePinsIn(){
	unordered_map<BaseManipulation*, int>::iterator it;
	//Warning!! If infinite while unset parent wrong
	while(m_parent.size()){
		it = m_parent.begin();
		mimmo::pin::removeAllPins(it->first, this);
	}
}

void
BaseManipulation::removePinsOut(){
	unordered_map<BaseManipulation*, int>::iterator it;
	//Warning!! If infinite while unset parent wrong
	while(m_child.size()){
		it = m_child.begin();
		mimmo::pin::removeAllPins(this, it->first);
	}
}

void
BaseManipulation::clearInput(){
	m_input.release();
}

void
BaseManipulation::clearResult(){
	m_result.release();
}

/*!It clears the object, by setting to zero/NULL each member/pointer in the object.
 */
void
BaseManipulation::clear(){
	unsetGeometry();
	removePins();
	clearInput();
	clearResult();
};

/*!Execution command. exec() runs the execution of output pins at the end of the execution.
 * execute is pure virtual and it has to be implemented in a derived class.
 * Temporary Inputs are cleared from class
 */
void
BaseManipulation::exec(){
	execute();
	for (int i=0; i<m_pinOut.size(); i++){
		if (m_pinOut[i]->getLink() != NULL){
			m_pinOut[i]->exec();
		}
	}
	clearInput();
}

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
	if(!m_child.count(child)){
		m_child.insert(pair<BaseManipulation*,int>(child,1)); //add new child with counter 1;
	}else{
		m_child[child]++; //just incrementing pre-existent child counter;
	}	
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
		if(m_parent[parent] <1)	m_parent.erase(parent);
	}
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
		if(m_child[child] <1) m_child.erase(child);
	}
};

vector<InOut*>
BaseManipulation::getPinsIn(){
	return (m_pinIn);
}

std::vector<InOut*>
BaseManipulation::getPinsOut(){
	return (m_pinOut);
}

int
BaseManipulation::findPinIn(InOut& pin){
	for (int i=0; i<m_pinIn.size(); i++){
		if (pin == *(m_pinIn[i])) return(i);
	}
	return(-1);
}

int
BaseManipulation::findPinOut(InOut& pin){
	for (int i=0; i<m_pinOut.size(); i++){
		if (pin == *(m_pinOut[i])) return(i);
	}
	return(-1);
}

void
BaseManipulation::removePinIn(int i){
	if (i<m_pinIn.size() && i != -1){
		delete m_pinIn[i];
		m_pinIn[i] = NULL;
		m_pinIn.erase(m_pinIn.begin()+i);
	}
}

void
BaseManipulation::removePinOut(int i){
	if (i<m_pinOut.size() && i != -1){
		delete m_pinOut[i];
		m_pinOut[i] = NULL;
		m_pinOut.erase(m_pinOut.begin()+i);
	}
}

