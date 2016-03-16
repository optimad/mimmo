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
};

/*!Custom constructor of BaseManipulation.
 * \param[in] geometry Pointer to target geometry to be linked.
 * \param[in] child Pointer to reference manipulator object to be linked (default value = NULL).
 */
BaseManipulation::BaseManipulation(MimmoObject* geometry, BaseManipulation* child){
	m_ndeg 			= 0;
	m_geometry 		= geometry;
	m_relInfo 		= false;
	m_displ.clear();
	m_parent.clear();
	if (child != NULL){
		addChild(child);
	}
};

/*!Custom constructor of BaseManipulation.
 * \param[in] child Pointer to reference manipulator object to be linked.
 */
BaseManipulation::BaseManipulation(BaseManipulation* child){
	m_ndeg 			= 0;
	m_geometry 		= NULL;
	m_relInfo 		= false;
	m_displ.clear();
	m_parent.clear();
	if (child != NULL){
		addChild(child);
	}
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

/*!It gets the number of degrees of freedom of a child of the manipulator object.
 * \param[in] i Index of target child.
 * \return #Degrees of freedom of ìtarget child.
 */
int
BaseManipulation::getNDegOut(int i){
	if (i>m_child.size()-1) return 0;
	return m_child[i]->getNDeg();
};

/*!It gets the displacement of the degree of freedom currently stored in the object.
 * \return Displacements of the degrees of freedom.
 */
dvecarr3E&
BaseManipulation::getDisplacements(){
	return m_displ;
};

/*!It gets the displacement of the degree of freedom currently stored in a child of the object.
 * \param[in] i Index of target child.
 * \return Displacements of the degrees of freedom of the target child.
 */
dvecarr3E*
BaseManipulation::getDisplacementsOut(int i){
	if (i>m_child.size()-1) return NULL;
	return &m_child[i]->getDisplacements();
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
	return m_parent[i];
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
	return m_child[i];
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

/*!It sets the number of degrees of freedom of a child of the manipulator object.
 * \param[in] i Index of target child.
 * \param[in] ndeg #Degrees of freedom to set in the target child.
 */
void
BaseManipulation::setNDegOut(int i, int ndeg){
	if (i>m_child.size()-1) return;
	m_child[i]->setNDeg(ndeg);
};

/*!It sets the displacement of the degree of freedom currently stored in the object.
 * \param[in] displacements Displacements of the degrees of freedom.
 */
void
BaseManipulation::setDisplacements(dvecarr3E displacements){
	m_displ = displacements;
	m_ndeg = m_displ.size();
};

/*!It sets the displacement of the degree of freedom of a child of the object.
 * \param[in] i Index of target child.
 * \param[in] displacements Displacements of the degrees of freedom to set in the target child.
 */
void
BaseManipulation::setDisplacementsOut(int i, dvecarr3E & displacements){
	if (i>m_child.size()-1) return;
	m_child[i]->setDisplacements(displacements);
};

/*!It adds a manipulator object linked by this object.
 * \param[in] parent Pointer to parent manipulator object.
 */
void
BaseManipulation::addParent(BaseManipulation* parent){
	m_parent.push_back(parent);
};

/*!It adds a child manipulator object to the children linked by this object.
 * \param[in] child Pointer to child manipulator object.
 */
void
BaseManipulation::addChild(BaseManipulation* child){
	int count = m_child.size();
	m_child.push_back(child);
	m_child[count]->addParent(this);
};

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

/*!It clears the manipulator objects linked by this object.
 * It sets to NULL the pointers to parent manipulator objects.
 */
void
BaseManipulation::unsetParent(){
	for (int i=0; i<m_parent.size(); i++) m_parent[i] = NULL;
	m_parent.clear();
};

/*!It clears the children objects linked by this object.
 * It sets to NULL the pointers to children manipulator object and clear the vector.
 */
void
BaseManipulation::unsetChild(){
	for (int i=0; i<m_child.size(); i++) m_child[i] = NULL;
	m_child.clear();
};

/*!It clears the displacement of the degree of freedom currently stored in the object.
 */
void
BaseManipulation::clearDisplacements(){
	m_displ.clear();
};

/*!It clears the displacements of the degrees of freedom of the children of the object.
 */
void
BaseManipulation::clearDisplacementsOut(){
	for (int i=0; m_child.size(); i++)	m_child[i]->clearDisplacements();
};

/*!It clears the displacements of the degrees of freedom of a child of the object.
 * \param[in] i Index of target child.
 */
void
BaseManipulation::clearDisplacementsOut(int i){
	if (i>m_child.size()-1) return;
	m_child[i]->clearDisplacements();
};

/*!It clears the pointer to the geometry linked by the object.
 */
void
BaseManipulation::unsetGeometry(){
	m_geometry = NULL;
};

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
}

//============================//
// OVERLOADED PINS METHODS
//============================//

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<double&(void)> getVal, function<void(double)> setVal){
	InOutD* pin = new InOutD();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<double(void)> getVal, function<void(double)> setVal){
	InOutD* pin = new InOutD();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(double)> setVal, function<double&(void)> getVal){
	InOutD* pin = new InOutD();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(double)> setVal, function<double(void)> getVal){
	InOutD* pin = new InOutD();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

//-----------------------------------------------------//

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<int&(void)> getVal, function<void(int)> setVal){
	InOutI* pin = new InOutI();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<int(void)> getVal, function<void(int)> setVal){
	InOutI* pin = new InOutI();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(int)> setVal, function<int&(void)> getVal){
	InOutI* pin = new InOutI();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(int)> setVal, function<int(void)> getVal){
	InOutI* pin = new InOutI();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

//-----------------------------------------------------//

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<dvector1D(void)> getVal, function<void(dvector1D)> setVal){
	InOutDV1* pin = new InOutDV1();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(dvector1D)> setVal, function<dvector1D(void)> getVal){
	InOutDV1* pin = new InOutDV1();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<dvector1D&(void)> getVal, function<void(dvector1D)> setVal){
	InOutDV1* pin = new InOutDV1();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(dvector1D)> setVal, function<dvector1D&(void)> getVal){
	InOutDV1* pin = new InOutDV1();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

//-----------------------------------------------------//

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<ivector1D(void)> getVal, function<void(ivector1D)> setVal){
	InOutIV1* pin = new InOutIV1();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(ivector1D)> setVal, function<ivector1D(void)> getVal){
	InOutIV1* pin = new InOutIV1();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<ivector1D&(void)> getVal, function<void(ivector1D)> setVal){
	InOutIV1* pin = new InOutIV1();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(ivector1D)> setVal, function<ivector1D&(void)> getVal){
	InOutIV1* pin = new InOutIV1();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

//-----------------------------------------------------//

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<darray3E(void)> getVal, function<void(darray3E)> setVal){
	InOutDA3* pin = new InOutDA3();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(darray3E)> setVal, function<darray3E(void)> getVal){
	InOutDA3* pin = new InOutDA3();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<darray3E&(void)> getVal, function<void(darray3E)> setVal){
	InOutDA3* pin = new InOutDA3();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(darray3E)> setVal, function<darray3E&(void)> getVal){
	InOutDA3* pin = new InOutDA3();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

//-----------------------------------------------------//

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<iarray3E(void)> getVal, function<void(iarray3E)> setVal){
	InOutIA3* pin = new InOutIA3();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(iarray3E)> setVal, function<iarray3E(void)> getVal){
	InOutIA3* pin = new InOutIA3();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<iarray3E&(void)> getVal, function<void(iarray3E)> setVal){
	InOutIA3* pin = new InOutIA3();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(iarray3E)> setVal, function<iarray3E&(void)> getVal){
	InOutIA3* pin = new InOutIA3();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

//-----------------------------------------------------//

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<dvecarr3E(void)> getVal, function<void(dvecarr3E)> setVal){
	InOutDVA3* pin = new InOutDVA3();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(dvecarr3E)> setVal, function<dvecarr3E(void)> getVal){
	InOutDVA3* pin = new InOutDVA3();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<dvecarr3E&(void)> getVal, function<void(dvecarr3E)> setVal){
	InOutDVA3* pin = new InOutDVA3();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(dvecarr3E)> setVal, function<dvecarr3E&(void)> getVal){
	InOutDVA3* pin = new InOutDVA3();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

//-----------------------------------------------------//

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<ivecarr3E(void)> getVal, function<void(ivecarr3E)> setVal){
	InOutIVA3* pin = new InOutIVA3();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(ivecarr3E)> setVal, function<ivecarr3E(void)> getVal){
	InOutIVA3* pin = new InOutIVA3();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

void
BaseManipulation::addPinIn(BaseManipulation* objIn, function<ivecarr3E&(void)> getVal, function<void(ivecarr3E)> setVal){
	InOutIVA3* pin = new InOutIVA3();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

void
BaseManipulation::addPinOut(BaseManipulation* objOut, function<void(ivecarr3E)> setVal, function<ivecarr3E&(void)> getVal){
	InOutIVA3* pin = new InOutIVA3();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};



//EXTERNAL METHODS



