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
#include "Chain.hpp"

using namespace std;

uint8_t Chain::sm_chaincounter(0);

/*!Default constructor of Chain.
 * It sets to zero/null each member/pointer.
 */
Chain::Chain(){
	m_id 			= sm_chaincounter;
	m_objcounter	= 0;
	m_objects.clear();
	m_idObjects.clear();
	sm_chaincounter++;
};


/*!Default destructor of Chain.
 */
Chain::~Chain(){
	clear();
	sm_chaincounter--;
};

/*!Copy constructor of Chain.
 */
Chain::Chain(const Chain & other){
	m_id 		= other.m_id;
	m_objects	= other.m_objects;
	m_idObjects	= other.m_idObjects;
	m_objcounter= other.m_objcounter;
};

/*!Assignement operator of Chain.
 */
Chain & Chain::operator=(const Chain & other){
	m_id 		= other.m_id;
	m_objects	= other.m_objects;
	m_idObjects	= other.m_idObjects;
	m_objcounter= other.m_objcounter;
};

/*!It clears the object, by setting to zero/NULL each member/pointer in the object.
 */
void
Chain::clear(){
	m_id 			= sm_chaincounter;
	m_objcounter	= 0;
	m_objects.clear();
	m_idObjects.clear();
};


/*!It gets the number of degrees of freedom of the manipulator object.
 * \return #Object in the chain.
 */
uint32_t
Chain::getNObjects(){
	return m_objects.size();
};

uint8_t
Chain::getID(){
	return m_id;
};

uint8_t
Chain::getNChains(){
	return sm_chaincounter;
};

/*!Return true if the chain is interrupted.
 *
 */
bool
Chain::deleteObject(int idobj){
	bool cut = false;
	vector<int>::iterator it = find(m_idObjects.begin(), m_idObjects.end(), idobj);
	if (it != m_idObjects.end()){
		int idx = distance(m_idObjects.begin(), it);
		if (m_objects[idx]->getParent() != NULL || m_objects[idx]->getNChild() != 0) cut = true;
		m_objects.erase(m_objects.begin()+idx);
		m_idObjects.erase(it);
	}
	return cut;
};

/*!Return false if the object cannot be inserted (wrong parent/child dependencies).
 *
 */
int
Chain::addObject(BaseManipulation* obj, int id_){
	int id;
	if (id_ == 0){
		id = m_objcounter;
		m_objcounter++;
	}
	if (m_objects.size() == 0){
		m_objects.push_back(obj);
		m_idObjects.push_back(id);
		return 0;
	}
	vector<BaseManipulation*>::iterator itchild = m_objects.end();
	vector<BaseManipulation*>::iterator itparent;
	int idxchild = m_objects.size();
	ivector1D idxparent(obj->getNParent(), -1);
	if (obj->getNChild()>0){
		for (int i=0; i<obj->getNChild(); i++){
			itchild = find(m_objects.begin(), m_objects.end(), obj->getChild(i));
			idxchild = min(idxchild, int(distance(m_objects.begin(), itchild)));
		}
	}
	if (obj->getNParent()>0){
		for (int i=0; i<obj->getNParent(); i++){
			if (obj->getParent(i) != NULL){
				itparent = find(m_objects.begin(), m_objects.end(), obj->getParent(i));
				if (itparent != m_objects.end()){
					idxparent[i] = distance(m_objects.begin(), itparent);
				}
			}
		}
	}
	for (int i=0; i<obj->getNParent(); i++){
		if (idxparent[i] == idxchild) return -1;
	}

	bool parentsMinor = true;
	for (int i=0; i<obj->getNParent(); i++){
		if (idxparent[i] > idxchild){
			parentsMinor = false;
		}
	}
	if (!parentsMinor){
		ivector1D idparent(obj->getNParent(), -1);
		vector<BaseManipulation*> parent(obj->getNParent());
		for (int i=0; i<obj->getNParent(); i++){
			if (idxparent[i] > idxchild){
				parent[i] = obj->getParent(i);
				itparent = find(m_objects.begin(), m_objects.end(), parent[i]);
				m_objects.erase(itparent);
				idparent[i] = m_idObjects[idxparent[i]];
				m_idObjects.erase(m_idObjects.begin()+idxparent[i]);
			}
		}
		int idx = addObject(obj, id);
		for (int i=0; i<obj->getNParent(); i++){
			if (idxparent[i] > idxchild){
				int idxp = addObject(parent[i], idparent[i]);
				parent[i] = NULL;
			}
		}
		return idx+1;
	}else{
		m_objects.insert(itchild, obj);
		m_idObjects.insert(m_idObjects.begin()+idxchild, id);
		return idxchild;
	}
	return 0;
}

void
Chain::exec(){
	vector<BaseManipulation*>::iterator it, itb = m_objects.begin();
	vector<BaseManipulation*>::iterator itend = m_objects.end();
	for (it = itb; it != itend; ++it){
		(*it)->exec();
	}
}

void
Chain::exec(int idobj){
	int idx = distance(m_idObjects.begin(), find(m_idObjects.begin(), m_idObjects.end(), idobj));
	cout << idx << "/" << m_objects.size() << endl;
	m_objects[idx]->exec();
	cout << idx << "/" << m_objects.size() << endl;
}




