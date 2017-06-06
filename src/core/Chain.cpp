/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
 *
 *  mimmo is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmo is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#include "Chain.hpp"

using namespace std;

namespace mimmo{

uint8_t Chain::sm_chaincounter(0);

/*!
 * Default constructor of Chain.
 * It sets to zero/null each member/pointer.
 */
Chain::Chain(){
    m_id 			= sm_chaincounter;
    m_objcounter	= 0;
    m_objects.clear();
    m_idObjects.clear();
    sm_chaincounter++;
};


/*!
 * Default destructor of Chain.
 */
Chain::~Chain(){
    clear();
    sm_chaincounter = 0;
};

/*!
 * Copy constructor of Chain.
 */
Chain::Chain(const Chain & other){
    *this = other;
};

/*!
 * Assignement operator of Chain.
 */
Chain & Chain::operator=(const Chain & other){
    m_id 		= other.m_id;
    m_objects	= other.m_objects;
    m_idObjects	= other.m_idObjects;
    m_objcounter= other.m_objcounter;
    return (*this);
};

/*!
 * It clears the object, by setting to zero/NULL each member/pointer in the object.
 */
void
Chain::clear(){
    m_objcounter	= 0;
    m_objects.clear();
    m_idObjects.clear();
};

/*!
 * It gets the number of manipulator objects in the chain.
 * \return Number of objects in the chain.
 */
uint32_t
Chain::getNObjects(){
    return m_objects.size();
};

/*!
 * It gets the ID of the chain.
 * \return ID of the chain.
 */
uint8_t
Chain::getID(){
    return m_id;
};

/*!
 * It gets the number of chains actually defined in the process.
 * \return Number of chains in the process.
 */
uint8_t
Chain::getNChains(){
    return sm_chaincounter;
};

/*!
 * It gets the ID of an object in the chain.
 * \param[in] i Index of the target object in the chain.
 * \return ID of the i-th object.
 */
int
Chain::getID(int i){
    return m_idObjects[i];
};

/*!
 * It gets the name of an object in the chain.
 * \param[in] i Index of the target object in the chain.
 * \return name of the i-th object.
 */
string
Chain::getName(int i){
    return m_objects[i]->getName();
};

/*!
 * It deletes a manipulator object in the chain.
 * \return True if some linking in the chain after the deletion are interrupted.
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

/*!
 * It adds a manipulator object in the chain.
 * The position of the objects after the insertion can be modified
 * in order to respect the parent/child dependencies during the execution of the chain.
 * \param[in] obj Pointer to object to be inserted.
 * \param[in] id_ Id of the object to be inserted (optional, if not passed or negative it is computed).
 * \return Index of insertion of the object in the chain. Return -1 if the object cannot be inserted
 * (wrong parent/child dependencies).
 */
int
Chain::addObject(BaseManipulation* obj, int id_){
    int id = id_;
    if (id_ < 0){
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
            if (itchild != m_objects.end()){
                idxchild = min(idxchild, int(distance(m_objects.begin(), itchild)));
            }
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
                int idxtemp = distance(m_objects.begin(), itparent);
                idparent[i] = m_idObjects[idxtemp];
                m_objects.erase(itparent);
                m_idObjects.erase(m_idObjects.begin()+idxtemp);
            }
        }
        int idx = addObject(obj, id);
        for (int i=0; i<obj->getNParent(); i++){
            if (idxparent[i] > idxchild){
                addObject(parent[i], idparent[i]);
                parent[i] = NULL;
            }
        }
        return idx+1;
    }else{
        m_objects.insert(m_objects.begin()+idxchild, obj);
        m_idObjects.insert(m_idObjects.begin()+idxchild, id);
        return idxchild;
    }
    return 0;
}

/*!
 * It executes the chain, i.e. it executes all the manipulator objects
 * contained in the chain following the correct order.
 * In the case that a loop exists in the chain the execution doesn't start and
 * the process ends with an error.
 * \param[in]	debug boolean to activate verbose execution mode.
 */
void
Chain::exec(bool debug){
    m_log = &bitpit::log::cout(MIMMO_LOG_FILE);
    m_log->setPriority(bitpit::log::DEBUG);
    vector<BaseManipulation*>::iterator it, itb = m_objects.begin();
    vector<BaseManipulation*>::iterator itend = m_objects.end();
    if(debug)
        m_log->setPriority(bitpit::log::NORMAL);
    (*m_log) << " " << endl;
    (*m_log) << "--------------------------------------------------" << std::endl;
    (*m_log) << " Execution of chain - "<< m_objcounter << " objects" << std::endl;
    (*m_log) << " " << std::endl;
    checkLoops();
    int i = 1;
    for (it = itb; it != itend; ++it){
        (*m_log) << " execution object " << i << "	: " << (*it)->getName() << std::endl;
        (*it)->exec();
        i++;
    }

    (*m_log) << " " << std::endl;
    (*m_log) << "--------------------------------------------------" << std::endl;
    (*m_log) << " " << std::endl;
    m_log->setPriority(bitpit::log::NORMAL);
}

/*!
 * It executes one manipulator object contained in the chain singularly.
 * \param[in] idobj ID of the target manipulator object.
 */
void
Chain::exec(int idobj){
    int idx = distance(m_idObjects.begin(), find(m_idObjects.begin(), m_idObjects.end(), idobj));
    if (idx <  (int)m_objects.size()) m_objects[idx]->exec();
}

/*!
 * It checks if a loop exists in the chain.
 * In the case that a loop exists the process ends with an error.
 */
void
Chain::checkLoops(){
    vector<BaseManipulation*>::iterator it, itb = m_objects.begin();
    vector<BaseManipulation*>::iterator itend = m_objects.end();
    int actualidx = 0;
    for (it = itb; it != itend; ++it){
        for (int i=0; i<(*it)->getNChild(); i++){
            int idxchild = distance(m_objects.begin(), find(m_objects.begin(), m_objects.end(), (*it)->getChild(i)));
            if (idxchild <= actualidx){
                (*m_log) << " error : loop in chain : "<< (*it)->getName() << " linked to " << (*it)->getChild(i)->getName() <<  std::endl;
                (*m_log) << " " << std::endl;
                throw std::runtime_error ("loop in chain : " + (*it)->getName() + " linked to " + (*it)->getChild(i)->getName());
            }
        }
        actualidx++;
    }
}

}

