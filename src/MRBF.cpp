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
using namespace bitpit;
using namespace mimmo;


/*! Default Constructor.*/
MRBF::MRBF(){
	m_name = "MiMMO.MRBF";
	setType(RBFType::INTERP);
	m_maxFields=-1;
	m_tol = 0.00001;
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
	*(static_cast<RBF * > (this)) = *(static_cast <const RBF*>(&other));
	*(static_cast<BaseManipulation * > (this)) = *(static_cast <const BaseManipulation * >(&other));
	m_tol = other.m_tol;
	return(*this);
};

/*!Adds a RBF point to the total control node list and activate it.
 * \param[in] node coordinates of control point.
 * \return RBF id.
 */
int MRBF::addNode(darray3E node){
	return(RBF::addNode(node));
};

/*!Adds a list of RBF points to the total control node list and activate them.
 * \param[in] nodes coordinates of control points.
 * \return Vector of RBF ids.
 */
std::vector<int> MRBF::addNode(dvecarr3E nodes){
	return(RBF::addNode(nodes));
};

/*!Adds a set of RBF points to the total control node list extracting
 * the vertices stored in a MimmoObject container. Return an unordered map < a,b > 
 * with key a equal to the global ID of the MimmoObject Vertex, and b equal to 
 * the RBF node id.
 * \param[in] geometry Pointer to MimmoObject that contains the geometry.
 * \return Vector of RBF ids.
 */
ivector1D MRBF::addNode(MimmoObject* geometry){
	dvecarr3E vertex = geometry->getVertex();
	return(RBF::addNode(vertex));
};


/*!Set a RBF point as unique control node and activate it.
 * \param[in] node coordinates of control point.
 */
void MRBF::setNode(darray3E node){
	removeAllNodes();
	RBF::addNode(node);
	return;
};

/*!Set a list of RBF points as control nodes and activate it.
 * \param[in] node coordinates of control points.
 */
void MRBF::setNode(dvecarr3E nodes){
	removeAllNodes();
	RBF::addNode(nodes);
	return;
};

/*!Set the RBF points as control nodes extracting
 * the vertices stored in a MimmoObject container.
 * \param[in] geometry Pointer to MimmoObject that contains the geometry.
 */
void MRBF::setNode(MimmoObject* geometry){
	removeAllNodes();
	dvecarr3E vertex = geometry->getVertex();
	RBF::addNode(vertex);
	return;
};

/*! Find all possible duplicated nodes within a prescribed distance tolerance.
 * Default tolerance value is 1.0E-12;
 * \param[in] tol distance tolerance
 * \return	list of duplicated nodes.
 */
ivector1D MRBF::checkDuplicatedNodes(double tol){
	ivector1D marked;
	int sizeEff = getTotalNodesCount();
	if( sizeEff == 0 ) return marked;
	
	bvector1D check(sizeEff, false);
	
	darray3E target = m_node[0];
	for(int i=1; i<sizeEff; ++i){
		double dist = norm2(m_node[i] - target);
		if(!check[i] && dist <= tol){
			marked.push_back(i);
			check[i] = true;
		}
	}
	return(marked);	
}

/*! Erase all nodes passed by their RBF id list. If no list is provided, the method find all 
 * possible duplicated nodes within a default tolerance of 1.0E-12 and erase them, if any.
 * \param[in] list pointer to a list of id's of RBF candidate nodes
 * \return	boolean, true if all duplicated nodes are erased, false if one or more of them are not.
 */
bool MRBF::removeDuplicatedNodes(ivector1D * list){
	ivector1D marked;
	if(list==NULL){
		marked = checkDuplicatedNodes();
		list = &marked;
	}
	return(removeNode(*list));
}

/*!It sets the tolerance for greedy algorithm.
 * \param[in] tol Target tolerance.
 */
void MRBF::setTol(double tol){
	m_tol = tol;
}

/*!
 * Define 3D displacements of your RBF parameterization. The method is not active in RBFType::INTERP mode
 * \param[in] displ list of nodal displacements
 */
void MRBF::setDisplacements(dvecarr3E displ){
	int size = displ.size();
	if(size != getTotalNodesCount()){
		std::cout << "MiMMO : ERROR : " << getName() << " displacements size (" << size << ") does not fit number of RBF nodes ("<< getTotalNodesCount() << ")" << std::endl;
		std::cout << "MiMMO :         " << getName() << " displacements resize to "<< getTotalNodesCount() << " (0 value for new displacements) " << std::endl;
		displ.resize(getTotalNodesCount(), {{0.0, 0.0, 0.0}});
		size = displ.size();
	}
	
	removeAllData();

	double maxdispl = 0.0;
	for(int i=0; i<size; ++i){
		maxdispl = std::max(maxdispl, norm2(displ[i]));
	}
	setSupportRadius(3*maxdispl);
	
	dvector1D temp(size);
	for(int loc=0; loc<3; ++loc){
		for(int i=0; i<size; ++i){
			temp[i] = displ[i][loc];
		}
		addData(temp);
	}
}

/*!
 * Define 3D displacements of your RBF parameterization, only on currently active Nodes. The method is not active in RBFType::INTERP mode
 * \param[in] displ list of nodal displacements on active nodes
 */
void MRBF::setActiveDisplacements(dvecarr3E displ){
	int size = displ.size();
	if(size != getActiveCount()){
		std::cout << "MiMMO : ERROR : " << getName() << " displacements size (" << size << ") does not fit number of active RBF nodes ("<< getActiveCount() << ")" << std::endl;
		std::cout << "MiMMO :         " << getName() << " displacements resize to "<< getActiveCount() << " (0 value for new displacements) " << std::endl;
		displ.resize(getActiveCount(), {{0.0, 0.0, 0.0}});
		size = displ.size();
	}
	
	removeAllData();
	
	dvector1D temp(size);
	for(int loc=0; loc<3; ++loc){
		
		for(int i=0; i<size; ++i){
			temp[i] = displ[i][loc];
		}
		dvector1D dummy = convertActiveToTotal(temp);
		addData(dummy);
	}
}

/*!Execution of RBF object. It evaluates the displacements (values) over the point of the
 * linked geometry, given as result of RBF technique implemented in bitpit::RBF base class.
 * The result is stored in the result member of BaseManipulation base class.
 *
 */
void MRBF::execute(){

	MimmoObject * container = getGeometry();
	if(container == NULL ) return;

	if (whichType() == RBFType::INTERP) greedy(m_tol);

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


/*! Expand vector of weights defined only for active RBF nodes
 * to a vector defined on all RBF nodes.
 * Zero weight is provided on inactive nodes
 */
dvector1D MRBF::convertActiveToTotal(dvector1D & target){
	
	dvector1D result(getTotalNodesCount(), 0.0);
	if(target.size() != getActiveCount()) return result;
	
	ivector1D list = getActiveSet();
	int size = list.size();
	for(int i=0; i<size; ++i){
		result[list[i]] = target[i];
	}
	return result;
}

