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
	m_maxFields=-1;
	m_tol = 0.00001;
	m_solver = MRBFSol::GREEDY;
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
	m_solver = other.m_solver;
	return(*this);
};


/*!It sets the geometry linked by the manipulator object (overloading of base class method).
 * \param[in] geometry Pointer to geometry to be deformed by the manipulator object.
 */
void
MRBF::setGeometry(MimmoObject* geometry){
	m_geometry = geometry;
};

/*! 
 * Return actual solver set for RBF data fields interpolation in MRBF::execute
 */
MRBFSol
MRBF::getSolver(){
	return m_solver;
};

/*!
 * Set type of solver set for RBF data fields interpolation in MRBF::execute
 * \param[in] solver type of MRBFSol enum; 
 */
void
MRBF::setSolver(MRBFSol solver){
	m_solver = solver;
};
/*!
 * Overloading of MRBF::setSolver(MRBFSol solver) with int input parameter
 * \param[in] int type of solver 1-WHOLE, 2-GREEDY, see MRBFSol enum; 
 */
void 
MRBF::setSolver(int type){
	m_solver = MRBFSol::GREEDY;
	if(type ==1){	m_solver = MRBFSol::WHOLE;}
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
 * Define 3D displacements of your RBF parameterization. 
 * The method temporary set your RBF supportRadius to 3*max norm value of vector field displ.
 * \param[in] displ list of nodal displacements
 */
void MRBF::setDisplacements(dvecarr3E displ){
	int size = displ.size();
	if(size != getTotalNodesCount()){
		std::cout << "MiMMO : WARNING : " << getName() << " sets displacements with size (" << size << ") that does not fit number of RBF nodes ("<< getTotalNodesCount() << ")" << std::endl;
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

/*!Execution of RBF object. It evaluates the displacements (values) over the point of the
 * linked geometry, given as result of RBF technique implemented in bitpit::RBF base class.
 * The result is stored in the result member of BaseManipulation base class.
 *
 */
void MRBF::execute(){

	MimmoObject * container = getGeometry();
	if(container == NULL ) return;

	int size;
	for (int i=0; i<getDataCount(); i++){
		size = m_value[i].size();
		if(size != getTotalNodesCount()){
			std::cout << "MiMMO : WARNING : " << getName() << " has displacements of " << i << " field with size (" << size << ") that does not fit number of RBF nodes ("<< getTotalNodesCount() << ")" << std::endl;
			fitDataToNodes(i);
		}
	}

	if (m_solver == MRBFSol::WHOLE) solve();
	else	greedy(m_tol);

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


