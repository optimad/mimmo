/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

// ========================================================================== //
//                         - SORTING ALGORITHMS -                             //
//                                                                            //
// Functions for data sorting.                                                //
// ========================================================================== //
// INFO                                                                       //
// Author    : Edoardo Lombardi                                               //
// Version   : v2.0                                                           //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

# include "MimmoNamespace.hpp"
# include "bitpit_patchkernel.hpp"

namespace mimmo{
namespace bvtree{

/*!
 \ingroup   SortAlgorithms
 \{
 */



/*!

    \class BvNode
    \brief class for bv-tree node.

    Store the informatino of a node in bv-tree data structure.

    Template parameters are:
    - T, container used to store node coordinates
    - T1, label associated to the node.

    Template parameters can be any type fulfilling the following requirements:
    1. operator* (dereferencing operator) must be defined for class T
    2. T1 must be a copy-constructible type.

 */

/*!
    Default constructor for class BvElement.
    Initialize an empty node in the bv-tree.
 */
BvElement::BvElement()
{
	m_object_ = NULL;
	return;
}


/*!
    Default destructor for class BvNode.
    Clear BvElement content and release memory.
 */
BvElement::~BvElement()
{
	m_object_ = NULL;
	return;
}



/*!

    \class BvNode
    \brief class for bv-tree node.

    Store the informatino of a node in bv-tree data structure.

    Template parameters are:
    - T, container used to store node coordinates
    - T1, label associated to the node.

    Template parameters can be any type fulfilling the following requirements:
    1. operator* (dereferencing operator) must be defined for class T
    2. T1 must be a copy-constructible type.

 */

/*!
    Default constructor for class BvNode.
    Initialize an empty node in the bv-tree.
 */
BvNode::BvNode()
{
	m_lchild = -1;
	m_rchild = -1;
	m_element[0] = -1;
	m_element[1] = -1;
	m_filled = false;
	return;
}


/*!
    Default destructor for class BvNode.
    Clear BvNode content and release memory.
 */
BvNode::~BvNode()
{
	m_lchild = -1;
	m_rchild = -1;
	m_element[0] = -1;
	m_element[1] = -1;
	m_filled = false;
	return;
}


// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS FOR BvTree                                        //
// ========================================================================== //

/*!

    \ingroup   SortAlgorithms
    \class BvTree
    \brief class for bv-tree data structure.

    Sort vertices in a d-dimensional Euclidean space into a bv-tree structure.

    Template parameters are:
    - d, number dimensions (i.e. number of coordinates for vertices)
    - T, container used to store vertex coordinates
    - T1, label type associated to each node in the bv-tree.

    Template parameters can be any type fulfilling the following requirements:
    1. operator* (dereferencing operator) must be defined for class T
    2. T1 must be a copy-constructible type.
    3. operator== must be defined for container T, which returns true
       if two objects of type T have the same content (i.e. two vertices have the
       same coordinates)

 */

/*!
    Default constructor for class BvTree.

    Initialize an empty bv-tree structure and reserve memory for the insertion
    of maxstack nodes.

    \param[in] maxstack memory reserve (in terms of number of elements)
 */
BvTree::BvTree(bitpit::PatchKernel *patch_)
{
	m_patch_ = patch_;
	m_dim = m_patch_->getDimension();
	m_MAXSTK = m_patch_->getVertexCount();
	m_nnodes = 0;
	m_nelements = 0;
	increaseStack();
	return;
}

/*!
    Default destructor for class BvTree.

    Clear bv-tree content and release memory.
 */
BvTree::~BvTree()
{
	m_patch_ = NULL;
	m_nelements = 0;
	m_nnodes = 0;
	m_MAXSTK = 0;
	m_nodes.clear();
	m_elements.clear();
	return;
}

void BvTree::initBoundingPoint(BvNode *nod_, bitpit::Cell *C_, int dim)
{
	int iV = 0;
	int nV = C_->getVertexCount();
	nod_->m_minP = m_patch_->getVertex(C_->getVertex(iV))[dim];
	nod_->m_maxP = m_patch_->getVertex(C_->getVertex(iV))[dim];
	for (iV=1; iV<nV; ++iV){
		nod_->m_minP = bitpit::min(nod_->m_minP, m_patch_->getVertex(C_->getVertex(iV))[dim]);
		nod_->m_maxP = bitpit::max(nod_->m_maxP, m_patch_->getVertex(C_->getVertex(iV))[dim]);
	}
	return;
}

void BvTree::updateBoundingPoint(BvNode *nod_, bitpit::Cell *C_, int dim)
{
	int iV, nV = C_->getVertexCount();
	for (iV=0; iV<nV; ++iV){
		nod_->m_minP = bitpit::min(nod_->m_minP, m_patch_->getVertex(C_->getVertex(iV))[dim]);
		nod_->m_maxP = bitpit::max(nod_->m_maxP, m_patch_->getVertex(C_->getVertex(iV))[dim]);
	}
	return;
}


///*!
//    Check whether a given vertex already exist in the bv-tree.
//    Check is performed via lexicographical comparison of vertex coordinates.
//
//    \param[in] P_ pointer to container storing the vertices of the triangle
//    to be checked.
//
//    \result on output returns the pointer to the element in the tree having the
//    same vertices of the input triangle. If no triangle is found in the bv-tree,
//    returns NULL.
// */
//template<int d, class T, class T3, class T2, class T1>
//BvNode<T2>* BvTree<d, T, T2, T1>::exist(
//		T      *P_
//) {
//
//	// Local variables
//	bool             check = false;
//	int              index = -1;
//	int              prev_ = -1, next_ = 0;
//	int              lev = 0, dim;
//
//	// ========================================================================== //
//	// EXIT FOR EMPTY TREE                                                        //
//	// ========================================================================== //
//	if (n_nodes == 0) { return(NULL); };
//
//	// ========================================================================== //
//	// MOVE ON TREE BRANCHES                                                      //
//	// ========================================================================== //
//
//	// Find node on leaf -------------------------------------------------------- //
//	while ((next_ >= 0) && (!check)) {
//		check = ((*P_) == (*(nodes[next_].object_)));
//		prev_ = next_;
//		dim = lev % d;
//		if ((*P_)[dim] <= (*(nodes[next_].object_))[dim]) {
//			next_ = nodes[next_].lchild_;
//		}
//		else {
//			next_ = nodes[next_].rchild_;
//		}
//		lev++;
//	} //next
//	if (check) {
//		index = prev_;
//	}
//
//	return(index); };
//
//// -------------------------------------------------------------------------- //
///*!
//    Check whether a given vertex already exist in the bv-tree.
//    Check is performed via lexicographical comparison of vertex coordinates.
//
//    \param[in] P_ pointer to container storing the coordinates of the vertex
//    to be checked.
//    \param[in,out] label on output stores the label associated with the BvNode
//    (if any) whose coordinates match those of the input vector
//
//    \result on output returns the index of the bv-node in the tree having the
//    same coordinates of the input vertex. If no node is found in the bv-tree,
//    returns -1.
// */
//template<int d, class T, class T3, class T2, class T1>
//int BvTree<d, T, T2, T1>::exist(
//		T               *P_,
//		T1              &label
//) {
//
//	// ========================================================================== //
//	// VARIABLES DECLARATION                                                      //
//	// ========================================================================== //
//
//	// Local variables
//	bool             check = false;
//	int              index = -1;
//	int              prev_ = -1, next_ = 0;
//	int              lev = 0, dim;
//
//	// Counters
//	// none
//
//	// ========================================================================== //
//	// EXIT FOR EMPTY TREE                                                        //
//	// ========================================================================== //
//	if (n_nodes == 0) { return(index); };
//
//	// ========================================================================== //
//	// MOVE ON TREE BRANCHES                                                      //
//	// ========================================================================== //
//
//	// Find node on leaf -------------------------------------------------------- //
//	while ((next_ >= 0) && (!check)) {
//		check = ((*P_) == (*(nodes[next_].object_)));
//		prev_ = next_;
//		dim = lev % d;
//		if ((*P_)[dim] <= (*(nodes[next_].object_))[dim]) {
//			next_ = nodes[next_].lchild_;
//		}
//		else {
//			next_ = nodes[next_].rchild_;
//		}
//		lev++;
//	} //next
//	if (check) {
//		index = prev_;
//		label = nodes[prev_].label;
//	}
//
//	return(index); };
//
//// -------------------------------------------------------------------------- //
///*!
//    Given an input vertex P, returns the index of the first node encountered
//    in the bv-tree which is in the 1-ball centered on P and having a radius of h.
//    The 1-ball is defined as:
//    B1(x; h) = {y: norm_1(y-x) <= h}
//    \param[in] P_ pointer to container storing the coordinates of P
//    \param[in] h 1-ball radius
//    \param[in] debug (default = false) flag for running the search in debug mode
//    \param[in] next_ (default = 0) index of element in the kd tree used as starting
//    point by the bv-tree search algorithm
//    \param[in] lev   (default = 0) level in bv-tree of node used as starting point
//    by the bv-tree search algorithm.
//
//    \result on output returns the index of the first bv-node encountered in the tree
//    which lies in the 1-ball centered on P. If no node is found, returns -1.
// */
//template<int d, class T, class T2, class T1>
//template< class T2>
//int BvTree<d, T, T2, T1>::hNeighbor(
//		T               *P_,
//		T2               h,
//		bool             debug,
//		int              next_,
//		int              lev
//) {
//
//	// ========================================================================== //
//	// VARIABLES DECLARATION                                                      //
//	// ========================================================================== //
//
//	// Local variables
//	int              index_l = -1, index_r = -1;
//	int              prev_ = next_;
//	int              dim;
//
//	// Counters
//	// none
//
//	// ========================================================================== //
//	// EXIT FOR EMPTY TREE                                                        //
//	// ========================================================================== //
//	if (n_nodes == 0) { return(-1); };
//
//	// ========================================================================== //
//	// MOVE ON TREE BRANCHES                                                      //
//	// ========================================================================== //
//
//	// Check if root is in the h-neighbor of P_ --------------------------------- //
//	//if (debug) { cout << "visiting: " << prev_ << endl; }
//	if (norm2((*(nodes[prev_].object_)) - (*P_)) <= h) { return(prev_); }
//
//	// Move on next branch ------------------------------------------------------ //
//	dim = lev % d;
//	if (((*(nodes[prev_].object_))[dim] >= (*P_)[dim] - h)
//			&& (nodes[prev_].lchild_ >= 0)) {
//		// if (nodes[prev_].lchild_ >= 0) {
//		next_ = nodes[prev_].lchild_;
//		index_l = hNeighbor(P_, h, debug, next_, lev+1);
//	}
//	if (((*(nodes[prev_].object_))[dim] <= (*P_)[dim] + h)
//			&& (nodes[prev_].rchild_ >= 0)) {
//		// if (nodes[prev_].rchild_ >= 0) {
//		next_ = nodes[prev_].rchild_;
//		index_r = hNeighbor(P_, h, debug, next_, lev+1);
//	}
//
//	//if (debug) { cout << "result is: " << max(index_l, index_r) << endl; }
//	return(std::max(index_l, index_r)); };
//
//// -------------------------------------------------------------------------- //
///*!
//    Given an input vertex P, returns the index of the first node encountered
//    in the bv-tree which is in the 1-ball centered on P and having a radius of h.
//    The 1-ball is defined as:
//    B1(x; h) = {y: norm_1(y-x) <= h}
//    \param[in] P_ pointer to container storing the coordinates of P
//    \param[in,out] label on output stores the label of the BvNode in the 1-ball
//    centered on P (if any).
//    \param[in] h 1-ball radius
//    \param[in] next_ (default = 0) index of element in the kd tree used as starting
//    point by the bv-tree search algorithm
//    \param[in] lev   (default = 0) level in bv-tree of node used as starting point
//    by the bv-tree search algorithm.
//
//    \result on output returns the index of the first bv-node encountered in the tree
//    which lies in the 1-ball centered on P. If no node is found, returns -1.
// */
//template<int d, class T, class T2, class T1 >
//template< class T2>
//int BvTree<d, T, T2, T1>::hNeighbor(
//		T               *P_,
//		T1              &label,
//		T2               h,
//		int              next_,
//		int              lev
//) {
//
//	// ========================================================================== //
//	// VARIABLES DECLARATION                                                      //
//	// ========================================================================== //
//
//	// Local variables
//	// bool             check = false;
//	int              index_l = -1, index_r = -1;
//	// int              dim;
//
//	// Counters
//	// none
//
//	// // ========================================================================== //
//	// // EXIT FOR EMPTY TREE                                                        //
//	// // ========================================================================== //
//	// if (n_nodes == 0) { return(index); };
//
//	// // ========================================================================== //
//	// // MOVE ON TREE BRANCHES                                                      //
//	// // ========================================================================== //
//
//	// // Find node on leaf -------------------------------------------------------- //
//	// while ((next_ >= 0) && (!check)) {
//	// check = (norm2((*P_) - (*(nodes[next_].object_))) <= h);
//	// prev_ = next_;
//	// dim = lev % d;
//	// if ((*P_)[dim] <= (*(nodes[next_].object_))[dim]) {
//	// next_ = nodes[next_].lchild_;
//	// }
//	// else {
//	// next_ = nodes[next_].rchild_;
//	// }
//	// lev++;
//	// } //next
//	// if (check) {
//	// index = prev_;
//	// label = nodes[prev_].label;
//	// }
//
//	return(std::max(index_l, index_r)); };
//
//// -------------------------------------------------------------------------- //
///*!
//    Given an input vertex P, returns the index of the first node encountered
//    in the bv-tree which is in the 1-ball centered on P and having a radius of h.
//    The 1-ball is defined as:
//    B1(x; h) = {y: norm_1(y-x) <= h}
//    \param[in] P_ pointer to container storing the coordinates of P
//    \param[in] h 1-ball radius
//    \param[in] debug (default = false) flag for running the search in debug mode
//    \param[in] next_ (default = 0) index of element in the kd tree used as starting
//    point by the bv-tree search algorithm
//    \param[in] lev   (default = 0) level in bv-tree of node used as starting point
//    by the bv-tree search algorithm.
//
//    \result on output returns the index of the first bv-node encountered in the tree
//    which lies in the 1-ball centered on P. If no node is found, returns -1.
// */
//template<int d, class T, class T2, class T1>
//template< class T2>
//int BvTree<d, T, T2, T1>::hNeighbor(
//		T               *P_,
//) {
//
//	// ========================================================================== //
//	// VARIABLES DECLARATION                                                      //
//	// ========================================================================== //
//
//	//Init variables
//	int              next_ = 0,
//	int              lev = 0,
//	T2               h = norm2((*(elements[next].object_)) - (*P_)),
//
//	// Local variables
//	int              index_l = -1, index_r = -1;
//	int              prev_ = next_;
//	int              dim;
//
//	// Counters
//	// none
//
//	// ========================================================================== //
//	// EXIT FOR EMPTY TREE                                                        //
//	// ========================================================================== //
//	if (n_nodes == 0) { return(-1); };
//
//	// ========================================================================== //
//	// MOVE ON TREE BRANCHES                                                      //
//	// ========================================================================== //
//
//	// Check if root is in the h-neighbor of P_ --------------------------------- //
//	//if (debug) { cout << "visiting: " << prev_ << endl; }
//	if (norm2((*(nodes[prev_].object_)) - (*P_)) <= h) { return(prev_); }
//
//	// Move on next branch ------------------------------------------------------ //
//	dim = lev % d;
//	if (((*(nodes[prev_].object_))[dim] >= (*P_)[dim] - h)
//			&& (nodes[prev_].lchild_ >= 0)) {
//		// if (nodes[prev_].lchild_ >= 0) {
//		next_ = nodes[prev_].lchild_;
//		index_l = hNeighbor(P_, h, debug, next_, lev+1);
//	}
//	if (((*(nodes[prev_].object_))[dim] <= (*P_)[dim] + h)
//			&& (nodes[prev_].rchild_ >= 0)) {
//		// if (nodes[prev_].rchild_ >= 0) {
//		next_ = nodes[prev_].rchild_;
//		index_r = hNeighbor(P_, h, debug, next_, lev+1);
//	}
//
//	//if (debug) { cout << "result is: " << max(index_l, index_r) << endl; }
//	return(std::max(index_l, index_r)); };
//

// -------------------------------------------------------------------------- //
///*!
//    Insert a new vertex in the bv-tree.
//
//    \param[in] P_ pointer to container storing element (simplex).
// */
//void BvTree::insert(bitpit::Cell *C_)
//{
//	// Local variables
//	bool            left = false;
//	bool			filled = true;
//	int             prev_ = -1, next_ = 0;
//	int             lev = 0, dim;
//	int				iV, nV = C_->getVertexCount();
//
//	// ========================================================================== //
//	// EXIT CONDITION FOR EMPTY TREE                                              //
//	// ========================================================================== //
//	if (m_nnodes == 0) {
//		m_elements[0].m_object_ = C_;
//		m_nodes[0].m_element[0] = 0;
//		m_nnodes++;
//		m_nelements++;
//		initBoundingPoint(&m_nodes[0], C_, lev);
//		return;
//	};
//	if (n_elements == 1) {
//		m_elements[1].m_object_ = C_;
//		m_nodes[0].m_element[1] = 1;
//		m_nnodes++;
//		m_nelements++;
//		m_nodes[0].m_filled = true;
//		updateBoundingPoint(&m_nodes[0], C_, lev);
//		return;
//	};
//
//	// ========================================================================== //
//	// MOVE ON TREE BRANCHES                                                      //
//	// ========================================================================== //
//
//	// Find node on leaf -------------------------------------------------------- //
//	while (next_ >= 0 && filled) {
//		prev_ = next_;
//		dim = lev % m_dim;
//		for (iV=0; iV<nV; ++iV){
//			if (m_patch_->getVertex(C_->getVertex(iV))[dim] > m_nodes[next_].m_minP && m_patch_->getVertex(C_->getVertex(iV))[dim] < m_nodes[next_].m_maxP) {
//				left = true;
//			}
//		}
//		if (left){
//			next_ = m_nodes[next_].m_lchild;
//		}
//		else {
//			next_ = m_nodes[next_].m_rchild;
//		}
//		if (next_ >=0 ) filled = m_nodes[next_].m_filled;
//		lev++;
//	} //next
//
//	// Insert new element ------------------------------------------------------- //
//
//	// Increase stack size
//	if (m_nnodes+1 > m_nodes.size() || m_nelements > m_elements.size()) {
//		increaseStack();
//	}
//
//	if (filled){
//		// Update parent
//		if (left) {
//			m_nodes[prev_].m_lchild = m_nnodes;
//		}
//		else {
//			m_nodes[prev_].m_rchild = m_nnodes;
//		}
//		// Insert children
//		m_elements[m_nelements].object_ = C_;
//		m_nodes[m_nnodes].element[0] = m_nelements;
//		dim = lev % m_dim;
//		initBoundingPoint(&m_nodes[m_nnodes], C_, dim);
//		m_nelements++;
//		m_nnodes++;
//	}
//	else{
//		// Insert element in not filled node
//		m_elements[n_elements].object_ = C_;
//		m_nodes[next_].m_element[1] = m_nelements;
//		m_nodes[next_].m_filled = true;
//		dim = lev % m_dim;
//		updateBoundingPoint(m_nodes[next_], C_, dim);
//		m_nelements++;
//	}
//
//	return;
//};
//

// -------------------------------------------------------------------------- //
/*!
    Insert a new vertex and the associated label into bv-tree.

    \param[in] P_ pointer to container storing vertex coordinates.
    \param[in] label label associated to P_
 */
void BvTree::insert(bitpit::Cell *C_, long label)
{

	// Local variables
	bool            left = false;
	bool			filled = true;
	int             prev_ = -1, next_ = 0;
	int             lev = 0, dim;
	int				iV, nV = C_->getVertexCount();

	// ========================================================================== //
	// EXIT CONDITION FOR EMPTY TREE                                              //
	// ========================================================================== //
	if (m_nnodes == 0) {
		m_elements[0].m_object_ = C_;
		m_nodes[0].m_element[0] = 0;
		m_nnodes++;
		m_nelements++;
		initBoundingPoint(&m_nodes[0], C_, lev);
		return;
	};
	if (m_nelements == 1) {
		m_elements[1].m_object_ = C_;
		m_nodes[0].m_element[1] = 1;
		m_nnodes++;
		m_nelements++;
		m_nodes[0].m_filled = true;
		updateBoundingPoint(&m_nodes[0], C_, lev);
		return;
	};

	// ========================================================================== //
	// MOVE ON TREE BRANCHES                                                      //
	// ========================================================================== //

	// Find node on leaf -------------------------------------------------------- //
	while (next_ >= 0 && filled) {
		prev_ = next_;
		dim = lev % m_dim;
		for (iV=0; iV<nV; ++iV){
			if (m_patch_->getVertex(C_->getVertex(iV))[dim] > m_nodes[next_].m_minP && m_patch_->getVertex(C_->getVertex(iV))[dim] < m_nodes[next_].m_maxP) {
				left = true;
			}
		}
		if (left){
			next_ = m_nodes[next_].m_lchild;
		}
		else {
			next_ = m_nodes[next_].m_rchild;
		}
		if (next_ >=0 ) filled = m_nodes[next_].m_filled;
		lev++;
	} //next

	// Insert new element ------------------------------------------------------- //

	// Increase stack size
	if (m_nnodes+1 > m_nodes.size() || m_nelements > m_elements.size()) {
		increaseStack();
	}

	if (filled){
		// Update parent
		if (left) {
			m_nodes[prev_].m_lchild = m_nnodes;
		}
		else {
			m_nodes[prev_].m_rchild = m_nnodes;
		}
		// Insert children
		m_elements[m_nelements].m_object_ = C_;
		m_elements[m_nelements].m_label = label;
		m_nodes[m_nnodes].m_element[0] = m_nelements;
		dim = lev % m_dim;
		initBoundingPoint(&m_nodes[m_nnodes], C_, dim);
		m_nelements++;
		m_nnodes++;
	}
	else{
		// Insert element in not filled node
		m_elements[m_nelements].m_object_ = C_;
		m_elements[m_nelements].m_label = label;
		m_nodes[next_].m_element[1] = m_nelements;
		m_nodes[next_].m_filled = true;
		dim = lev % m_dim;
		updateBoundingPoint(&m_nodes[next_], C_, dim);
		m_nelements++;
	}

	return;
};

// -------------------------------------------------------------------------- //
/*!
    Whenever the bv-tree reaches its full capacity (i.e. the number of BvNode
    stored in the tree is equal to the memory reserved), increase the memory
    reserve by maxstack. The parameters maxstack is set at construction.
 */
void BvTree::increaseStack()
{
	// ========================================================================== //
	// INCREASE STACK SIZE                                                        //
	// ========================================================================== //
	m_nodes.resize(m_nodes.size() + m_MAXSTK);
	m_elements.resize(m_elements.size() + m_MAXSTK);
	return;
};

// -------------------------------------------------------------------------- //
/*!
    Decrease the memory reserved for bv-tree by maxstack.
    The parameters maxstack is set at construction.
 */
void BvTree::decreaseStack()
{
	// ========================================================================== //
	// INCREASE STACK SIZE                                                        //
	// ========================================================================== //
	m_nodes.resize(bitpit::max(m_MAXSTK, int(m_nodes.size()) - m_MAXSTK));
	m_elements.resize(bitpit::max(m_MAXSTK, int(m_elements.size()) - m_MAXSTK));
	return;
};


//// -------------------------------------------------------------------------- //
///*!
//    Given an input vertex P, returns the index of the first node encountered
//    in the bv-tree which is in the 1-ball centered on P and having a radius of h.
//    The 1-ball is defined as:
//    B1(x; h) = {y: norm_1(y-x) <= h}
//    \param[in] P_ pointer to container storing the coordinates of P
//    \param[in] h 1-ball radius
//    \param[out] L_ pointer to container filled with labels of the
//     bv-nodes encountered in the tree which lies in the 1-ball centered on P.
//    \param[in] EX_ pointers to vector with labels to be excluded from the results
//    \param[in] next_ (default = 0) index of element in the kd tree used as starting
//    point by the bv-tree search algorithm
//    \param[in] lev   (default = 0) level in bv-tree of node used as starting point
//    by the bv-tree search algorithm.
//
// */
//template<int d, class T, class T2, class T1 >
//template< class T2>
//void BvTree<d, T, T2, T1>::hNeighbor(
//		T               *P_,
//		T2               h,
//		std::vector<T1> *L_,
//		std::vector<T1> *EX_,
//		int              next_,
//		int              lev
//) {
//
//	// ========================================================================== //
//	// VARIABLES DECLARATION                                                      //
//	// ========================================================================== //
//
//	// Local variables
//	int              index_l = -1, index_r = -1;
//	int              prev_ = next_;
//	int              dim;
//
//	// Counters
//	// none
//
//	// ========================================================================== //
//	// EXIT FOR EMPTY TREE                                                        //
//	// ========================================================================== //
//	if (n_nodes == 0) { return; };
//
//	// ========================================================================== //
//	// MOVE ON TREE BRANCHES                                                      //
//	// ========================================================================== //
//
//	// Check if root is in the h-neighbor of P_ --------------------------------- //
//	if (norm2((*(nodes[prev_].object_)) - (*P_)) <= h) {
//
//		if (EX_ != NULL){
//			if (std::find(EX_->begin(), EX_->end(), nodes[prev_].label) == EX_->end() ) {
//				L_->push_back(nodes[prev_].label);
//			}
//		}
//		else{
//			L_->push_back(nodes[prev_].label);
//		}
//	}
//
//
//
//	// Move on next branch ------------------------------------------------------ //
//	dim = lev % d;
//	if (((*(nodes[prev_].object_))[dim] >= (*P_)[dim] - h)
//			&& (nodes[prev_].lchild_ >= 0)) {
//		next_ = nodes[prev_].lchild_;
//		hNeighbor(P_, h, L_, EX_, next_, lev+1);
//	}
//	if (((*(nodes[prev_].object_))[dim] <= (*P_)[dim] + h)
//			&& (nodes[prev_].rchild_ >= 0)) {
//		next_ = nodes[prev_].rchild_;
//		hNeighbor(P_, h, L_, EX_, next_, lev+1);
//	}
//
//	return ;
//
//};



/*!
 \}
 */

}
}
