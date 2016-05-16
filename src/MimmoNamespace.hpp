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
#ifndef __MIMMONAMESPACE_HPP__
#define __MIMMONAMESPACE_HPP__

#include <functional>
#include <map>

#include "MimmoObject.hpp"

namespace mimmo{

namespace pin{

enum PinsType{BOTH, BACKWARD, FORWARD}; 	/**< Type of pins of the object: bidirectional,
												only input or only output.*/

template<typename OO, typename G, typename OI, typename S, typename VAL>
bool addPin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL));

template<typename OO, typename G, typename OI, typename S, typename VAL>
bool addPin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL));

template<typename OO, typename G, typename OI, typename S, typename VAL>
bool addPin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL));

template<typename OO, typename G, typename OI, typename S, typename VAL>
bool addPin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL*));

template<typename OO, typename G, typename OI, typename S, typename VAL>
bool addPin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL*));

template<typename OO, typename G, typename OI, typename S, typename VAL>
bool addPin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL*));

template<typename T, typename U, typename VAL>
std::function<VAL(void)> pinGet(VAL (T::*fget) (), U* obj);

template<typename T, typename U, typename VAL>
std::function<VAL&(void)> pinGetR(VAL& (T::*fget) (), U* obj);

template<typename T, typename U, typename VAL>
std::function<VAL*(void)> pinGetP(VAL* (T::*fget) (), U* obj);

template<typename T, typename U, typename VAL>
std::function<void(VAL)> pinSet(void (T::*fset) (VAL), U* obj);

template<typename T, typename U, typename VAL>
std::function<void(VAL*)> pinSetP(void (T::*fset) (VAL*), U* obj);

template<typename OO, typename OI>
void removeAllPins(OO* objSend, OI* objRec);

template<typename OO, typename G, typename OI, typename S, typename VAL>
void findPin(OO* objSend, OI* objRec);

template<typename OO, typename G, typename OI, typename S, typename VAL>
void removePin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL));

template<typename OO, typename G, typename OI, typename S, typename VAL>
void removePin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL));

template<typename OO, typename G, typename OI, typename S, typename VAL>
void removePin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL));

template<typename OO, typename G, typename OI, typename S, typename VAL>
void removePin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL*));

template<typename OO, typename G, typename OI, typename S, typename VAL>
void removePin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL*));

template<typename OO, typename G, typename OI, typename S, typename VAL>
void removePin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL*));

};


namespace bvtree{

// BvTree for bitpit::patch------------------------------------------------------------------- //
class BvElement {
public:
	bitpit::Cell	*m_object_;                   /**< pointer to object (bitpit::Cell) */
	long 			m_label;                     /**< label of the object (ID of the Cell)*/

public:
	BvElement();
	~BvElement();
};

class BvNode {
public:
	int					m_lchild;            /**< index of left child couple */
	int					m_rchild;            /**< index of right child couple */
	std::array<int,2>	m_element;            /**< index of elements of the couple */
	bool				m_filled;				/**< if two elements in node filled = true. */
	double				m_minP;
	double				m_maxP;

public:
	BvNode();
	~BvNode();
};

class BvTree {
public:
	int						m_dim;
	bitpit::PatchKernel     *m_patch_;                           /**< Patch*/
	int                     m_MAXSTK;                            /**< max stack size */
	int                    	m_nnodes;                            /**< number of nodes */
	std::vector<BvNode>		m_nodes;                             /**< bv-tree nodes */
	int                     m_nelements;                         /**< number of elements */
	std::vector<BvElement>	m_elements;                          /**< bv-tree elements */

public:
	BvTree(bitpit::PatchKernel *patch_);
	~BvTree();

	void initBoundingPoint(BvNode *nod_, bitpit::Cell *C_, int dim);
	void updateBoundingPoint(BvNode *nod_, bitpit::Cell *C_, int dim);
	void buildBvTree();

	/*	int exist(                                                               // Check if element exist in the bv-tree
			  T           *                                                  // (input) pointer to element to be tested
			  );
	int exist(                                                               // Check if element exist in the bv-tree
			  T           *,                                                 // (input) pointer to element to be tested
			  T1          &                                                  // (input/output) label of the bv node matching test object
			  );

	template <class T2>
	int hNeighbor(                                                           // Check if a bv-node exists in the h-neighborhood of a given item
				  T           *,                                                        // (input) pointer to element to be tested
				  T2           ,                                                        // (input) radius of ball
				  bool         ,
				  int         n = 0,                                                    // (input/optional) root for binary search algorithm
				  int         l = 0                                                     // (input/optional) level of root on binary tree
				  );

	template <class T2>
	void hNeighbor(                                                           // Check if a bv-node exists in the h-neighborhood of a given item
				   T           	*,                                                        // (input) pointer to element to be tested
				   T2           	,                                                        // (input) radius of ball
				   std::vector<T1>	*,                                                        // (output) pointer to container of labels of h-neighbors
				   std::vector<T1>	*,
				   int         	n = 0,                                                    // (input/optional) root for binary search algorithm
				   int         	l = 0                                                     // (input/optional) level of root on binary tree
				   );

	template <class T2>
	int hNeighbor(                                                           // Check if a bv-node exists in the h-neighborhood of a given item
				  T           *,                                                        // (input) pointer to element to be tested
				  T1          &,                                                      // (input/output) label of the bv node matching test object
				  T2           ,                                                        // (input) radius of ball
				  int         n = 0,                                                    // (input/optional) root for binary search algorithm
				  int         l = 0                                                     // (input/optional) level of root on binary tree
				  );
	 */


	//	template <class T2>
	//	int minhNeighbor(                                                           // Check if a bv-node exists in the h-neighborhood of a given item
	//				  T           *,                                                        // (input) pointer to element to be tested
	//				  );

	//	void insert(bitpit::Cell *C);
	void insert(bitpit::Cell *C, long id = -1);

private:
	void increaseStack();
	void decreaseStack();

};

}//end namespace bvtree

};//end namespace mimmo

#include "MimmoNamespace.tpp"


#endif
