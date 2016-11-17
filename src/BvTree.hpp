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
# ifndef __BVTREE_HPP__
# define __BVTREE_HPP__

# include "bitpit_patchkernel.hpp"
# include "bitpit_surfunstructured.hpp"


namespace mimmo{

/*!
 *	\date			20/may/2016
 *	\authors		Edoardo Lombardi
 *
 *	\brief Bv-Element is the class to manage the elements linked by a Bv-Tree.
 *
 */
class BvElement {
public:
	long 					m_label;	/**< Label of the object (ID of the Cell). */
	std::array<double,3>	m_centroid;	/**< Centroid of the cell. */

public:
	BvElement(long label = -1);
	~BvElement();
	BvElement(const BvElement & other);
	BvElement & operator=(const BvElement & other);
};

/*!
 *	\date			20/may/2016
 *	\authors		Edoardo Lombardi
 *
 *	\brief Bv-Node is the class of a node of a Bv-Tree.
 *
 */
class BvNode {
public:
	std::array<double,3>	m_minPoint;	/**<Minimum coordinates of the bounding box of the node. */
	std::array<double,3>	m_maxPoint;	/**<Maximum coordinates of the bounding box of the node. */
	std::array<int,2>	 	m_element;	/**<Range index of elements of the node. */
	int					 	m_nrange;	/**<Number of elements in the bounding box of the node. */
	int						m_lchild;	/**<Index of left child. */
	int						m_rchild;	/**<Index of right child. */
	bool					m_leaf;		/**<True if the node is a leaf node. */

public:
	BvNode();
	~BvNode();
	BvNode(const BvNode & other);
	BvNode & operator=(const BvNode & other);
};

/*!
 *	\date			20/may/2016
 *	\authors		Edoardo Lombardi
 *
 *	\brief Bv-Tree is the class to manage a Bounding Volume Hierarchy tree of a bitpit patch.
 *
 * A Bv-Tree is composed by its nodes and the elements related to each node.
 * A node has two possible child given by splitting in two sub-volume its bounding box.
 * The bounding box of each node is computed as the bounding volume of all the simplex
 * (elements) of the geometry contained in the box.
 * In order to split the box the average of the centroids of the elements
 * is computede and then the volume is divided by using the plane normal to
 * the direction of maximum span of the centroids.
 * Each node has the information about many elements are contained in its
 * bounding volume and the first and last index of them as stored in the structure
 * elements. The elements are ordered in function of the coordinates of their
 * centroid during the construction of the tree.
 * The leaf nodes can be have a number of elements in their bounding volume
 * greater then one.
 *
 */
class BvTree {
public:
	int							m_dim;			/**<Dimension of the space (currently only 3d space admitted.). */
	bitpit::PatchKernel			*m_patch_;		/**< Patch linked by the bv-tree. */
	int                     	m_nelements;	/**< Number of elements of the patch. */
	std::vector<BvElement>		m_elements;		/**< Bv-tree elements structure. */
	int                    		m_nnodes;		/**< Number of nodes in the tree. */
	std::vector<BvNode>			m_nodes;		/**< Bv-tree nodes structure. */
	int                     	m_MAXSTK;		/**< Max stack size. */

	int							m_nleaf;		/**<Number of leaf nodes in the bv-tree. */
	int							m_maxsize;		/**<Maximum number of elements for leaf nodes. */

private:
	double 						m_tol;			/**<Internal tolerance.*/

public:
	BvTree(bitpit::PatchKernel *patch_ = NULL);
	~BvTree();
	BvTree(const BvTree & other);
	BvTree & operator=(const BvTree & other);

	void setPatch(bitpit::PatchKernel *patch_);
	void setMaxLeafSize(int maxsize);
    bool inBoundingBox(std::array<double,3> *P_, BvNode *nod_, double r = 0.0);
    bool SphereBoundingBox(std::array<double,3> *P_, BvNode *nod_, double r = 0.0);

	void clean();
	void setup();
	void buildTree();


private:
	void fillTree(int iparent);
	std::array<double,3> computeMeanPoint(int istart, int iend);
	int findFirstGreater(int inode, std::array<double,3> meanC, int dir);
	int pseudoSort(std::vector<BvElement>::iterator itstart,
			std::vector<BvElement>::iterator itend, std::array<double,3> meanC, int dir);
	void computeBoundingBox(int inode);
	void increaseStack();
	void decreaseStack();

};

/*!
 *	\date			20/may/2016
 *	\authors		Edoardo Lombardi
 *
 *	\brief Utilities employing bvTree.
 */
namespace bvTreeUtils{
    double signedDistance(std::array<double,3> *P_, BvTree *bvtree_, long &id, std::array<double,3>  &n, double &r, int method = 1, bitpit::SurfUnstructured *spatch_ = NULL, int next = 0, double h = 1.0e+18);
	double distance(std::array<double,3> *P_, BvTree* bvtree_, long &id, double &r, int method = 1, int next = 0, double h = 1.0e+18);
	std::array<double,3> projectPoint(std::array<double,3> *P_, BvTree *bvtree_, double r_ = 1.0e+18);

	std::vector<double> signedDistance(std::vector<std::array<double,3> > *P_, BvTree *bvtree_, std::vector<long> &id, std::array<double,3>  &n, double r_ = 1.0e+18, int method = 1 );
	std::vector<double> distance(std::vector<std::array<double,3> > *P_, BvTree *bvtree_, std::vector<long> &id, double r_ = 1.0e+18, int method = 1 );
	std::vector<std::array<double,3> > projectPoint(std::vector<std::array<double,3> > *P_, BvTree *bvtree_, double r_ = 1.0e+18);

	std::vector<long> selectByPatch(BvTree *selection, BvTree *target, double tol = 1.0e-04);
	void extractTarget(BvTree *target, std::vector<BvNode*> leafSelection, std::vector<long> &extracted, double tol, int next = 0);
	
	
}//end namespace bvTreeUtils

}//end namespace mimmo

#endif
