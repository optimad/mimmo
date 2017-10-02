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

# include "SkdTreeUtils.hpp"
# include "CG.hpp"
# include "mimmoTypeDef.hpp"
# include <cmath>

namespace mimmo{
//
///*!
// * Default constructor for class BvElement.
// * Initialize an empty element of the bv-tree.
// */
//BvElement::BvElement(long label)
//{
//    m_label        = label;
//    m_centroid    = {{0.0, 0.0, 0.0}};
//}
//
///*!
// * Default destructor for class BvElement.
// * Clear BvElement content and release memory.
// */
//BvElement::~BvElement(){}
//
///*!
// * Assignement operator of BvElement.
// */
//BvElement & BvElement::operator=(const BvElement & other)
//{
//    m_label        = other.m_label;
//    m_centroid    = other.m_centroid;
//    return *this;
//}
//
///*!
// * \class compareElements
// * \ingroup core
// * \brief compareElements is an ad-hoc class used to sort the elements of a Bv-Tree by their
// * centroid coordinates.
// *
// */
//class compareElements{
//public:
//    int dir; /**< target centroid coordinate (0,1,2) to compare */
//public:
//    /*!Custom constructor for class compareElements.
//     * \param[in] dir_ Direction used to compare the coordinates of the elements centroid.
//     */
//    compareElements(int dir_) : dir(dir_){}
//
//    /*!Custom operator () for class compareElements.
//     * \param[in] a First element to be compared.
//     * \param[in] b Secon element to be compared.
//     */
//    bool operator()(BvElement a, BvElement b)
//    {
//        return (a.m_centroid[dir] < b.m_centroid[dir] );
//    }
//};
//
///*!
// * Default constructor for class BvNode.
// * Initialize an empty node in the bv-tree.
// */
//BvNode::BvNode()
//{
//    m_lchild        = -1;
//    m_rchild         = -1;
//    m_element[0]    = -1;
//    m_element[1]    = -1;
//    m_nrange        = 0;
//    m_leaf            = false;
//    m_minPoint        = {{m_MAXLIMIT, m_MAXLIMIT, m_MAXLIMIT}};
//    m_maxPoint        = {{-m_MAXLIMIT, -m_MAXLIMIT, -m_MAXLIMIT}};
//}
//
///*!
// * Default destructor for class BvNode.
// * Clear BvNode content and release memory.
// */
//BvNode::~BvNode(){}
//
///*! Assignement operator for class BvNode.
// */
//BvNode & BvNode::operator=(const BvNode & other)
//{
//    m_lchild        = other.m_lchild;
//    m_rchild         = other.m_rchild;
//    m_element        = other.m_element;
//    m_nrange        = other.m_nrange;
//    m_leaf            = other.m_leaf;
//    m_minPoint        = other.m_minPoint;
//    m_maxPoint        = other.m_maxPoint;
//    return *this;
//}
//
///*!
// * Default constructor for class BvTree.
// * Initialize an empty bv-tree structure and reserve memory for the insertion
// * of maxstack nodes.
// * If a linked patch is passed by argument, the contained information are used
// * to initialize the member of the tree.
// *
// *  \param[in] patch_ Pointer to the target bitpit::PatchKernel.
// *
// */
//BvTree::BvTree(bitpit::PatchKernel *patch_)
//{
//    m_patch    = patch_;
//    m_dim        = 3;
//    m_nnodes    = 0;
//    m_nleaf        = 0;
//    if (patch_ != NULL)
//    {
//        m_maxstk    = std::max(10, int(m_patch->getVertexCount()));
//        m_nelements = m_patch->getCellCount();
//    }
//    else
//    {
//        m_maxstk    = 10;
//        m_nelements = 0;
//    }
//    m_elements.resize(m_nelements);
//    increaseStack();
//
//    m_tol             = 1.0e-08;
//    m_maxsize    = 1;
//}
//
///*!
// * Default destructor for class BvTree.
// * Clear bv-tree content and release memory.
// */
//BvTree::~BvTree()
//{
//    m_patch = NULL;
//    clean();
//}
//
///*!
// * Assignement operator for class BvTree.
// */
//BvTree & BvTree::operator=(const BvTree & other)
//{
//    clean();
//    m_patch        = other.m_patch;
//    setup();
//    m_dim             = other.m_dim;
//    m_maxstk        = other.m_maxstk;
//    m_nodes            = other.m_nodes;
//    m_nnodes        = other.m_nnodes;
//    m_nleaf            = other.m_nleaf;
//    m_maxsize        = other.m_maxsize;
//    m_nelements        = other.m_nelements;
//    m_elements        = other.m_elements;
//    m_tol             = other.m_tol;
//    return *this;
//}
//
///*!
// * It sets the bitpit::PatchKernel linked by the tree. The Simplex of the patch
// * are ordered during the construction of the bv-tree.
// * \param[in] patch_ Pointer to the target bitpit::PatchKernel.
// */
//void BvTree::setPatch(bitpit::PatchKernel *patch_)
//{
//    m_patch     = NULL;
//    m_patch     = patch_;
//    m_dim         = 3;
//    if (patch_ != NULL) m_maxstk = std::max(10, int(m_patch->getVertexCount()));
//    m_nelements = 0;
//    m_nnodes     = 0;
//    m_nleaf        = 0;
//    m_nodes.clear();
//    m_elements.clear();
//    increaseStack();
//}
//
///*!
// * It sets the maximum number of elements in the leaf nodes of the tree.
// *  \param[in] maxsize Maximum number of elements in a leaf node (default =1).
// */
//void BvTree::setMaxLeafSize(int maxsize){
//    m_maxsize = maxsize;
//}
//
///*!
// * It builds the bv-tree starting from the root node (its bounding volume is the
// * AABB bounding volume of the whole geometry contained in the patch).
// */
//void BvTree::buildTree()
//{
//    if (m_patch == NULL) return;
//
//    if (m_nelements == 0 || m_nodes.size() == 0) setup();
//    if (m_nelements == 0) return;
//
//    int iel = 0;
//    long id;
//
//    //fill ids and centroid
//    for ( auto & cell : m_patch->getCells() )
//    {
//        id = cell.getId();
//        m_elements[iel].m_label = id;
//        m_elements[iel].m_centroid = m_patch->evalCellCentroid(id);
//        iel++;
//    }
//
//    //node root
//    BvNode node;
//    node.m_element[0] = 0;
//    node.m_element[1] = m_nelements;
//    node.m_nrange = m_nelements;
//    m_nodes[m_nnodes] = node;
//    if ( m_nelements == 1 )
//    {
//        m_nodes[m_nnodes].m_leaf = true;
//        computeBoundingBox(0);
//        m_nleaf++;
//        m_nnodes++;
//    }
//    else
//    {
//        m_nnodes++;
//        fillTree(0);
//    }
//    decreaseStack();
//}
//
///*!
// * It fills the bv-tree starting from an existing node.
// * \param[in] iparent Index of the node used to start the fill procedure.
// */
//void BvTree::fillTree(int iparent)
//{
//
//
//    //compute Bounding Box of node
//    computeBoundingBox(iparent);
//
//
//    if ( m_nodes[iparent].m_nrange > m_maxsize )
//    {
//
//        //split parent by mean of centroid in node
//        darray3E meanC = computeMeanPoint(m_nodes[iparent].m_element[0], m_nodes[iparent].m_element[1]);
//
//        //split by plane normal to maximum dimension of the centroid box
//        darray3E minPC{{m_MAXLIMIT, m_MAXLIMIT, m_MAXLIMIT}};
//        darray3E maxPC{{-m_MAXLIMIT, -m_MAXLIMIT, -m_MAXLIMIT}};
//        int ics = m_nodes[iparent].m_element[0], ice = m_nodes[iparent].m_element[1];
//        for (int ic=ics; ic<ice; ++ic )
//        {
//            for (int j=0; j<m_dim; ++j )
//            {
//                minPC[j] = std::min(minPC[j], m_elements[ic].m_centroid[j]);
//                maxPC[j] = std::max(maxPC[j], m_elements[ic].m_centroid[j]);
//            }
//        }
//
//        int dir = 0;
//        double maxRange = maxPC[dir] - minPC[dir];
//        for ( int idim=1; idim<m_dim; ++idim )
//        {
//            if ((maxPC[idim] - minPC[idim]) > maxRange)
//            {
//                dir = idim;
//                maxRange = maxPC[dir] - minPC[dir];
//            }
//        }
//
//        //Sort elements by coordinate dir.
//        std::vector<BvElement>::iterator itstart = m_elements.begin()+m_nodes[iparent].m_element[0];
//        std::vector<BvElement>::iterator itend = m_elements.begin()+(m_nodes[iparent].m_element[1]);
//
//        //Old algorithm
//        //        sort( itstart, itend, compareElements(dir) );
//        //        //find first element with coords[dir] > of meanC[dir]
//        //        int firstRight = findFirstGreater(iparent, meanC, dir);
//
//
//        int firstRight = pseudoSort( itstart, itend, meanC, dir );
//
//        //insert lchild
//        if ( firstRight > m_nodes[iparent].m_element[0]  )
//        {
//            m_nodes[iparent].m_lchild = m_nnodes;
//            if ( (int)m_nodes.size() <= m_nnodes+1) increaseStack();
//            BvNode node;
//            node.m_element[0] = m_nodes[iparent].m_element[0];
//            node.m_element[1] = firstRight;
//            node.m_nrange = node.m_element[1] - node.m_element[0];
//            m_nodes[m_nnodes] = node;
//            if ( node.m_nrange == 1 )
//            {
//                m_nodes[m_nnodes].m_leaf = true;
//                computeBoundingBox(m_nnodes);
//                m_nnodes++;
//                m_nleaf++;
//            }
//            else
//            {
//                m_nnodes++;
//                fillTree(m_nnodes-1);
//            }
//        }
//
//        //insert rchild
//        if ( firstRight < m_nodes[iparent].m_element[1] )
//        {
//            m_nodes[iparent].m_rchild = m_nnodes;
//            if ( (int)m_nodes.size() <= m_nnodes+1) increaseStack();
//            BvNode node;
//            node.m_element[0] = firstRight;
//            node.m_element[1] = m_nodes[iparent].m_element[1];
//            node.m_nrange = node.m_element[1] - node.m_element[0];
//            m_nodes[m_nnodes] = node;
//            if (node.m_nrange == 1){
//                m_nodes[m_nnodes].m_leaf = true;
//                computeBoundingBox(m_nnodes);
//                m_nnodes++;
//                m_nleaf++;
//            }
//            else{
//                m_nnodes++;
//                fillTree(m_nnodes-1);
//            }
//        }
//
//    }
//    else
//    {
//        //insert this node as leaf
//        m_nodes[iparent].m_leaf = true;
//        m_nleaf++;
//
//    }
//
//}
//
///*!
// * It computes the mean centroid of a set of elements given the starting index and
// * the final index of the elements in the m_elements structure.
// * \param[in] istart Index of the first element to be used.
// * \param[in] iend Index of the last element to be used.
// * \return Coordinates of the mean centroid.
// */
//darray3E BvTree::computeMeanPoint(int istart, int iend)
//{
//
//    long id;
//    darray3E meanC({{0.0, 0.0, 0.0}});
//    for ( int i=istart; i<iend; ++i )
//    {
//        id = m_elements[i].m_label;
//        meanC += m_patch->evalCellCentroid(id);
//    }
//    meanC /= double(iend-istart);
//    return meanC;
//}
//
///*!
// * Find the index of the first BvElement in a given list, whose "dir" centroid coordinate is
// * immediately at right of a a target mean centroid "dir "coordinate.
// * \param[in] itstart begin iterator of the target list
// * \param[in] itend end iterator of a target list
// * \param[in] meanC mean centroid
// * \param[in] dir centroid coordinate (0-x,1-y, 2-z)
// */
//int    BvTree::pseudoSort(std::vector<BvElement>::iterator itstart,
//        std::vector<BvElement>::iterator itend, darray3E meanC, int dir){
//
//    double thres = meanC[dir];
//
//    itend--;
//    while(itstart != itend){
//
//        while (itstart->m_centroid[dir] <= thres && itstart != itend){
//            itstart++;
//        }
//        while (itend->m_centroid[dir] > thres && itend != itstart){
//            itend--;
//        }
//
//        if (itstart != itend){
//            std::iter_swap(itstart, itend);
//        }
//    }
//
//    return (std::max(0, int(distance(m_elements.begin(), itend))));
//}
//
///*!
// * It finds the first element, contained din a target node, with a coordinate
// * of its centroid greater than a reference point (the mean centroid in the
// * bv-tree construction).
// * \param[in] inode Node index used to search.
// * \param[in] meanC Reference point coordinates (mean centroid in bv-tree fill procedure).
// * \param[in] dir Coordinate used to search the first greater point.
// * \return Index of the first element in the bounding volume
// * of the inode-th node with dir-th coordinate greater than meanC[dir].
// *
// */
//int BvTree::findFirstGreater(int inode, darray3E meanC, int dir)
//{
//    double mean = meanC[dir];
//    int guess = (m_nodes[inode].m_element[1] - m_nodes[inode].m_element[0])/2;
//    int step = guess;
//
//    while( step > 0 )
//    {
//        step /= 2;
//        guess += step * int(sign(mean - m_elements[guess].m_centroid[dir]));
//        guess = std::max(m_nodes[inode].m_element[0], guess);
//        guess = std::min(m_nodes[inode].m_element[1]-1, guess);
//    }
//
//    while ( mean < m_elements[guess].m_centroid[dir] && guess > m_nodes[inode].m_element[0] )
//    {
//        guess--;
//    }
//    while ( mean > m_elements[guess].m_centroid[dir] && guess < m_nodes[inode].m_element[1]-1 )
//    {
//        guess++;
//    }
//    return guess;
//}
//
///*!
// * It computes the bounding box coordinate of a given node.
// * \param[in] inode Index of the target node.
// *
// * After the call of the method the member m_minPpoint and m_maxPoint
// * of the inode.th node are filled by the coordinates of the extrema
// * points of the bounding volume that contains the elements
// * of the node.
// */
//void BvTree::computeBoundingBox(int inode)
//{
//    BvNode &node(m_nodes[inode]);
//
//    int is = node.m_element[0];
//    int ie = node.m_element[1];
//
//    int nV;
//
//    node.m_minPoint = {{m_MAXLIMIT, m_MAXLIMIT, m_MAXLIMIT}};
//    node.m_maxPoint = {{-m_MAXLIMIT, -m_MAXLIMIT, -m_MAXLIMIT}};
//
//    bitpit::Cell     cell;
//    darray3E         coords;
//
//    for (int i=is; i<ie; i++ )
//    {
//        cell = m_patch->getCell(m_elements[i].m_label);
//        nV = cell.getVertexCount();
//        for (int iV=0; iV<nV; ++iV )
//        {
//            coords = m_patch->getVertexCoords(cell.getVertex(iV));
//            for (int j=0; j<m_dim; ++j )
//            {
//                node.m_minPoint[j] = std::min(node.m_minPoint[j], coords[j]);
//                node.m_maxPoint[j] = std::max(node.m_maxPoint[j], coords[j]);
//            }
//        }
//    }
//    node.m_minPoint -= m_tol;
//    node.m_maxPoint += m_tol;
//}
//
///*!
// * It evaluates if a point is inside the bounding box
// * of a node broaden by a given quantity.
// * \param[in] P_ Pointer to the input point coordinates.
// * \param[in] node_ Pointer to the target node.
// * \param[in] r Length used to grow the bounding volume (each span of the volume
// * is expanded by r in each direction).
// * \return True if the point is inside the expanded bounding volume.
// *
// */
//bool BvTree::inBoundingBox(darray3E *P_, BvNode *node_, double r)
//{
//    bool in = true;
//    for ( int i=0; i<m_dim; i++ )
//    {
//        if ( ( (*P_)[i] < ( node_->m_minPoint[i] - r) ) || ( (*P_)[i] > ( node_->m_maxPoint[i] + r ) ) ){
//            in = false;
//            break;
//        }
//    }
//    return in;
//}
//
///*!
// * It evaluates if a sphere centered in a point intersects the bounding box
// * of a node broaden by a given quantity.
// * \param[in] P_ Pointer to the input point coordinates.
// * \param[in] node_ Pointer to the target node.
// * \param[in] r Radius of the sphere.
// * \return True if there is intersection between sphere and bounding box (or
// * the sphere is completely inside the bounding box).
// */
//bool BvTree::SphereBoundingBox(darray3E *P_, BvNode *node_, double r)
//{
//
//    darray3E closestP;
//    for ( int i=0; i<m_dim; i++ )
//    {
//        closestP[i] = std::min(std::max((*P_)[i], node_->m_minPoint[i]), node_->m_maxPoint[i]);
//    }
//    double dist = dotProduct(closestP-(*P_), closestP-(*P_));
//
//    return (dist < r*r);
//}
//
///*!
// * It cleans the bv-tree and release memory.
// */
//void BvTree::clean(){
//
//    m_dim        = 3;
//    m_nelements = 0;
//    m_nnodes    = 0;
//    m_maxstk    = 10;
//    m_nodes.clear();
//    m_elements.clear();
//    m_nleaf        = 0;
//}
//
///*!
// * It sets the bv-tree parameters by using the linked patch.
// */
//void BvTree::setup(){
//
//    m_dim        = 3;
//    m_nnodes    = 0;
//    m_nleaf        = 0;
//    if (m_patch != NULL)
//    {
//        m_maxstk    = std::max(10, int(m_patch->getVertexCount()));
//        m_nelements = m_patch->getCellCount();
//    }
//    else
//    {
//        m_maxstk    = 10;
//        m_nelements = 0;
//    }
//    m_elements.resize(m_nelements);
//    increaseStack();
//}
//
///*!
// * Whenever the bv-tree reaches its full capacity (i.e. the number of BvNode
// * stored in the tree is equal to the memory reserved), increase the memory
// * reserve by maxstack.
// */
//void BvTree::increaseStack()
//{
//    m_nodes.resize(m_nodes.size() + m_maxstk);
//}
//
///*!
// * Decrease the memory reserved for bv-tree by maxstack.
// * Resize the m_nodes to effective number of nodes.
// */
//void BvTree::decreaseStack()
//{
//    m_nodes.resize(m_nnodes);
//}



namespace skdTreeUtils{

/*!
 * It computes the unsigned distance of a point to a geometry linked in a BvTree
 * object. The geometry has to be a surface mesh, in particular an object of type
 * bitpit::SurfUnstructured (a static cast is hardly coded in the method).
 * It searches the element with minimum distance in a box of length 2*r, if none
 * is found return a default value of distance equal to 1.0e+18;
 * \param[in] P_ Pointer to coordinates of input point.
 * \param[in] bvtree_ Pointer to Boundary Volume Hierarchy tree that stores the geometry.
 * \param[out] id Label of the element found as minimum distance element in the bv-tree.
 * \param[in] r Length of the side of the box or radius of the sphere used to search. (The algorithm checks
 * every element encountered inside the box/sphere).
 * \return Unsigned distance of the input point from the patch in the bv-tree.
 */
double distance(std::array<double,3> *P_, bitpit::PatchSkdTree *bvtree_, long &id, double &r)
{

    double h;
    static_cast<bitpit::SurfaceSkdTree*>(bvtree_)->findPointClosestCell(*P_, r, &id, &h);
    return h;

}


/*!
 * It computes the signed distance of a point to a geometry linked in a BvTree
 * object. The geometry has to be a surface mesh, in particular an object of type
 * bitpit::SurfUnstructured (a static cast is hardly coded in the method).
 * The sign of the distance is provided by the normal to the geometry locally
 * computed. A positive distance means that the point is located, in respect to the
 * surface mesh, on the same side of the outer normal vector.
 * It searches the element with minimum distance in a box of length 2*r, if none
 * is found return a default value of distance equal to 1.0e+18;
 * \param[in] P_ Pointer to coordinates of input point.
 * \param[in] bvtree_ Pointer to Boundary Volume Hierarchy tree that stores the geometry.
 * \param[out] id Label of the element found as minimum distance element in the bv-tree.
 * \param[out] n Pseudo-normal of the element (i.e. unit vector with direction (P-xP),
 * where P is the input point and xP the nearest point on the element (simplex) to
 * the projection of P on the plane of the simplex.
 * \param[in] r Length of the side of the box or radius of the sphere used to search. (The algorithm checks
 * every element encountered inside the box/sphere).
 * \return Signed distance of the input point from the patch in the bv-tree.
 */
double signedDistance(std::array<double,3> *P_, bitpit::PatchSkdTree *bvtree_, long &id, std::array<double,3> &n, double &r)
{

    double h;
    static_cast<bitpit::SurfaceSkdTree*>(bvtree_)->findPointClosestCell(*P_, r, &id, &h);
    const bitpit::PatchKernel &patch = (bvtree_->getPatch());
    const bitpit::SurfUnstructured &spatch = static_cast<const bitpit::SurfUnstructured &>(patch);

    //signed distance only for patch with triangles or segments
    bitpit::Cell cell = spatch.getCell(id);
    int nV = cell.getVertexCount();
    long * connCC_ = cell.getConnect();
    dvecarr3E VS(nV);
    for (int iV=0; iV<nV; ++iV)
    {
        VS[iV] = spatch.getVertexCoords(connCC_[iV]);
    }
    if ( nV == 3 )
    {
        darray3E lambda ;
        darray3E xP, normal;
        h = bitpit::CGElem::distancePointTriangle((*P_), VS[0], VS[1], VS[2],lambda);
        normal  = lambda[0] * spatch.evalVertexNormal(id,0) ;
        normal += lambda[1] * spatch.evalVertexNormal(id,1) ;
        normal += lambda[2] * spatch.evalVertexNormal(id,2) ;
        xP = lambda[0]*VS[0] + lambda[1]*VS[1] + lambda[2]*VS[2];
        double s =  sign( dotProduct(normal, (*P_) - xP) );
        h = s * h;
        //pseudo-normal (direction P and xP closest point on triangle)
        n = s * ((*P_) - xP);
        double normX = norm2(n);
        if(normX < 1.E-15){
            n = normal/norm2(normal);
        }else{
            n /= norm2(n);
        }
    }
    else if ( nV == 2 )
    {
        darray2E lambda ;
        darray3E xP, normal;
        h = bitpit::CGElem::distancePointSegment((*P_), VS[0], VS[1], lambda);
        normal  = lambda[0] * spatch.evalVertexNormal(id,0) ;
        normal += lambda[1] * spatch.evalVertexNormal(id,1) ;
        xP = lambda[0]*VS[0] + lambda[1]*VS[1];
        double s = sign( dotProduct(normal, (*P_) - xP) );
        h = s * h;
        //pseudo-normal (direction P and xP closest point on triangle)
        n = s * ((*P_) - xP);
        double normX = norm2(n);
        if(normX < 1.E-15){
            n = normal/norm2(normal);
        }else{
            n /= norm2(n);
        }
    }

    return h;

}

/*!
 * It selects the elements of a geometry stored in a bv-tree by a distance criterion
 * in respect to an other geometry stored in a different bv-tree.
 * \param[in] selection Pointer to bv-tree used as selection patch.
 * \param[in] target Pointer to bv-tree that store the target geometry.
 * \param[in] tol Distance threshold used to select the elements of target.
 * \return Vector of the label of all the elements of the target bv-tree placed
 * at a distance <= tol from the bounding boxes of the leaf nodes of the bv-tree
 * selection.
 */
std::vector<long> selectByPatch(bitpit::PatchSkdTree *selection, bitpit::PatchSkdTree *target, double tol){

    int nleafs = selection->getLeafCount();
    int nnodess = selection->getNodeCount();
    std::vector<const bitpit::SkdNode*> leafSelection(nleafs);
    int count = 0;
    for (int i=0; i<nnodess; i++){
        if (selection->getNode(i).isLeaf()){
            if (bitpit::CGElem::intersectBoxBox(selection->getNode(i).getBoxMin()-tol,
                    selection->getNode(i).getBoxMax()+tol,
                    target->getNode(0).getBoxMin(),
                    target->getNode(0).getBoxMax() ) ){
                leafSelection[count] = &(selection->getNode(i));
                count++;

            }
        }
    }
    leafSelection.resize(count);


    std::vector<long> extracted;
    extractTarget(target, leafSelection, extracted, tol);

    return extracted;

}

/*!
 * It extracts the elements of a leaf node of geometry stored in a bv-tree
 * by a distance criterion in respect to an other geometry stored
 * in a different bv-tree. It is a recursive method used in selectByPatch method.
 * \param[in] target Pointer to bv-tree that store the target geometry.
 * \param[in,out] leafSelection Vector of pointers to the leaf nodes currently interesting
 * for the selection procedure.
 * \param[in,out] extracted of the label of all the elements of the target bv-tree,
 * currently found placed at a distance <= tol from the bounding boxes of the
 * leaf nodes in leafSelection.
 * \param[in] tol Distance threshold used to select the elements of target.
 * \param[in] next Index of the node of the bv-tree target to be checked. If
 * the next-th node is not a leaf node the method is recursively called.
 *
 * The leafSelection is (potentially) reduced during the procedure by deleting
 * from the vector the nodes that have a distance greater than tol from the bounding
 * volume of the next-th node.
 *
 */
void extractTarget(bitpit::PatchSkdTree *target, std::vector<const bitpit::SkdNode*> leafSelection, std::vector<long> &extracted, double tol, int next){

    bool check = false;
    std::vector<const bitpit::SkdNode*> tocheck;
    for (int i=0; i<(int)leafSelection.size(); i++){
        if (bitpit::CGElem::intersectBoxBox(leafSelection[i]->getBoxMin()-tol,
                leafSelection[i]->getBoxMax()+tol,
                target->getNode(next).getBoxMin(),
                target->getNode(next).getBoxMax() ) ){
            check = true;
            tocheck.push_back(leafSelection[i]);
        }
    }

    leafSelection.clear();
    leafSelection = tocheck;

    if (check){
        bitpit::SkdNode node = target->getNode(next);
        if (node.isLeaf()){
            std::vector<long> cellids = node.getCells();
            extracted.insert(extracted.end(),cellids.begin(), cellids.end());
        }
        else{
            for (int i = bitpit::SkdNode::CHILD_BEGIN; i != bitpit::SkdNode::CHILD_END; ++i) {
                bitpit::SkdNode::ChildLocation childLocation = static_cast<bitpit::SkdNode::ChildLocation>(i);
                if(node.getChildId(childLocation) != bitpit::SkdNode::NULL_ID)   {
                    extractTarget(target, leafSelection, extracted, tol, node.getChildId(childLocation));
                }
            }
        }
    }

}

/*!
 * It computes the projection of a point on a geometry linked in a BvTree
 * object. The geometry has to be a surface mesh, in particular an object of type
 * bitpit::SurfUnstructured (a static cast is hardly coded in the method).
 * It searches the elements of the geometry with minimum distance
 * recursively in a sphere of radius r, by increasing the size r at each step
 * until at least one element is found.
 * \param[in] P_ Pointer to coordinates of input point.
 * \param[in] bvtree_ Pointer to Boundary Volume Hierarchy tree that stores the geometry.
 * \param[in] r Initial length of the sphere radius used to search. (The algorithm checks
 * every element encountered inside the sphere).
 * \return Coordinates of the projected point.
 */
darray3E projectPoint( std::array<double,3> *P_, bitpit::PatchSkdTree *bvtree_, double r )
{

    long         id;
    darray3E     normal;

    double         dist = 1.0e+18;
    while (std::abs(dist) >= 1.0e+18){
        //use method sphere by default
        dist = signedDistance(P_, bvtree_, id, normal, r);
        r *= 1.5;
    }

    darray3E     projP = (*P_) - dist*normal;

    return (projP);
}




///*!
// * It computes the signed distance of a set of points to a geometry linked in a BvTree
// * object. The geometry has to be a surface mesh, in particular an object of type
// * bitpit::SurfUnstructured (a static cast is hardly coded in the method).
// * The sign of the distance is provided by the normal to the geometry locally
// * computed. A positive distance means that a point is located, in respect to the
// * surface mesh, on the same side of the outer normal vector.
// * It searches the elements with minimum distance in a box of length 2*r or sphere of radius r alternatively. If none
// * is found return a default value of distance equal to 1.0e+18;
// * \param[in] P_ Pointer to a vector with the coordinates of input points.
// * \param[in] bvtree_ Pointer to Boundary Volume Hierarchy tree that stores the geometry.
// * \param[out] id Vector of labels of the elements found as minimum distance elements in the bv-tree.
// * \param[out] n Vector of pseudo-normals of the elements (i.e. unit vectors with direction (P-xP),
// * where P is an input point and xP the nearest point on the element (simplex) to
// * the projection of P on the plane of the simplex.
// * \param[in] r_ Length of the side of the box or radius of the sphere used to search. (The algorithm checks
// * every element encountered inside the box/sphere).
// * \param[in] method Method used to search the element (0=bounding box, 1=sphere).
// * \return Vector with signed distances of the input points from the patch in the bv-tree.
// */
//dvector1D signedDistance(std::vector<std::array<double,3> > *P_, BvTree *bvtree_, std::vector<long> &id, std::vector<std::array<double,3> > &n, double r_, int method)
//{
//    bitpit::SurfUnstructured *spatch_ = static_cast<bitpit::SurfUnstructured*>(bvtree_->m_patch);
//
//    double         r = r_;
//    int            nP = P_->size();
//    dvector1D     dist(nP);
//
//    for (int i = 0; i < nP; ++i){
//        dist[i] = signedDistance(&((*P_)[i]), bvtree_, id[i], n[i], r, method, spatch_);
//    }
//    return dist;
//}
//
///*!
// * It computes the unsigned distance of a set of points to a geometry linked in a BvTree
// * object. The geometry has to be a surface mesh, in particular an object of type
// * bitpit::SurfUnstructured (a static cast is hardly coded in the method).
// * It searches the elements with minimum distance in a box of length 2*r or sphere of radius r alternatively. If none
// * is found return a default value of distance equal to 1.0e+18;
// * \param[in] P_ Pointer to vector with the coordinates of input points.
// * \param[in] bvtree_ Pointer to Boundary Volume Hierarchy tree that stores the geometry.
// * \param[out] id Vector with labels of the elements found as minimum distance elements in the bv-tree.
//  * \param[in] r_ Length of the side of the box or radius of the sphere used to search. (The algorithm checks
// * every element encountered inside the box/sphere).
// * \param[in] method Method used to search the element (0=bounding box, 1=sphere).
// * \return Vector with unsigned distances of the input points from the patch in the bv-tree.
// */
//dvector1D distance(std::vector<std::array<double,3> > *P_, BvTree *bvtree_, std::vector<long> &id, double r_, int method)
//{
//
//    double         r = r_;
//    int            nP = P_->size();
//    dvector1D     dist(nP);
//
//    for (int i = 0; i < nP; ++i){
//        dist[i] = distance(&((*P_)[i]), bvtree_, id[i], r, method);
//    }
//    return dist;
//}
//
///*!
// * It computes the projection of a set of points on a geometry linked in a BvTree
// * object. The geometry has to be a surface mesh, in particular an object of type
// * bitpit::SurfUnstructured (a static cast is hardly coded in the method).
// * It searches the elements of the geometry with minimum distance
// * recursively in a box of length 2*r or sphere of radius r, by increasing the size r at each step
// * until at least one element is found.
// * \param[in] P_ Pointer to vector with coordinates of input points.
// * \param[in] bvtree_ Pointer to Boundary Volume Hierarchy tree that stores the geometry.
// * \param[in] r_ Initial length of the sphere radius used to search. (The algorithm checks
// * every element encountered inside the sphere).
// * \return Vector with coordinates of the projected points.
// */
//dvecarr3E projectPoint(std::vector<std::array<double,3> > *P_, BvTree *bvtree_, double r_)
//{
//
//    int            nP = P_->size();
//    dvector1D     dist(nP);
//    dvecarr3E     projPoint(nP);
//
//    for (int i = 0; i < nP; ++i){
//        projPoint[i] = projectPoint(&(*P_)[i], bvtree_, r_);
//    }
//    return projPoint;
//}
//
///*!
// * It selects the elements of a geometry stored in a bv-tree by a distance criterion
// * in respect to an other geometry stored in a different bv-tree.
// * \param[in] selection Pointer to bv-tree used as selection patch.
// * \param[in] target Pointer to bv-tree that store the target geometry.
// * \param[in] tol Distance threshold used to select the elements of target.
// * \return Vector of the label of all the elements of the target bv-tree placed
// * at a distance <= tol from the bounding boxes of the leaf nodes of the bv-tree
// * selection.
// */
//std::vector<long> selectByPatch(BvTree *selection, BvTree *target, double tol){
//
//    int nleafs = selection->m_nleaf;
//    int nnodess = selection->m_nnodes;
//    std::vector<BvNode*> leafSelection(nleafs);
//    int count = 0;
//    for (int i=0; i<nnodess; i++){
//        if (selection->m_nodes[i].m_leaf){
//            if (bitpit::CGElem::intersectBoxBox(selection->m_nodes[i].m_minPoint-tol,
//                    selection->m_nodes[i].m_maxPoint+tol,
//                    target->m_nodes[0].m_minPoint,
//                    target->m_nodes[0].m_maxPoint ) ){
//                leafSelection[count] = &(selection->m_nodes[i]);
//                count++;
//
//            }
//        }
//    }
//    leafSelection.resize(count);
//
//
//    std::vector<long> extracted;
//    extractTarget(target, leafSelection, extracted, tol);
//
//    return extracted;
//
//
//}

}

}; // end namespace mimmo
