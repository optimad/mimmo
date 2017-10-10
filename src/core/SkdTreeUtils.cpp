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
 * It computes the signed distance of a point to a geometry linked in a SkdTree
 * object. The geometry must be a surface mesh, in particular an object of type
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

    //signed distance only for 2D element patches(quads, pixels, triangles or segments)
    const bitpit::Cell & cell = spatch.getCell(id);
    bitpit::ConstProxyVector<long> vertIds = cell.getVertexIds();
    dvecarr3E VS(vertIds.size());
    int count = 0;
    for (const auto & iV: vertIds){
        VS[count] = spatch.getVertexCoords(iV);
        ++count;
    }
 
    darray3E xP = {{0.0,0.0,0.0}};
    darray3E normal= {{0.0,0.0,0.0}};

    if ( vertIds.size() == 3 ){ //TRIANGLE
        darray3E lambda;
        h = bitpit::CGElem::distancePointTriangle((*P_), VS[0], VS[1], VS[2],lambda);
        int count = 0;
        for(const auto &val: lambda){
            normal += val * spatch.evalVertexNormal(id,count) ;
            xP += val * VS[count];
            ++count;
        }
    }else if ( vertIds.size() == 2 ){ //LINE/SEGMENT
        darray2E lambda;
        h = bitpit::CGElem::distancePointSegment((*P_), VS[0], VS[1], lambda);
        int count = 0;
        for(const auto &val: lambda){
            normal += val * spatch.evalVertexNormal(id,count) ;
            xP += val * VS[count];
            ++count;
        }
    }else{ //GENERAL POLYGON
        std::vector<double> lambda;
        h = bitpit::CGElem::distancePointPolygon((*P_), VS,lambda);
        int count = 0;
        for(const auto &val: lambda){
            normal += val * spatch.evalVertexNormal(id,count) ;
            xP += val * VS[count];
            ++count;
        }
    }

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
 * object. The geometry must be a surface mesh, in particular an object of type
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
