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
# include "MimmoCGUtils.hpp"
# include <bitpit_surfunstructured.hpp>
# include <surface_skd_tree.hpp>
# include <CG.hpp>

namespace mimmo{

namespace skdTreeUtils{

/*!
 * It computes the unsigned distance of a point to a geometry linked in a SkdTree
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

    double h = 1.E+18;
    if(!bvtree_ ){
        throw std::runtime_error("Invalid use of skdTreeUtils::distance method: a void tree is detected.");
    }
    if(!dynamic_cast<const bitpit::SurfUnstructured*>(&(bvtree_->getPatch()))){
        throw std::runtime_error("Invalid use of skdTreeUtils::distance method: a not surface patch tree is detected.");
    }

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

    double h = 1.E+18;

    if(!bvtree_ ){
    	throw std::runtime_error("Invalid use of skdTreeUtils::signedDistance method: a void tree is detected.");
    }
    const bitpit::SurfUnstructured *spatch = dynamic_cast<const bitpit::SurfUnstructured*>(&(bvtree_->getPatch()));
    if(!spatch){
    	throw std::runtime_error("Invalid use of skdTreeUtils::signedDistance method: a not surface patch tree is detected.");
    }

    static_cast<bitpit::SurfaceSkdTree*>(bvtree_)->findPointClosestCell(*P_, r, &id, &h);

    //signed distance only for 2D element patches(quads, pixels, triangles or segments)
    if (id != bitpit::Cell::NULL_ID){
    	const bitpit::Cell & cell = spatch->getCell(id);
    	bitpit::ConstProxyVector<long> vertIds = cell.getVertexIds();
    	dvecarr3E VS(vertIds.size());
    	int count = 0;
    	for (const auto & iV: vertIds){
    		VS[count] = spatch->getVertexCoords(iV);
    		++count;
    	}

    	darray3E xP = {{0.0,0.0,0.0}};
    	darray3E normal= {{0.0,0.0,0.0}};

    	if ( vertIds.size() == 3 ){ //TRIANGLE
    		darray3E lambda;
    		h = bitpit::CGElem::distancePointTriangle((*P_), VS[0], VS[1], VS[2],lambda);
    		int count = 0;
    		for(const auto &val: lambda){
    			normal += val * spatch->evalVertexNormal(id,count) ;
    			xP += val * VS[count];
    			++count;
    		}
    	}else if ( vertIds.size() == 2 ){ //LINE/SEGMENT
    		darray2E lambda;
    		h = bitpit::CGElem::distancePointSegment((*P_), VS[0], VS[1], lambda);
    		int count = 0;
    		for(const auto &val: lambda){
    			normal += val * spatch->evalVertexNormal(id,count) ;
    			xP += val * VS[count];
    			++count;
    		}
    	}else{ //GENERAL POLYGON
    		std::vector<double> lambda;
    		h = bitpit::CGElem::distancePointPolygon((*P_), VS,lambda);
    		int count = 0;
    		for(const auto &val: lambda){
    			normal += val * spatch->evalVertexNormal(id,count) ;
    			xP += val * VS[count];
    			++count;
    		}
    	}

    	double s =  sign( dotProduct(normal, (*P_) - xP) );
    	if(s == 0.0)    s =1.0;
    	h = s * h;
    	//pseudo-normal (direction P and xP closest point on triangle)
    	n = s * ((*P_) - xP);
    	double normX = norm2(n);
    	if(normX < 1.E-15){
    		n = normal/norm2(normal);
    	}else{
    		n /= norm2(n);
    	}
    }//end if not id null
    return h;

}

/*!
 * It selects the elements of a geometry stored in a skdtree by a distance criterion
 * in respect to an other geometry stored in a different skdtree.
 * \param[in] selection Pointer to bv-tree used as selection patch.
 * \param[in] target Pointer to bv-tree that store the target geometry.
 * \param[in] tol Distance threshold used to select the elements of target.
 * \return Vector of the label of all the elements of the target skdtree placed
 * at a distance <= tol from the bounding boxes of the leaf nodes of the skdtree
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
 * It extracts the elements of a leaf node of geometry stored in a skdtree
 * by a distance criterion in respect to an other geometry stored
 * in a different skdtree. It is a method used in selectByPatch method.
 * \param[in] target Pointer to skdtree that store the target geometry.
 * \param[in,out] leafSelection Vector of pointers to the leaf nodes currently interesting
 * for the selection procedure.
 * \param[in,out] extracted of the label of all the elements of the target skdtree,
 * currently found placed at a distance <= tol from the bounding boxes of the
 * leaf nodes in leafSelection.
 * \param[in] tol Distance threshold used to select the elements of target.
 * the next-th node is not a leaf node the method is recursively called.
 *
 *
 */
void extractTarget(bitpit::PatchSkdTree *target, const std::vector<const bitpit::SkdNode*> & leafSelection, std::vector<long> &extracted, double tol){

    if(leafSelection.empty())   return;
    std::size_t rootId  =0;

    std::vector<std::size_t> candidates;
    std::vector<std::pair< std::size_t, std::vector<const bitpit::SkdNode*> > > nodeStack;

    std::vector<const bitpit::SkdNode*> tocheck;
    for (int i=0; i<(int)leafSelection.size(); i++){
        if (bitpit::CGElem::intersectBoxBox(leafSelection[i]->getBoxMin()-tol,
                                            leafSelection[i]->getBoxMax()+tol,
                                            target->getNode(rootId).getBoxMin(),
                                            target->getNode(rootId).getBoxMax() ) )
        {
            tocheck.push_back(leafSelection[i]);
        }
    }
    nodeStack.push_back(std::make_pair(rootId, tocheck) );

    while(!nodeStack.empty()){

        std::pair<std::size_t,  std::vector<const bitpit::SkdNode*> >  touple = nodeStack.back();
        const bitpit::SkdNode & node = target->getNode(touple.first);
        nodeStack.pop_back();


        bool isLeaf = true;
        for (int i = bitpit::SkdNode::CHILD_BEGIN; i != bitpit::SkdNode::CHILD_END; ++i) {
            bitpit::SkdNode::ChildLocation childLocation = static_cast<bitpit::SkdNode::ChildLocation>(i);
            std::size_t childId = node.getChildId(childLocation);
            if (childId != bitpit::SkdNode::NULL_ID) {
                isLeaf = false;
                tocheck.clear();
                for (int i=0; i<(int)touple.second.size(); i++){
                    if (bitpit::CGElem::intersectBoxBox(touple.second[i]->getBoxMin()-tol,
                                                        touple.second[i]->getBoxMax()+tol,
                                                        target->getNode(childId).getBoxMin(),
                                                        target->getNode(childId).getBoxMax() ) )
                    {
                        tocheck.push_back(touple.second[i]);
                    }
                }
                if(!tocheck.empty())    nodeStack.push_back(std::make_pair(childId, tocheck) );
            }
        }

        if (isLeaf) {
            candidates.push_back(touple.first);
        }
    }

    std::unordered_set<long> cellExtracted;
    for(const auto & nodeId: candidates){
        const bitpit::SkdNode & node = target->getNode(nodeId);
        std::vector<long> cellids = node.getCells();
        cellExtracted.insert(cellids.begin(), cellids.end());
    }

    extracted.insert(extracted.end(), cellExtracted.begin(), cellExtracted.end());
}

/*!
 * It computes the projection of a point on a geometry linked in a skdtree
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

    r = std::max(r, bvtree_->getPatch().getTol());

    double         dist = 1.0e+18;
    while (std::abs(dist) >= 1.0e+18){
        //use method sphere by default
        dist = signedDistance(P_, bvtree_, id, normal, r);
        r *= 1.5;
    }

    darray3E     projP = (*P_) - dist*normal;

    return (projP);
}

/*!
 * Given the specified point find the cell of a surface patch it is into.
 * The method works only with trees generated with bitpit::SurfUnstructured mesh.
 *
 * \param[in] point is the point
 * \param[in] tree reference to SkdTree relative to the target surface geometry.
 * \return id of the geometry cell the point is into. Return bitpit::Cell::NULL_ID if no cell is found.
 */
long locatePointOnPatch(const std::array<double, 3> &point, bitpit::PatchSkdTree &tree)
{
    if(!dynamic_cast<const bitpit::SurfUnstructured*>(&(tree.getPatch()))){
        throw std::runtime_error("Invalid use of skdTreeUtils::locatePointOnPatch method: a non surface patch or void patch was detected.");
    }
    // Initialize the cell id
    long id = bitpit::Cell::NULL_ID;

    // Initialize the distance with an estimate
    //
    // The real distance will be lesser than or equal to the estimate.
    std::size_t rootId = 0;
    const bitpit::SkdNode &root = tree.getNode(rootId);
    double distance = root.evalPointMaxDistance(point);

    // Get a list of candidates nodes
    //
    std::vector<std::size_t>    m_candidateIds;

    std::vector<std::size_t> nodeStack;
    nodeStack.push_back(rootId);
    while (!nodeStack.empty()) {
        std::size_t nodeId = nodeStack.back();
        const bitpit::SkdNode &node = tree.getNode(nodeId);
        nodeStack.pop_back();

        // Do not consider nodes with a minimum distance greater than
        // the distance estimate
        double nodeMinDistance = node.evalPointMinDistance(point);
        if (nodeMinDistance > distance) {
            continue;
        }

        // Update the distance estimate
        //
        // The real distance will be lesser than or equal to the
        // estimate.
        double nodeMaxDistance = node.evalPointMaxDistance(point);
        distance = std::min(nodeMaxDistance, distance);

        // If the node is a leaf add it to the candidates, otherwise
        // add its children to the stack.
        bool isLeaf = true;
        for (int i = bitpit::SkdNode::CHILD_BEGIN; i != bitpit::SkdNode::CHILD_END; ++i) {
            bitpit::SkdNode::ChildLocation childLocation = static_cast<bitpit::SkdNode::ChildLocation>(i);
            std::size_t childId = node.getChildId(childLocation);
            if (childId != bitpit::SkdNode::NULL_ID) {
                isLeaf = false;
                nodeStack.push_back(childId);
            }
        }

        if (isLeaf) {
            m_candidateIds.push_back(nodeId);
        }
    }

    // Process the candidates and find which cell contains the point
    for (std::size_t k = 0; k < m_candidateIds.size(); ++k) {

        // Evaluate the min distance between all cells in the candidate node
        std::size_t nodeId = m_candidateIds[k];
        const bitpit::SkdNode &node = tree.getNode(nodeId);

        bitpit::ElementType eltype ;
        bool checkBelong;
        for(const auto & cellId : node.getCells() ){
            eltype = tree.getPatch().getCell(cellId).getType();
            auto vList = tree.getPatch().getCell(cellId).getVertexIds();
            dvecarr3E vertCoords(vList.size());
            int counterList = 0;
            for(const auto & idV : vList){
                vertCoords[counterList] = tree.getPatch().getVertexCoords(idV);
                ++counterList;
            }
            switch(eltype){
                case bitpit::ElementType::LINE:
                    checkBelong = mimmoCGUtils::isPointInsideSegment(point, vertCoords[0], vertCoords[1]);
                    break;
                case bitpit::ElementType::TRIANGLE:
                    checkBelong = mimmoCGUtils::isPointInsideTriangle(point, vertCoords[0], vertCoords[1], vertCoords[2]);
                    break;
                case bitpit::ElementType::PIXEL:
                case bitpit::ElementType::QUAD:
                case bitpit::ElementType::POLYGON:
                    checkBelong = mimmoCGUtils::isPointInsidePolygon(point, vertCoords);
                    break;
                default:
                    throw std::runtime_error("Not supported cell type");
                    break;
            }
            if (checkBelong) {
                id       = cellId;
            }
        }
    }
    return id;
}

/*!
 * Try to return the closest cell of a target mesh to a prescribed orientation point
 * it works with all those meshes who support a SkdTree.
 *
 * \param[in] point is the target point
 * \param[in] tree reference to SkdTree relative to the target mesh (can be a surface or a volume).
 * \return id of the geometry cell the point is close. Return bitpit::Cell::NULL_ID if no cell is found.
 */
long closestCellToPoint(const std::array<double, 3> &point, bitpit::PatchSkdTree &tree)
{
    // Initialize the distance with an estimate
    //
    // The real distance will be lesser than or equal to the estimate.
    std::size_t rootId = 0;
    const bitpit::SkdNode &root = tree.getNode(rootId);
    double distance = root.evalPointMaxDistance(point);

    // Get a list of candidates nodes
    //
    std::vector<std::size_t>    m_candidateIds;

    std::vector<std::size_t> nodeStack;
    nodeStack.push_back(rootId);
    while (!nodeStack.empty()) {
        std::size_t nodeId = nodeStack.back();
        const bitpit::SkdNode &node = tree.getNode(nodeId);
        nodeStack.pop_back();

        // Do not consider nodes with a minimum distance greater than
        // the distance estimate
        double nodeMinDistance = node.evalPointMinDistance(point);
        if (nodeMinDistance > distance) {
            continue;
        }

        // Update the distance estimate
        //
        // The real distance will be lesser than or equal to the
        // estimate.
        double nodeMaxDistance = node.evalPointMaxDistance(point);
        distance = std::min(nodeMaxDistance, distance);

        // If the node is a leaf add it to the candidates, otherwise
        // add its children to the stack.
        bool isLeaf = true;
        for (int i = bitpit::SkdNode::CHILD_BEGIN; i != bitpit::SkdNode::CHILD_END; ++i) {
            bitpit::SkdNode::ChildLocation childLocation = static_cast<bitpit::SkdNode::ChildLocation>(i);
            std::size_t childId = node.getChildId(childLocation);
            if (childId != bitpit::SkdNode::NULL_ID) {
                isLeaf = false;
                nodeStack.push_back(childId);
            }
        }

        if (isLeaf) {
            m_candidateIds.push_back(nodeId);
        }
    }

    // Process the candidates and find which cell is the nearest to the target point
    double distanceMin = std::numeric_limits<double>::max();
    distance = 1.0E18;
    // Initialize the cell id
    long id = bitpit::Cell::NULL_ID;
    long idwork;

    for (std::size_t cand : m_candidateIds) {
        tree.getNode(cand).findPointClosestCell(point, &idwork, &distance);
        if(distance < distanceMin){
            distanceMin = distance;
            id = idwork;
        }
    }
    return id;
}



}

}; // end namespace mimmo
