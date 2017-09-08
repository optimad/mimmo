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
 \ *---------------------------------------------------------------------------*/
#include "BasicShapes.hpp"
#include "customOperators.hpp"
#include "SortAlgorithms.hpp"
#include <algorithm>

using namespace bitpit;

namespace mimmo{


/*! 
 * Basic Constructor 
 */
BasicShape::BasicShape(){

	m_shape = ShapeType::CUBE;
	m_origin.fill(0.0);
	m_span.fill(0.0);
	m_infLimits.fill(0.0);
	for(int i=0; i<3; ++i){
		m_sdr[i].fill(0.0);
	}
	m_sdr[0][0] = m_sdr[1][1] = m_sdr[2][2] = 1.0;
	m_typeCoord.fill(CoordType::CLAMPED);
	m_scaling.fill(1.0);
};

/*! Basic Destructor */
BasicShape::~BasicShape(){};

/*!
 * Swap method, assign data of the current class to an external 
 * one of the same type and vice-versa
 * 
 * \param[in] x BasicShape object to be swapped;
 */
void BasicShape::swap(BasicShape & x) noexcept
{
    std::swap(m_shape, x.m_shape);
    std::swap(m_origin, x.m_origin);
    std::swap(m_span, x.m_span);
    std::swap(m_infLimits, x.m_infLimits);
    std::swap(m_sdr, x.m_sdr);
    std::swap(m_scaling, x.m_scaling);
    std::swap(m_typeCoord, x.m_typeCoord);
}

/*! 
 * Set origin of your shape. The origin is meant as the baricenter of your shape in 
 * absolute reference system.
 * \param[in] origin new origin point
 */
void BasicShape::setOrigin(darray3E origin){
	m_origin = origin;
}

/*! 
 * Set span of your shape, according to its local reference system 
 * \param[in] s0 first coordinate span
 * \param[in] s1 second coordinate span
 * \param[in] s2 third coordinate span
 */
void BasicShape::setSpan(double s0, double s1, double s2){
	checkSpan(s0,s1,s2);
	setScaling(s0,s1,s2);
}

/*! 
 * Set span of your shape, according to its local reference system
 * \param[in] span coordinate span
 */
void BasicShape::setSpan(darray3E span){
	checkSpan(span[0],span[1],span[2]);
	setScaling(span[0],span[1],span[2]);
}

/*! 
 * Set inferior limits of your shape, in its local reference system.
 *  The method is useful when drawing part of your original shape, as in the 
 *  case of cylinders or sphere portions, by manipulating the adimensional 
 *  curvilinear coordinate.
 *     
 * \param[in] orig first coordinate origin
 * \param[in] dir 0,1,2 int flag identifying coordinate
 */
void BasicShape::setInfLimits(double orig, int dir){
	if(dir<0 || dir>2) return;
	
	darray3E span = getSpan();
	bool check = checkInfLimits(orig,dir);
	if(check){
		m_infLimits[dir] = orig; 
		setSpan(span[0],span[1],span[2]);
	}
}

/*! 
 * Set inferior limits your shape, in its local reference system.
 *  The method is useful when drawing part of your original shape, as in the 
 *  case of cylinders or sphere portions, by manipulating the adimensional 
 *  curvilinear coordinate.
 *     
 * \param[in] val lower value origin for all three coordinates
 */
void BasicShape::setInfLimits(darray3E val){
	for (int dir=0; dir<3; dir++){
		setInfLimits(val[dir], dir);
	}
}

/*! 
 * Set new axis orientation of the local reference system
 * \param[in] axis0 first axis
 * \param[in] axis1 second axis
 * \param[in] axis2 third axis
 * 
 * if chosen axes are not orthogonal, doing nothing
 */
void BasicShape::setRefSystem(darray3E axis0, darray3E axis1, darray3E axis2){
	
	axis0 = axis0/norm2(axis0);
	axis1 = axis1/norm2(axis1);
	axis2 = axis2/norm2(axis2);
	
	double tol = 1.0e-12;
	double check = std::abs(dotProduct(axis0,axis1)) + std::abs(dotProduct(axis1,axis2)) + std::abs(dotProduct(axis0,axis2));
	if(check > tol) return;
	m_sdr[0] = axis0;
	m_sdr[1] = axis1;
	m_sdr[2] = axis2;
}

/*! 
 * Set new axis orientation of the local reference system
 * \param[in] label 0,1,2 identify local x,y,z axis of the primitive shape
 * \param[in] axis new direction of selected local axis.
 */
void BasicShape::setRefSystem(int label, darray3E axis){
	
	if(label <0 || label >2 ) return;
	
	m_sdr[label] = axis/norm2(axis);
	dvecarr3E point_mat(3,darray3E{0,0,0});
	point_mat[0][0] = point_mat[1][1]= point_mat[2][2]=1.0;
	
	int next_label = (label + 1)%3;
	int fin_label = (label + 2)%3;
	
	double pj = dotProduct(point_mat[next_label], m_sdr[label]);
	m_sdr[next_label] = point_mat[next_label] - pj*m_sdr[label];
	m_sdr[next_label] = m_sdr[next_label]/norm2(m_sdr[next_label]);
	
	m_sdr[fin_label] = crossProduct(m_sdr[label],m_sdr[next_label]);
	m_sdr[fin_label] = m_sdr[fin_label]/norm2(m_sdr[fin_label]);
}

/*! 
 * Set new axes orientation of the local reference system
 * \param[in] axes new direction of local axes.
 */
void BasicShape::setRefSystem(dmatrix33E axes){
	for (int i=0; i<3; i++){
		m_sdr[i] = axes[i]/norm2(axes[i]);
	}
}

/*! 
 * Set type to treat your shape coordinates. 
 * \param[in] type CoordType enum.
 * \param[in] dir  0,1,2 int flag identifying coordinate
 */
void BasicShape::setCoordinateType(CoordType type, int dir){
	m_typeCoord[dir] = type;
}

/*!
 * \return current origin of your shape
 */
darray3E BasicShape::getOrigin(){
	return(m_origin);
}

/*! 
 * \return current span of your shape
 * 
 */
darray3E BasicShape::getSpan(){
	darray3E result = getLocalSpan();
	darray3E scale = getScaling();
	for(int i=0; i<3; ++i){
		result[i] = scale[i]*result[i];
	}
	return(result);
}

/*! 
 * \return current inferior limits of your shape, in local coord reference system
 */
darray3E BasicShape::getInfLimits(){
	return(m_infLimits);
}

/*! 
 * \return actual axis of global relative sdr
 */
dmatrix33E BasicShape::getRefSystem(){
	return(m_sdr);
}

/*! 
 * \return type of your current shape coordinate "dir". See CoordType enum 
 * \param[in] dir   0,1,2 int flag identifying coordinate
 */
CoordType BasicShape::getCoordinateType(int dir){
	return(m_typeCoord[dir]);
}
/*!
 * \return current type of shape instantiated
 */
ShapeType BasicShape::getShapeType(){
	return(m_shape);
};

/*! 
 * \return current type of shape instantiated. Const method overloading
 */
ShapeType BasicShape::getShapeType() const {
	return(m_shape);
};

/*! 
 * \return current scaling w.r.t the local sdr elemental shape
 */
darray3E BasicShape::getScaling(){
	return(m_scaling);
}

/*!
 * \return span of your elementary shape, in local coord system
 */
darray3E BasicShape::getLocalSpan(){
	return(m_span);
}

/*!
 * Given a geometry by MimmoObject class, return cell identifiers of those simplex inside the volume of
 * the BasicShape object. The methods implicitly use search algorithms based on the BvTree 
 * of the class MimmoObject.  
 * \param[in] geo target tessellation
 * \return list-by-ids of simplicies included in the volumetric patch
 */
livector1D BasicShape::includeGeometry(mimmo::MimmoObject * geo ){
	if(geo == NULL)	return livector1D(0);
	if(!(geo->isBvTreeSupported()))	return livector1D(0);
	
	//create BvTree and fill it w/ cell list
	if(!(geo->isBvTreeBuilt()))	geo->buildBvTree();
	//get recursively all the list element in the shape
	livector1D elements(geo->getNCells());
	int countElements = 0;
	searchBvTreeMatches(*(geo->getBvTree()), geo->getPatch(), 0, elements, countElements);
	elements.resize(countElements);
	
	return(elements);
};

/*! 
 * Given a geometry by MimmoObject class, return cell identifiers of those simplex outside the volume of
 * the BasicShape object. The methods implicitly use search algorithms based on the BvTree 
 * of the class MimmoObject.
 * \param[in] geo target tesselation
 * \return list-by-ids of simplicies outside the volumetric patch
 */
livector1D BasicShape::excludeGeometry(mimmo::MimmoObject * geo){
	
	if(geo == NULL)	return livector1D(0);
	livector1D internals = includeGeometry(geo);
	std::sort(internals.begin(), internals.end());
	
	livector1D originals = geo->getCellsIds();
	std::sort(originals.begin(), originals.end());
	
	livector1D result(originals.size() - internals.size());
	int counter = 0;
	auto internal_itr = internals.begin();
	auto original_itr = originals.begin();
	while (original_itr != originals.end()) {
		long origval = *original_itr;
		if (internal_itr == internals.end() || origval != *internal_itr) {
			result[counter] = origval;
			++counter;
		} else {
			++internal_itr;
		}
		++original_itr;
	}
	return result;
};

/*! 
 * Given a bitpit class bitpit::PatchKernel tessellation, return cell identifiers of those simplex inside the volume of
 * the BasicShape object. The method searches simplex one-by-one on the whole tesselation.
 * \param[in] tri target tessellation
 * \return list-by-ids of simplicies included in the volumetric patch
 */
livector1D BasicShape::includeGeometry(bitpit::PatchKernel * tri ){
	
	if(tri == NULL)	return livector1D(0);
	livector1D result(tri->getCellCount());
	long id;
	int counter = 0;
	for(auto & cell : tri->getCells()){
		id = cell.getId();
		if(isSimplexIncluded(tri,id)){
			result[counter] = id;
			++counter;
		}
	}	
	result.resize(counter);
	return(result);
 
};

/*!
 * Given a bitpit class bitpit::PatchKernel tessellation, return cell identifiers of those simplex outside the volume of
 * the BasicShape object.The method searches simplex one-by-one on the whole tesselation. 
 * \param[in] tri target tesselation
 * \return list-by-ids of simplicies outside the volumetric patch
 */
livector1D BasicShape::excludeGeometry(bitpit::PatchKernel * tri){
	
	if(tri == NULL)	return livector1D(0);
	livector1D result(tri->getCellCount());
	long id;
	int counter = 0;
	for(auto & cell : tri->getCells()){
		id = cell.getId();
		if(!isSimplexIncluded(tri,id)){
			result[counter] = id;
			++counter;
		}
	}	
	result.resize(counter);
	return(result);
	
};

/*! 
 * Given a list of vertices of a point cloud, return indices of those vertices included into 
 * the volume of the object. The method searches vertex one-by-one on the whole point cloud.  
 * \param[in] list list of cloud points
 * \return list-by-indices of vertices included in the volumetric patch
 */
livector1D BasicShape::includeCloudPoints(dvecarr3E & list){

	if(list.empty())	return livector1D(0);
	livector1D result(list.size());
	int counter = 0;
	long real = 0;
	for(auto & vert : list){
		if(isPointIncluded(vert)){
			result[counter] = real;
			++counter;
		}
		++real;	
	}	
	result.resize(counter);
	return(result);
};

/*! 
 * Given a list of vertices of a point cloud, return indices of those vertices outside 
 * the volume of BasicShape object 
 * \param[in] list list of cloud points
 * \return list-by-indices of vertices outside the volumetric patch
 */
livector1D BasicShape::excludeCloudPoints(dvecarr3E & list){
	if(list.empty())	return livector1D(0);
	livector1D result(list.size());
	long real = 0;
	int counter = 0;
	for(auto & vert : list){
		if(!isPointIncluded(vert)){
			result[counter] = real;
			++counter;
		}
		++real;	
	}	
	result.resize(counter);
	return(result);
	
};

/*! 
 * Given a bitpit class bitpit::PatchKernel point cloud, return identifiers of those points inside the volume of
 * the BasicShape object. The method force a one-by-one vertex search.   
 * \param[in] tri pointer to bitpit::PatchKernel object retaining the cloud point
 * \return list-by-ids of vertices included in the volumetric patch
 */
livector1D BasicShape::includeCloudPoints(bitpit::PatchKernel * tri){
	
	if(tri == NULL)		return livector1D(0);
	livector1D result(tri->getVertexCount());
	int counter = 0;
	long id;
	for(auto & vert : tri->getVertices()){
		id = vert.getId();
		if(isPointIncluded(tri, id)){
			result[counter] = id;
			++counter;
		}	
	}	
	result.resize(counter);
	return(result);
};

/*! 
 * Given a bitpit class bitpit::PatchKernel point cloud, return identifiers of those points outside the volume of
 * the BasicShape object.The method force a one-by-one vertex search.     
 * \param[in] tri pointer to bitpit::PatchKernel object retaining the cloud point
 * \return list-by-ids of vertices outside the volumetric patch
 */
livector1D BasicShape::excludeCloudPoints(bitpit::PatchKernel * tri){

	if(tri == NULL)		return livector1D(0);
	livector1D result(tri->getVertexCount());
	int counter = 0;
	long id;
	for(auto & vert : tri->getVertices()){
		id = vert.getId();
		if(!isPointIncluded(tri, id)){
			result[counter] = id;
		++counter;
		}
	}	
	result.resize(counter);
	return(result);
	
};

/*! 
 * Given a geometry by MimmoObject class, return vertex identifiers of those vertices inside the volume of
 * the BasicShape object. The methods implicitly use search algorithms based on the bitpit::KdTree 
 * of the class MimmoObject.  
 * \param[in] geo target geometry
 * \return list-by-ids of vertices included in the volumetric patch
 */
livector1D BasicShape::includeCloudPoints(mimmo::MimmoObject * geo ){
	if(geo == NULL)	return livector1D(0);
	
	livector1D elements(geo->getNVertex()); 
	//create BvTree and fill it w/ cell list
	if(!(geo->isKdTreeBuilt()))	geo->buildKdTree();
	
	getTempBBox();
	int countVertex = 0;
	//get recursively all the list element in the shape
	searchKdTreeMatches(*(geo->getKdTree()), 0, 0, elements, countVertex);
	elements.resize(countVertex);
	
	return(elements);
};

/*! Given a geometry by MimmoObject class, return vertex identifiers of those vertices outside the volume of
 * the BasicShape object. The methods implicitly use search algorithms based on the bitpit::KdTree 
 * of the class MimmoObject.
 * \param[in] geo target geometry
 * \return list-by-ids of vertices outside the volumetric patch
 */
livector1D BasicShape::excludeCloudPoints(mimmo::MimmoObject * geo){
	
	if(geo == NULL)	return livector1D(0);
	livector1D internals = includeCloudPoints(geo);
	std::sort(internals.begin(), internals.end());
	
	livector1D originals = geo->getCellsIds();
	std::sort(originals.begin(), originals.end());
	
	livector1D result(originals.size() - internals.size());
	int counter = 0;
	auto internal_itr = internals.begin();
	auto original_itr = originals.begin();
	while (original_itr != originals.end()) {
		long origval = *original_itr;
		if (internal_itr == internals.end() || origval != *internal_itr) {
 			result[counter] = origval;
			++counter;
		} else {
			++internal_itr;
		}
		++original_itr;
	}
	return result;
};

/*! 
 * \return true if all vertices of a given simplex are included in the volume of the shape
 * \param[in] simplexVert 3 vertices of the given Triangle
 */   
bool BasicShape::isSimplexIncluded(dvecarr3E & simplexVert){
  
  bool check = true;
  for(auto & val : simplexVert){
    check = check && isPointIncluded(val); 
  }
  return(check);
};

/*!
 * \return true if if all vertices of a given simplex are included in the volume of the shape
 * \param[in] tri pointer to a Class_SurfTri tesselation 
 * \param[in] indexT triangle index of tri.
 */ 
bool BasicShape::isSimplexIncluded(bitpit::PatchKernel * tri, long int indexT){

  Cell & cell = tri->getCell(indexT);
  int nVertices = cell.getVertexCount();
  bool check = true;
  for(int i=0; i<nVertices; ++i){ 
	//recover vertex index
	check = check && isPointIncluded(tri, cell.getVertex(i)); 
  }
  return(check);
};

/*!
 * \return true if the given point is included in the volume of the patch
 * \param[in] point given vertex
 */
bool BasicShape::isPointIncluded(darray3E point){
	
	bool check = true;
	darray3E temp = toLocalCoord(point);
	darray3E temp2 = localToBasic(temp);
	double tol = 1.0E-12;
	
	for(int i=0; i<3; ++i){  
		check = check && ((temp2[i] >= -1.0*tol) && (temp2[i]<=(1.0+tol)));
	}
	
	return(check);  
};

/*! 
 * \return true if the given point is included in the volume of the patch
 * \param[in] tri pointer to a bitpit::Patch tesselation / point cloud
 * \param[in] indexV id of a vertex belonging to tri;
 */
bool BasicShape::isPointIncluded(bitpit::PatchKernel * tri, long int indexV){
	
	darray3E coords = tri->getVertex(indexV).getCoords();
	return(isPointIncluded(coords));  
};

/*!
 * Check if current shape totally contains the given Axis Aligned Bounding Box
 * \param[in]	bMin	inferior extremal point of the AABB	
 * \param[in]	bMax	superior extremal point of the AABB	
 * \return true if completely contains AABB, false otherwise
 */
bool BasicShape::containShapeAABBox(darray3E bMin,darray3E bMax){
	bool check = true;
	dvecarr3E points(8);
	points[0] = bMin;
	points[1] = bMax;
	points[2] = bMax; points[2][2] = bMin[2];
	points[3] = bMin; points[3][2] = bMax[2];
	points[4] = bMin; points[4][0] = bMax[0];
	points[5] = bMax; points[5][0] = bMin[0];
	points[6] = bMin; points[6][1] = bMax[1];
	points[7] = bMax; points[7][1] = bMin[1];
	
	int counter=0;
	while (check && counter<8)	{
		check = check && isPointIncluded(points[counter]);	
		++counter;
	}	
	return check;
};


/*!
 * \return the nearest point belonging to an Axis Aligned Bounding Box, given a target vertex.
 * If the target is internal to or on surface of the AABB, return the target itself.
 * \param[in]	point	target vertex
 * \param[in]	bMin	inferior extremal point of the AABB	
 * \param[in]	bMax	superior extremal point of the AABB
 */
darray3E BasicShape::checkNearestPointToAABBox(darray3E point, darray3E bMin, darray3E bMax){
	
	darray3E result = {{0,0,0}};
	int counter=0;
	for (auto && val: point){
		result[counter] = std::fmin(bMax[counter], std::fmax(val,bMin[counter]));
		++counter;
	}
	return result;
};

/*!
 * Visit recursively KdTree relative to a cloud points and extract possible vertex candidates included in the current shape.
 * Identifiers of extracted matches are collected in result structure
 *\param[in] tree			KdTree of cloud points
 *\param[in] indexKdNode	KdNode index of tree, which start seaching from  
 *\param[in] level 			leaf level of indexKdNode in the tree
 *\param[in,out] result		list of KdNode labels, which are included in the shape.
 *\param[in,out] counter	last filled position of result vector. 
 * 
 */
void	BasicShape::searchKdTreeMatches(bitpit::KdTree<3,bitpit::Vertex,long> & tree, int indexKdNode, int level, livector1D & result, int &counter ){
	
	//check indexKdNode admissible.
	if(indexKdNode <0 || indexKdNode > tree.n_nodes)	return;
	
	//1st step get data
	bitpit::KdNode<bitpit::Vertex, long> & target = tree.nodes[indexKdNode];
	//check target point is in the shape
	if (bitpit::CGElem::intersectPointBox(target.object_->getCoords(), m_bbox[0], m_bbox[1], 3)){
		if(isPointIncluded(target.object_->getCoords())){
			result[counter] = target.label;
			++counter;
		}
		searchKdTreeMatches(tree, target.lchild_, level+1, result, counter);
		searchKdTreeMatches(tree, target.rchild_, level+1, result, counter);
	}else{
		
		uint32_t check = intersectShapePlane(level%3, target.object_->getCoords());
		bool right = (check > 0);
		bool left = (check%2 == 0);
		if (left && target.lchild_ >= 0){
			searchKdTreeMatches(tree, target.lchild_, level+1, result, counter);
		}
		if (right && target.rchild_ >= 0){
			searchKdTreeMatches(tree, target.rchild_, level+1, result, counter);
		}
	}	
	
	return;
};

/*!
 * Visit recursively BvTree relative to a PatchKernel structure and extract possible simplex candidates included in the current shape.
 * Identifiers of extracted matches are collected in result structure
 *\param[in] tree			BvTree of PatchKernel simplicies
 *\param[in] geo            pointer to tessellation the tree refers to. 
 *\param[in] indexBvNode	BvTree node index of tree, which start seaching from  
 *\param[in,out] result		list of PatchKernel ids , which are included in the shape.
 *\param[in,out] counter	last filled position of result vector.
 * 
 */
void	BasicShape::searchBvTreeMatches(bitpit::PatchSkdTree & tree,  bitpit::PatchKernel * geo, int indexBvNode, livector1D & result, int &counter ){
	
	//check indexBvNode admissible.
	if(indexBvNode <0 || indexBvNode > (int)tree.getNodeCount())	return;

    bitpit::SkdNode node = tree.getNode(indexBvNode);

	//1st step get data and control if it's leaf or if shape contain bounding box
	if(containShapeAABBox(node.getBoxMin(), node.getBoxMax())){
	    std::vector<long> cellids = node.getCells();
        for(std::vector<long>::iterator it = cellids.begin(); it != cellids.end(); ++it){
			result[counter] = *it;
			++counter;
		}	
		return;
	};
	
	if(node.isLeaf()){
        std::vector<long> cellids = node.getCells();
        for(std::vector<long>::iterator it = cellids.begin(); it != cellids.end(); ++it){
	        if(isSimplexIncluded(geo, *it)){
	            result[counter] = *it;
	            ++counter;
	        }
	    }
	    return;
	}

    for (int i = SkdNode::CHILD_BEGIN; i != SkdNode::CHILD_END; ++i) {
        SkdNode::ChildLocation childLocation = static_cast<SkdNode::ChildLocation>(i);
        if(node.getChildId(childLocation) != SkdNode::NULL_ID)   {
            if( intersectShapeAABBox(tree.getNode(node.getChildId(childLocation)).getBoxMin(), tree.getNode(node.getChildId(childLocation)).getBoxMax()) )
                searchBvTreeMatches(tree, geo, node.getChildId(childLocation), result, counter);
        }
    }

	return;
};


/*!
 * Given a 3D point, Check if the Axis Aligned Bounding box of your current shape intersects 
 * a given fundamental plane, x=a, y=b, z=c passing from such point. 
 * 
 * Return unsigned integer 0,1,2 with the following meaning:
 * 	- 0 : no intersection , shape on the left/bottom/before the plane
 * 	- 1 : no intersection , shape on the right/top/after the plane
 * 	- 2 : intersection occurs
 * 
 * \param[in] level current level of kdTree. 0-> xplane, 1-> yplane, 2->zplane, no other value are allowed.
 * \param[in] target 3D point
 * \return	 unsigned integer flag between 0 and 2
 */
uint32_t		BasicShape::intersectShapePlane(int level, darray3E target) {
	
	if(target[level] < m_bbox[0][level])	return 1; //shape is on the right
	if(target[level] > m_bbox[1][level])	return 0; //shape is on the left
	return 2;	// in the middle of something
	
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Cube IMPLEMENTATION 

/*! 
 * Basic Constructor
 */
Cube::Cube(){
	m_shape=ShapeType::CUBE;
	m_span = {{1.0, 1.0, 1.0}};
};

 /*! 
  * Custom Constructor. Set shape origin and its span,
  * ordered as width, height and depth. 
  * \param[in] origin point origin in global reference system
  * \param[in] span span in each shape local coordinate x,y,z;
  */
 Cube::Cube(const darray3E &origin, const darray3E & span): Cube(){
	      
	setOrigin(origin);
	setSpan(span[0], span[1], span[2]);
}; 

/*! Basic Destructor */
Cube::~Cube(){};


/*! 
 * Copy Constructor 
 * \param[in] other Cube object where copy from
 */
Cube::Cube(const Cube & other):BasicShape(other){};


/*! 
 * Transform point from local reference system of the shape,
 * to world reference system. 
 * \param[in] point target
 * \return transformed point
 */
darray3E	Cube::toWorldCoord(darray3E  point){
	
	darray3E work, work2;
	//unscale your local point
	for(int i =0; i<3; ++i){
		work[i] = point[i]*m_scaling[i];
	}
	
	//return to local xyz system
	// -> cube, doing nothing
	
	//unapply change to local sdr transformation
	linearalgebra::matmul(work, m_sdr, work2);
	
	//unapply origin translation
	work = work2 + m_origin;
	return(work);
};

/*!
 * Transform point from world coordinate system, to local reference system 
 * of the shape.
 * \param[in] point target
 * \return transformed point
 */
darray3E	Cube::toLocalCoord(darray3E  point){

	darray3E work, work2;

	//unapply origin translation
	work = point - m_origin;

	//apply change to local sdr transformation
	dmatrix33E transp = linearalgebra::transpose(m_sdr);
	linearalgebra::matmul(work, transp, work2);
	
	//get to proper local system
	// -> cube, doing nothing
	
	//scale your local point
	for(int i =0; i<3; ++i){
		work[i] = work2[i]/m_scaling[i];
	}
	
	return(work);
	
};

/*! 
 * \return local origin of your primitive shape
 */
darray3E	Cube::getLocalOrigin(){
	return(darray3E{-0.5,-0.5,-0.5});
};

/*!
 * Transform point from unitary cube reference system, to local reference system 
 * of the shape.
 * \param[in] point target
 * \return transformed point
 */
darray3E	Cube::basicToLocal(darray3E point){
	return(point + getLocalOrigin());
};

/*!
 * Transform point from local reference system of the shape,
 * to unitary cube reference system. 
 * \param[in] point target
 * \return transformed point
 */
darray3E	Cube::localToBasic(darray3E point){
	return(point - getLocalOrigin());
};

/*! 
 * Check if your new span values fit your current shape set up
 * and eventually return correct values.
 * \param[in] s0 first span dimension
 * \param[in] s1 second span dimension
 * \param[in] s2 third span dimension
 */
void 		Cube::checkSpan(double &s0, double &s1, double &s2){
			s0 = std::abs(s0);
			s1 = std::abs(s1);
			s2 = std::abs(s2);
};

/*! 
 * Check if your inf limits coords values fit your current shape set up
 * and eventually return correct values.
 * \param[in] o0 inf limits origin array
 * \param[in] dir coordinate index
 * \return true, if valid new value is set.
 */
bool 		Cube::checkInfLimits(double &o0, int &dir){
	BITPIT_UNUSED(o0);
	BITPIT_UNUSED(dir);
	//really doing nothing here.
	//whatever origin for whatever coordinate must return always 0 for cubic/cuboidal shape
	return(false);
};

/*! 
 * Set local span & scaling vectors of your object
 * \param[in] s0 first span dimension
 * \param[in] s1 second span dimension
 * \param[in] s2 third span dimension
 */
void 		Cube::setScaling( double &s0, double &s1, double &s2){
			m_span.fill(1.0);
			m_scaling[0] = s0;
			m_scaling[1] = s1;
			m_scaling[2] = s2;
};

/*!
 * Evaluate temporary Axis Aligned Bounding Box of the shape and 
 * store value in m_bbox protected member.
 */
void Cube::getTempBBox(){
	
	m_bbox.clear();
	m_bbox.resize(2);
	m_bbox[0].fill(1.e18);
	m_bbox[1].fill(-1.e18);
	
	dvecarr3E locals(8,darray3E{{0.0,0.0,0.0}});
	darray3E temp;
	locals[1].fill(1.0);
	locals[2][0]= locals[2][1]=1.0;
	locals[3][2]=1.0;
	locals[4][0]= 1.0;
	locals[5][1]= locals[5][2] = 1.0;
	locals[6][1] = 1.0;
	locals[7][0] = locals[7][2] = 1.0;
	
	for(auto &vv : locals){
		temp = toWorldCoord(basicToLocal(vv));
		for(int i=0; i<3; ++i)	{
			m_bbox[0][i] = std::fmin(m_bbox[0][i], temp[i]);
			m_bbox[1][i] = std::fmax(m_bbox[1][i], temp[i]);
		}	
	}
}

/*!
 * Check if current shape intersects or is totally contained into the given Axis Aligned Bounding Box
 * \param[in]	bMin	inferior extremal point of the AABB	
 * \param[in]	bMax	superior extremal point of the AABB	
 * \return true if intersects, false otherwise
 * 
 * TODO NEED TO BE OPTIMIZED!
 */
bool Cube::intersectShapeAABBox(darray3E bMin,darray3E bMax){
	
	dvecarr3E points(8, bMin);
	points[1][0] = points[2][0] = points[5][0]=points[6][0]= bMax[0];
	points[2][1] = points[3][1] = points[6][1]=points[7][1]= bMax[1];
	points[4][2] = points[5][2] = points[6][2]=points[7][2]= bMax[2];
	
	for(auto &val:points)	val += -1.0*m_origin;
	
	for(int j=0; j<3; ++j){
		double ref1 = -0.5*m_scaling[j];
		double ref2 =  0.5*m_scaling[j];
	
		double tmin, tmax, t;
		t= dotProduct(points[0],m_sdr[j]);
		tmax=tmin=t;
		
		for(int i=1; i<8; ++i){
			t= dotProduct(points[i],m_sdr[j]);
			if(t<tmin){
				tmin=t;
			}else if(t>tmax){
				tmax=t;
			}
		}
	
	//check height dimension if overlaps;
	if(tmin>ref2 || tmax<ref1)	return false;
	}
	
	return true;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Cylinder IMPLEMENTATION 

/*! 
 * Basic Constructor
 */
Cylinder::Cylinder(){
	m_shape=ShapeType::CYLINDER;
	m_span = {{1.0, 2*M_PI, 1.0}};
	setCoordinateType(CoordType::PERIODIC, 1);
};

/*! 
 * Custom Constructor. Set shape origin and its dimensions, 
 * ordered as basis radius, azimuthal/tangential coordinate and height. 
 * \param[in] origin point origin in global reference system
 * \param[in] span  characteristic dimension of your cylinder;
 */
Cylinder::Cylinder(const darray3E &origin, const darray3E & span): Cylinder(){
	
	setOrigin(origin);
	setSpan(span[0], span[1], span[2]);
}; 

/*! Basic Destructor */
Cylinder::~Cylinder(){};

/*! 
 * Copy Constructor 
 * \param[in] other Cylinder object where copy from
 */
Cylinder::Cylinder(const Cylinder & other):BasicShape(other){};

/*!
 * Transform point from local reference system of the shape,
 * to world reference system. 
 * \param[in] point target
 * \return transformed point
 */
darray3E	Cylinder::toWorldCoord(darray3E  point){
	
	darray3E work, work2;
	//unscale your local point
	for(int i =0; i<3; ++i){
		work[i] = point[i]*m_scaling[i];
	}
	
	//return to local xyz system
	work2[0] = work[0]*std::cos(work[1] + m_infLimits[1]); 
	work2[1] = work[0]*std::sin(work[1] + m_infLimits[1]); 
	work2[2] = work[2];
	
	//unapply change to local sdr transformation
	work.fill(0.0);
	linearalgebra::matmul(work2, m_sdr, work);
	
	//unapply origin translation
	work2 = work + m_origin;
	
	return(work2);
};

/*! 
 * Transform point from world coordinate system, to local reference system 
 * of the shape.
 * \param[in] point target
 * \return transformed point
 */
darray3E	Cylinder::toLocalCoord(darray3E  point){
	darray3E work, work2;
	
	//unapply origin translation
	work = point - m_origin;
	
	//apply change to local sdr transformation
	dmatrix33E transp = linearalgebra::transpose(m_sdr);
	linearalgebra::matmul(work, transp, work2);
	
	//get to proper local system
	if(work2[0] ==0.0 && work2[1] ==0.0){work[0] = 0.0; work[1] = 0.0;}
	else{
		work[0] = pow(work2[0]*work2[0] + work2[1]*work2[1],0.5);
		double pdum = std::atan2(work2[1],work2[0]);
		work[1] = pdum - 4.0*(getSign(pdum)-1.0)*std::atan(1.0); 
	}
	//get to the correct m_thetaOrigin mark
	double param = 8*std::atan(1.0);
	work[1] = work[1] - m_infLimits[1];
	if(work[1] < 0) 		work[1] = param + work[1];
	if(work[1] > param) 	work[1] = work[1] - param;
	
	work[2] = work2[2];
	
	//scale your local point
	for(int i =0; i<3; ++i){
		work2[i] = work[i]/m_scaling[i];
	}
	return(work2);
};

/*! 
 *\return local origin of your primitive shape
 */
darray3E	Cylinder::getLocalOrigin(){
	return(darray3E{0.0,0.0,-0.5});
};

/*!
 * Transform point from unitary cube reference system, to local reference system 
 * of the shape.
 * \param[in] point target
 * \return transformed point
 */
darray3E	Cylinder::basicToLocal(darray3E  point){
	point[1] = point[1]*m_span[1];
	point[2] = point[2]-0.5;
	return(point);
};

/*!
 * Transform point from local reference system of the shape,
 * to unitary cube reference system.
 * \param[in] point target
 * \return transformed point
 */
darray3E	Cylinder::localToBasic(darray3E  point){
	point[1] = point[1]/m_span[1];
	point[2] = point[2]+0.5;
	return(point);
};

/*!
 * Check if your new span values fit your current shape set up
 * and eventually return correct values.
 * \param[in] s0 first span dimension
 * \param[in] s1 second span dimension
 * \param[in] s2 third span dimension
 */
void 		Cylinder::checkSpan(double &s0, double &s1, double &s2){
	s0 = std::abs(s0);
	s1 = std::abs(s1);
	s2 = std::abs(s2);
	
	double thetalim = 8.0* std::atan(1.0);
	s1 = std::min(s1, thetalim);
	//check closedLoops;
	if(!(s1 < thetalim)){setCoordinateType(CoordType::PERIODIC,1);}
};

/*!
 * Check if your inf limits origin values fit your current shape set up
 * and eventually return correct values. 
 * \param[in] orig inf limits origin
 * \param[in] dir coordinate index to check
 * \return true, if valid new value is set.
 */
bool 	Cylinder::checkInfLimits(double &orig, int &dir){
	double thetalim = 8.0* std::atan(1.0);
	bool check = false;
	switch(dir){
		case 1: 
			orig = std::min(thetalim, std::max(0.0, orig));
			check = true;
			break;
		default:	// doing nothing
			break;
	}
	return(check);
};

/*! 
 * Set span and scaling vectors of your object
 * \param[in] s0 first span dimension
 * \param[in] s1 second span dimension
 * \param[in] s2 third span dimension
 */
void 		Cylinder::setScaling(double &s0, double &s1, double &s2){
		m_span.fill(1.0);
		m_scaling.fill(1.0);
		m_span[1] = s1;
		m_scaling[0] = s0;
		m_scaling[2] = s2;
};


/*!
 * Evaluate temporary Axis Aligned Bounding Box of the shape and 
 * store value in m_bbox protected member.
 */
void Cylinder::getTempBBox(){
	
	m_bbox.clear();
	m_bbox.resize(2);
	m_bbox[0].fill(1.e18);
	m_bbox[1].fill(-1.e18);
	
	dvecarr3E locals(20,darray3E{{0.0,0.0,0.0}});
	darray3E temp;
	int counter = 0;
	for(int i=0; i<2; ++i){
		for(int j=0; j<5; ++j){
			for(int k=0; k<2; ++k){
				locals[counter][0] = double(i)*1.0;
				locals[counter][1] = double(j)*0.25;
				locals[counter][2] = double(k)*1.0;
				++counter;
			}
		}
	}
	for(auto &vv : locals){
		temp = toWorldCoord(basicToLocal(vv));
		for(int i=0; i<3; ++i)	{
			m_bbox[0][i] = std::fmin(m_bbox[0][i], temp[i]);
			m_bbox[1][i] = std::fmax(m_bbox[1][i], temp[i]);
		}	
	}
}

/*!
 * Check if current shape intersects or is totally contained into the given Axis Aligned Bounding Box
 * \param[in]	bMin	inferior extremal point of the AABB	
 * \param[in]	bMax	superior extremal point of the AABB	
 * \return true if intersects, false otherwise
 * 
 * TODO NEED TO BE OPTIMIZED!
 */
bool Cylinder::intersectShapeAABBox(darray3E bMin,darray3E bMax){
	
	dvecarr3E points(8, bMin);
	points[1][0] = points[2][0] = points[5][0]=points[6][0]= bMax[0];
	points[2][1] = points[3][1] = points[6][1]=points[7][1]= bMax[1];
	points[4][2] = points[5][2] = points[6][2]=points[7][2]= bMax[2];
	
	for(auto &val:points)	val += -1.0*m_origin;
	
	double ref1, ref2;
	ref1 = -0.5*m_scaling[2];
	ref2 =  0.5*m_scaling[2];
	
	double tmin, tmax, t;
	t= dotProduct(points[0],m_sdr[2]);
	tmax=tmin=t;
	points[0] += -1.0*t*m_sdr[2]; //project point on plane
	
	for(int i=1; i<8; ++i){
		t= dotProduct(points[i],m_sdr[2]);
		points[i] += -1.0*t*m_sdr[2];
		
		if(t<tmin){
			tmin=t;
		}else if(t>tmax){
			tmax=t;
		}
		
	}
	
	//check height dimension if overlaps;
	if(tmin>ref2 || tmax<ref1)	return false;
	darray3E bMin2, bMax2;
	bMin2 = points[0]; bMax2=points[0];
	for(int i=1; i<8; ++i){
		for(int j=0; j<3; ++j){
			bMin2[j] = std::fmin(bMin2[j], points[i][j]);
			bMax2[j] = std::fmax(bMax2[j], points[i][j]);
		}
	}
	
	return(isPointIncluded(checkNearestPointToAABBox({{0.0,0.0,0.0}},bMin2,bMax2) + m_origin));
};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Sphere IMPLEMENTATION 

/*! 
 * Basic Constructor
 */
Sphere::Sphere(){
	m_shape=ShapeType::SPHERE;
	m_span = {{1.0, 2*M_PI, M_PI}};
	setCoordinateType(CoordType::PERIODIC, 1);
	setCoordinateType(CoordType::CLAMPED, 2);
};

/*! 
 * Custom Constructor. Set shape originand its dimensions, 
 * ordered as overall radius, azimuthal/tangential coordinate, polar coordinate. 
 * \param[in] origin point origin in global reference system
 * \param[in] span  characteristic dimension of your sphere/ portion of;
 */
Sphere::Sphere(const darray3E &origin, const darray3E & span): Sphere(){
	
	setOrigin(origin);
	setSpan(span[0], span[1], span[2]);
}; 

/*! Basic Destructor */
Sphere::~Sphere(){};

/*! 
 * Copy Constructor 
 * \param[in] other Sphere object where copy from
 */
Sphere::Sphere(const Sphere & other):BasicShape(other){};

/*! 
 * Transform point from local reference system of the shape,
 * to world reference system. 
 * \param[in] point target
 * \return transformed point
 */
darray3E	Sphere::toWorldCoord(darray3E  point){
	
	darray3E work, work2;
	//unscale your local point
	for(int i =0; i<3; ++i){
		work[i] = point[i]*m_scaling[i];
	}
	
	//return to local xyz system
	work2[0] = work[0]*std::cos(work[1] + m_infLimits[1])*std::sin(work[2] + m_infLimits[2]); 
	work2[1] = work[0]*std::sin(work[1] + m_infLimits[1])*std::sin(work[2] + m_infLimits[2]); 
	work2[2] = work[0]*std::cos(work[2] + m_infLimits[2]);
	
	//unapply change to local sdr transformation
	work.fill(0.0);
	linearalgebra::matmul(work2, m_sdr, work);
	
	//unapply origin translation
	work2 = work + m_origin;
	return(work2);
};

/*!
 * Transform point from world coordinate system, to local reference system 
 * of the shape.
 * \param[in] point target
 * \return transformed point
 */
darray3E	Sphere::toLocalCoord(darray3E  point){
	
	darray3E work, work2;
	//unapply origin translation
	work = point - m_origin;
	
	//apply change to local sdr transformation
	dmatrix33E transp = linearalgebra::transpose(m_sdr);
	linearalgebra::matmul(work, transp, work2);
	
	//get to proper local system
	work[0] = norm2(work2);
	
	if(work[0]>0.0){
		if(work2[0] ==0.0 && work2[1] ==0.0){
			work[1] = 0.0;
		}else{
			double pdum = std::atan2(work2[1],work2[0]);
			work[1] = pdum - 4.0*(getSign(pdum)-1.0)*std::atan(1.0); 
		}
		//get to the correct m_thetaOrigin mark
		double param = 8*std::atan(1.0);
		work[1] = work[1] - m_infLimits[1];
		if(work[1] < 0) 		work[1] = param + work[1];
		if(work[1] > param) 	work[1] = work[1] - param;
	
		work[2] = std::acos(work2[2]/work[0]);
		work[2] = work[2] - m_infLimits[2];
	}
	
	//scale your local point
	for(int i =0; i<3; ++i){
		work2[i] = work[i]/m_scaling[i];
	}
	return(work2);
};

/*! 
 * \return local origin of your primitive shape
 */
darray3E	Sphere::getLocalOrigin(){
	return(darray3E{0.0,0.0,0.0});
};

/*! 
 * Transform point from unitary cube reference system, to local reference system 
 * of the shape.
 * \param[in] point target
 * \return transformed point
 */
darray3E	Sphere::basicToLocal(darray3E  point){
	point[1] = point[1]*m_span[1];
	point[2] = point[2]*m_span[2];
	return(point);
};

/*! 
 * Transform point from local reference system of the shape,
 * to unitary cube reference system.
 * \param[in] point target
 * \return transformed point
 */
darray3E	Sphere::localToBasic(darray3E  point){
	point[1] = point[1]/m_span[1];
	point[2] = point[2]/m_span[2];
	return(point);
};

/*! 
 * Check if your new span values fit your current shape set up
 * and eventually return correct values.
 * \param[in] s0 first span dimension
 * \param[in] s1 second span dimension
 * \param[in] s2 third span dimension
 */
void 		Sphere::checkSpan(double &s0, double &s1, double &s2){
	s0 = std::abs(s0);
	s1 = std::abs(s1);
	s2 = std::abs(s2);
	
	double thetalim = 8.0* std::atan(1.0);
	s1 = std::min(s1, thetalim);
	
	double maxS2 = 0.5*thetalim - m_infLimits[2];
	s2 = std::min(s2, maxS2);
	
	//check closedLoops;
	if(!(s1 < thetalim)){setCoordinateType(CoordType::PERIODIC,1);}
	
};

/*! 
 * Check if inf limits origin values fit your current shape set up
 * and eventually return correct values.
 * \param[in] orig  inf limits origin
 * \param[in] dir   coordinate index to check
 * \return true, if valid new value is set
 * 
 */
bool 		Sphere::checkInfLimits(double &orig, int &dir){
	
	double thetalim = 8.0* std::atan(1.0);
	double tol = 1.e-12;
	bool check = false;
	switch(dir){
		case 1: 
			orig = std::min(thetalim, std::max(0.0, orig));
			check = true;
			break;
		case 2: 
			orig = std::min((0.5-tol)*thetalim, std::max(0.0, orig));
			check = true;
			break;
		default:	// doing nothing
			break;
	}
	return(check);
};

/*! 
 * Set scaling vector of your object 
 * \param[in] s0 first span dimension
 * \param[in] s1 second span dimension
 * \param[in] s2 third span dimension
 */
void 		Sphere::setScaling(double &s0, double &s1, double &s2){
	m_span.fill(1.0);
	m_scaling.fill(1.0);
	
	m_scaling[0] = s0;
	m_span[1] = s1;
	m_span[2] = s2;
};

/*!
 * Evaluate temporary Axis Aligned Bounding Box of the shape and 
 * store value in m_bbox protected member.
 */
void Sphere::getTempBBox(){
	
	m_bbox.clear();
	m_bbox.resize(2);
	m_bbox[0].fill(1.e18);
	m_bbox[1].fill(-1.e18);
	
	dvecarr3E locals(30,darray3E{{0.0,0.0,0.0}});
	darray3E temp;
	int counter = 0;
	for(int i=0; i<2; ++i){
		for(int j=0; j<5; ++j){
			for(int k=0; k<3; ++k){
				locals[counter][0] = double(i)*1.0;
				locals[counter][1] = double(j)*0.25;
				locals[counter][2] = double(k)*0.5;
				++counter;
			}
		}
	}
	
	for(auto &vv : locals){
		temp = toWorldCoord(basicToLocal(vv));
		for(int i=0; i<3; ++i)	{
			m_bbox[0][i] = std::fmin(m_bbox[0][i], temp[i]);
			m_bbox[1][i] = std::fmax(m_bbox[1][i], temp[i]);
		}	
	}
}

/*!
 * Check if current shape intersects or is totally contained into the given Axis Aligned Bounding Box
 * \param[in]	bMin	inferior extremal point of the AABB	
 * \param[in]	bMax	superior extremal point of the AABB	
 * \return true if intersects, false otherwise
 */
bool Sphere::intersectShapeAABBox(darray3E bMin,darray3E bMax){
	return isPointIncluded(checkNearestPointToAABBox(m_origin, bMin, bMax));
};

}
