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
#include "BasicShapes.hpp"
#include "customOperators.hpp"
#include <chrono>

using namespace bitpit;
using namespace mimmo;
using namespace std::chrono;

/*
 *	\date			03/01/2015
 *  \authors		Edoardo Lombardi
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *
 *	\brief Abstract Interface class for Elementary Shape Representation
 *
 *	Interface class for Volumetric Core Element, suitable for interaction with Data Structure stored in a MimmoObject class.
 *  Object orientation in 3D space can be externally manipulated with dedicated transformation blocks. Class 
 *  internally implement transformation to/from local sdr to/from world sdr, that can be used in derived objects from it.
 *	Class works with three reference systems:
 * 	1) Global Absolute SDR: is the external World reference system
 *  2) Local Relative SDR: is the local reference system, not affected by Rigid Transformations as RotoTranslations or Scalings 
 *  3) basic SDR: local system remapping to unitary cube, not accounting of the shape type.  
 *   
 */
 
/*! Basic Constructor */
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

/*! Set origin of your shape. The origin is meant as the baricenter of your shape in absolute r.s.
 * \param[in] origin new origin point
 */
void BasicShape::setOrigin(darray3E origin){
	m_origin = origin;
}

/*! Set span of your shape, according to its local reference system 
 * \param[in] s0 first coordinate span
 * \param[in] s1 second coordinate span
 * \param[in] s2 third coordinate span
 */
void BasicShape::setSpan(double s0, double s1, double s2){
	checkSpan(s0,s1,s2);
	setScaling(s0,s1,s2);
}

/*! Set inferior limits your shape, in its local reference system.
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

/*! Set new axis orientation of the local reference system
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

/*! Set new axis orientation of the local reference system
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

/*! Set new axes orientation of the local reference system
 * \param[in] axes new direction of local axes.
 */
void BasicShape::setRefSystem(dmatrix33E axes){
	for (int i=0; i<3; i++){
		m_sdr[i] = axes[i]/norm2(axes[i]);
	}
}

/*! Set type to treat your shape coordinates. 
 * \param[in] type CoordType enum.
 * \param[in] dir  0,1,2 int flag identifying coordinate
 */
void BasicShape::setCoordinateType(CoordType type, int dir){
	m_typeCoord[dir] = type;
}

/*! Return current origin of your shape
 * \return origin
 */
darray3E BasicShape::getOrigin(){
	return(m_origin);
}

/*! Return current span of your shape
 * \return span
 */
darray3E BasicShape::getSpan(){
	darray3E result = getLocalSpan();
	darray3E scale = getScaling();
	for(int i; i<3; ++i){
		result[i] = scale[i]*result[i];
	}
	return(result);
}

/*! Return current inferior limits of your shape, in local coord reference system
 * \return coords origin
 */
darray3E BasicShape::getInfLimits(){
	return(m_infLimits);
}

/*! Return actual axis of global relative sdr
 * \return relative sdr
 */
dmatrix33E BasicShape::getRefSystem(){
	return(m_sdr);
}

/*! Return type of your current shape coordinate "dir". See CoordType enum 
 * \param[in] dir   0,1,2 int flag identifying coordinate
 */
CoordType BasicShape::getCoordinateType(int dir){
	return(m_typeCoord[dir]);
}
/*! Get current type of shape instantiated
 * \return ShapeType enum
 */
ShapeType BasicShape::getShapeType(){
	return(m_shape);
};

/*! Get current type of shape instantiated. Const method overloading
 * \return const ShapeType enum
 */
const ShapeType BasicShape::getShapeType() const {
	return(m_shape);
};

/*! Get current scaling w.r.t the local sdr elemental shape*/
darray3E BasicShape::getScaling(){
	return(m_scaling);
}

/*! Return span of your elementary shape, in local coord system
 * \return span
 */
darray3E BasicShape::getLocalSpan(){
	return(m_span);
}

/*! Given a geometry by MimmoObject class, return cell identifiers of those simplex inside the volume of
 * the BasicShape object. The methods implicitly use search algorithms based on the bvTree 
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
	livector1D elements;
	searchBvTreeMatches(*(geo->getBvTree()), geo->getPatch(), 0, elements);
	
	return(elements);
};

/*! Given a geometry by MimmoObject class, return cell identifiers of those simplex outside the volume of
 * the BasicShape object. The methods implicitly use search algorithms based on the bvTree 
 * of the class MimmoObject.
 * \param[in] geo target tesselation
 * \return list-by-ids of simplicies outside the volumetric patch
 */
livector1D BasicShape::excludeGeometry(mimmo::MimmoObject * geo){
	
	if(geo == NULL)	return livector1D(0);
	if(!(geo->isBvTreeSupported()))	return livector1D(0);

	livector1D interiors, result;
	interiors = includeGeometry(geo);
	result = geo->getCells().getIds();
	
	std::unordered_set<long> map;
	map.insert(result.begin(), result.end());
	result.clear();
	
	for(auto && val: interiors)	map.erase(val);
	
	result.insert(result.end(), map.begin(), map.end());
	return result;
};



/*! Given a bitpit class bitpit::Patch tessellation, return cell identifiers of those simplex inside the volume of
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

/*! Given a bitpit class bitpit::Patch tessellation, return cell identifiers of those simplex outside the volume of
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

/*! Given a list of vertices of a point cloud, return indices of those vertices included into 
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

/*! Given a list of vertices of a point cloud, return indices of those vertices outside 
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

/*! Given a bitpit class bitpit::Patch point cloud, return identifiers of those points inside the volume of
 * the BasicShape object. The method force a build up of kdTree on point cloud, and use it to refine search.   
 * \param[in] tri pointer to bitpit::PatchKernel object retaining the cloud point
 * \return list-by-ids of vertices included in the volumetric patch
 */
livector1D BasicShape::includeCloudPoints(bitpit::PatchKernel * tri){
	
	if(tri == NULL)		return livector1D(0);
	livector1D elements; 
	//create kdTree and fill it w/ vertex list
	bitpit::KdTree<3,bitpit::Vertex,long> tree;
	int counter=0;
	long label;
	
	for(auto & val : tri->getVertices()){
		label = val.getId();
		tree.insert(&val, label);
		++counter;
	}
	getShapeTempPoints();
	//get recursively all the list element in the shape
	searchKdTreeMatches(tree, 0, 0, elements);
	return(elements);
};

/*! Given a bitpit class bitpit::Patch point cloud, return identifiers of those points outside the volume of
 * the BasicShape object.The method force a build up of kdTree on point cloud, and use it to refine search.  
 * \param[in] tri pointer to bitpit::PatchKernel object retaining the cloud point
 * \return list-by-ids of vertices outside the volumetric patch
 */
livector1D BasicShape::excludeCloudPoints(bitpit::PatchKernel * tri){

	if(tri == NULL)		return livector1D(0);
	
	livector1D interiors, result;
	interiors = includeCloudPoints(tri);
	
	std::set<long> map;
	long size = tri->getVertexCount();
	for (long i=0; i<size; ++i)	map.insert(i);
	for(auto && val: interiors)	map.erase(val);
	
	result.insert(result.end(), map.begin(), map.end());
	return result;
};

/*! Given a geometry by MimmoObject class, return vertex identifiers of those vertices inside the volume of
 * the BasicShape object. The methods implicitly use search algorithms based on the kdTree 
 * of the class MimmoObject.  
 * \param[in] geo target geometry
 * \return list-by-ids of vertices included in the volumetric patch
 */
livector1D BasicShape::includeCloudPoints(mimmo::MimmoObject * geo ){
	if(geo == NULL)	return livector1D(0);
	
	livector1D elements; 
	//create BvTree and fill it w/ cell list
	if(!(geo->isKdTreeBuilt()))	geo->buildKdTree();
	
	getShapeTempPoints();
	//get recursively all the list element in the shape
	searchKdTreeMatches(*(geo->getKdTree()), 0, 0, elements);
	return(elements);
};

/*! Given a geometry by MimmoObject class, return vertex identifiers of those vertices outside the volume of
 * the BasicShape object. The methods implicitly use search algorithms based on the kdTree 
 * of the class MimmoObject.
 * \param[in] geo target geometry
 * \return list-by-ids of vertices outside the volumetric patch
 */
livector1D BasicShape::excludeCloudPoints(mimmo::MimmoObject * geo){
	
	if(geo == NULL)	return livector1D(0);
	livector1D interiors, result;
	interiors = includeCloudPoints(geo);
	result = geo->getVertices().getIds();
	
	std::unordered_set<long> map;
	map.insert(result.begin(), result.end());
	result.clear();
	
	for(auto && val: interiors)	map.erase(val);
	
	result.insert(result.end(), map.begin(), map.end());
	return result;
};



/*! Return True if all vertices of a given simplex are included in the volume of the shape
 * \param[in] simplexVert 3 vertices of the given Triangle
 * \return boolean
 */   
bool BasicShape::isSimplexIncluded(dvecarr3E & simplexVert){
  
  bool check = true;
  for(auto && val : simplexVert){
   check = check && isPointIncluded(val); 
  }
  return(check);
};

/*! Return True if if all vertices of a given simplex are included in the volume of the shape
 * \param[in] tri pointer to a Class_SurfTri tesselation 
 * \param[in] indexT triangle index of tri.
 * \return boolean
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

/*! Return True if the given point is included in the volume of the patch
 * \param[in] point given vertex
 * \return boolean
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

/*! Return True if the given point is included in the volume of the patch
 * \param[in] tri pointer to a bitpit::Patch tesselation / point cloud
 * \param[in] indexV id of a vertex belonging to tri;
 * \return boolean
 */
bool BasicShape::isPointIncluded(bitpit::PatchKernel * tri, long int indexV){
	
	bool check = true;
	darray3E coords = tri->getVertex(indexV).getCoords();
	return(isPointIncluded(coords));  
};

/*!
 * Check if current shape intersects or is totally contained into the given Axis Aligned Bounding Box
 * \param[in]	bMin	inferior extremal point of the AABB	
 * \param[in]	bMax	superior extremal point of the AABB	
 * \return true if intersects, false otherwise
 */
bool BasicShape::intersectShapeAABBox(darray3E bMin,darray3E bMax){
	return isPointIncluded(checkNearestPointToAABBox(m_origin, bMin, bMax));
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
 * Return the nearest point belonging to an Axis Aligned Bounding Box, given a target vertex.
 * If the target is internal to or on surface of the AABB, return the target itself.
 * \param[in]	point	target vertex
 * \param[in]	bMin	inferior extremal point of the AABB	
 * \param[in]	bMax	superior extremal point of the AABB
 * \return		the nearest point of AABB wrt to target vertex
 */
darray3E BasicShape::checkNearestPointToAABBox(darray3E point, darray3E bMin, darray3E bMax){
	
	darray3E result;
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
 * 
 */
void	BasicShape::searchKdTreeMatches(bitpit::KdTree<3,bitpit::Vertex,long> & tree, int indexKdNode, int level, livector1D & result ){
	
	//check indexKdNode admissible.
	if(indexKdNode <0 || indexKdNode > tree.n_nodes)	return;
	
	//1st step get data
	bitpit::KdNode<bitpit::Vertex, long> & target = tree.nodes[indexKdNode];
	//check target point is in the shape
	if(isPointIncluded(target.object_->getCoords())){
		result.push_back(target.label);
	}
	
	switch(intersectShapePlane(getKdPlane(level, target.object_->getCoords()))){
		case 0: 
			searchKdTreeMatches(tree, target.lchild_, level+1, result);
			break;
		case 1: 		
			searchKdTreeMatches(tree, target.rchild_, level+1, result);
			break;
		default:
			searchKdTreeMatches(tree, target.lchild_, level+1, result);
			searchKdTreeMatches(tree, target.rchild_, level+1, result);
			break;
	}		
	return;
};

/*!
 * Visit recursively BvTree relative to a PatchKernel structure and extract possible simplex candidates included in the current shape.
 * Identifiers of extracted matches are collected in result structure
 *\param[in] tree			BvTree of PatchKernel simplicies
 *\param[in] geo            pointer to tessellation the tree refers to. 
 *\param[in] indexKdNode	BvTree node index of tree, which start seaching from  
 *\param[in,out] result		list of PatchKernel ids , which are included in the shape.
 * 
 */
void	BasicShape::searchBvTreeMatches(mimmo::BvTree & tree,  bitpit::PatchKernel * geo, int indexBvNode, livector1D & result ){
	
	//check indexBvNode admissible.
	if(indexBvNode <0 || indexBvNode > tree.m_nnodes)	return;
	   
	//1st step get data and control if it's leaf or if shape contain bounding box
	if(containShapeAABBox(tree.m_nodes[indexBvNode].m_minPoint, tree.m_nodes[indexBvNode].m_maxPoint)){
	
		for(int i = tree.m_nodes[indexBvNode].m_element[0]; i<tree.m_nodes[indexBvNode].m_element[1]; ++i){
			result.push_back(tree.m_elements[i].m_label);		
		}	
		return;
	};
	
	if(tree.m_nodes[indexBvNode].m_leaf){
		for(int i = tree.m_nodes[indexBvNode].m_element[0]; i<tree.m_nodes[indexBvNode].m_element[1]; ++i){
			if(isSimplexIncluded(geo, tree.m_elements[i].m_label))
								result.push_back(tree.m_elements[i].m_label);		
		}
		return;
	}
	
	if(tree.m_nodes[indexBvNode].m_lchild >= 0)	{
		if( intersectShapeAABBox(tree.m_nodes[tree.m_nodes[indexBvNode].m_lchild].m_minPoint, tree.m_nodes[tree.m_nodes[indexBvNode].m_lchild].m_maxPoint) )
			searchBvTreeMatches(tree, geo, tree.m_nodes[indexBvNode].m_lchild, result);
	}
	
	if(tree.m_nodes[indexBvNode].m_rchild >= 0)	{
		if( intersectShapeAABBox(tree.m_nodes[tree.m_nodes[indexBvNode].m_rchild].m_minPoint, tree.m_nodes[tree.m_nodes[indexBvNode].m_rchild].m_maxPoint) )
			searchBvTreeMatches(tree, geo, tree.m_nodes[indexBvNode].m_rchild, result);
	}
	
	return;
};


/*!
 * Return x,y,z planes passing by a given point
 * \param[in]	level int counter. level%3 select x-plane(0), yplane(1), zplane(2)
 * \param[in]	point point belonging to request plane.
 * \return		plane  
 */
darray4E 	BasicShape::getKdPlane(int level, darray3E point){
	darray4E plane({{0.0,0.0,0.0,0.0}});
	
	int type = level%3;
	plane[type] = 1.0;
	plane[3] = -1.0*point[type];
	return plane;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Cube IMPLEMENTATION 
/*
 *	\date			03/01/2015
 *	\authors		Edoardo Lombardi
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *
 *	\brief Elementary Shape Representation of a Cube
 *
 *	Volumetric Core Element, shaped as a cube, directly derived from BasicShape class.   
 */
   
/*! Basic Constructor */
Cube::Cube(){
	m_shape=ShapeType::CUBE;
	m_span = {{1.0, 1.0, 1.0}};
};

 /*! Custom Constructor. Set shape origin and its span,
  * ordered as width, height and depth. 
   * \param[in] origin point origin in global reference system
   * \param[in] span span in each shape local coordinate x,y,z;
   */
 Cube::Cube(darray3E &origin, darray3E & span): Cube(){
	      
	setOrigin(origin);
	setSpan(span[0], span[1], span[2]);
}; 

/*! Basic Destructor */
Cube::~Cube(){};


/*! Copy Constructor 
 * \param[in] other Cube object where copy from
 */
Cube::Cube(const Cube & other){
	*this = other;
};

/*! Copy Operator 
 * \param[in] other Cube object where copy from
 */
Cube & Cube::operator=(const Cube & other){
	
	m_shape = other.m_shape;
	m_origin = other.m_origin;
	m_span = other.m_span;
	m_infLimits = other.m_infLimits;
	m_sdr = other.m_sdr;
	m_typeCoord = other.m_typeCoord;
	m_scaling = other.m_scaling;
	return(*this);
};


/*! Transform point from local reference system of the shape,
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
/*! Transform point from world coordinate system, to local reference system 
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

/*! Return local origin of your primitive shape*/
darray3E	Cube::getLocalOrigin(){
	return(darray3E{-0.5,-0.5,-0.5});
};

/*! Transform point from unitary cube reference system, to local reference system 
 * of the shape.
 * \param[in] point target
 * \return transformed point
 */
darray3E	Cube::basicToLocal(darray3E point){
	return(point + getLocalOrigin());
};

/*! Transform point from local reference system of the shape,
 * to unitary cube reference system. 
 * \param[in] point target
 * \return transformed point
 */
darray3E	Cube::localToBasic(darray3E point){
	return(point - getLocalOrigin());
};

/*! Check if your new span values fit your current shape set up
 * and eventually return correct values.
 */
void 		Cube::checkSpan(double &s0, double &s1, double &s2){
			s0 = std::abs(s0);
			s1 = std::abs(s1);
			s2 = std::abs(s2);
};

/*! Check if your coords origin values fit your current shape set up
 * and eventually return correct values. Return true, if valid new value is set.
 */
bool 		Cube::checkInfLimits(double &o0, int &dir){
			//really doing nothing here.
			//whatever origin for whatever coordinate must return always 0 for cubic/cuboidal shape
			return(false);
};

/*! set local span & scaling vectors of your object */
void 		Cube::setScaling( double &s0, double &s1, double &s2){
			m_span.fill(1.0);
			m_scaling[0] = s0;
			m_scaling[1] = s1;
			m_scaling[2] = s2;
};

/*!
 * Check if you current shape intersects a given plane, in its implicit form implicit plane a*x + b*y + c*z + d = 0.
 * Return unsigned integer 0,1,2 with the following meaning:
 * 	- 0 : no intersection , shape on the left/bottom/before the plane
 * 	- 1 : no intersection , shape on the right/top/after the plane
 * 	- 2 : intersection occurs
 * \param[in] plane coefficient of plane
 * \return	 unsigned integer flag between 0 and 2
 */
uint32_t		Cube::intersectShapePlane(darray4E plane) {
	//get eight vertex of the cube.
	bool check = false;
	double distToCheck,  dist = plane[3];

	int counter = 1;
	for(int i=0; i<3; ++i)	dist +=  plane[i]*m_tempPoints[0][i];
	if(std::abs(dist) < 1.e-12)	return 2;

	while(!check && counter<8){
		distToCheck = plane[3];
		for(int i=0; i<3; ++i)	distToCheck +=  plane[i]*m_tempPoints[counter][i];
		check = (dist*distToCheck <= 0.0 || std::abs(distToCheck) < 1.E-12 );
		++counter;
	}
	
	if(check) 	return 2;
	else{
		dist = plane[3];
		for(int i=0; i<3; ++i)	dist += m_origin[i]*plane[i];
		if(dist > 0)	return 1;
		return 0;
	}
}

void Cube::getShapeTempPoints(){
	
	m_tempPoints.clear();
	m_tempPoints.resize(8);
	
	dvecarr3E locals(8,darray3E{{0.0,0.0,0.0}});
	locals[1].fill(1.0);
	locals[2][0]= locals[2][1]=1.0;
	locals[3][2]=1.0;
	locals[4][0]= 1.0;
	locals[5][1]= locals[5][2] = 1.0;
	locals[6][1] = 1.0;
	locals[7][0] = locals[7][2] = 1.0;
	
	int counter = 0;
	for(auto &vv : locals){
		m_tempPoints[counter] = toWorldCoord(basicToLocal(vv));
		++counter;
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Cylinder IMPLEMENTATION 
/*
 *	\date			03/01/2015
 *	\authors		Edoardo Lombardi
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *
 *	\brief Elementary Shape Representation of a Cylinder or portion of it
 *
 *	Volumetric Core Element, shaped as a cylinder, directly derived from BasicShape class.   
 */

/*! Basic Constructor */
Cylinder::Cylinder(){
	m_shape=ShapeType::CYLINDER;
	m_span = {{1.0, 2*M_PI, 1.0}};
	setCoordinateType(CoordType::PERIODIC, 1);
};

/*! Custom Constructor. Set shape origin and its dimensions, 
 * ordered as basis radius, azimuthal/tangential coordinate and height. 
 * \param[in] origin point origin in global reference system
 * \param[in] span  characteristic dimension of your cylinder;
 */
Cylinder::Cylinder(darray3E &origin, darray3E & span): Cylinder(){
	
	setOrigin(origin);
	setSpan(span[0], span[1], span[2]);
}; 

/*! Basic Destructor */
Cylinder::~Cylinder(){};

/*! Copy Constructor 
 * \param[in] other Cylinder object where copy from
 */
Cylinder::Cylinder(const Cylinder & other){
	*this = other;
};

/*! Copy Operator 
 * \param[in] other Cylinder object where copy from
 */
Cylinder & Cylinder::operator=(const Cylinder & other){
	
	m_shape = other.m_shape;
	m_origin = other.m_origin;
	m_span = other.m_span;
	m_infLimits = other.m_infLimits;
	m_sdr = other.m_sdr;
	m_typeCoord = other.m_typeCoord;
	m_scaling = other.m_scaling;
	return(*this);
};


/*! Transform point from local reference system of the shape,
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
/*! Transform point from world coordinate system, to local reference system 
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

/*! Return local origin of your primitive shape*/
darray3E	Cylinder::getLocalOrigin(){
	return(darray3E{0.0,0.0,-0.5});
};

/*! Transform point from unitary cube reference system, to local reference system 
 * of the shape.
 * \param[in] point target
 * \return transformed point
 */
darray3E	Cylinder::basicToLocal(darray3E  point){
	point[1] = point[1]*m_span[1];
	point[2] = point[2]-0.5;
	return(point);
};

/*! Transform point from local reference system of the shape,
 * to unitary cube reference system.
 * \param[in] point target
 * \return transformed point
 */
darray3E	Cylinder::localToBasic(darray3E  point){
	point[1] = point[1]/m_span[1];
	point[2] = point[2]+0.5;
	return(point);
};

/*! Check if your new span values fit your current shape set up
 * and eventually return correct values.
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

/*! Check if your coords origin values fit your current shape set up
 * and eventually return correct values. Return true, if valid new value is set.
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

/*! set span and scaling vectors of your object */
void 		Cylinder::setScaling(double &s0, double &s1, double &s2){
		m_span.fill(1.0);
		m_scaling.fill(1.0);
		m_span[1] = s1;
		m_scaling[0] = s0;
		m_scaling[2] = s2;
};

/*!
 * Check if you current shape intersects a given plane, in its implicit form implicit plane a*x + b*y + c*z + d = 0.
 * Return unsigned integer 0,1,2 with the following meaning:
 * 	- 0 : no intersection , shape on the left/bottom/before the plane
 * 	- 1 : no intersection , shape on the right/top/after the plane
 * 	- 2 : intersection occurs
 * \param[in] plane coefficient of plane
 * \return	 unsigned integer flag between 0 and 2
 */
uint32_t	Cylinder::intersectShapePlane(darray4E plane) {
	//get center of circular basis
	int pcand = 0;
	std::array<double,2> dist;
	dist.fill(plane[3]);
	
	for(int i=0; i<3; ++i){
		dist[0] += plane[i]*m_tempPoints[0][i];
		dist[1] += plane[i]*m_tempPoints[1][i];
	}	
	
	if(std::abs(dist[0]) < 1.e-12 || std::abs(dist[1]) < 1.e-12 || dist[0]*dist[1] < 0 )	return 2;
	
	if(std::abs(dist[1]) < std::abs(dist[0]))	pcand = 1;
	darray3E pointOnPlane;
	for(int i=0; i<3; ++i)	pointOnPlane[i] = m_tempPoints[pcand][i] - plane[i] * dist[pcand]; 

	if(isPointIncluded(pointOnPlane))	return 2;
	else{
		dist[0] = plane[3];
		for(int i=0; i<3; ++i)	dist[0] +=m_origin[i]*plane[i];
		if(dist[0]> 0)	return 1;
		return 0;
	}
}

void Cylinder::getShapeTempPoints(){
	
	m_tempPoints.clear();
	m_tempPoints.resize(2);
	
	dvecarr3E locals(2,darray3E{{0.0,0.0,0.0}});
	locals[1][2]=1.0;
	
	int counter = 0;
	for(auto &vv : locals){
		m_tempPoints[counter] = toWorldCoord(basicToLocal(vv));
		++counter;
	}
	
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Sphere IMPLEMENTATION 
/*
 *	\date			03/01/2015
 *	\authors		Edoardo Lombardi
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *
 *	\brief Elementary Shape Representation of a sphere or portion of it
 *
 *	Volumetric Core Element, shaped as a sphere, directly derived from BasicShape class.   
 */

/*! Basic Constructor */
Sphere::Sphere(){
	m_shape=ShapeType::SPHERE;
	m_span = {{1.0, 2*M_PI, M_PI}};
	setCoordinateType(CoordType::PERIODIC, 1);
	setCoordinateType(CoordType::CLAMPED, 2);
};

/*! Custom Constructor. Set shape originand its dimensions, 
 * ordered as overall radius, azimuthal/tangential coordinate, polar coordinate. 
 * \param[in] origin point origin in global reference system
 * \param[in] span  characteristic dimension of your sphere/ portion of;
 */
Sphere::Sphere(darray3E &origin, darray3E & span): Sphere(){
	
	setOrigin(origin);
	setSpan(span[0], span[1], span[2]);
}; 

/*! Basic Destructor */
Sphere::~Sphere(){};

/*! Copy Constructor 
 * \param[in] other Sphere object where copy from
 */
Sphere::Sphere(const Sphere & other){
	*this = other;
};

/*! Copy Operator 
 * \param[in] other Sphere object where copy from
 */
Sphere & Sphere::operator=(const Sphere & other){
	
	m_shape = other.m_shape;
	m_origin = other.m_origin;
	m_span = other.m_span;
	m_infLimits = other.m_infLimits;
	m_sdr = other.m_sdr;
	m_typeCoord = other.m_typeCoord;
	m_scaling = other.m_scaling;
	return(*this);
};

/*! Transform point from local reference system of the shape,
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
/*! Transform point from world coordinate system, to local reference system 
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

/*! Return local origin of your primitive shape*/
darray3E	Sphere::getLocalOrigin(){
	return(darray3E{0.0,0.0,0.0});
};

/*! Transform point from unitary cube reference system, to local reference system 
 * of the shape.
 * \param[in] point target
 * \return transformed point
 */
darray3E	Sphere::basicToLocal(darray3E  point){
	point[1] = point[1]*m_span[1];
	point[2] = point[2]*m_span[2];
	return(point);
};

/*! Transform point from local reference system of the shape,
 * to unitary cube reference system.
 * \param[in] point target
 * \return transformed point
 */
darray3E	Sphere::localToBasic(darray3E  point){
	point[1] = point[1]/m_span[1];
	point[2] = point[2]/m_span[2];
	return(point);
};

/*! Check if your new span values fit your current shape set up
 * and eventually return correct values.
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

/*! Check if your coords origin values fit your current shape set up
 * and eventually return correct values.
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

/*! set scaling vector of your object */
void 		Sphere::setScaling(double &s0, double &s1, double &s2){
	m_span.fill(1.0);
	m_scaling.fill(1.0);
	
	m_scaling[0] = s0;
	m_span[1] = s1;
	m_span[2] = s2;
};

/*!
 * Check if you current shape intersects a given plane, in its implicit form implicit plane a*x + b*y + c*z + d = 0.
 * Return unsigned integer 0,1,2 with the following meaning:
 * 	- 0 : no intersection , shape on the left/bottom/before the plane
 * 	- 1 : no intersection , shape on the right/top/after the plane
 * 	- 2 : intersection occurs
 * \param[in] plane coefficient of plane
 * \return	 unsigned integer flag between 0 and 2
 */
uint32_t	Sphere::intersectShapePlane(darray4E plane) {
	//get center of sphere
	
	darray3E p = toWorldCoord(basicToLocal({{0.0,0.0,0.0}}));
	double dist = plane[3];
	for(int i=0; i<3; ++i)	dist += plane[i]*p[i];
	if(dist <= getSpan()[0])	return 2;
	else if(dist> 0)	return 1;
	else return 0;
}

void Sphere::getShapeTempPoints(){
	
	m_tempPoints.clear();
}

