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

using namespace bitpit;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//BasicShape IMPLEMENTATION 
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
 *  internally implement transformation to/from local sdr to/from world sdr, that can be used in derived objects from it   
 */
 
/*! Basic Constructor */
BasicShape::BasicShape(){

	for(int i=0; i<3; ++i){
		m_limitSpan[i].fill(0.0);
		m_sdr[i].fill(0.0);
	}
	m_sdr[0][0] = m_sdr[1][1] = m_sdr[2][2] = 1.0;
	m_closedLoops.resize(3, false);
	m_shape = ShapeType::CUBE;
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
	m_limitSpan[0][1] = m_limitSpan[0][0] + s0;
	m_limitSpan[1][1] = m_limitSpan[1][0] + s1;
	m_limitSpan[2][1] = m_limitSpan[2][0] + s2;
	setScaling();
}

/*! Set coordinate's origin of your shape, according to its local reference system  
 * \param[in] o0 first coordinate origin
 * \param[in] o1 second coordinate origin
 * \param[in] o2 third coordinate origin
 */
void BasicShape::setInfLimits(double o0, double o1, double o2){
	
	darray3E span = getSpan();
	
	checkInfLimits(o0,o1,o2);
	m_limitSpan[0][0] =o0;
	m_limitSpan[1][0] =o1;
	m_limitSpan[2][0] =o2;
	
	setSpan(span[0],span[1],span[2]);
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
	double check = dotProduct(axis0,axis1) + dotProduct(axis1,axis2) + dotProduct(axis0,axis2);
	if(check > tol) return;
	m_sdr[0] = axis0;
	m_sdr[1] = axis1;
	m_sdr[2] = axis2;
}

/*! Set new axis orientation of the local reference system
 * \param[in] axis target axis
 * \param[in] point external point, not belonging to the axis, to recover the three axis system
 */
void BasicShape::setRefSystem(darray3E axis, darray3E point){

	axis = axis/norm2(axis);
	
	double pj = dotProduct(point, axis);
	darray3E axisN = point - pj*axis;
	
	axisN = axisN/norm2(axisN);

	darray3E axisT = crossProduct(axis,axisN);

	setRefSystem(axis, axisN, axisT);
}

/*! Set booleans to treat your shape coordinates as periodic (true) or regular (false)
 * \param[in] per0 first coordinate
 * \param[in] per1 second coordinate
 * \param[in] per2 third coordinate
 */
void BasicShape::setClosedLoops(bool per0, bool per1, bool per2){
	m_closedLoops[0] = per0;
	m_closedLoops[1] = per1;
	m_closedLoops[2] = per2;
}

/*! Return current origin of your shape
 * \param[out] result origin
 */
darray3E BasicShape::getOrigin(){
	return(m_origin);
}

/*! Return current span of your shape
 * \param[out] result span
 */
darray3E BasicShape::getSpan(){
	darray3E result;
	for(int i=0; i<3; ++i) {
		result[i] = m_limitSpan[i][1] - m_limitSpan[i][0];
	}
	return(result);
}

/*! Return current coordinates' origin of your shape
 * \param[out] result coords origin
 */
darray3E BasicShape::getInfLimits(){
	darray3E result;
	for(int i=0; i<3; ++i) {
		result[i] = m_limitSpan[i][0];
	}
	return(result);
}

/*! Return actual axis of global relative sdr
 * \param[out] result relative sdr
 */
dmatrix33E BasicShape::getRefSystem(){
	return(m_sdr);
}

/*! Return if your current shape coordinates are set as periodic or not
 * \param[in] result boolean vector
 */
bvector1D BasicShape::areClosedLoops(){
	return(m_closedLoops);
}
/*! Get current type of shape instantiated
 * \param[out] result BasicShape::ShapeType enum
 */
BasicShape::ShapeType BasicShape::getShapeType(){
	return(m_shape);
};

/*! Get current type of shape instantiated. Const method version
 * \param[out] result const BasicShape::ShapeType enum
 */
const BasicShape::ShapeType BasicShape::getShapeType() const {
	return(m_shape);
};

/*! Get current scaling w.r.t the primitive unitary shape*/
darray3E BasicShape::getScaling(){
	return(m_scaling);
}

/*! Given a bitpit class bitpit::Patch tessellation, return cell identifiers of those simplex inside the volume of
 * the BasicShape object
 * \param[in] tri target tessellation
 * \param[out] result list-by-ids of simplicies included in the volumetric patch
 */
ivector1D BasicShape::includeGeometry(bitpit::Patch * tri ){
  
  int nCells = tri->getCellCount();	
  ivector1D result(nCells); 
  int counter=0;
  
  for(auto &cell : tri->cells()){
   
    if(isSimplexIncluded(tri, cell.get_id())){
      result[counter] = cell.get_id();
      ++counter;
    }
  }
  result.resize(counter);
  return(result);
};

/*! Given a bitpit class bitpit::Patch tessellation, return cell identifiers of those simplex outside the volume of
 * the BasicShape object 
 * \param[in] tri target tesselation
 * \param[out] result list-by-ids of simplicies outside the volumetric patch
 */
ivector1D BasicShape::excludeGeometry(bitpit::Patch * tri){
  
	int nCells = tri->getCellCount();	
	ivector1D result(nCells); 
	int counter=0;
	
	for(auto &cell : tri->cells()){
		
		if(!isSimplexIncluded(tri, cell.get_id())){
			result[counter] = cell.get_id();
			++counter;
		}
	}
	result.resize(counter);
	return(result);
};

/*! Given a list of vertices of a point cloud, return indices of those vertices included into 
 * the volume of BaseSelPatch object 
 * \param[in] list list of cloud points
 * \param[out] result list-by-indices of vertices included in the volumetric patch
 */
ivector1D BasicShape::includeCloudPoints(dvecarr3E & list){

  int size = list.size();
  ivector1D result(size); 
  int counter=0;
  
  for(int i=0; i<size; ++i){
   
    if(isPointIncluded(list[i])){
      result[counter] = i;
      ++counter;
    }
  }
  result.resize(counter);
  return(result);
};

/*! Given a list of vertices of a point cloud, return indices of those vertices outside 
 * the volume of BasicShape object 
 * \param[in] list list of cloud points
 * \param[out] result list-by-indices of vertices outside the volumetric patch
 */
ivector1D BasicShape::excludeCloudPoints(dvecarr3E & list){
  
  int size = list.size();
  ivector1D result(size); 
  int counter=0;
  
  for(int i=0; i<size; ++i){
   
    if(!isPointIncluded(list[i])){
      result[counter] = i;
      ++counter;
    }
  }
  result.resize(counter);
  return(result);
  
};

/*! Given a bitpit class bitpit::Patch point cloud, return identifiers of those points inside the volume of
 * the BasicShape object  
 * \param[in] list list of cloud points
 * \param[out] result list-by-indices of vertices included in the volumetric patch
 */
ivector1D BasicShape::includeCloudPoints(bitpit::Patch * tri){
	
	int nVert = tri->getVertexCount();	
	ivector1D result(nVert); 
	int counter=0;
	
	for(auto &vertex : tri->vertices()){
		
		if(isPointIncluded(tri, vertex.get_id())){
			result[counter] = vertex.get_id();
			++counter;
		}
	}
	result.resize(counter);
	return(result);
};

/*! Given a bitpit class bitpit::Patch point cloud, return identifiers of those points outside the volume of
 * the BasicShape object  
 * \param[in] list list of cloud points
 * \param[out] result list-by-indices of vertices outside the volumetric patch
 */
ivector1D BasicShape::excludeCloudPoints(bitpit::Patch * tri){
	
	int nVert = tri->getVertexCount();	
	ivector1D result(nVert); 
	int counter=0;
	
	for(auto &vertex : tri->vertices()){
		
		if(!isPointIncluded(tri, vertex.get_id())){
			result[counter] = vertex.get_id();
			++counter;
		}
	}
	result.resize(counter);
	return(result);	
};

/*! Return True if at least one vertex of a given triangle is included in the volume of the shape
 * \param[in] simplexVert 3 vertices of the given Triangle
 * \param[out] result boolean
 */   
bool BasicShape::isSimplexIncluded(dvecarr3E & simplexVert){
  
  bool check = false;
  for(int i=0; i<simplexVert.size(); ++i){
   check = check || isPointIncluded(simplexVert[i]); 
  }
  return(check);
};

/*! Return True if at least one vertex of a given triangle is included in the volume of the shape
 * \param[in] tri pointer to a Class_SurfTri tesselation 
 * \param[in] indexT triangle index of tri.
 * \param[out] result boolean
 */ 
bool BasicShape::isSimplexIncluded(bitpit::Patch * tri, int indexT){

  Cell cell = tri->getCell(indexT);
  long * conn = cell.getConnect();
  int nVertices = cell.getVertexCount();
  bool check = false;
  for(int i=0; i<nVertices; ++i){ 
	//recover vertex index
	check = check || isPointIncluded(tri, conn[i]); 
  }
  return(check);
};

/*! Return True if the given point is included in the volume of the patch
 * \param[in] point given vertex
 * \param[out] result boolean
 */
bool BasicShape::isPointIncluded(darray3E point){
	
	bool check = true;
	darray3E temp = toLocalCoord(point);
	darray3E temp2 = localToBasic(temp);
	
	for(int i=0; i<3; ++i){  
		check = check && ((temp2[i] >= 0.0) && (temp2[i]<=1.0));
	}
	
	return(check);  
};

/*! Return True if the given point is included in the volume of the patch
 * \param[in] tri pointer to a bitpit::Patch tesselation / point cloud
 * \param[in] indexV id of a vertex belonging to tri;
 * \param[out] result boolean
 */
bool BasicShape::isPointIncluded(bitpit::Patch * tri, int indexV){
	
	bool check = true;
	darray3E coords = tri->getVertex(indexV).getCoords();
	return(isPointIncluded(coords));  
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
};

 /*! Custom Constructor. Set shape origin, inferior/superior limits of the cube,
  * ordered as width, height, span. 
   * \param[in] origin_ point origin in global reference system
   * \param[in] limits inf/sup limit for each shape coordinate;
   */
Cube::Cube(darray3E &origin, dmatrix32E & limits): Cube(){
	      
	setOrigin(origin);
	darray3E span;
	darray3E infLim;
	for(int i=0; i<3; ++i){
		span[i] = limits[i][1] - limits[i][0];
		infLim[i] = limits[i][0];
	}
	setInfLimits(infLim[0], infLim[1], infLim[2]);
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
	m_limitSpan = other.m_limitSpan;
	m_sdr = other.m_sdr;
	m_closedLoops = other.m_closedLoops;
	m_scaling = other.m_scaling;
	return(*this);
};

/*! Return current span of your shape
 * \param[out] result span
 */
darray3E Cube::getLocalSpan(){
	darray3E result{1.0,1.0,1.0};
	return(result);
}

/*! Return current coordinates' origin of your shape
 * \param[out] result coords origin
 */
darray3E Cube::getLocalInfLimits(){
	darray3E result{-0.5,-0.5,-0.5};
	return(result);
}

/*! Transform point from local reference system of the shape,
 * to world reference system. 
 * \param[in] point target
 * \param[out] result transformed point
 */
darray3E	Cube::toWorldCoord(darray3E & point){
	
	darray3E work, work2;
	//unscale your local point
	for(int i =0; i<3; ++i){
		work[i] = point[i]*m_scaling[i];
	}
	
	//return to local xyz system
	// -> cube, doing nothing
	
	//unapply change to local sdr transformation
	dmatrix33E transp = linearalgebra::transpose(m_sdr);
	linearalgebra::matmul(work, transp, work2);
	
	//unapply origin translation
	work = work2 + m_origin;
	return(work);
};
/*! Transform point from world coordinate system, to local reference system 
 * of the shape.
 * \param[in] point target
 * \param[out] result transformed point
 */
darray3E	Cube::toLocalCoord(darray3E & point){
	darray3E work, work2;

	//unapply origin translation
	work = point - m_origin;

	//apply change to local sdr transformation
	linearalgebra::matmul(work, m_sdr, work2);

	//get to proper local system
	// -> cube, doing nothing
	
	//scale your local point
	for(int i =0; i<3; ++i){
		work[i] = work2[i]/m_scaling[i];
	}
	
	return(work);
	
};

/*! Transform point from unitary cube reference system, to local reference system 
 * of the shape.
 * \param[in] point target
 * \param[out] result transformed point
 */
darray3E	Cube::basicToLocal(darray3E & point){
	darray3E basicOr{-0.5,-0.5,-0.5};
	return(point + basicOr);
};

/*! Transform point from local reference system of the shape,
 * to unitary cube reference system.
 * \param[in] point target
 * \param[out] result transformed point
 */
darray3E	Cube::localToBasic(darray3E & point){
	darray3E basicOr{-0.5,-0.5,-0.5};
	return(point - basicOr);
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
 * and eventually return correct values.
 */
void 		Cube::checkInfLimits(double &o0, double &o1, double &o2){
			//really doing nothing here.
};

/*! set scaling vector of your object */
void 		Cube::setScaling(){
			for(int i=0; i<3; ++i){
				m_scaling[i] = m_limitSpan[i][1] - m_limitSpan[i][0];
			}
};


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
};

/*! Custom Constructor. Set shape origin, inferior/superior limits of the cylinder
 * ordered as basis radius, azimuthal/tangential coordinate, height. 
 * \param[in] origin_ point origin in global reference system
 * \param[in] limits inf/sup limit for each shape coordinate;
 */
Cylinder::Cylinder(darray3E &origin, dmatrix32E & limits): Cylinder(){
	
	setOrigin(origin);
	darray3E span;
	darray3E infLim;
	for(int i=0; i<3; ++i){
		span[i] = limits[i][1] - limits[i][0];
		infLim[i] = limits[i][0];
	}
	setInfLimits(infLim[0], infLim[1], infLim[2]);
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
	m_limitSpan = other.m_limitSpan;
	m_sdr = other.m_sdr;
	m_closedLoops = other.m_closedLoops;
	m_scaling = other.m_scaling;
	return(*this);
};

/*! Return current span of your shape
 * \param[out] result span
 */
darray3E Cylinder::getLocalSpan(){
	darray3E result{1.0,1.0,1.0};
	result[1] = m_limitSpan[1][1] - m_limitSpan[1][0];
	return(result);
}

/*! Return current coordinates' origin of your shape
 * \param[out] result coords origin
 */
darray3E Cylinder::getLocalInfLimits(){
	darray3E result{0.0,0.0,-0.5};
	return(result);
}



/*! Transform point from local reference system of the shape,
 * to world reference system. 
 * \param[in] point target
 * \param[out] result transformed point
 */
darray3E	Cylinder::toWorldCoord(darray3E & point){
	
	darray3E work, work2, work3;
	//unscale your local point
	for(int i =0; i<3; ++i){
		work[i] = point[i]*m_scaling[i];
	}
	
	//return to local xyz system
	work2[0] = work[0]*std::cos(work[1] + m_limitSpan[1][0]); 
	work2[1] = work[0]*std::sin(work[1] + m_limitSpan[1][0]); 
	work2[2] = work[2];
	
	//unapply change to local sdr transformation
	dmatrix33E transp = linearalgebra::transpose(m_sdr);
	linearalgebra::matmul(work2, transp, work3);
	
	//unapply origin translation
	work = work3 + m_origin;
	return(work);
};
/*! Transform point from world coordinate system, to local reference system 
 * of the shape.
 * \param[in] point target
 * \param[out] result transformed point
 */
darray3E	Cylinder::toLocalCoord(darray3E & point){
	darray3E work, work2,work3;
	
	//unapply origin translation
	work = point - m_origin;
	
	//apply change to local sdr transformation
	linearalgebra::matmul(work, m_sdr, work2);
	
	//get to proper local system
	if(work2[0] ==0.0 && work2[1] ==0.0){work3[0] = 0.0; work3[1] = 0.0;}
	else{
		work3[0] = pow(work2[0]*work2[0] + work2[1]*work2[1],0.5);
		double pdum = std::atan2(work2[1],work2[0]);
		work3[1] = pdum - 4.0*(getSign(pdum)-1.0)*std::atan(1.0); 
	}
		//get to the correct m_thetaOrigin mark
		double param = 8*std::atan(1.0);
		work3[1] = work3[1] - m_limitSpan[1][0];
		if(work3[1] < 0) 		work3[1] = param + work3[1];
		if(work3[1] > param) 	work3[1] = work3[1] - param;
		work3[2] = work2[2];
	
	//scale your local point
	for(int i =0; i<3; ++i){
		work[i] = work3[i]/m_scaling[i];
	}
	return(work);
};

/*! Transform point from unitary cube reference system, to local reference system 
 * of the shape.
 * \param[in] point target
 * \param[out] result transformed point
 */
darray3E	Cylinder::basicToLocal(darray3E & point){
	darray3E basicOr{0.0,0.0,-0.5};
	point[1] = point[1]*(m_limitSpan[1][1] - m_limitSpan[1][0]);
	return(point + basicOr);
};

/*! Transform point from local reference system of the shape,
 * to unitary cube reference system.
 * \param[in] point target
 * \param[out] result transformed point
 */
darray3E	Cylinder::localToBasic(darray3E & point){
	darray3E basicOr{0.0,0.0,-0.5};
	point[1] = point[1]/(m_limitSpan[1][1] - m_limitSpan[1][0]);
	return(point - basicOr);
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
	if(s1 <= 0.0) {s1 = thetalim;}
	//check closedLoops;
	m_closedLoops[1] = !(s1 < thetalim);
	
};

/*! Check if your coords origin values fit your current shape set up
 * and eventually return correct values.
 */
void 		Cylinder::checkInfLimits(double &o0, double &o1, double &o2){
	
	o0 = std::max(o0,0.0);
	double thetalim = 8.0* std::atan(1.0);
	o1 = std::min(thetalim, std::max(0.0, o1));
};

/*! set scaling vector of your object */
void 		Cylinder::setScaling(){
	m_scaling[0] = m_limitSpan[0][1] - m_limitSpan[0][0];
	m_scaling[1] = 1.0;
	m_scaling[2] = m_limitSpan[2][1] - m_limitSpan[2][0];
};

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
};

/*! Custom Constructor. Set shape origin, inferior/superior limits of the sphere
 * ordered as overall radius, azimuthal/tangential coordinate, polar coordinate. 
 * \param[in] origin_ point origin in global reference system
 * \param[in] limits inf/sup limit for each shape coordinate;
 */
Sphere::Sphere(darray3E &origin, dmatrix32E & limits): Sphere(){
	
	setOrigin(origin);
	darray3E span;
	darray3E infLim;
	for(int i=0; i<3; ++i){
		span[i] = limits[i][1] - limits[i][0];
		infLim[i] = limits[i][0];
	}
	setInfLimits(infLim[0], infLim[1], infLim[2]);
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
	m_limitSpan = other.m_limitSpan;
	m_sdr = other.m_sdr;
	m_closedLoops = other.m_closedLoops;
	m_scaling = other.m_scaling;
	return(*this);
};

/*! Return current span of your shape
 * \param[out] result span
 */
darray3E Sphere::getLocalSpan(){
	darray3E result{1.0,1.0,1.0};
	result[1] = m_limitSpan[1][1] - m_limitSpan[1][0];
	result[2] = m_limitSpan[2][1] - m_limitSpan[2][0];
	return(result);
}

/*! Return current coordinates' origin of your shape
 * \param[out] result coords origin
 */
darray3E Sphere::getLocalInfLimits(){
	darray3E result{0.0,0.0,0.0};
	return(result);
}

/*! Transform point from local reference system of the shape,
 * to world reference system. 
 * \param[in] point target
 * \param[out] result transformed point
 */
darray3E	Sphere::toWorldCoord(darray3E & point){
	
	darray3E work, work2, work3;
	//unscale your local point
	for(int i =0; i<3; ++i){
		work[i] = point[i]*m_scaling[i];
	}
	
	//return to local xyz system
	work2[0] = work[0]*std::cos(work[1] + m_limitSpan[1][0])*std::sin(work[2] + m_limitSpan[2][0]); 
	work2[1] = work[0]*std::sin(work[1] + m_limitSpan[1][0])*std::sin(work[2] + m_limitSpan[2][0]); 
	work2[2] = work[0]*std::cos(work[2] + m_limitSpan[2][0]);
	
	//unapply change to local sdr transformation
	dmatrix33E transp = linearalgebra::transpose(m_sdr);
	linearalgebra::matmul(work2, transp, work3);
	
	//unapply origin translation
	work = work3 + m_origin;
	return(work);
};
/*! Transform point from world coordinate system, to local reference system 
 * of the shape.
 * \param[in] point target
 * \param[out] result transformed point
 */
darray3E	Sphere::toLocalCoord(darray3E & point){
	
	darray3E work, work2, work3{0,0,0};
	//unapply origin translation
	work = point - m_origin;
	
	//apply change to local sdr transformation
	linearalgebra::matmul(work, m_sdr, work2);
	
	//get to proper local system
	work3[0] = norm2(work2);
	
	if(work3[0]>0.0){
		if(work2[0] ==0.0 && work2[1] ==0.0){
			work3[1] = 0.0;
		}else{
			double pdum = std::atan2(work2[1],work2[0]);
			work3[1] = pdum - 4.0*(getSign(pdum)-1.0)*std::atan(1.0); 
		}
		//get to the correct m_thetaOrigin mark
		double param = 8*std::atan(1.0);
		work3[1] = work3[1] - m_limitSpan[1][0];
		if(work3[1] < 0) 		work3[1] = param + work3[1];
		if(work3[1] > param) 	work3[1] = work3[1] - param;
	
		work3[2] = std::acos(work2[2]/work3[0]);
		work3[2] = work3[2] - m_limitSpan[2][0];
	}
	
	//scale your local point
	for(int i =0; i<3; ++i){
		work[i] = work3[i]/m_scaling[i];
	}
	return(work);
};

/*! Transform point from unitary cube reference system, to local reference system 
 * of the shape.
 * \param[in] point target
 * \param[out] result transformed point
 */
darray3E	Sphere::basicToLocal(darray3E & point){
	darray3E  basicOr{0.0,0.0,0.0};
	
	point[1] = point[1]*(m_limitSpan[1][1] - m_limitSpan[1][0]);
	point[2] = point[2]*(m_limitSpan[2][1] - m_limitSpan[2][0]);
	return(point + basicOr);
};

/*! Transform point from local reference system of the shape,
 * to unitary cube reference system.
 * \param[in] point target
 * \param[out] result transformed point
 */
darray3E	Sphere::localToBasic(darray3E & point){
	darray3E  basicOr{0.0,0.0,0.0};

	point[1] = point[1]/(m_limitSpan[1][1] - m_limitSpan[1][0]);
	point[2] = point[2]/(m_limitSpan[2][1] - m_limitSpan[2][0]);
	return(point - basicOr);
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
	if(s1<=0.0){s1 = thetalim;}
	
	double maxS2 = 0.5*thetalim - m_limitSpan[2][0];
	
	s2 = std::min(s2, maxS2);
	if(s2<=0.0){s2 = maxS2;};
	
	//check closedLoops;
	m_closedLoops[1] = !(s1 < thetalim);
	
};

/*! Check if your coords origin values fit your current shape set up
 * and eventually return correct values.
 */
void 		Sphere::checkInfLimits(double &o0, double &o1, double &o2){
	
	double tol = 1.e-12;
	o0 = std::max(o0,0.0);
	double thetalim = 8.0* std::atan(1.0);
	double polarlim = (0.5 - tol)*thetalim;
	o1 = std::min(thetalim, std::max(0.0, o1));
	o2 = std::min(polarlim, std::max(0.0, o2));
};

/*! set scaling vector of your object */
void 		Sphere::setScaling(){
	m_scaling[0] = m_limitSpan[0][1] - m_limitSpan[0][0];
	m_scaling[1] = 1.0;
	m_scaling[2] = 1.0;
};
