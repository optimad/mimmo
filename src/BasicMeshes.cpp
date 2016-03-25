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

#include "BasicMeshes.hpp"
#include "Operators.hpp"
#include "customOperators.hpp"
#include "MiMMO_VTKInterfaces.hpp"

using namespace std;

//*****************************************************************************************************************************
// USTRUCTMESH IMPLEMENTATION 
/*
*	\date			31/12/2015
*	\authors		Edoardo Lombardi
*	\authors		Rocco Arpa
*	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
*
*	\brief Class for 3D uniform structured mesh --> based on ucartmesh scheme of Bitpit 
*
*	Interface class for uniform structured grids, suitable for derivation of meshes on pure cartesian, cylindrical or spherical
*	coordinates system. The class retains internal members of Class BasicShape who determine the shape of your current grid.
*  The mesh works and its nodal structures are defined in the its local reference frame, that is, the local reference frame 
*	retained by its core BasicShape object. 
*/


/*! Basic Constructor. Default Shape CUBE */
UStructMesh::UStructMesh(){
	m_shape1 = NULL;
	m_setmesh = false;
	m_setorigin = false;
	m_setspan = false;
	m_origin = {{0.0,0.0,0.0}};
	m_span = {{0.0,0.0,0.0}};
	m_nx = 0; m_ny=0; m_nz=0;
	m_dx = 0; m_dy=0; m_dz=0;
	setShape(BasicShape::ShapeType::CUBE);
};

///*! Custom Constructor. Set your mesh, according to the following input parameters
// * \param[in] origin 3D point origin of your mesh
// * \param[in] span span for each coordinate defining your mesh
// * \param[in] type   shape of your mesh, based on BasicShape::ShapeType enum.(option available are: CUBE(default), CYLINDER, SPHERE)
// * \param[in] dimensions number of mesh points for each coordinate.
// */
//UStructMesh::UStructMesh(darray3E & origin, darray3E & span, BasicShape::ShapeType type, ivector1D & dimensions):
//			UStructMesh()
//{
//	setMesh(origin,span,type,dimensions);
//};
//
///*! Custom Constructor. Set your mesh, according to the following input parameters
// * \param[in] origin 3D point baricenter of your mesh
// * \param[in] span span for each coordinate defining your mesh
// * \param[in] type   shape of your mesh, based on BasicShape::ShapeType enum.(option available are: CUBE(default), CYLINDER, SPHERE)
// * \param[in] spacing fixed spacing for each coordinate
// */
//UStructMesh::UStructMesh(darray3E & origin, darray3E &span, BasicShape::ShapeType type, dvector1D &spacing):
//		UStructMesh()
//{
//	setMesh(origin,span,type,spacing);
//};
//
///*! Custom Constructor. Set your mesh, according to the following input parameters
// * \param[in] shape pointer to an external allocated BasicShape object
// * \param[in] dimensions number of mesh points for each coordinate.
// */
//UStructMesh::UStructMesh(BasicShape * shape, ivector1D & dimensions):UStructMesh(){
//	setMesh(shape, dimensions);
//};
//
///*! Custom Constructor. Set your mesh, according to the following input parameters
// * \param[in] shape pointer to an external allocated BasicShape object
// * \param[in] spacing fixed spacing for each coordinate
// */
//UStructMesh::UStructMesh(BasicShape * shape, dvector1D & spacing):UStructMesh(){
//	setMesh(shape, spacing);
//};

/*! Basic destructor */
UStructMesh::~UStructMesh(){
	m_xnode.clear(); m_ynode.clear(); m_znode.clear();
	m_xedge.clear(); m_yedge.clear(); m_zedge.clear();
	m_shape1 = NULL;
}

/*! Copy Constructor.
 * \param[in] other UStructMesh object where copy from
 */
UStructMesh::UStructMesh(const UStructMesh & other){
	*this = other;
};

/*! Copy Operator.
 * \param[in] other UStructMesh object where copy from
 */
UStructMesh & UStructMesh::operator=(const UStructMesh & other){
	
	// Number of cells
	m_nx = other.m_nx;
	m_ny = other.m_ny;
	m_nz = other.m_nz;
	
	// Cell Spacing
	m_dx = other.m_dx;
	m_dy = other.m_dy;
	m_dz = other.m_dz;

	// Copy cell edges and cell centers ----------------------------------------- //
	m_xnode = other.m_xnode;
	m_ynode = other.m_ynode;
	m_znode = other.m_znode;
	m_xedge = other.m_xedge;
	m_yedge = other.m_yedge;
	m_zedge = other.m_zedge;

	m_setorigin = other.m_setorigin;
	m_setspan = other.m_setspan;
	m_setmesh = other.m_setmesh;
	m_shape1 = other.m_shape1;

	if(other.m_shape2){
		switch(other.getShape()->getShapeType()){
			case BasicShape::ShapeType::CUBE :
				m_shape2 = std::unique_ptr<BasicShape>(new Cube(*(dynamic_cast<const Cube *> (&other))));
				break;
			case BasicShape::ShapeType::CYLINDER :
				m_shape2 = std::unique_ptr<BasicShape>(new Cylinder(*(dynamic_cast<const Cylinder*> (&other))));
				break;
			case BasicShape::ShapeType::SPHERE :
				m_shape2 = std::unique_ptr<BasicShape>(new Sphere(*(dynamic_cast<const Sphere*> (&other))));
				break;
			default:
				//never been reached
				break;
		}
	}
	return(*this); 
};

/*! Return a const pointer to the inner BasicShape object the current mesh is built on. Const method 
 * \param[out] result BasicShape of the mesh
 */
const BasicShape * UStructMesh::getShape() const {
	
	BasicShape * result;
	if(m_shape2){
		result = m_shape2.get();
	}else{
		result = m_shape1;
	}
	return(result);
}

/*! Return current origin of BasicShape core of the mesh*/
darray3E UStructMesh::getOrigin(){
	if (getShape() == NULL) return(darray3E({{-9999, -9999, -9999}}));
	return(getShape()->getOrigin());
}

/*! Return current span of BasicShape core of the mesh*/
darray3E UStructMesh::getSpan(){
	return(getShape()->getSpan());
}

/*! Return current lower limits of coordinates in BasicShape core of the mesh*/
darray3E UStructMesh::getInfLimits(){
	return(getShape()->getInfLimits());
}

/*! Return current local Reference System af axes*/
dmatrix33E UStructMesh::getRefSystem(){
	return(getShape()->getRefSystem());
}

/*! Return actual scaling to primitive shape used in BasicShape core of the mesh*/
darray3E UStructMesh::getScaling(){
	return(getShape()->getScaling());
}

/*! Return local span of the primitive shape associated to BasicShape core of the mesh*/
darray3E UStructMesh::getLocalSpan(){
	return(getShape()->getLocalSpan());
}

/*! Return type of shape associated to mesh core. See BasicShape::ShapeType enum */
BasicShape::ShapeType UStructMesh::getShapeType(){
	return(getShape()->getShapeType());
}

/*! Return coordinate type of component 0 a BasicShape mesh core.
 * See BasicShape::CoordType enum.
 */
BasicShape::CoordType UStructMesh::getCoordTypex(){
	return(getShape()->getCoordinateType(0));
}

/*! Return coordinate type of a component of a BasicShape mesh core.
 * See BasicShape::CoordType enum.
 * \param[in] i index of component.
 */
BasicShape::CoordType UStructMesh::getCoordType(int i){
	return(getShape()->getCoordinateType(i));
}

/*! Return coordinate type of component 1 a BasicShape mesh core.
 * See BasicShape::CoordType enum.
 */
BasicShape::CoordType UStructMesh::getCoordTypey(){
	return(getShape()->getCoordinateType(1));
}

/*! Return coordinate type of component 2 a BasicShape mesh core.
 * See BasicShape::CoordType enum.
 */
BasicShape::CoordType UStructMesh::getCoordTypez(){
	return(getShape()->getCoordinateType(2));
}

/*! Return coordinates type of a BasicShape mesh core. See BasicShape::CoordType enum.
 * \return Type of all the cooordinates.
 */
array<BasicShape::CoordType, 3> UStructMesh::getCoordType(){
	array<BasicShape::CoordType, 3> types;
	for (int i=0; i<3; i++){
		types[i] = getShape()->getCoordinateType(i);
	}
	return(types);
}

/*! Return current mesh spacing */
darray3E UStructMesh::getSpacing(){
	darray3E res;
	darray3E scale = getScaling();
	res[0] = m_dx*scale[0]; res[1] =m_dy*scale[1]; res[2] = m_dz*scale[2];
	return(res); 
};

/*! Return current dimension of the mesh (number of mesh nodes in each direction) */
ivector1D UStructMesh::getDimension(){
	
	ivector1D res(3,0);
	res[0] = m_nx + 1; 
	res[1] = m_ny + 1; 
	res[2] = m_nz + 1;
	return(res); 
};

/*! Get n-th center cell coordinates in local mesh reference frame, 
 * given its global cell index on the mesh.  
 * \param[in] index cell index in the global nodal list.
 */
darray3E UStructMesh::getLocalCCell(int index){
	
	int i,j,k;
	accessCellIndex(index, i,j,k);
	darray3E res = getLocalCCell(i,j,k);
	
	return(res);
};
/*! Get n-th center cell coordinates in local mesh reference frame, 
 * given its cartesian cell indices on the mesh. 
 * \param[in] i_ x cartesian index.
 * \param[in] j_ y cartesian index.
 * \param[in] k_ z cartesian index.
 */
darray3E UStructMesh::getLocalCCell(int i_, int j_, int k_){
	
	darray3E res;
	res[0] = m_xnode[i_];
	res[1] = m_ynode[j_];
	res[2] = m_znode[k_];
	
	return(res);
};

/*! Get n-th nodal vertex coordinates in local mesh reference frame, 
 * given its global point index on the mesh. 
 * \param[in] index point index in the global nodal list.
 */
darray3E UStructMesh::getLocalPoint(int index){
	
	int i,j,k;
	accessPointIndex(index, i,j,k);
	darray3E res = getLocalPoint(i,j,k);
	
	return(res);
};
/*! Get n-th nodal vertex coordinates in local reference frame, 
 * given its cartesian point indices on the mesh. 
 * \param[in] i_ x cartesian index.
 * \param[in] j_ y cartesian index.
 * \param[in] k_ z cartesian index.
 */
darray3E UStructMesh::getLocalPoint(int i_, int j_, int k_){
	
	darray3E res;
	res[0] = m_xedge[i_];
	res[1] = m_yedge[j_];
	res[2] = m_zedge[k_];
	
	return(res);
};    

/*! Get n-th center cell coordinates in global absolute reference frame, 
 * given its global cell index on the mesh. 
 * \param[in] index cell index in the global nodal list.
 */
darray3E UStructMesh::getGlobalCCell(int index){
	
	darray3E res = getLocalCCell(index);
	return(transfToGlobal(res));
};
/*! Get n-th center cell coordinates in global absolute reference frame, 
 * given its cartesian cell indices on the mesh. 
 * \param[in] i_ x cartesian index.
 * \param[in] j_ y cartesian index.
 * \param[in] k_ z cartesian index.
 */
darray3E UStructMesh::getGlobalCCell(int i_, int j_, int k_){
	
	darray3E res = getLocalCCell(i_,j_,k_);
	return(transfToGlobal(res));
};

/*! Get n-th nodal vertex coordinates in global absolute reference frame, 
 * given its global point index on the mesh. 
 * \param[in] index point index in the global nodal list.
 */
darray3E UStructMesh::getGlobalPoint(int index){
	darray3E res = getLocalPoint(index);
	return(transfToGlobal(res));
};
/*! Get n-th nodal vertex coordinates in global absolute reference frame, 
 * given its cartesian point indices on the mesh. 
 * \param[in] i_ x cartesian index.
 * \param[in] j_ y cartesian index.
 * \param[in] k_ z cartesian index.
 */
darray3E UStructMesh::getGlobalPoint(int i_, int j_, int k_){
	
	darray3E res = getLocalPoint(i_,j_,k_);
	return(transfToGlobal(res));
};    

/*! Get neighbor vertices (by their global indices and in VTK-hexahedra ordered) of a given cell
 * \param[in] i_ x coordinate cell index
 * \param[in] j_ y coordinate cell index
 * \param[in] k_ z coordinate cell index
 */ 
ivector1D UStructMesh::getCellNeighs(int i_, int j_, int k_){
	
	ivector1D result(8);
	result[0] = accessPointIndex(i_, j_, k_);
	result[1] = accessPointIndex(i_+1, j_, k_);
	result[2] = accessPointIndex(i_+1, j_+1, k_);
	result[3] = accessPointIndex(i_, j_+1, k_);
	result[4] = accessPointIndex(i_, j_, k_+1);
	result[5] = accessPointIndex(i_+1, j_, k_+1);
	result[6] = accessPointIndex(i_+1, j_+1, k_+1);
	result[7] = accessPointIndex(i_, j_+1, k_+1);
	return(result);  
}

/*! Get neighbor vertices (by their global indices and in VTK-hexahedra ordered) of a given cell
 * \param[in] index cell global index
 */ 
ivector1D UStructMesh::getCellNeighs(int index){
	
	ivector1D pp(3,0);
	accessCellIndex(index, pp[0],pp[1], pp[2]);
	return(getCellNeighs(pp[0],pp[1],pp[2]));
}


dvecarr3E UStructMesh::getLocalCoords(){
	int np = (m_nx+1)*(m_ny+1)*(m_nz+1);
	dvecarr3E coords(np);
	for (int i=0; i<np; i++){
		coords[i] = getLocalPoint(i);
	}
	return coords;
};

dvecarr3E UStructMesh::getGlobalCoords(){
	int np = (m_nx+1)*(m_ny+1)*(m_nz+1);
	dvecarr3E coords(np);
	for (int i=0; i<np; i++){
		coords[i] = getGlobalPoint(i);
	}
	return coords;
};


/*! Set origin of your shape. The origin is meant as the baricenter of your shape in absolute r.s.
 * \param[in] origin new origin point
 */
void UStructMesh::setOrigin(darray3E origin){
	if (getShape() == NULL){
		std::cout << "none Shape set ---> exit" << std::endl;
		exit(2);
	}
	m_origin = origin;
	m_setorigin = true;
	getShape()->setOrigin(origin);
}

/*! Set span of your shape, according to its local reference system 
 * \param[in] s0 first coordinate span
 * \param[in] s1 second coordinate span
 * \param[in] s2 third coordinate span
 * \param[in] flag if true, mesh is rebuilt according to the new input.TRUE is default; 
 */
void UStructMesh::setSpan(double s0, double s1, double s2, bool flag){
	getShape()->setSpan( s0, s1,s2);
	m_span[0] = s0;
	m_span[1] = s1;
	m_span[2] = s2;
	m_setspan = true;
	if(flag){
		rebaseMesh();
	}
}

/*! Set span of your shape, according to its local reference system. After set the span
 * the mesh i rebuilt according to the new input.
 * \param[in] s coordinates span
 */
void UStructMesh::setSpan(darray3E s){
	getShape()->setSpan( s[0], s[1], s[2]);
	m_span = s;
	m_setspan = true;
	rebaseMesh();
}

/*! Set coordinate's origin of your shape, according to its local reference system  
 * \param[in] orig first coordinate origin
 * \param[in] dir 0,1,2 int flag identifying coordinate
 * \param[in] flag if true mesh is rebuilt according to the new input.TRUE is default
 */
void UStructMesh::setInfLimits(double orig, int dir, bool flag){
	getShape()->setInfLimits(orig, dir);
	m_origin[dir] = orig;
	m_setorigin = true;
	if(flag){
		rebaseMesh();
	}
}

/*! Set coordinates' origin of your shape, according to its local reference system.
 * Mesh is rebuilt according to the new input.
 * \param[in] orig coordinates' origin
 */
void UStructMesh::setInfLimits(darray3E orig){
	getShape()->setInfLimits(orig[0], 0);
	getShape()->setInfLimits(orig[1], 1);
	getShape()->setInfLimits(orig[2], 2);
	m_origin = orig;
	m_setorigin = true;
	rebaseMesh();
}

/*! Set new axis orientation of the local reference system of your mesh core shape
 * \param[in] axis0 first axis
 * \param[in] axis1 second axis
 * \param[in] axis2 third axis
 * 
 * if chosen axes are not orthogonal, doing nothing
 */
void UStructMesh::setRefSystem(darray3E axis0, darray3E axis1, darray3E axis2){
	getShape()->setRefSystem(axis0, axis1, axis2);
}

/*! Set new axis orientation of the local reference system of your mesh core shape
 * \param[in] int 0,1,2 identify local x,y,z axis of the primitive shape
 * \param[in] axis new direction of selected local axis.
 */
void UStructMesh::setRefSystem(int label, darray3E axis){
	getShape()->setRefSystem(label,axis);
}

/*! Set new axis orientation of the local reference system of your mesh core shape
 * \param[in] axis new direction of all local axes.
 */
void UStructMesh::setRefSystem(dmatrix33E axes){
	getShape()->setRefSystem(axes);
}

/*! Set the dimensions of the mesh (number of mesh nodes in each direction) */
 void UStructMesh::setDimension(ivector1D dim){
	 m_nx = dim[0] - 1;
	 m_ny = dim[1] - 1;
	 m_nz = dim[2] - 1;
	 rebaseMesh();
};

 /*! Set the dimensions of the mesh (number of mesh nodes in each direction) */
  void UStructMesh::setDimension(iarray3E dim){
 	 m_nx = dim[0] - 1;
 	 m_ny = dim[1] - 1;
 	 m_nz = dim[2] - 1;
 	 rebaseMesh();
 };

  /*! Set your shape, according to the following input parameters and the already saved/default parmaters.
   * It rebuilds the mesh after set the shape.
   * \param[in] type shape of your mesh, casted to enum.(option available are: 0-CUBE(default), 1-CYLINDER, 2-SPHERE)
   */
  void UStructMesh::setShape(int itype){
	  BasicShape::ShapeType type = static_cast<BasicShape::ShapeType>(itype);
	  UStructMesh::setShape(type);
  }

  /*! Set your shape, according to the following input parameters and the already saved/default parmaters.
   * It rebuilds the mesh after set the shape.
   * \param[in] type shape of your mesh, based on BasicShape::ShapeType enum.(option available are: CUBE(default), CYLINDER, SPHERE)
   */
  void UStructMesh::setShape(BasicShape::ShapeType type){
	  ivector1D dimLimit(3,2);
	  darray3E span, origin;
	  //create internal shape using unique_ptr member.
	  // unlink external shape eventually
  	m_shape1 = NULL;
  	if(m_shape2){m_shape2.release();}

  	//default values of origin (0,0,0) and span (1,1,1)cube (1,2*pi,pi)sphere (1,2*pi,1)cylinder

  	if (!m_setspan && !m_setorigin){
  		switch(type){
  		case BasicShape::ShapeType::CYLINDER :
  			m_shape2 = std::unique_ptr<BasicShape>(new Cylinder());
  			dimLimit[1] = 5;
  			break;
  		case BasicShape::ShapeType::SPHERE :
  			m_shape2 = std::unique_ptr<BasicShape>(new Sphere());
  			dimLimit[1] = 5; dimLimit[2] = 3;
  			break;
  		default://CUBE
  			m_shape2 = std::unique_ptr<BasicShape>(new Cube());
  			break;
  		}
  	}else if (!m_setspan && m_setorigin){
  		switch(type){
  		case BasicShape::ShapeType::CYLINDER :
  			span = {{1, 1, 1}};
  			m_shape2 = std::unique_ptr<BasicShape>(new Cylinder(m_origin, span));
  			dimLimit[1] = 5;
  			break;
  		case BasicShape::ShapeType::SPHERE :
  			span = {{1, 2*M_PI, M_PI}};
  			m_shape2 = std::unique_ptr<BasicShape>(new Sphere(m_origin, span));
  			dimLimit[1] = 5; dimLimit[2] = 3;
  			break;
  		default://CUBE
  			span = {{1, 2*M_PI, 1}};
  			m_shape2 = std::unique_ptr<BasicShape>(new Cube(m_origin, span));
  			break;
  		}
  	}else if (m_setspan && !m_setorigin){
  		switch(type){
  		case BasicShape::ShapeType::CYLINDER :
  			origin = {{0.0, 0.0, 0.0}};
  			m_shape2 = std::unique_ptr<BasicShape>(new Cylinder(origin, m_span));
  			dimLimit[1] = 5;
  			break;
  		case BasicShape::ShapeType::SPHERE :
  			origin = {{0.0, 0.0, 0.0}};
  			m_shape2 = std::unique_ptr<BasicShape>(new Sphere(origin, m_span));
  			dimLimit[1] = 5; dimLimit[2] = 3;
  			break;
  		default://CUBE
  			origin = {{0.0, 0.0, 0.0}};
  			m_shape2 = std::unique_ptr<BasicShape>(new Cube(origin, m_span));
  			break;
  		}
  	}else{
  		switch(type){
  		case BasicShape::ShapeType::CYLINDER :
  			m_shape2 = std::unique_ptr<BasicShape>(new Cylinder(m_origin, m_span));
  			dimLimit[1] = 5;
  			break;
  		case BasicShape::ShapeType::SPHERE :
  			m_shape2 = std::unique_ptr<BasicShape>(new Sphere(m_origin, m_span));
  			dimLimit[1] = 5; dimLimit[2] = 3;
  			break;
  		default://CUBE
  			m_shape2 = std::unique_ptr<BasicShape>(new Cube(m_origin, m_span));
  			break;
  		}
  	}

  	//check on dimensions and eventual closed loops on coordinates.
  	m_nx = std::max(m_nx, dimLimit[0])-1;
  	m_ny = std::max(m_ny, dimLimit[1])-1;
  	m_nz = std::max(m_nz, dimLimit[2])-1;

  	rebaseMesh();
  	m_setmesh = true;
  }


/*! Set mesh shape, according to the following input parameters.
 * It sets to default values of number of point in each direction [ (2,2,2) for a cube,
 * (2,5,2) for a cylinder, (2,5,3) for a sphere ].
 * It rebuilds a mesh after the set of the shape.
 * \param[in] shape pointer to an external allocated BasicShape object
 */
void UStructMesh::setShape(BasicShape * shape){

	ivector1D dimLimit(3,2);
	//create internal shape using unique_ptr member.
	// unlink external shape eventually
	m_shape1 = NULL;
	if(m_shape2){m_shape2.release();}

	m_shape1 = shape;

	switch(shape->getShapeType()){
		case BasicShape::ShapeType::CYLINDER :
			dimLimit[1] = 5;
			break;
		case BasicShape::ShapeType::SPHERE :
			dimLimit[1] = 5; dimLimit[2] = 3;
			break;
		default://CUBE
			break;
	}

	//check on dimensions and eventual closed loops on coordinates.
	m_nx = std::max(m_nx, dimLimit[0])-1;
	m_ny = std::max(m_ny, dimLimit[1])-1;
	m_nz = std::max(m_nz, dimLimit[2])-1;

	rebaseMesh();
	m_setmesh = true;

};

/*! Set your mesh, according to the following input parameters
 * \param[in] origin 3D point baricenter of your mesh 
 * \param[in] span span for each coordinate defining your mesh
 * \param[in] type   shape of your mesh, based on BasicShape::ShapeType enum.(option available are: CUBE(default), CYLINDER, SPHERE)
 * \param[in] dimensions number of mesh points for each coordinate.
 */
void UStructMesh::setMesh(darray3E & origin, darray3E &span, BasicShape::ShapeType type, ivector1D & dimensions){
	
	ivector1D dimLimit(3,2);
	//create internal shape using unique_ptr member.
	// unlink external shape eventually
	m_shape1 = NULL;
	if(m_shape2){m_shape2.release();}
	
	switch(type){
		case BasicShape::ShapeType::CYLINDER :
			m_shape2 = std::unique_ptr<BasicShape>(new Cylinder(origin, span));
			dimLimit[1] = 5;
			break;
		case BasicShape::ShapeType::SPHERE :
			m_shape2 = std::unique_ptr<BasicShape>(new Sphere(origin, span));
			dimLimit[1] = 5; dimLimit[2] = 3;
			break;
		default://CUBE
			m_shape2 = std::unique_ptr<BasicShape>(new Cube(origin, span));
		break;
	}
	
	//check on dimensions and eventual closed loops on coordinates.
	m_nx = std::max(dimensions[0], dimLimit[0])-1;
	m_ny = std::max(dimensions[1], dimLimit[1])-1;
	m_nz = std::max(dimensions[2], dimLimit[2])-1;
	
	m_origin = origin;
	m_setorigin = true;
	m_span = span;
	m_setspan = true;

	rebaseMesh();
	m_setmesh = true;
};

/*! Set your mesh, according to the following input parameters
 * \param[in] origin 3D point baricenter of your mesh 
 * \param[in] span span for each coordinate defining your mesh
 * \param[in] type   shape of your mesh, based on BasicShape::ShapeType enum.(option available are: CUBE(default), CYLINDER, SPHERE)
 * \param[in] spacing fixed spacing for each coordinate
 */
void UStructMesh::setMesh(darray3E & origin, darray3E &span, BasicShape::ShapeType type, dvector1D & spacing){

	ivector1D dimLimit(3,2);
	//create internal shape using unique_ptr member.
	// unlink external shape eventually
	m_shape1 = NULL;
	if(m_shape2){m_shape2.release();}
	
	switch(type){
		case BasicShape::ShapeType::CYLINDER :
			m_shape2 = std::unique_ptr<BasicShape>(new Cylinder(origin, span));
			dimLimit[1] = 5;
			break;
		case BasicShape::ShapeType::SPHERE :
			m_shape2 = std::unique_ptr<BasicShape>(new Sphere(origin, span));
			dimLimit[1] = 5; dimLimit[2] = 3;
			break;
		default://CUBE
			m_shape2 = std::unique_ptr<BasicShape>(new Cube(origin, span));
			break;
	}
	
	darray3E span2 = getSpan();
	ivector1D dim(3,0);
	
	for(int i=0; i<3; ++i){
		if(spacing[i] != 0.0) {
			dim[i] = (int) std::floor(span2[i]/spacing[i] +0.5) + 1;
		}else{
			dim[i] = dimLimit[i];
		}
	}
	
	//check on dimensions and eventual closed loops on coordinates.
	m_nx = std::max(dim[0], dimLimit[0])-1;
	m_ny = std::max(dim[1], dimLimit[1])-1;
	m_nz = std::max(dim[2], dimLimit[2])-1;
	
	m_origin = origin;
	m_setorigin = true;
	m_span = span;
	m_setspan = true;

	rebaseMesh();
	m_setmesh = true;
};

/*! Set your mesh, according to the following input parameters
 * \param[in] shape pointer to an external allocated BasicShape object
 * \param[in] dimensions number of mesh points for each coordinate.
 */
void UStructMesh::setMesh(BasicShape * shape, ivector1D & dimensions){

	ivector1D dimLimit(3,2);
	//create internal shape using unique_ptr member.
	// unlink external shape eventually
	m_shape1 = NULL;
	if(m_shape2){m_shape2.release();}

	m_shape1 = shape;
	
	switch(shape->getShapeType()){
		case BasicShape::ShapeType::CYLINDER :
			dimLimit[1] = 5;
			break;
		case BasicShape::ShapeType::SPHERE :
			dimLimit[1] = 5; dimLimit[2] = 3;
			break;
		default://CUBE
			break;
	}
	
	//check on dimensions and eventual closed loops on coordinates.
	m_nx = std::max(dimensions[0], dimLimit[0])-1;
	m_ny = std::max(dimensions[1], dimLimit[1])-1;
	m_nz = std::max(dimensions[2], dimLimit[2])-1;
	
	m_origin = {{0.0,0.0,0.0}};
	m_setorigin = false;
	m_span = {{0.0,0.0,0.0}};
	m_setspan = false;

	rebaseMesh();
	m_setmesh = true;
};

/*!Set your mesh, according to the following input parameters
 * \param[in] shape pointer to an external allocated BasicShape object
 * \param[in] spacing fixed spacing for each coordinate
 */
void UStructMesh::setMesh(BasicShape * shape, dvector1D & spacing){
	
	ivector1D dimLimit(3,2);
	//create internal shape using unique_ptr member.
	// unlink external shape eventually
	m_shape1 = NULL;
	if(m_shape2){m_shape2.release();}
	
	m_shape1 = shape;
	
	switch(shape->getShapeType()){
		case BasicShape::ShapeType::CYLINDER :
			dimLimit[1] = 5;
			break;
		case BasicShape::ShapeType::SPHERE :
			dimLimit[1] = 5; dimLimit[2] = 3;
			break;
		default://CUBE
			break;
	}
	
	darray3E span2 = shape->getSpan();
	ivector1D dim(3,0);
	
	for(int i=0; i<3; ++i){
		
		if(spacing[i] != 0.0) {
			dim[i] = (int) std::floor(span2[i]/spacing[i] +0.5) + 1;
		}else{
			dim[i] = dimLimit[i];
		}
	}
	
	//check on dimensions and eventual closed loops on coordinates.
	m_nx = std::max(dim[0], dimLimit[0]) -1;
	m_ny = std::max(dim[1], dimLimit[1]) -1;
	m_nz = std::max(dim[2], dimLimit[2]) -1;
	
	m_origin = {{0.0,0.0,0.0}};
	m_setorigin = false;
	m_span = {{0.0,0.0,0.0}};
	m_setspan = false;

	rebaseMesh();
	m_setmesh = true;
};


/*!Clear the Mesh structure. Unlink external shapes or destroy internal shapes, destroy nodal structure.*/
void UStructMesh::clearMesh(){
	
	m_shape1 = NULL;
	m_shape2.release();
	m_nx=0; m_ny=0; m_nz=0;
	m_dx=0.0; m_dy=0.0; m_dz=0.0;
	m_origin = {{0.0,0.0,0.0}};
	m_setorigin = false;
	m_span = {{0.0,0.0,0.0}};
	m_setspan = false;

	destroyNodalStructure();
	m_setmesh = false;
};  


/*! Return cartesian indices of the cell containing the target point in global reference frame
 * \param[in] point 3D coordinate of target point
 * \param[out] i x cell index
 * \param[out] j y cell index
 * \param[out] k z cell index 
 */ 
void UStructMesh::locateCellByPoint(darray3E & point, int &i, int &j, int &k){
	
	darray3E P = transfToLocal(point);
	darray3E locOr = getShape()->getLocalOrigin();

	i = min(m_nx-1, max(0, (int) floor((P[0]-locOr[0])/m_dx)));
	j = min(m_ny-1, max(0, (int) floor((P[1]-locOr[1])/m_dy)));
	k = min(m_nz-1, max(0, (int) floor((P[2]-locOr[2])/m_dz)));
};

/*! Return cartesian indices of the cell containing the target point in global reference framne
 * \param[in] point 3D coordinate of target point
 * \param[out] i x cell index
 * \param[out] j y cell index
 * \param[out] k z cell index 
 */ 
void UStructMesh::locateCellByPoint(dvector1D & point, int &i, int &j, int &k){
		darray3E temp = conArray<double,3>(point);
		locateCellByPoint(temp,i,j,k);
};

/*! Return global index of the cell given its cartesian indices. Follows the ordering sequences z-y-x
 * \param[in] i x cartesian index
 *\param[in] j y cartesian index
 *\param[in] k z cartesian index
 *\param[out] result global index 
 */
int UStructMesh::accessCellIndex(int i, int j, int k){
	
	int index = m_ny * m_nz * i + m_nz * j + k;
	return(index);
};

/*! Return cartesian indices of the cell given its global index. Follows the ordering sequences z-y-x
 * \param[in] N_ global index 
 *\param[out] i x cartesian index
 *\param[out] j y cartesian index
 *\param[out] k z cartesian index
 */
void UStructMesh::accessCellIndex(int N_, int & i, int & j, int & k){
	k = N_ % m_nz;
	int index = N_ / m_nz;
	j = index % m_ny;
	i = index / m_ny; 
};

/*! Return cartesian indices of the point given its global index. Follows the ordering sequences z-y-x
 * \param[in] N_ global index 
 *\param[out] i x cartesian index
 *\param[out] j y cartesian index
 *\param[out] k z cartesian index
 */
void UStructMesh::accessPointIndex(int N_, int &i, int &j, int &k){
	
	k = N_ % (m_nz+1);
	int index = N_ / (m_nz+1);
	j = index % (m_ny+1);
	i = index / (m_ny+1); 
};

/*! Transform point from Local to Global ref system*/
darray3E	UStructMesh::transfToGlobal( darray3E & point){
	return(getShape()->toWorldCoord(point));
};
/*! Transform point from Local to Global ref system*/
dvector1D	UStructMesh::transfToGlobal( dvector1D & point){
	darray3E temp = conArray<double,3>(point);
	darray3E temp2 = getShape()->toWorldCoord(temp);
	return(conVect(temp2));
};
/*! Transform list of points from Local to Global ref system*/
dvecarr3E	UStructMesh::transfToGlobal( dvecarr3E & list_points){
	int size = list_points.size();
	dvecarr3E result(size);
	for(int i=0; i<size; ++i){
		result[i] = transfToGlobal(list_points[i]); 
	}
	return(result);
};    

/*! Transform point from Global to Local ref system*/
darray3E 	UStructMesh::transfToLocal( darray3E & point){
	return(getShape()->toLocalCoord(point));
};
/*! Transform point from Global to Local ref system*/
dvector1D 	UStructMesh::transfToLocal( dvector1D & point){
	darray3E temp = conArray<double,3>(point);
	darray3E temp2 = getShape()->toLocalCoord(temp);
	return(conVect(temp2));
	
};
/*! Transform point from Global to Local ref system*/
dvecarr3E 	UStructMesh::transfToLocal( dvecarr3E & list_points){
	int size = list_points.size();
	dvecarr3E result(size);
	for(int i=0; i<size; ++i){
		result[i] = transfToLocal(list_points[i]); 
	}
	return(result);
	
};  


/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on centercells
 * \param[out] interpResult interpolated value
 */
double UStructMesh::interpolateCellData(darray3E & point, dvector1D & celldata){
	
	int i0, j0, k0, ip, jp, kp;
	double wx0,wx1,wy0,wy1,wz0,wz1;
	
	locateCellByPoint(point, i0, j0, k0);
	darray3E P = transfToLocal(point);
	if (P[0] > m_xnode[i0]) 	{ ip = min(i0+1, m_nx-1);   }
	else                  	{ ip = max(0, i0-1);      }
	if (P[1] > m_ynode[j0]) 	{ jp = min(j0+1, m_ny-1);   }
	else                  	{ jp = max(0, j0-1);      }
	if (P[2] > m_znode[k0]) 	{ kp = min(k0+1, m_nz-1);   }
	else                  	{ kp = max(0, k0-1);      }
	
	// Interpolation weights
	wx1 = max(0.0, min(1.0, abs((P[0] - m_xnode[i0])/m_dx)));     wx0 = 1.0 - wx1;
	wy1 = max(0.0, min(1.0, abs((P[1] - m_ynode[j0])/m_dy)));     wy0 = 1.0 - wy1;
	wz1 = max(0.0, min(1.0, abs((P[2] - m_znode[k0])/m_dz)));     wz0 = 1.0 - wz1;
	
	// Interpolation
	double result  = wz0 * wx0 * wy0 * celldata[accessCellIndex(i0,j0,k0)]
	+ wz0 * wx0 * wy1 * celldata[accessCellIndex(i0,jp,k0)]
	+ wz0 * wx1 * wy0 * celldata[accessCellIndex(ip,j0,k0)]
	+ wz0 * wx1 * wy1 * celldata[accessCellIndex(ip,jp,k0)]
	+ wz1 * wx0 * wy0 * celldata[accessCellIndex(i0,j0,kp)]
	+ wz1 * wx0 * wy1 * celldata[accessCellIndex(i0,jp,kp)]
	+ wz1 * wx1 * wy0 * celldata[accessCellIndex(ip,j0,kp)]
	+ wz1 * wx1 * wy1 * celldata[accessCellIndex(ip,jp,kp)];
	
	return(result);      
};

/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on centercells
 * \param[out] interpResult interpolated value
 */
int UStructMesh::interpolateCellData(darray3E & point, ivector1D & celldata){
	
	int i0, j0, k0, ip, jp, kp;
	double wx0,wx1,wy0,wy1,wz0,wz1;
	
	locateCellByPoint(point, i0, j0, k0);
	darray3E P = transfToLocal(point);
	if (P[0] > m_xnode[i0]) 	{ ip = min(i0+1, m_nx-1);   }
	else                  	{ ip = max(0, i0-1);      }
	if (P[1] > m_ynode[j0]) 	{ jp = min(j0+1, m_ny-1);   }
	else                  	{ jp = max(0, j0-1);      }
	if (P[2] > m_znode[k0]) 	{ kp = min(k0+1, m_nz-1);   }
	else                  	{ kp = max(0, k0-1);      }
	
	// Interpolation weights
	wx1 = max(0.0, min(1.0, abs((P[0] - m_xnode[i0])/m_dx)));     wx0 = 1.0 - wx1;
	wy1 = max(0.0, min(1.0, abs((P[1] - m_ynode[j0])/m_dy)));     wy0 = 1.0 - wy1;
	wz1 = max(0.0, min(1.0, abs((P[2] - m_znode[k0])/m_dz)));     wz0 = 1.0 - wz1;
	
	// Interpolation
	double result  = 	wz0 * wx0 * wy0 * (double)celldata[accessCellIndex(i0,j0,k0)]
	+ wz0 * wx0 * wy1 * (double)celldata[accessCellIndex(i0,jp,k0)]
	+ wz0 * wx1 * wy0 *(double) celldata[accessCellIndex(ip,j0,k0)]
	+ wz0 * wx1 * wy1 *(double) celldata[accessCellIndex(ip,jp,k0)]
	+ wz1 * wx0 * wy0 *(double) celldata[accessCellIndex(i0,j0,kp)]
	+ wz1 * wx0 * wy1 *(double) celldata[accessCellIndex(i0,jp,kp)]
	+ wz1 * wx1 * wy0 *(double) celldata[accessCellIndex(ip,j0,kp)]
	+ wz1 * wx1 * wy1 *(double) celldata[accessCellIndex(ip,jp,kp)]; 
	
	int result2 =std::floor(result+0.5); 		
	return(result2);
};

/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on centercells
 * \param[out] interpResult interpolated value
 */
darray3E UStructMesh::interpolateCellData(darray3E & point, dvecarr3E & celldata){
	
	int i0, j0, k0, ip, jp, kp;
	double wx0,wx1,wy0,wy1,wz0,wz1;
	
	locateCellByPoint(point, i0, j0, k0);
	darray3E P = transfToLocal(point);
	
	if (P[0] > m_xnode[i0]) 	{ ip = min(i0+1, m_nx-1);   }
	else                  	{ ip = max(0, i0-1);      }
	if (P[1] > m_ynode[j0]) 	{ jp = min(j0+1, m_ny-1);   }
	else                  	{ jp = max(0, j0-1);      }
	if (P[2] > m_znode[k0]) 	{ kp = min(k0+1, m_nz-1);   }
	else                  	{ kp = max(0, k0-1);      }
	
	// Interpolation weights
	wx1 = max(0.0, min(1.0, abs((P[0] - m_xnode[i0])/m_dx)));     wx0 = 1.0 - wx1;
	wy1 = max(0.0, min(1.0, abs((P[1] - m_ynode[j0])/m_dy)));     wy0 = 1.0 - wy1;
	wz1 = max(0.0, min(1.0, abs((P[2] - m_znode[k0])/m_dz)));     wz0 = 1.0 - wz1;
	
	// Interpolation
	darray3E result = wz0 * wx0 * wy0 * celldata[accessCellIndex(i0,j0,k0)]
	+ wz0 * wx0 * wy1 * celldata[accessCellIndex(i0,jp,k0)]
	+ wz0 * wx1 * wy0 * celldata[accessCellIndex(ip,j0,k0)]
	+ wz0 * wx1 * wy1 * celldata[accessCellIndex(ip,jp,k0)]
	+ wz1 * wx0 * wy0 * celldata[accessCellIndex(i0,j0,kp)]
	+ wz1 * wx0 * wy1 * celldata[accessCellIndex(i0,jp,kp)]
	+ wz1 * wx1 * wy0 * celldata[accessCellIndex(ip,j0,kp)]
	+ wz1 * wx1 * wy1 * celldata[accessCellIndex(ip,jp,kp)];
	
	return(result);	      
};

/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] pointdata data field defined on grid points
 * \param[out] interpResult interpolated value
 */
double UStructMesh::interpolatePointData(darray3E & point, dvector1D & pointdata){
	
	int i0, j0, k0, ip, jp, kp;
	double wx0,wx1,wy0,wy1,wz0,wz1;
	
	darray3E P = transfToLocal(point);
	darray3E locOr = getShape()->getLocalOrigin();
	i0 = max(0, min(m_nx, (int) floor((P[0]-locOr[0])/m_dx)));
	j0 = max(0, min(m_ny, (int) floor((P[1]-locOr[1])/m_dy)));
	k0 = max(0, min(m_nz, (int) floor((P[2]-locOr[2])/m_dz)));
	
	if (P[0] >= m_xedge[i0]) 	{ ip = min(i0+1, m_nx);   }
	else                  	{ ip = max(0, i0-1);      }
	if (P[1] >= m_yedge[j0]) 	{ jp = min(j0+1, m_ny);   }
	else                  	{ jp = max(0, j0-1);      }
	if (P[2] >= m_zedge[k0]) 	{ kp = min(k0+1, m_nz);   }
	else                  	{ kp = max(0, k0-1);      }
	
	// Interpolation weights
	wx1 = max(0.0, min(1.0, abs((P[0] - m_xedge[i0])/m_dx)));     wx0 = 1.0 - wx1;
	wy1 = max(0.0, min(1.0, abs((P[1] - m_yedge[j0])/m_dy)));     wy0 = 1.0 - wy1;
	wz1 = max(0.0, min(1.0, abs((P[2] - m_zedge[k0])/m_dz)));     wz0 = 1.0 - wz1;
	
	// Interpolation
	double result =  wz0 * wx0 * wy0 * pointdata[accessPointIndex(i0,j0,k0)]
	+ wz0 * wx0 * wy1 * pointdata[accessPointIndex(i0,jp,k0)]
	+ wz0 * wx1 * wy0 * pointdata[accessPointIndex(ip,j0,k0)]
	+ wz0 * wx1 * wy1 * pointdata[accessPointIndex(ip,jp,k0)]
	+ wz1 * wx0 * wy0 * pointdata[accessPointIndex(i0,j0,kp)]
	+ wz1 * wx0 * wy1 * pointdata[accessPointIndex(i0,jp,kp)]
	+ wz1 * wx1 * wy0 * pointdata[accessPointIndex(ip,j0,kp)]
	+ wz1 * wx1 * wy1 * pointdata[accessPointIndex(ip,jp,kp)];
	return(result);
};

/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] pointdata data field defined on grid points
 * \param[out] interpResult interpolated value
 */
int UStructMesh::interpolatePointData(darray3E & point, ivector1D & pointdata){
	
	int i0, j0, k0, ip, jp, kp;
	double wx0,wx1,wy0,wy1,wz0,wz1;
	
	darray3E P = transfToLocal(point);
	darray3E locOr = getShape()->getLocalOrigin();
	i0 = max(0, min(m_nx, (int) floor((P[0]-locOr[0])/m_dx)));
	j0 = max(0, min(m_ny, (int) floor((P[1]-locOr[1])/m_dy)));
	k0 = max(0, min(m_nz, (int) floor((P[2]-locOr[2])/m_dz)));
	
	if (P[0] >= m_xedge[i0]) 	{ ip = min(i0+1, m_nx);   }
	else                  	{ ip = max(0, i0-1);      }
	if (P[1] >= m_yedge[j0]) 	{ jp = min(j0+1, m_ny);   }
	else                  	{ jp = max(0, j0-1);      }
	if (P[2] >= m_zedge[k0]) 	{ kp = min(k0+1, m_nz);   }
	else                  	{ kp = max(0, k0-1);      }
	
	// Interpolation weights
	wx1 = max(0.0, min(1.0, abs((P[0] - m_xedge[i0])/m_dx)));     wx0 = 1.0 - wx1;
	wy1 = max(0.0, min(1.0, abs((P[1] - m_yedge[j0])/m_dy)));     wy0 = 1.0 - wy1;
	wz1 = max(0.0, min(1.0, abs((P[2] - m_zedge[k0])/m_dz)));     wz0 = 1.0 - wz1;
	
	// Interpolation
	double result = wz0 * wx0 * wy0 * (double)pointdata[accessPointIndex(i0,j0,k0)]
	+ wz0 * wx0 * wy1 * (double)pointdata[accessPointIndex(i0,jp,k0)]
	+ wz0 * wx1 * wy0 *(double) pointdata[accessPointIndex(ip,j0,k0)]
	+ wz0 * wx1 * wy1 *(double) pointdata[accessPointIndex(ip,jp,k0)]
	+ wz1 * wx0 * wy0 *(double) pointdata[accessPointIndex(i0,j0,kp)]
	+ wz1 * wx0 * wy1 *(double) pointdata[accessPointIndex(i0,jp,kp)]
	+ wz1 * wx1 * wy0 *(double) pointdata[accessPointIndex(ip,j0,kp)]
	+ wz1 * wx1 * wy1 *(double) pointdata[accessPointIndex(ip,jp,kp)];
	
	int result2 = std::floor(result+0.5);			
	return(result2);
};

/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] pointdata data field defined on grid points
 * \param[out] interpResult interpolated value
 */
darray3E UStructMesh::interpolatePointData(darray3E & point, dvecarr3E & pointdata){
	
	int i0, j0, k0, ip, jp, kp;
	double wx0,wx1,wy0,wy1,wz0,wz1;
	
	darray3E P = transfToLocal(point);
	darray3E locOr = getShape()->getLocalOrigin();
	i0 = max(0, min(m_nx, (int) floor((P[0]-locOr[0])/m_dx)));
	j0 = max(0, min(m_ny, (int) floor((P[1]-locOr[1])/m_dy)));
	k0 = max(0, min(m_nz, (int) floor((P[2]-locOr[2])/m_dz)));
	
	
	if (P[0] >= m_xedge[i0]) 	{ ip = min(i0+1, m_nx);   }
	else                  	{ ip = max(0, i0-1);      }
	if (P[1] >= m_yedge[j0]) 	{ jp = min(j0+1, m_ny);   }
	else                  	{ jp = max(0, j0-1);      }
	if (P[2] >= m_zedge[k0]) 	{ kp = min(k0+1, m_nz);   }
	else                  	{ kp = max(0, k0-1);      }
	
	// Interpolation weights
	wx1 = max(0.0, min(1.0, abs((P[0] - m_xedge[i0])/m_dx)));     wx0 = 1.0 - wx1;
	wy1 = max(0.0, min(1.0, abs((P[1] - m_yedge[j0])/m_dy)));     wy0 = 1.0 - wy1;
	wz1 = max(0.0, min(1.0, abs((P[2] - m_zedge[k0])/m_dz)));     wz0 = 1.0 - wz1;
	
	// Interpolation
	darray3E result = wz0 * wx0 * wy0 * pointdata[accessPointIndex(i0,j0,k0)]
	+ wz0 * wx0 * wy1 * pointdata[accessPointIndex(i0,jp,k0)]
	+ wz0 * wx1 * wy0 * pointdata[accessPointIndex(ip,j0,k0)]
	+ wz0 * wx1 * wy1 * pointdata[accessPointIndex(ip,jp,k0)]
	+ wz1 * wx0 * wy0 * pointdata[accessPointIndex(i0,j0,kp)]
	+ wz1 * wx0 * wy1 * pointdata[accessPointIndex(i0,jp,kp)]
	+ wz1 * wx1 * wy0 * pointdata[accessPointIndex(ip,j0,kp)]
	+ wz1 * wx1 * wy1 * pointdata[accessPointIndex(ip,jp,kp)];
	
	return(result);
};


/*! Write your grid as a point cloud in a *.vtu file. 
 * \param[in] folder output directory path
 * \param[in] outfile output namefile w/out tag
 * \param[in] counterFile integer to mark output with a counter number
 * \param[in] codexFlag boolean to distinguish between "ascii" false, and "appended" true
 * \param[in] extPoints OPTIONAL. If defined, use the current vertexList and write them as point cloud, provided that coherent dimensions with the mesh are set. 
 *			If Not, write the current mesh vertices. 
 * */
void UStructMesh::plotCloud( std::string & folder , std::string outfile, int counterFile, bool codexFlag, dvecarr3E * extPoints){
	
	bitpit::VTKFormat codex = bitpit::VTKFormat::ASCII;
	if(codexFlag){codex=bitpit::VTKFormat::APPENDED;}
	
	ivector1D dim = getDimension();
	int sizeTot = dim[0]*dim[1]*dim[2];
	
	VTK_BASICCLOUD handle_vtk_output(folder, outfile, codex, sizeTot);
	
	dvecarr3E activeP;
	if(extPoints != NULL && extPoints->size() == sizeTot){activeP = *extPoints; }
	else{
		activeP.resize(sizeTot);
		for(int i=0; i<sizeTot; i++){
			activeP[i] = getGlobalPoint(i);
		}
	}
	
	if(counterFile>=0){handle_vtk_output.setCounter(counterFile);}
	handle_vtk_output.setGeomTypes(bitpit::VTKDataType::Float64,bitpit::VTKDataType::Int32, bitpit::VTKDataType::Int32, bitpit::VTKDataType::Int32);
	handle_vtk_output.linkData(activeP);
	handle_vtk_output.write();
};

/*! Write your grid as a point cloud in a *.vtu file. 
 * \param[in] folder output directory path
 * \param[in] outfile output namefile w/out tag
 * \param[in] counterFile integer to mark output with a counter number
 * \param[in] codexFlag boolean to distinguish between "ascii" false, and "appended" true
 * \param[in] vertexList list of global indices of selected vertices that the user wants to write on file. 
 * \param[in] extPoints OPTIONAL. If defined, use the current vertexList and write them as point cloud, provided that coherent dimensions with the mesh are set. 
 *			If Not, write the current mesh vertices. 
 * */
void UStructMesh::plotCloud( std::string & folder , std::string outfile, int counterFile, bool codexFlag, ivector1D & vertexList, dvecarr3E * extPoints){
	
	bitpit::VTKFormat codex = bitpit::VTKFormat::ASCII;
	if(codexFlag){codex=bitpit::VTKFormat::APPENDED;}
	
	ivector1D dim = getDimension();
	int sizeTot = dim[0]*dim[1]*dim[2];
	int sizeMap = std::min((int)vertexList.size(), sizeTot);
	
	VTK_BASICCLOUD handle_vtk_output(folder, outfile, codex, sizeMap);
	
	dvecarr3E activeP(sizeMap);
	if(extPoints != NULL && extPoints->size() == sizeTot){
		for(int i=0; i<sizeMap; i++){
			activeP[i] = (*extPoints)[vertexList[i]];
		}  
	}
	else{
		for(int i=0; i<sizeMap; i++){
			activeP[i] = getGlobalPoint(vertexList[i]);
		}
	}
	
	if(counterFile>=0){handle_vtk_output.setCounter(counterFile);}
	handle_vtk_output.setGeomTypes(bitpit::VTKDataType::Float64,bitpit::VTKDataType::Int32, bitpit::VTKDataType::Int32, bitpit::VTKDataType::Int32);
	handle_vtk_output.linkData(activeP);
	handle_vtk_output.write();
	
};

/*! Write your grid as a hexahedrical one in a *.vtu file. 
 * \param[in] folder output directory path
 * \param[in] outfile output namefile w/out tag
 * \param[in] counterFile integer to mark output with a counter number
 * \param[in] codexFlag boolean to distinguish between "ascii" false, and "appended" true
 * \param[in] extPoints OPTIONAL. If defined, use the current vertexList and write them as point cloud, provided that coherent dimensions with the mesh are set. 
 *			If Not, write the current mesh vertices. 
 * */
void UStructMesh::plotGrid(std::string & folder, std::string outfile , int counterFile, bool codexFlag, dvecarr3E * extPoints){
	
	bitpit::VTKFormat codex = bitpit::VTKFormat::ASCII;
	if(codexFlag){codex=bitpit::VTKFormat::APPENDED;}
	
	ivector1D dim = getDimension();
	int sizePt = dim[0]*dim[1]*dim[2];
	int sizeCl = (dim[0]-1)*(dim[1]-1)*(dim[2]-1);
	VTK_BASICMESH handle_vtk_output(folder, outfile, codex, sizePt, sizeCl, 8*sizeCl);
	
	dvecarr3E activeP(sizePt);
	ivector2D activeConn(sizeCl, ivector1D(8,0));
	
	if(extPoints != NULL && extPoints->size() == sizePt){activeP = *extPoints;}
	else{
		for(int i=0; i<sizePt; i++){
			activeP[i] = getGlobalPoint(i);
		}
	}
	
	for(int i=0; i<sizeCl; ++i){
		activeConn[i] = getCellNeighs(i); 
	}
	
	if(counterFile>=0){handle_vtk_output.setCounter(counterFile);}
	handle_vtk_output.setGeomTypes(bitpit::VTKDataType::Float64,bitpit::VTKDataType::Int32, bitpit::VTKDataType::Int32, bitpit::VTKDataType::Int32);
	handle_vtk_output.linkData(activeP,activeConn);
	handle_vtk_output.write();
};

/*! Write your grid as a hexahedrical one in a *.vtu file. 
 * \param[in] folder output directory path
 * \param[in] outfile output namefile w/out tag
 * \param[in] counterFile integer to mark output with a counter number
 * \param[in] codexFlag boolean to distinguish between "ascii" false, and "appended" true
 * \param[in] vertexList list of global indices of selected cells that the user wants to write on file. 
 * \param[in] extPoints OPTIONAL. If defined, use the current vertexList and write them as point cloud, provided that coherent dimensions with the mesh are set. 
 *			If Not, write the current mesh vertices. 
 * */

void UStructMesh::plotGrid(std::string & folder, std::string outfile, int counterFile, bool codexFlag, ivector1D & cellList, dvecarr3E * extPoints){
	
	bitpit::VTKFormat codex = bitpit::VTKFormat::ASCII;
    if(codexFlag){codex=bitpit::VTKFormat::APPENDED;}
	
	ivector1D dim = getDimension();
	int sizePt = dim[0]*dim[1]*dim[2];
	int sizeCl = cellList.size();
	
	if(sizeCl > (dim[0]-1)*(dim[1]-1)*(dim[2]-1)){return;}
	
	std::map<int,int> mapPoints;
	ivector2D activeConn(sizeCl, ivector1D(8,0));
	
	for(int i=0; i<sizeCl; ++i){
		activeConn[i] = getCellNeighs(i);
		for(int j=0; j<activeConn[i].size(); ++j){
			mapPoints[activeConn[i][j]] = activeConn[i][j];
		}
	}
	
	dvecarr3E activeP(mapPoints.size());
	ivector1D listP(mapPoints.size());
	int counter = 0;
	std::map<int,int>::iterator itF;
	if(extPoints != NULL && extPoints->size() == sizePt){
		for(itF = mapPoints.begin(); itF != mapPoints.end(); ++itF){
			int pos = itF->second;
			activeP[counter] = (*extPoints)[pos];
			listP[counter] = pos;
			++counter;
		}
	}
	else{
		for(itF = mapPoints.begin(); itF != mapPoints.end(); ++itF){
			int pos = itF->second;
			activeP[counter] = getGlobalPoint(pos);
			listP[counter]=  pos;
			++counter;
		}
	}
	
	//update connectivity w/ local indexing of vertex;
	
	for(int i=0; i<sizeCl; ++i){
		for(int j=0; j<activeConn[i].size(); ++j){
			int posV =posVectorFind(listP, activeConn[i][j]);
			//update local connectivity
			activeConn[i][j] = posV;
		}
	}
	
	VTK_BASICMESH handle_vtk_output(folder, outfile, codex, counter, sizeCl, 8*sizeCl);
	if(counterFile>=0){handle_vtk_output.setCounter(counterFile);}
	handle_vtk_output.setGeomTypes(bitpit::VTKDataType::Float64,bitpit::VTKDataType::Int32, bitpit::VTKDataType::Int32, bitpit::VTKDataType::Int32);
	handle_vtk_output.linkData(activeP,activeConn);
	handle_vtk_output.write();
};

/*! Return a pointer to the inner BasicShape object the current mesh is built on 
 * \param[out] result BasicShape of the mesh
 */
BasicShape * UStructMesh::getShape(){
	
	BasicShape * result;
	if(m_shape2){
		result = m_shape2.get();
	}else{
		result = m_shape1;
	}
	return(result);
}

/*! Destroy the all nodal structures of the mesh. */
void UStructMesh::destroyNodalStructure(){
	m_xnode.clear();
	m_ynode.clear();
	m_znode.clear();
	m_xedge.clear();
	m_yedge.clear();
	m_zedge.clear();  
};

/*! Destroy the all nodal structures of the mesh, and reinitialize them to current mesh dimensions. */
void UStructMesh::reshapeNodalStructure(){
	destroyNodalStructure();
	resizeMesh();
};

/*! Resize nodal structure to current dimension set */
void UStructMesh::resizeMesh(){
	// Cell centers
	m_xnode.resize(m_nx, 0.0);
	m_ynode.resize(m_ny, 0.0);
	m_znode.resize(m_nz, 0.0);

	// Points
	m_xedge.resize(m_nx+1, 0.0);
	m_yedge.resize(m_ny+1, 0.0);
	m_zedge.resize(m_nz+1, 0.0);
};

/*! Check current mesh dimensions and refresh nodal structures*/
void UStructMesh::rebaseMesh(){
	
	if(getShape() == NULL){return;}
	darray3E spanEff = getLocalSpan();
	
	reshapeNodalStructure();
	
	m_dx = spanEff[0]/m_nx;
	m_dy = spanEff[1]/m_ny;
	m_dz = spanEff[2]/m_nz;
	darray3E locOr = getShape()->getLocalOrigin();
	// get point distro;
	for (int i = 0; i < m_nx+1; i++) {m_xedge[i] = locOr[0] + ((double) i) * m_dx;} 
	for (int i = 0; i < m_ny+1; i++) {m_yedge[i] = locOr[1] + ((double) i) * m_dy;}
	for (int i = 0; i < m_nz+1; i++) {m_zedge[i] = locOr[2] + ((double) i) * m_dz;}
	// get cell distro
	for (int i = 0; i < m_nx; i++) {m_xnode[i] = m_xedge[i] + 0.5 * m_dx;}
	for (int i = 0; i < m_ny; i++) {m_ynode[i] = m_yedge[i] + 0.5 * m_dy;}
	for (int i = 0; i < m_nz; i++) {m_znode[i] = m_zedge[i] + 0.5 * m_dz;}
	
};

void UStructMesh::execute(){
	rebaseMesh();
};













