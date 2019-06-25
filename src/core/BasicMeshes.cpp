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

#include "BasicMeshes.hpp"
#include "Operators.hpp"
#include "customOperators.hpp"

using namespace std;
namespace mimmo{


/*!
 * Basic Constructor. Default Shape CUBE
 */
UStructMesh::UStructMesh(){
	m_nx = 0; m_ny=0; m_nz=0;
	m_dx = 0; m_dy=0; m_dz=0;

	m_setorigin = false;
	m_setspan = false;
	m_setInfLimits = false;
	m_setRefSys = false;
	m_isBuild = false;
	m_origin_temp = {{0.0,0.0,0.0}};
	m_span_temp = {{1.0,1.0,1.0}};
	m_inflimits_temp = {{0.0,0.0,0.0}};
	for(int i=0; i<3; ++i){m_refsystem_temp[i].fill(0.0); m_refsystem_temp[i][i] = 1.0;}
	m_shapetype_temp = ShapeType::CUBE;
};

/*!
 * Basic destructor
 */
UStructMesh::~UStructMesh(){
	m_xnode.clear(); m_ynode.clear(); m_znode.clear();
	m_xedge.clear(); m_yedge.clear(); m_zedge.clear();
}

/*!
 * Copy Constructor.
 * \param[in] other UStructMesh object where copy from
 */
UStructMesh::UStructMesh(const UStructMesh & other){

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

    if(other.m_shape){
        switch(other.m_shape->getShapeType()){
            case ShapeType::WEDGE :
            {
                const Wedge * pp = dynamic_cast<const Wedge*>(other.m_shape.get());
                if (pp != NULL)  m_shape = std::unique_ptr<BasicShape>(new Wedge(*(pp)));
            }
            break;
            case ShapeType::CYLINDER :
                {
                    const Cylinder * pp = dynamic_cast<const Cylinder*>(other.m_shape.get());
                    if (pp != NULL)  m_shape = std::unique_ptr<BasicShape>(new Cylinder(*(pp)));
                }
                break;
            case ShapeType::SPHERE :
                {
                    const Sphere * pp = dynamic_cast<const Sphere*>(other.m_shape.get());
                    if (pp != NULL)  m_shape = std::unique_ptr<BasicShape>(new Sphere(*(pp)));
                }
                break;
            default://CUBE
                {
                    const Cube * pp = dynamic_cast<const Cube*>(other.m_shape.get());
                    if (pp != NULL)  m_shape = std::unique_ptr<BasicShape>(new Cube(*(pp)));
                }
                break;
        }
    }

    //copy temporary members
    m_origin_temp = other.m_origin_temp;
    m_span_temp = other.m_span_temp;
    m_inflimits_temp = other.m_inflimits_temp;
    m_refsystem_temp = other.m_refsystem_temp;
    m_shapetype_temp = other.m_shapetype_temp;
    m_setorigin = other.m_setorigin;
    m_setspan = other.m_setspan;
    m_setInfLimits = other.m_setInfLimits;
    m_setRefSys = other.m_setRefSys;
    m_isBuild = other.m_isBuild;

};

/*! Assignement Operator.
 * \param[in] other UStructMesh object where copy from
 */
UStructMesh & UStructMesh::operator=(UStructMesh other){
    swap(other);
    return  *this;
};

/*!
 * Swap method, exchange data between the current class and another of the same type.
 * Data of one class will be assigned to the other and vice versa.
 *
 * \param[in] x UstructMesh object to be swapped.
 */
void UStructMesh::swap(UStructMesh & x) noexcept
{
        // Number of cells
        std::swap(m_nx, x.m_nx);
        std::swap(m_ny, x.m_ny);
        std::swap(m_nz, x.m_nz);

        // Cell Spacing
        std::swap(m_dx, x.m_dx);
        std::swap(m_dy, x.m_dy);
        std::swap(m_dz, x.m_dz);

        // Copy cell edges and cell centers ----------------------------------------- //
        std::swap(m_xnode, x.m_xnode);
        std::swap(m_ynode, x.m_ynode);
        std::swap(m_znode, x.m_znode);
        std::swap(m_xedge, x.m_xedge);
        std::swap(m_yedge, x.m_yedge);
        std::swap(m_zedge, x.m_zedge);

        std::swap(m_origin_temp,x.m_origin_temp);
        std::swap(m_span_temp,x.m_span_temp);
        std::swap(m_inflimits_temp,x.m_inflimits_temp);
        std::swap(m_refsystem_temp,x.m_refsystem_temp);
        std::swap(m_shapetype_temp,x.m_shapetype_temp);

        std::swap(m_setorigin,x.m_setorigin);
        std::swap(m_setspan,x.m_setspan);
        std::swap(m_setInfLimits,x.m_setInfLimits);
        std::swap(m_setRefSys,x.m_setRefSys);
        std::swap(m_isBuild,x.m_isBuild);

        if(x.m_shape){
            std::unique_ptr<BasicShape> tempshape;
            switch(x.m_shape->getShapeType()){
                case ShapeType::WEDGE :
                {
                    const Wedge * pp = dynamic_cast<const Wedge*>(x.m_shape.get());
                    if (pp != NULL)  tempshape = std::unique_ptr<BasicShape>(new Wedge(*(pp)));
                }
                break;
                case ShapeType::CYLINDER :
                {
                    const Cylinder * pp = dynamic_cast<const Cylinder*>(x.m_shape.get());
                    if (pp != NULL)  tempshape = std::unique_ptr<BasicShape>(new Cylinder(*(pp)));
                }
                break;
                case ShapeType::SPHERE :
                {
                    const Sphere * pp = dynamic_cast<const Sphere*>(x.m_shape.get());
                    if (pp != NULL)  tempshape = std::unique_ptr<BasicShape>(new Sphere(*(pp)));
                }
                break;
                default://CUBE
                {
                    const Cube * pp = dynamic_cast<const Cube*>(x.m_shape.get());
                    if (pp != NULL)  tempshape = std::unique_ptr<BasicShape>(new Cube(*(pp)));
                }
                break;
            }

            x.m_shape = std::move(m_shape);
            m_shape = std::move(tempshape);
        }
}


/*!
 * \return a const pointer to the inner BasicShape object the current mesh is built on. Const method
 */
const BasicShape * UStructMesh::getShape() const {
	return(m_shape.get());
}

/*!
 * \return a pointer to the inner BasicShape object the current mesh is built on
 */
BasicShape * UStructMesh::getShape(){
	return(m_shape.get());
}

/*!
 * \return current origin of BasicShape core of the mesh
 */
darray3E UStructMesh::getOrigin(){
	if (getShape() == NULL) return(m_origin_temp);
	return(getShape()->getOrigin());
}

/*!
 * \return current span of BasicShape core of the mesh
 */
darray3E UStructMesh::getSpan(){
	if (getShape() == NULL) return(m_span_temp);
	return(getShape()->getSpan());
}

/*!
 * \return current lower limits of coordinates in BasicShape core of the mesh
 */
darray3E UStructMesh::getInfLimits(){
	if (getShape() == NULL) return(m_inflimits_temp);
	return(getShape()->getInfLimits());
}

/*!
 *\return current local Reference System af axes
 */
dmatrix33E UStructMesh::getRefSystem(){
	if (getShape() == NULL) return(m_refsystem_temp);
	return(getShape()->getRefSystem());
}

/*!
 *\return actual scaling to primitive shape used in BasicShape core of the mesh
 */
darray3E UStructMesh::getScaling(){
	if (getShape() == NULL) return(darray3E{{1,1,1}});
	return(getShape()->getScaling());
}

/*!
 *\return local span of the primitive shape associated to BasicShape core of the mesh
 */
darray3E UStructMesh::getLocalSpan(){
	if (getShape() == NULL) return(darray3E{{1,1,1}});
	return(getShape()->getLocalSpan());
}

/*!
 *\return type of shape associated to mesh core. See ShapeType enum
 */
ShapeType UStructMesh::getShapeType(){
	if (getShape() == NULL) return(m_shapetype_temp);
	return(getShape()->getShapeType());
}

/*!
 * \return coordinate type of component 0 a BasicShape mesh core.
 * See CoordType enum.
 */
CoordType UStructMesh::getCoordTypex(){
	if (getShape() == NULL) return(CoordType::CLAMPED);
	return(getShape()->getCoordinateType(0));
}

/*! \return coordinate type of a component of a BasicShape mesh core.
 * See CoordType enum.
 * \param[in] i index of component.
 */
CoordType UStructMesh::getCoordType(int i){
	if (getShape() == NULL) return(CoordType::CLAMPED);
	return(getShape()->getCoordinateType(i));
}

/*!
 * \return coordinate type of component 1 a BasicShape mesh core.
 * See CoordType enum.
 */
CoordType UStructMesh::getCoordTypey(){
	if (getShape() == NULL) return(CoordType::CLAMPED);
	return(getShape()->getCoordinateType(1));
}

/*!
 * \return coordinate type of component 2 a BasicShape mesh core.
 * See CoordType enum.
 */
CoordType UStructMesh::getCoordTypez(){
	if (getShape() == NULL) return(CoordType::CLAMPED);
	return(getShape()->getCoordinateType(2));
}

/*!
 * \return coordinates type of a BasicShape mesh core. See CoordType enum.
 * \return Type of all the cooordinates.
 */
array<CoordType, 3> UStructMesh::getCoordType(){
	array<CoordType, 3> types;
	for (int i=0; i<3; i++){
		types[i] = getShape()->getCoordinateType(i);
	}
	return(types);
}

/*!
 * \return current mesh spacing in all three coordinates
 */
darray3E UStructMesh::getSpacing(){
	darray3E res;
	darray3E scale = getScaling();
	res[0] = m_dx*scale[0]; res[1] =m_dy*scale[1]; res[2] = m_dz*scale[2];
	return(res);
};

/*!
 * \return current dimensions of the mesh (number of mesh nodes in each direction)
 */
iarray3E UStructMesh::getDimension(){

	iarray3E res;
	res[0] = m_nx + 1;
	res[1] = m_ny + 1;
	res[2] = m_nz + 1;
	return(res);
};

/*!
 * \return the n-th center-cell coordinates in local mesh reference frame,
 * given its global cell index on the mesh.
 * \param[in] index cell index in the global nodal list.
 */
darray3E UStructMesh::getLocalCCell(int index){

	int i,j,k;
	accessCellIndex(index, i,j,k);
	darray3E res = getLocalCCell(i,j,k);

	return(res);
};

/*!
 * \return the n-th center-cell coordinates in local mesh reference frame,
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

/*!
 * \return teh n-th nodal vertex coordinates in local mesh reference frame,
 * given its global point index on the mesh.
 * \param[in] index point index in the global nodal list.
 */
darray3E UStructMesh::getLocalPoint(int index){

	int i,j,k;
	accessPointIndex(index, i,j,k);
	darray3E res = getLocalPoint(i,j,k);

	return(res);
};
/*!
 * \return the n-th nodal vertex coordinates in local reference frame,
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

/*!
 * \return the n-th center-cell coordinates in global absolute reference frame,
 * given its global cell index on the mesh.
 * \param[in] index cell index in the global nodal list.
 */
darray3E UStructMesh::getGlobalCCell(int index){

	darray3E res = getLocalCCell(index);
	return(transfToGlobal(res));
};
/*!
 * \return the n-th center cell coordinates in global absolute reference frame,
 * given its cartesian cell indices on the mesh.
 * \param[in] i_ x cartesian index.
 * \param[in] j_ y cartesian index.
 * \param[in] k_ z cartesian index.
 */
darray3E UStructMesh::getGlobalCCell(int i_, int j_, int k_){

	darray3E res = getLocalCCell(i_,j_,k_);
	return(transfToGlobal(res));
};

/*!
 * \return the n-th nodal vertex coordinates in global absolute reference frame,
 * given its global point index on the mesh.
 * \param[in] index point index in the global nodal list.
 */
darray3E UStructMesh::getGlobalPoint(int index){
	darray3E res = getLocalPoint(index);
	return(transfToGlobal(res));
};
/*!
 * \return the n-th nodal vertex coordinates in global absolute reference frame,
 * given its cartesian point indices on the mesh.
 * \param[in] i_ x cartesian index.
 * \param[in] j_ y cartesian index.
 * \param[in] k_ z cartesian index.
 */
darray3E UStructMesh::getGlobalPoint(int i_, int j_, int k_){

	darray3E res = getLocalPoint(i_,j_,k_);
	return(transfToGlobal(res));
};

/*!
 * \return the list of neighbor vertices (by their global indices and in VTK-hexahedron ordering)
 * of a given cell
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

/*!
 * \return the list of neighbor vertices (by their global indices and in VTK-hexahedron ordering)
 * of a given cell
 * \param[in] index cell global index
 */
ivector1D UStructMesh::getCellNeighs(int index){

	ivector1D pp(3,0);
	accessCellIndex(index, pp[0],pp[1], pp[2]);
	return(getCellNeighs(pp[0],pp[1],pp[2]));
}

/*!
 *\return the complete list of mesh nodes, in local shape reference system
 */
dvecarr3E UStructMesh::getLocalCoords(){

	int np = (m_nx+1)*(m_ny+1)*(m_nz+1);
	dvecarr3E coords(np);
	for (int i=0; i<np; i++){
		coords[i] = getLocalPoint(i);
	}
	return coords;
};

/*!
 * \return teh complete list of mesh nodes, in global absolute reference system
 */
dvecarr3E UStructMesh::getGlobalCoords(){
	int np = (m_nx+1)*(m_ny+1)*(m_nz+1);
	dvecarr3E coords(np);
	for (int i=0; i<np; i++){
		coords[i] = getGlobalPoint(i);
	}
	return coords;
};

/*!
 * \return complete list of mesh cell centers, in local shape reference system
 */
dvecarr3E UStructMesh::getLocalCellCentroids(){

	int np = (m_nx)*(m_ny)*(m_nz);
	dvecarr3E coords(np);
	for (int i=0; i<np; i++){
		coords[i] = getLocalCCell(i);
	}
	return coords;
};

/*!
 * \return complete list of mesh cell centers, in global absolute reference system
 */
dvecarr3E UStructMesh::getGlobalCellCentroids(){
	int np = (m_nx)*(m_ny)*(m_nz);
	dvecarr3E coords(np);
	for (int i=0; i<np; i++){
		coords[i] = getGlobalCCell(i);
	}
	return coords;
};

/*!
 * Set the origin of your shape. The origin is meant as the baricenter of your shape in absolute r.s.
 * Info is just passed and stored in memory, but no modifications are applied to your current mesh.
 * To apply current modifications use UStructMesh::execute()/build() method.
 * \param[in] origin new origin point
 */
void UStructMesh::setOrigin(darray3E origin){
	if (getShape() == NULL){
		m_origin_temp = origin;
		m_setorigin = true;
	}else{
		getShape()->setOrigin(origin);
	}
}

/*!
 * Set the span of your shape, according to its local reference system.
 * Info is just passed and stored in memory, but no modifications are applied to your current mesh.
 * To apply current modifications use UStructMesh::execute()/build() method.
 * \param[in] s0 first coordinate span
 * \param[in] s1 second coordinate span
 * \param[in] s2 third coordinate span
 */
void UStructMesh::setSpan(double s0, double s1, double s2){

	if (getShape() == NULL){
		m_span_temp[0] = s0;
		m_span_temp[1] = s1;
		m_span_temp[2] = s2;
		m_setspan = true;
	}else{
		getShape()->setSpan( s0, s1,s2);
	}
	m_isBuild = false;
}

/*!
 * Set the span of your shape, according to its local reference system.
 *  Info is just passed and stored in memory, but no modifications are applied to your current mesh.
 *  To apply current modifications use UStructMesh::execute()/build() method.
 * \param[in] s coordinates span
 */
void UStructMesh::setSpan(darray3E s){
	setSpan( s[0], s[1], s[2]);
}

/*!
 * Set inferior limits of your shape, according to its local reference system.
 *  Info is just passed and stored in memory, but no modifications are a*pplied to your current mesh.
 *  To apply current modifications use UStructMesh::execute()/build() method.
 * \param[in] inflim coordinate inferior limit
 * \param[in] dir 0,1,2 int flag identifying coordinate
 */
void UStructMesh::setInfLimits(double inflim, int dir){
	if (getShape() == NULL){
		m_inflimits_temp[dir] = inflim;
		m_setInfLimits = true;
	}else{
		getShape()->setInfLimits(inflim, dir);
	}
}

/*!
 * Set coordinates' inferior limits of your shape, according to its local reference system.
 *  Info is just passed and stored in memory, but no modifications are applied to your current mesh.
 *  To apply current modifications use UStructMesh::execute()/build() method.
 * \param[in] inflim coordinates inferior limits
 */
void UStructMesh::setInfLimits(darray3E inflim){
	setInfLimits(inflim[0], 0);
	setInfLimits(inflim[1], 1);
	setInfLimits(inflim[2], 2);
}

/*!
 * Set new axis orientation of the local reference system of your mesh core shape.
 *  Info is just passed and stored in memory, but no modifications are applied to your current mesh.
 *  To apply current modifications use UStructMesh::execute()/build() method.
 * \param[in] axis0 first axis
 * \param[in] axis1 second axis
 * \param[in] axis2 third axis
 *
 * if chosen axes are not orthogonal, doing nothing
 */
void UStructMesh::setRefSystem(darray3E axis0, darray3E axis1, darray3E axis2){

	if (getShape() == NULL){
		m_refsystem_temp[0] = axis0;
		m_refsystem_temp[1] = axis1;
		m_refsystem_temp[2] = axis2;
		m_setRefSys = true;
	}else{
		getShape()->setRefSystem(axis0, axis1, axis2);
	}
}

/*!
 * Set new axis orientation of the local reference system of your mesh core shape
 *  Info is just passed and stored in memory, but no modifications are applied to your current mesh.
 *  To apply current modifications use UStructMesh::execute()/build() method.
 * \param[in] label 0,1,2 identify local x,y,z axis of the primitive shape
 * \param[in] axis new direction of selected local axis.
 */
void UStructMesh::setRefSystem(int label, darray3E axis){

	if(getShape() == NULL){
		BasicShape * temp = new Cube();
		temp->setRefSystem(label, axis);
		m_refsystem_temp = temp->getRefSystem();
		delete temp;
		temp = NULL;
		m_setRefSys = true;
	}else{
		getShape()->setRefSystem(label,axis);
	}
}

/*!
 * Set new axis orientation of the local reference system of your mesh core shape
 *  Info is just passed and stored in memory, but no modifications are applied to your current mesh.
 *  To apply current modifications use UStructMesh::execute()/build() method.
 * \param[in] axes new direction of all local axes.
 */
void UStructMesh::setRefSystem(dmatrix33E axes){
	if (getShape() == NULL){

    	m_refsystem_temp[0] = axes[0];
		m_refsystem_temp[1] = axes[1];
		m_refsystem_temp[2] = axes[2];
		m_setRefSys = true;
	}else{
        getShape()->setRefSystem(axes);
	}
}

/*!
 * Set the dimensions of the mesh.
 *  Info is just passed and stored in memory, but no modifications are applied to your current mesh.
 *  To apply current modifications use UStructMesh::execute()/build() method.
 *  \param[in] dim  number of mesh nodes in each direction
 */
 void UStructMesh::setDimension(ivector1D dim){
	 m_nx = dim[0] - 1;
	 m_ny = dim[1] - 1;
	 m_nz = dim[2] - 1;

	 m_isBuild = false;
};

/*!
 * Set the dimensions of the mesh.
 *  Info is just passed and stored in memory, but no modifications are applied to your current mesh.
 *  To apply current modifications use UStructMesh::execute()/build() method.
 *  \param[in] dim  number of mesh nodes in each direction
 */
  void UStructMesh::setDimension(iarray3E dim){
 	 m_nx = dim[0] - 1;
 	 m_ny = dim[1] - 1;
 	 m_nz = dim[2] - 1;

	 m_isBuild = false;
 };

  /*!
   * Set your shape, according to the following input parameters and the already saved/default parmaters.
   * Mesh is still not build. use UStructMesh::execute()/build() to build the mesh.
   * \param[in] itype shape of your mesh, casted to enum.(option available are: 0-CUBE(default), 1-CYLINDER, 2-SPHERE)
   */
  void UStructMesh::setShape(int itype){
	  ShapeType type = static_cast<ShapeType>(itype);
	  UStructMesh::setShape(type);
  }

  /*!
   * Set your shape, according to the following input parameters and the already saved/default parmaters.
   * Mesh is still not build. use UStructMesh::execute()/build() to build the mesh.
   * \param[in] type shape of your mesh, based on ShapeType enum.(option available are: CUBE(default), CYLINDER, SPHERE)
   */
  void UStructMesh::setShape(ShapeType type){
	  //create internal shape using unique_ptr member.
	  // unlink external shape eventually
		if(m_shape){

			m_origin_temp = m_shape.get()->getOrigin();
			m_span_temp = m_shape.get()->getSpan();
			m_inflimits_temp = m_shape.get()->getInfLimits();
			m_refsystem_temp = m_shape.get()->getRefSystem();
			m_setorigin = true;
			m_setspan = true;
			m_setInfLimits = true;
			m_setRefSys = true;

			m_shape.reset(nullptr);
		}

		darray3E origin{{0,0,0}};
		dmatrix33E spanMat;
		spanMat[0].fill(1);	spanMat[2].fill(1);	spanMat[2].fill(1);
		spanMat[1][1] = spanMat[2][1] = 2*M_PI;
		spanMat[2][2] = M_PI;

		if(m_setorigin){origin = m_origin_temp; m_setorigin= false;}
		if(m_setspan){
			for(int i=0; i<3; ++i) spanMat[i] = m_span_temp;
			m_setspan = false;
		}

  		switch(type){
            case ShapeType::WEDGE :
                m_shape = std::unique_ptr<BasicShape>(new Wedge(origin,spanMat[0]));
                break;
            case ShapeType::CYLINDER :
                m_shape = std::unique_ptr<BasicShape>(new Cylinder(origin,spanMat[1]));
                break;
            case ShapeType::SPHERE :
                m_shape = std::unique_ptr<BasicShape>(new Sphere(origin, spanMat[2]));
                break;
            default://CUBE
                m_shape = std::unique_ptr<BasicShape>(new Cube(origin, spanMat[0]));
                break;
		}

		if(m_setInfLimits){
			m_shape.get()->setInfLimits(m_inflimits_temp[0],0);
			m_shape.get()->setInfLimits(m_inflimits_temp[1],1);
			m_shape.get()->setInfLimits(m_inflimits_temp[2],2);
			m_setInfLimits = false;
		}

		if(m_setRefSys){
			m_shape.get()->setRefSystem(m_refsystem_temp);
			m_setRefSys = false;
		}

		m_isBuild = false;
}

/*!
 * Set mesh shape, copying an external BasicShape object .
 * Mesh is still not build. use UStructMesh::execute()/build() to build the mesh.
 * \param[in] shape pointer to an external allocated BasicShape object
 */
void UStructMesh::setShape(const BasicShape * shape){

	if(shape == NULL) return;
	m_shape.reset(nullptr);
	m_setorigin = false;
	m_setspan = false;
	m_setInfLimits = false;
	m_setRefSys = false;

	switch(shape->getShapeType()){
        case ShapeType::WEDGE :
        {
            const Wedge * pp = dynamic_cast<const Wedge*>(shape);
            if (pp != NULL)  m_shape = std::unique_ptr<BasicShape>(new Wedge(*(pp)));
        }
        break;
        case ShapeType::CYLINDER :
            {
                const Cylinder * pp = dynamic_cast<const Cylinder*>(shape);
                if (pp != NULL)  m_shape = std::unique_ptr<BasicShape>(new Cylinder(*(pp)));
            }
            break;
		case ShapeType::SPHERE :
            {
                const Sphere * pp = dynamic_cast<const Sphere*>(shape);
                if (pp != NULL)  m_shape = std::unique_ptr<BasicShape>(new Sphere(*(pp)));
            }
            break;
		default://CUBE
            {
                const Cube * pp = dynamic_cast<const Cube*>(shape);
                if (pp != NULL)  m_shape = std::unique_ptr<BasicShape>(new Cube(*(pp)));
            }
            break;
	}

	m_isBuild = false;
};


/*!
 * Set coordinate type of component 0 of a BasicShape mesh core.
 * See CoordType enum.
 */
void UStructMesh::setCoordTypex(CoordType type){
	if (getShape() == NULL) return;
	getShape()->setCoordinateType(type,0);
}

/*!
 * Set coordinate type of a component of a BasicShape mesh core.
 * See CoordType enum.
 * \param[in] type CoordType enum input
 * \param[in] i index of component.
 */
void UStructMesh::setCoordType(CoordType type, int i){
	if (getShape() == NULL) return;
	getShape()->setCoordinateType(type,i);
}

/*!
 * Set coordinate type of component 1 a BasicShape mesh core.
 * See CoordType enum.
 * * \param[in] type CoordType enum input
 */
void UStructMesh::setCoordTypey(CoordType type){
	if (getShape() == NULL) return;
	getShape()->setCoordinateType(type,1);
}

/*!
 * Set coordinate type of component 2 a BasicShape mesh core.
 * See CoordType enum.
 * * \param[in] type CoordType enum input
 */
void UStructMesh::setCoordTypez(CoordType type){
	if (getShape() == NULL) return;
	getShape()->setCoordinateType(type,2);
}

/*!
 * Set coordinates type of a BasicShape mesh core. See CoordType enum.
 * \param[in] types array of CoordType enum for all the cooordinates.
 */
void UStructMesh::setCoordType(std::array<CoordType, 3> types){
	for (int i=0; i<3; i++){
		getShape()->setCoordinateType(types[i],i);
	}
}

/*!
 * Set your mesh, according to the following input parameters
 * \param[in] origin 3D point baricenter of your mesh
 * \param[in] span span for each coordinate defining your mesh
 * \param[in] type   shape of your mesh, based on ShapeType enum.(option available are: CUBE(default), CYLINDER, SPHERE, WEDGE)
 * \param[in] dimensions number of mesh points for each coordinate.
 */
void UStructMesh::setMesh(darray3E & origin, darray3E &span, ShapeType type, iarray3E & dimensions){

	if(m_shape){m_shape.reset(nullptr);}

	setShape(type);
	setOrigin(origin);
	setSpan(span);
	setDimension(dimensions);
	build();
};


/*!
 * Set your mesh, according to the following input parameters
 * \param[in] origin 3D point baricenter of your mesh
 * \param[in] span span for each coordinate defining your mesh
 * \param[in] type   shape of your mesh, based on ShapeType enum.(option available are: CUBE(default), CYLINDER, SPHERE,WEDGE)
 * \param[in] spacing fixed spacing for each coordinate
 */
void UStructMesh::setMesh(darray3E & origin, darray3E &span, ShapeType type, dvector1D & spacing){

	ivector1D dimLimit(3,2);
	//create internal shape using unique_ptr member.
	if(m_shape){m_shape.reset(nullptr);}

	switch(type){
		case ShapeType::CYLINDER :
			dimLimit[1] = 4;
			break;
		case ShapeType::SPHERE :
			dimLimit[1] = 5; dimLimit[2] = 3;
			break;
		default://CUBE //WEDGE
			break;
	}

	setShape(type);
	setOrigin(origin);
	setSpan(span);

	darray3E span2 = getSpan();
	iarray3E dim;

	for(int i=0; i<3; ++i){
		if(spacing[i] != 0.0) {
			dim[i] = (int) std::floor(span2[i]/spacing[i] +0.5) + 1;
		}else{
			dim[i] = dimLimit[i];
		}
	}

	setDimension(dim);
	build();
};

/*!
 * Set your mesh, according to the following input parameters
 * \param[in] shape pointer to an external allocated BasicShape object
 * \param[in] dimensions number of mesh points for each coordinate.
 */
void UStructMesh::setMesh(BasicShape * shape, iarray3E & dimensions){

	if(m_shape){m_shape.reset(nullptr);}

	setShape(shape);
	setDimension(dimensions);
	build();
};

/*!
 * Set your mesh, according to the following input parameters
 * \param[in] shape pointer to an external allocated BasicShape object
 * \param[in] spacing fixed spacing for each coordinate
 */
void UStructMesh::setMesh(BasicShape * shape, dvector1D & spacing){

	ivector1D dimLimit(3,2);
	//create internal shape using unique_ptr member.
	if(m_shape){m_shape.reset(nullptr);}

	switch(shape->getShapeType()){
		case ShapeType::CYLINDER :
			dimLimit[1] = 4;
			break;
		case ShapeType::SPHERE :
			dimLimit[1] = 5; dimLimit[2] = 3;
			break;
		default://CUBE //WEDGE
			break;
	}

	setShape(shape);

	darray3E span2 = getSpan();
	iarray3E dim;

	for(int i=0; i<3; ++i){
		if(spacing[i] != 0.0) {
			dim[i] = (int) std::floor(span2[i]/spacing[i] +0.5) + 1;
		}else{
			dim[i] = dimLimit[i];
		}
	}

	setDimension(dim);
	build();
};


/*!
 *Clear the Mesh structure. Unlink external shapes
 * or destroy internal shapes. Destroy nodal structure.
 */
void UStructMesh::clearMesh(){

	m_shape.reset(nullptr);
	m_nx=0; m_ny=0; m_nz=0;
	m_dx=0.0; m_dy=0.0; m_dz=0.0;


	m_setorigin = false;
	m_setspan = false;
	m_setInfLimits = false;
	m_setRefSys = false;
	m_isBuild = false;

	m_origin_temp = {{0.0,0.0,0.0}};
	m_span_temp = {{1.0,1.0,1.0}};
	m_inflimits_temp = {{0.0,0.0,0.0}};
	for(int i=0; i<3; ++i){m_refsystem_temp[i].fill(0.0); m_refsystem_temp[i][i] = 1.0;}
	m_shapetype_temp = ShapeType::CUBE;

	destroyNodalStructure();
};


/*!
 * \return cartesian indices of the cell containing
 * the target point in global reference frame.
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

/*!
 * \return cartesian indices of the cell containing the target point
 * in global reference framne
 * \param[in] point 3D coordinate of target point
 * \param[out] i x cell index
 * \param[out] j y cell index
 * \param[out] k z cell index
 */
void UStructMesh::locateCellByPoint(dvector1D & point, int &i, int &j, int &k){
		darray3E temp = conArray<double,3>(point);
		locateCellByPoint(temp,i,j,k);
};

/*!
 * \return global index of the cell given its cartesian indices.
 * Follows the ordering sequences z-y-x
 * \param[in] i x cartesian index
 *\param[in] j y cartesian index
 *\param[in] k z cartesian index
 *\return global index
 */
int UStructMesh::accessCellIndex(int i, int j, int k){

	int index = m_ny * m_nz * i + m_nz * j + k;
	return(index);
};

/*!
 *\return cartesian indices of the cell given its global index.
 * Follows the ordering sequences z-y-x
 *\param[in] N_ global index
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

/*!
 * \return cartesian indices of the point given its global index.
 * Follows the ordering sequences z-y-x
 * \param[in] N_ global index
 * \param[out] i x cartesian index
 * \param[out] j y cartesian index
 * \param[out] k z cartesian index
 */
void UStructMesh::accessPointIndex(int N_, int &i, int &j, int &k){

	k = N_ % (m_nz+1);
	int index = N_ / (m_nz+1);
	j = index % (m_ny+1);
	i = index / (m_ny+1);
};

/*!
 * \return transformed point from Local to Global reference system
 * \param[in] point in local reference system
 */
darray3E	UStructMesh::transfToGlobal( darray3E & point){
	return(getShape()->toWorldCoord(point));
};

/*!
 * \return transformed point from Local to Global reference system
 * \param[in] point in local reference system
 */
dvector1D	UStructMesh::transfToGlobal( dvector1D & point){
	darray3E temp = conArray<double,3>(point);
	darray3E temp2 = getShape()->toWorldCoord(temp);
	return(conVect(temp2));
};

/*!
 * \return list of transformed point from Local to Global reference system
 * \param[in] list_points in local reference system
 */
dvecarr3E	UStructMesh::transfToGlobal( dvecarr3E & list_points){
	int size = list_points.size();
	dvecarr3E result(size);
	for(int i=0; i<size; ++i){
		result[i] = transfToGlobal(list_points[i]);
	}
	return(result);
};

/*!
 * \return transformed point from Global to Local reference system
 * \param[in] point in global reference system
 */
darray3E 	UStructMesh::transfToLocal( darray3E & point){
	return(getShape()->toLocalCoord(point));
};

/*! Transform point from Global to Local ref system*/
dvector1D 	UStructMesh::transfToLocal( dvector1D & point){
	darray3E temp = conArray<double,3>(point);
	darray3E temp2 = getShape()->toLocalCoord(temp);
	return(conVect(temp2));

};

/*!
 * \return list of transformed point from Global to Local reference system
 * \param[in] list_points in global reference system
 */
dvecarr3E 	UStructMesh::transfToLocal( dvecarr3E & list_points){
	int size = list_points.size();
	dvecarr3E result(size);
	for(int i=0; i<size; ++i){
		result[i] = transfToLocal(list_points[i]);
	}
	return(result);

};

/*!
 * \return interpolated value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on centercells
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

/*!
 * \return interpolated value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on centercells
 *
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

/*!
 * \return interpolated value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on centercells
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

/*!
 * \return interpolated value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] pointdata data field defined on grid points
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

/*!
 * \return interpolated value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] pointdata data field defined on grid points
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

/*!
 * \return interpolated value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] pointdata data field defined on grid points
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


/*!
 * Write your grid as a point cloud in a *.vtu file.
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

	iarray3E dim = getDimension();
	int sizeTot = dim[0]*dim[1]*dim[2];

	dvecarr3E activeP;
	if(extPoints != NULL && (int)extPoints->size() == sizeTot){
		activeP = *extPoints;
	}else{
		activeP.resize(sizeTot);
		for(int i=0; i<sizeTot; i++){
			activeP[i] = getGlobalPoint(i);
		}
	}

	ivector1D conn(activeP.size());
	{
		int counter = 0;
		for(auto & val: conn){
			val = counter;
			++counter;
		}
	}
	bitpit::VTKUnstructuredGrid vtk(folder, outfile, bitpit::VTKElementType::VERTEX);
    vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, activeP) ;
    vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, conn) ;
	vtk.setDimensions(conn.size(), activeP.size());
	vtk.setCodex(codex);
	if(counterFile>=0){vtk.setCounter(counterFile);}

	vtk.write();

};

/*!
 * Write your grid as a point cloud in a *.vtu file.
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

	iarray3E dim = getDimension();
	int sizeTot = dim[0]*dim[1]*dim[2];
	int sizeMap = std::min((int)vertexList.size(), sizeTot);

	dvecarr3E activeP(sizeMap);
	if(extPoints != NULL && (int)extPoints->size() == sizeTot){
		for(int i=0; i<sizeMap; i++){
			activeP[i] = (*extPoints)[vertexList[i]];
		}
	}
	else{
		for(int i=0; i<sizeMap; i++){
			activeP[i] = getGlobalPoint(vertexList[i]);
		}
	}

	ivector1D conn(activeP.size());
	{
		int counter = 0;
		for(auto & val: conn){
			val = counter;
			++counter;
		}
	}

	bitpit::VTKUnstructuredGrid vtk(folder, outfile, bitpit::VTKElementType::VERTEX);
    vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, activeP) ;
    vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, conn) ;
	vtk.setDimensions(conn.size(), activeP.size());
	vtk.setCodex(codex);
	if(counterFile>=0){vtk.setCounter(counterFile);}

	vtk.write();

};

/*!
 * Write your grid as a hexahedrical one in a *.vtu file.
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

	iarray3E dim = getDimension();
	int sizePt = dim[0]*dim[1]*dim[2];
	int sizeCl = (dim[0]-1)*(dim[1]-1)*(dim[2]-1);

	dvecarr3E activeP(sizePt);
	ivector2D activeConn(sizeCl, ivector1D(8,0));

	if(extPoints != NULL && (int)extPoints->size() == sizePt){
		activeP = *extPoints;
	}
	else{
		for(int i=0; i<sizePt; i++){
			activeP[i] = getGlobalPoint(i);
		}
	}

	for(int i=0; i<sizeCl; ++i){
		activeConn[i] = getCellNeighs(i);
	}

	bitpit::VTKElementType elDM = bitpit::VTKElementType::HEXAHEDRON;
	bitpit::VTKUnstructuredGrid vtk(folder, outfile, elDM);
    vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, activeP) ;
    vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, activeConn) ;
	vtk.setDimensions(sizeCl, sizePt);
	vtk.setCodex(codex);
	if(counterFile>=0){vtk.setCounter(counterFile);}

	vtk.write();

};

/*!
 * Write your grid as a hexahedrical one in a *.vtu file.
 * \param[in] folder output directory path
 * \param[in] outfile output namefile w/out tag
 * \param[in] counterFile integer to mark output with a counter number
 * \param[in] codexFlag boolean to distinguish between "ascii" false, and "appended" true
 * \param[in] cellList list of global indices of selected cells that the user wants to write on file.
 * \param[in] extPoints OPTIONAL. If defined, use the current vertexList and write them as point cloud, provided that coherent dimensions with the mesh are set.
 *			If Not, write the current mesh vertices.
 * */
void UStructMesh::plotGrid(std::string & folder, std::string outfile, int counterFile, bool codexFlag, ivector1D & cellList, dvecarr3E * extPoints){

	bitpit::VTKFormat codex = bitpit::VTKFormat::ASCII;
    if(codexFlag){codex=bitpit::VTKFormat::APPENDED;}

	iarray3E dim = getDimension();
	int sizePt = dim[0]*dim[1]*dim[2];
	int sizeCl = cellList.size();

	if(sizeCl > (dim[0]-1)*(dim[1]-1)*(dim[2]-1)){return;}

	std::map<int,int> mapPoints;
	ivector2D activeConn(sizeCl, ivector1D(8,0));

	for(int i=0; i<sizeCl; ++i){
		activeConn[i] = getCellNeighs(i);
		for(int j=0; j<(int)activeConn[i].size(); ++j){
			mapPoints[activeConn[i][j]] = activeConn[i][j];
		}
	}

	dvecarr3E activeP(mapPoints.size());
	ivector1D listP(mapPoints.size());
	int counter = 0;
	std::map<int,int>::iterator itF;
	if(extPoints != NULL && (int)extPoints->size() == sizePt){
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
		for(int j=0; j<(int)activeConn[i].size(); ++j){
			int posV =posVectorFind(listP, activeConn[i][j]);
			//update local connectivity
			activeConn[i][j] = posV;
		}
	}

	bitpit::VTKElementType elDM = bitpit::VTKElementType::HEXAHEDRON;
	bitpit::VTKUnstructuredGrid vtk(folder, outfile, elDM);
    vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, activeP) ;
    vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, activeConn) ;
	vtk.setDimensions(sizeCl, sizePt);
	vtk.setCodex(codex);
	if(counterFile>=0){vtk.setCounter(counterFile);}

	vtk.write();

};

/*!
 * Write your grid as a hexahedrical one in a *.vtu file.
 * \param[in] folder output directory path
 * \param[in] outfile output namefile w/out tag
 * \param[in] counterFile integer to mark output with a counter number
 * \param[in] codexFlag boolean to distinguish between "ascii" false, and "appended" true
 * \param[in] data Scalar field
 * */
void UStructMesh::plotGridScalar(std::string folder, std::string outfile , int counterFile, bool codexFlag, dvector1D & data){

	bitpit::VTKFormat codex = bitpit::VTKFormat::ASCII;
	if(codexFlag){codex=bitpit::VTKFormat::APPENDED;}

	iarray3E dim = getDimension();
	int sizePt = dim[0]*dim[1]*dim[2];
	int sizeCl = (dim[0]-1)*(dim[1]-1)*(dim[2]-1);

	dvecarr3E activeP(sizePt);
	ivector2D activeConn(sizeCl, ivector1D(8,0));

	for(int i=0; i<sizePt; i++){
		activeP[i] = getGlobalPoint(i);
	}

	for(int i=0; i<sizeCl; ++i){
		activeConn[i] = getCellNeighs(i);
	}

	bitpit::VTKElementType  elDM = bitpit::VTKElementType::HEXAHEDRON;
	bitpit::VTKUnstructuredGrid vtk(folder, outfile, elDM);
    vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, activeP) ;
    vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, activeConn) ;
	vtk.setDimensions(sizeCl, sizePt);

	vtk.setCodex(codex);
	if(counterFile>=0){vtk.setCounter(counterFile);}

	bitpit::VTKLocation loc = bitpit::VTKLocation::POINT;
	if ((int)data.size() == sizeCl) loc = bitpit::VTKLocation::CELL;

	vtk.addData("field",bitpit::VTKFieldType::SCALAR, loc, data) ;
	vtk.write();

};


/*!
 * Destroy all nodal structures of the mesh.
 */
void UStructMesh::destroyNodalStructure(){
	m_xnode.clear();
	m_ynode.clear();
	m_znode.clear();
	m_xedge.clear();
	m_yedge.clear();
	m_zedge.clear();
};

/*!
 * Destroy all nodal structures of the mesh,
 * and reinitialize them to current mesh dimensions.
 */
void UStructMesh::reshapeNodalStructure(){
	destroyNodalStructure();
	resizeMesh();
};

/*! Resize nodal structures to current
 * dimensions set
 */
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

/*!
 * Apply stored mesh information and build the mesh. If no core shape is
 *  build, does nothing and exit.
 */
bool UStructMesh::build(){

	if(getShape() == NULL){return false;}

	ivector1D dimLimit(3,2);
	//create internal shape using unique_ptr member.
	// unlink external shape eventually
	switch(getShapeType()){
		case ShapeType::CYLINDER :
			dimLimit[1] = 4;
			break;
		case ShapeType::SPHERE :
			dimLimit[1] = 5; dimLimit[2] = 3;
			break;
		default://CUBE //WEDGE
			break;
	}

	//check on dimensions and eventual closed loops on coordinates.
	m_nx = std::max(m_nx, dimLimit[0]-1);
	m_ny = std::max(m_ny, dimLimit[1]-1);
	m_nz = std::max(m_nz, dimLimit[2]-1);

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

	m_isBuild = true;
    return true;
};


/*!
 * Execute the object, that is build the mesh with current stored parameters.
 */
void UStructMesh::execute(){
	build();
}

/*!
 * \return true if mesh is built according to the currently set parameters,
 */
bool UStructMesh::isBuilt(){
	return(m_isBuild);
}

};
