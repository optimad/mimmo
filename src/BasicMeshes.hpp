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

#ifndef __MiMMO_BASICMESHES_HH
#define __MiMMO_BASICMESHES_HH

// local libs
#include <memory>
#include <array>
#include "BasicShapes.hpp"
#include "MiMMO_TypeDef.hpp"

/*!
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

class UStructMesh{

private:
	BasicShape * 				m_shape1; 	/**!< generic pointer to BasicShape core of the mesh, for External USE*/
	std::unique_ptr<BasicShape>	m_shape2;	/**!< unique pointer to BasicShape core of the mesh, for Internal USE*/
							
protected:
	double				m_dx, m_dy, m_dz;	/**< Mesh spacing in each direction */
	int					m_nx, m_ny, m_nz;	/**< Mesh number of cells in each direction */
	dvector1D 	m_xnode, m_ynode, m_znode; 	/**< Lists holding the center cells coordinates of the mesh, in local reference sistem */
	dvector1D 	m_xedge, m_yedge, m_zedge; 	/**< Lists holding the point coordinates of the mesh, in local reference system */
	//TODO m_setmesh useful?
	bool						m_setmesh; 	/**< check true if mesh is actually set */

public:	     
	//Building stuffs	    
	UStructMesh();
//	UStructMesh(darray3E & origin, darray3E & span, BasicShape::ShapeType, ivector1D & dimensions);
//	UStructMesh(darray3E & origin, darray3E & span, BasicShape::ShapeType, dvector1D & spacing);
//	UStructMesh(BasicShape *, ivector1D & dimensions);
//	UStructMesh(BasicShape *, dvector1D & spacing);
	virtual ~UStructMesh();
	
	// Copy constructor & assignment operators
	UStructMesh(const UStructMesh & other);
	UStructMesh & operator=(const UStructMesh & other);
	
	
	//get-set methods  
	const BasicShape*		getShape() const;
	BasicShape*				getShape();
	darray3E				getOrigin();
	darray3E				getSpan();
	darray3E				getInfLimits();
	dmatrix33E				getRefSystem();
	darray3E				getScaling();
	darray3E				getLocalSpan();
	BasicShape::ShapeType	getShapeType();
	BasicShape::CoordType	getCoordType(int);
	BasicShape::CoordType	getCoordTypex();
	BasicShape::CoordType	getCoordTypey();
	BasicShape::CoordType	getCoordTypez();

	std::array<BasicShape::CoordType,3>	getCoordType();
	
	darray3E 				getSpacing();
	ivector1D				getDimension();

	//TODO Useful for pin in/out interface? (Now they can't be used with pins)
	darray3E 				getLocalCCell(int);
	darray3E 				getLocalCCell(int, int, int);
	darray3E 				getLocalPoint(int);
	darray3E 				getLocalPoint(int, int, int);
	darray3E 				getGlobalCCell(int);
	darray3E 				getGlobalCCell(int, int, int);
	darray3E 				getGlobalPoint(int);
	darray3E 				getGlobalPoint(int, int, int);
	
	ivector1D 				getCellNeighs(int);
	ivector1D 				getCellNeighs(int, int, int);

	//TODO each set method does a rebuild of the mesh, can be found an alternative strategy?
	void	setOrigin(darray3E origin);
	void	setSpan(double, double, double, bool flag = true);
	void	setSpan(darray3E span);
	void	setInfLimits(double val, int dir, bool flag = true);
	void	setInfLimits(darray3E val);

	void	setRefSystem(darray3E, darray3E, darray3E);
	void	setRefSystem(int, darray3E);
	void	setRefSystem(dmatrix33E);

	void	setDimension(ivector1D dim);
	void	setDimension(iarray3E dim);

	void 	setShape(BasicShape::ShapeType type);
	void 	setShape(BasicShape *);

	void 	setMesh(darray3E & origin, darray3E & span, BasicShape::ShapeType, ivector1D & dimensions);
	void 	setMesh(darray3E & origin, darray3E & span, BasicShape::ShapeType, dvector1D & spacing);
	void 	setMesh(BasicShape *, ivector1D & dimensions);
	void 	setMesh(BasicShape *, dvector1D & spacing);

	//generic manteinance of the mesh
	void	clearMesh();
	bool	isBuilt();
	
	//functionalities
	void 	locateCellByPoint(darray3E & point, int &i, int &j, int &k);
	void 	locateCellByPoint(dvector1D & point, int &i, int &j, int &k);
	int  	accessCellIndex(int i, int j, int k);
	void 	accessCellIndex(int N_, int &i, int &j, int &k);
	int  	accessPointIndex(int i, int j, int k);
	void 	accessPointIndex(int N_, int &i, int &j, int &k);
	
	darray3E	transfToGlobal( darray3E & point);
	dvector1D	transfToGlobal( dvector1D & point);
	dvecarr3E	transfToGlobal( dvecarr3E & list_points);    
	darray3E 	transfToLocal( darray3E & point);
	dvector1D 	transfToLocal( dvector1D & point);
	dvecarr3E 	transfToLocal( dvecarr3E & list_points);    
	
	// interpolators  
	double 		interpolateCellData(darray3E & point, dvector1D & celldata);
	int 		interpolateCellData(darray3E & point, ivector1D & celldata);
	darray3E	interpolateCellData(darray3E & point, dvecarr3E & celldata);
	
	double		interpolatePointData(darray3E & point, dvector1D & pointdata);
	int 		interpolatePointData(darray3E & point, ivector1D & pointdata);
	darray3E	interpolatePointData(darray3E & point, dvecarr3E & pointdata);

	//plotting
	void 		plotCloud( std::string & , std::string, int , bool,  dvecarr3E * extPoints=NULL);
	void 		plotCloud( std::string & , std::string, int , bool,  ivector1D & vertexList, dvecarr3E * extPoints=NULL);
	void 		plotGrid(std::string &, std::string , int, bool, dvecarr3E * extPoints=NULL);
	void 		plotGrid(std::string &, std::string , int, bool, ivector1D & cellList, dvecarr3E * extPoints=NULL);

protected:
	//internal maintenance of the mesh
//	BasicShape*	getShape();
	void 		resizeMesh();
	void 		destroyNodalStructure();
	void 		reshapeNodalStructure();
	void 		rebaseMesh();
	
};



/*! Return global index of the point given its cartesian indices. Follows the ordering sequences z-y-x
 * \param[in] i x cartesian index
 *\param[in] j y cartesian index
 *\param[in] k z cartesian index
 *\param[out] result global index
 */
inline int  UStructMesh::accessPointIndex(int i, int j, int k){
	int index = (m_ny+1) * (m_nz+1) * i + (m_nz+1) * j + k;
	return(index);
};



#endif //__MiMMO_BASICMESHES_HH
