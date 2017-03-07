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
#ifndef __SPECULARPOINTS_HPP__
#define __SPECULARPOINTS_HPP__

#include "ProjectCloud.hpp"

namespace mimmo{

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief SpecularPoints is a class that mirrors a point cloud w.r.t. a reference plane, on a target geometry if any
 *
 *	SpecularPoints is derived from ProjectCloud class. Given a certain number of points and a reference plane,
 *  the class mirrors such points with respect to this plane. If a geometry is linked, project mirrored points on 
 *  this target geometry. Any data attached, as scalar/vector float data format, are mirrored as well. 
 *
 *	=========================================================
 * ~~~
 *	|-----------------------------------------------------------------------|
 *	|                 Port Input                       		 		 		|
 *	|-------|---------------|-------------------|---------------------------|
 *	|PortID | PortType 		| variable/function | DataType 		 		 	|
 *	|-------|---------------|-------------------|---------------------------|
 *	| 0     | M_COORDS 		| setPoints         | (VECARR3, FLOAT)		 	|
 *	| 10    | M_DISPLS 		| setVectorData     | (VECARR3, FLOAT) 	     	|
 *	| 19    | M_SCALARFIELD | setScalarData     | (VECTOR, FLOAT)		 	|
 *	| 29    | M_PLANE		| setPlane	     	| (ARRAY4, FLOAT)		 	|
 *  | 32    | M_VALUEB      | setInsideOut      | (SCALAR, BOOL)            |
 *  | 140   | M_VALUEB      | setForce          | (SCALAR, BOOL)            |
 *	| 99    | M_GEOM 		| m_inside          | (SCALAR, MIMMO_)			|
 *	|-------|---------------|-------------------|---------------------------|
 *
 *
 *	|--------------------------------------------------------------------|
 *	|            Port Output     	                     		 		 |
 *	|-------|---------------|-------------------|------------------------|
 *	|PortID | PortType 		| variable/function | DataType 		 		 |
 *	|-------|---------------|-------------------|------------------------|
 *	| 0     | M_COORDS 		| getCloudResult	| (VECARR3, FLOAT)		 |
 *	| 10    | M_DISPLS 		| getCloudVectorData| (VECARR3, FLOAT) 	     |
 *	| 19    | M_SCALARFIELD | getCloudScalarData| (VECTOR, FLOAT)		 |
 *	|-------|---------------|-------------------|------------------------|
 * ~~~
 *	=========================================================
 *
 */
class SpecularPoints: public ProjectCloud{
private:
    bool        m_insideout; /**< plane direction for mirroring */
    bool        m_force;    /**< if true force the mirroring of points that belong to the symmetry plane. */
	darray4E 	m_plane;  /**< reference plane */
	dvector1D	m_scalar; /**< float scalar data attached to original points*/
	dvecarr3E   m_vector; /**< float vector data attached to original points*/
	dvector1D	m_scalarMirrored; /**< resulting float scalar data after mirroring*/
	dvecarr3E   m_vectorMirrored; /**< resulting float scalar data after mirroring*/
	
public:
	SpecularPoints();
	SpecularPoints(const bitpit::Config::Section & rootXML);
	~SpecularPoints();

	SpecularPoints(const SpecularPoints & other);
	SpecularPoints & operator=(const SpecularPoints & other);

	void	buildPorts();

	dvector1D getOriginalScalarData();
	dvecarr3E getOriginalVectorData();

	dvector1D getCloudScalarData();
	dvecarr3E getCloudVectorData();
	
	darray4E  getPlane();
    bool      isInsideOut();
    bool      isForce();
	
	void	setVectorData(dvecarr3E data);
	void	setScalarData(dvector1D data);
	void	setPlane(darray4E plane);
	void	setPlane(darray3E origin, darray3E normal);
    void    setInsideOut(bool flag);
    void    setForce(bool flag);
	
	void 	execute();
	
	void clear();
	
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name= "");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");
	
protected:
	virtual void plotOptionalResults();
};

REGISTER(BaseManipulation, SpecularPoints, "MiMMO.SpecularPoints")
};

#endif /* __SPECULARPOINTS_HPP__ */
