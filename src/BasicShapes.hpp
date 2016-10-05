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

#ifndef __MIMMOBASICSHAPES_HH
#define __MIMMOBASICSHAPES_HH

#include "bitpit.hpp"
#include "MiMMO_TypeDef.hpp"
#include "MimmoObject.hpp"
#include "SortAlgorithms.hpp"
#include <unordered_set>

namespace mimmo{

/*!
 * @enum ShapeType
 * @brief Identifies the type of elemental shape supported by BasicShape class
 */
enum class ShapeType{
	CUBE /**< cubic or generic voxel-shaped objects */,
	CYLINDER /**< cylindric objects */, 
	SPHERE /**< spherical objects */
	
};	

/*!
 * @enum CoordType
 * @brief Specify type of conditions to distribute NURBS node in a given coordinate of the shape
 */
enum class CoordType{
	UNCLAMPED /**< free clamping conditions */, 
	CLAMPED /**< force nurbs to pass on the extremal point of the interval (clamping) */, 
	PERIODIC /**< provide periodic condition on the interval extrema */, 
	SYMMETRIC /**< provide simmetry condition on the interval extrema */
}; 

/*!
 *	\date			03/01/2015
 *	\authors		Edoardo Lombardi
 *	\authors		Arpa Rocco
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
class BasicShape {

protected:
	ShapeType	m_shape; 		/**< shape identifier, see BasicShape::ShapeType enum */
	darray3E	m_origin; 		/**< origin of your shape. May change according to the shape type*/
	darray3E	m_span;			/**< coordinate span of your current shape, in its local reference system*/
	darray3E    m_infLimits;	/**< inferior limits of coordinate of your current shape */
	dmatrix33E	m_sdr;			/**< axis position of the local reference system w.r.t absolute one*/
	darray3E 	m_scaling; 		/**< scaling vector of dimensional coordinates */

	std::array<CoordType,3>		m_typeCoord;	/**< identifiers for coordinate type definition.DEFAULT is clamped*/

	dvecarr3E m_bbox;	//point of bounding box of the shape, temporary not copied build with get bounding box !!!
	
public:
	
	//Costructors, Destructor, Copy/Assignment operators
    BasicShape();
    virtual ~BasicShape();
 	
	//set-get methods
	void		setOrigin(darray3E);
	void		setSpan(double, double, double);
	void		setSpan(darray3E);
	void		setInfLimits(double val, int dir);
	void		setInfLimits(darray3E val);
	void		setRefSystem(darray3E, darray3E, darray3E);
	void		setRefSystem(int, darray3E);
	void		setRefSystem(dmatrix33E axes);
	void		setCoordinateType(CoordType, int dir);
	
	darray3E			getOrigin();
	darray3E			getSpan();
	darray3E			getInfLimits();
	dmatrix33E			getRefSystem();
	CoordType			getCoordinateType(int dir);
	ShapeType   		getShapeType();
	const ShapeType		getShapeType() const;
	
	darray3E	getScaling();
	darray3E	getLocalSpan();
	
	
	//functionalities
	livector1D	includeGeometry(mimmo::MimmoObject * ); 
	livector1D	excludeGeometry(mimmo::MimmoObject * ); 
	
	livector1D	includeGeometry(bitpit::PatchKernel * ); 
	livector1D	excludeGeometry(bitpit::PatchKernel * ); 
	
	livector1D	includeCloudPoints(dvecarr3E &); 
	livector1D	excludeCloudPoints(dvecarr3E &); 
	
	livector1D	includeCloudPoints(bitpit::PatchKernel * ); 
	livector1D	excludeCloudPoints(bitpit::PatchKernel * ); 

	livector1D	includeCloudPoints(mimmo::MimmoObject * ); 
	livector1D	excludeCloudPoints(mimmo::MimmoObject * ); 
	
	bool		isSimplexIncluded(dvecarr3E &);
	bool		isSimplexIncluded(bitpit::PatchKernel * , long int indexT);
    bool		isPointIncluded(darray3E);
	bool		isPointIncluded(bitpit::PatchKernel * , long int indexV);

	virtual bool	intersectShapeAABBox(darray3E bMin, darray3E bMax)=0;
	bool		containShapeAABBox(darray3E bMin, darray3E bMax);
	

	/*! pure virtual method to convert from Local to World Coordinates*/
	virtual	darray3E	toWorldCoord(darray3E  point)=0;
	/*! pure virtual method to convert from World to Local Coordinates*/
	virtual	darray3E	toLocalCoord(darray3E  point)=0;
	/*! pure virtual method to get local Coordinate inferior limits of primitive shape*/
	virtual darray3E	getLocalOrigin()=0;  

protected:
	uint32_t 	intersectShapePlane(int level, darray3E target); 
	darray3E	checkNearestPointToAABBox(darray3E point, darray3E bMin, darray3E bMax);
private:	
	virtual	darray3E	basicToLocal(darray3E  point)=0;
	virtual	darray3E	localToBasic(darray3E  point)=0;
	virtual void 		checkSpan(double &, double &, double &)=0;
	virtual bool 		checkInfLimits(double &, int &dir)=0;
	virtual void		setScaling(double &, double &, double &)=0;
	void				searchKdTreeMatches(bitpit::KdTree<3,bitpit::Vertex,long> & tree,  int indexKdNode, int level, livector1D & result, int &counter );
	void				searchBvTreeMatches(mimmo::BvTree & tree, bitpit::PatchKernel * geo, int indexBvNode, livector1D & result, int &counter);
	void				searchBvTreeNotMatches(mimmo::BvTree & tree, bitpit::PatchKernel * geo, int indexBvNode, livector1D & result, int &counter);	
	virtual void		getTempBBox()=0;
};

/*!
 *	\date			03/01/2015
 *	\authors		Edoardo Lombardi
 *	\authors		Arpa Rocco
 *
 *	\brief Elementary Shape Representation of a Cube
 *
 *	Volumetric Core Element, shaped as a cube, directly derived from BasicShape class.   
 */

class Cube: public BasicShape {
	
public:
	
	//Costructors, Destructor, Copy/Assignment operators
	Cube();
	Cube(darray3E &origin, darray3E &span);
	~Cube();

	Cube(const Cube &);
	Cube & operator=(const Cube &);
	
	//reimplementing pure virtuals	

	darray3E	toWorldCoord(darray3E  point);
	darray3E	toLocalCoord(darray3E  point);
	darray3E	getLocalOrigin();
	bool	intersectShapeAABBox(darray3E bMin, darray3E bMax);
	
private:	
	darray3E	basicToLocal(darray3E  point);
	darray3E	localToBasic(darray3E  point);
	void 		checkSpan(double &, double &, double &);
	bool 		checkInfLimits(double &, int & dir);
	void 		setScaling(double &, double &, double &);
	void		getTempBBox();
};

/*!
 *	\date			03/01/2015
 *	\authors		Edoardo Lombardi
 *	\authors		Arpa Rocco
 *
 *	\brief Elementary Shape Representation of a Cylinder or portion of it
 *
 *	Volumetric Core Element, shaped as a cylinder, directly derived from BasicShape class.   
 */

class Cylinder: public BasicShape {
	
public:
	//Costructors, Destructor, Copy/Assignment operators
	Cylinder();
	Cylinder(darray3E &origin, darray3E &span);
	~Cylinder();

	Cylinder(const Cylinder &);
	Cylinder & operator=(const Cylinder &);
	
	
	//reimplementing pure virtuals
	darray3E	toWorldCoord(darray3E  point);
	darray3E	toLocalCoord(darray3E  point);
	darray3E	getLocalOrigin();
	bool		intersectShapeAABBox(darray3E bMin, darray3E bMax);
	
private:	
	darray3E	basicToLocal(darray3E  point);
	darray3E	localToBasic(darray3E  point);
	void 		checkSpan(double &, double &, double &);
	bool 		checkInfLimits(double &, int &);
	void 		setScaling(double &, double &, double &);
	void		getTempBBox();
};


/*!
 *	\date			03/01/2015
 *	\authors		Edoardo Lombardi
 *	\authors		Arpa Rocco
 *
 *	\brief Elementary Shape Representation of a Sphere or portion of it
 *
 *	Volumetric Core Element, shaped as a sphere, directly derived from BasicShape class.   
 */

class Sphere: public BasicShape {
	
public:
	//Costructors, Destructor, Copy/Assignment operators
	Sphere();
	Sphere(darray3E &origin, darray3E &span);
	~Sphere();
	
	Sphere(const Sphere &);
	Sphere & operator=(const Sphere &);
	
	
	//reimplementing pure virtuals
	darray3E	toWorldCoord(darray3E  point);
	darray3E	toLocalCoord(darray3E  point);
	darray3E	getLocalOrigin();
	bool		intersectShapeAABBox(darray3E bMin, darray3E bMax);

private:	
	darray3E	basicToLocal(darray3E  point);
	darray3E	localToBasic(darray3E  point);
	void 		checkSpan(double &, double &, double &);
	bool 		checkInfLimits(double &, int &);
	void 		setScaling(double &, double &, double &);
	void		getTempBBox();
};

}

#endif //  __MIMMOBASICSHAPES_HH

