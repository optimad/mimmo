#ifndef CAMILOBASICELESHAPES_HH
#define CAMILOBASICELESHAPES_HH

// local libs
#include "Operators.hpp"
#include "customOperators.hpp"
#include "LinearAlgebra.hpp"
#include "Class_SurfTri.hpp"
#include "Class_VolTri.hpp"

/*!
*	\date			31/12/2015
*	\authors		Gallizio Federico
* 	\authors 		Ratzinger Joseph
*	\authors		Arpa Rocco
*	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
*	\par			License:\n
*	This class version is released under .
*
*	\brief Virtual class for Basic Elemental Shapes Representation 
*
*	Base class for elemental shapes definition,from cubes and hexahedra to cylinders and spheres
*/

class BASE_ElementalShapes
{

protected:
	darray3E shapeOrigin; /**< Shape origin*/
	std::string eleShapeType; /**< string for Shape identification */  
	Class_SurfTri * surface; /**< Storage of the superficial tesselation */
	Class_VolTri * volume; /**< Storage of the volume unstructured Mesh */
	short int surfaceTessType; /**<check type of your tessellation.*/
public:	     
//Building stuffs	    
	BASE_ElementalShapes();
	BASE_ElementalShapes(darray3E & origin_);
	~BASE_ElementalShapes();

// Copy constructor & assignment operators
	BASE_ElementalShapes(const BASE_ElementalShapes & other);
	BASE_ElementalShapes & operator=(const BASE_ElementalShapes & other);

//generic mantenaince - get/set stuffs  
	virtual void clearShape();
	virtual void setShapeOrigin(double, double, double);
	darray3E getShapeOrigin();
	std::string getClassType();

	Class_SurfTri * getSurfaceTess();
	Class_VolTri * getVolumeTess();

	short int checkTessellation();

	virtual void tessellate(int )=0;
	virtual void triangulate(int )=0;
	virtual void tetrahedralize(int )=0;

	void exportSurfaceStl(std::string, std::string);
	void exportSurface(std::string);
	void exportVolume(std::string);
};

/*!
*	\date			31/12/2015
*	\authors		Gallizio Federico
* 	\authors 		Ratzinger Joseph
*	\authors		Arpa Rocco
*	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
*	\par			License:\n
*	This class version is released under .
*
*	\brief Cube Class
*
*	Elemental shape Cube class, derived from BASE_ElementalShapes
*/

class ElementalCube: public BASE_ElementalShapes 
{

protected:
	double width; /**< Width of the shape*/
	double height; /**< Height of the shape*/
	double depth; /**< Depth of the shape*/

public:	     
//Building stuffs	    
	ElementalCube();
	ElementalCube(darray3E origin_, double w_, double h_, double d_);
	~ElementalCube();
// Copy constructor & assignment operators
	ElementalCube(const ElementalCube & other);
	ElementalCube & operator=(const ElementalCube & other);

//generic mantenaince - get/set stuffs  
	void clearShape();
	void setShapeSpan(double, double, double);
	dvector1D getShapeSpan();

	virtual void tessellate(int );
	virtual void triangulate(int );
	virtual void tetrahedralize(int );
};

/*!
*	\date			31/12/2015
*	\authors		Gallizio Federico
* 	\authors 		Ratzinger Joseph
*	\authors		Arpa Rocco
*	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
*	\par			License:\n
*	This class version is released under .
*
*	\brief Cylinder Class
*
*	Elemental shape Cylinder class, derived from BASE_ElementalShapes
*/

class ElementalCylinder:  public BASE_ElementalShapes {

protected:
	double radius; /**< Base Radius of the shape*/
	double height; /**< Height of the shape*/

public:	     
//Building stuffs	    
	ElementalCylinder();
	ElementalCylinder(darray3E origin_, double r_, double h_);
	~ElementalCylinder();

// Copy constructor & assignment operators
	ElementalCylinder(const ElementalCylinder & other);
	ElementalCylinder & operator=(const ElementalCylinder & other);

//generic mantenaince - get/set stuffs  
	void clearShape();
	void setShapeSpan(double, double);
	dvector1D getShapeSpan();

	virtual void tessellate(int );
	virtual void triangulate(int );
	virtual void tetrahedralize(int );
};

/*!
*	\date			31/12/2015
*	\authors		Gallizio Federico
* 	\authors 		Ratzinger Joseph
*	\authors		Arpa Rocco
*	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
*	\par			License:\n
*	This class version is released under .
*
*	\brief Sphere Class
*
*	Elemental shape Sphere class, derived from BASE_ElementalShapes
*/

class ElementalSphere :  public BASE_ElementalShapes{

protected:
	double radius; /**< Radius of the shape*/

public:	     
//Building stuffs	    
	ElementalSphere();
	ElementalSphere(darray3E origin_, double r_);
	~ElementalSphere();

// Copy constructor & assignment operators
	ElementalSphere(const ElementalSphere & other);	
	ElementalSphere & operator=(const ElementalSphere & other);

//generic mantenaince - get/set stuffs  
	void clearShape();
	void setShapeSpan(double);
	double getShapeSpan();

	virtual void tessellate(int );
	virtual void triangulate(int );
	virtual void tetrahedralize(int );
};  

/*!
*	\date			31/12/2015
*	\authors		Gallizio Federico
* 	\authors 		Ratzinger Joseph
*	\authors		Arpa Rocco
*	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
*	\par			License:\n
*	This class version is released under .
*
*	\brief Generic HexaHedron Class
*
*	Elemental shape HexaHedron class, derived from BASE_ElementalShapes
*/

class ElementalHexaHedron :  public BASE_ElementalShapes{

protected:
	dvecarr3E vertex; /**< list of 8 vertices composing the hexahedron in x-y-z ordering*/
	dvecarr4E planes; /**< list of inplicit planes coeffs for each hexahedron face*/
	ivector2D face_conn; /**< Connectivity matrix of faces/vertices */
	ivector2D vert_conn; /**< Connectivity matrix of vertices faces */
	
public:	     
//Building stuffs	    
	ElementalHexaHedron();
	ElementalHexaHedron(dvecarr3E & vertList);
	~ElementalHexaHedron();

// Copy constructor & assignment operators
	ElementalHexaHedron(const ElementalHexaHedron & other);
	ElementalHexaHedron & operator=(const ElementalHexaHedron & other);

//generic mantenaince - get/set stuffs  
	void clearShape();
	void resetHexaHedron();
	void setShapeOrigin(double, double, double);

	void setBaseVertex(dvecarr3E * list = NULL);
	void changeVertex(int index, darray3E &, bool);
	void changeVertices(ivector1D & indices, dvecarr3E & vert);

	dvecarr3E getHexaVertices();
	dvecarr4E getHexaPlanes();
	ivector2D getFaceConnectivity();
	ivector2D getVertConnectivity();

	darray3E getHexaVertex(int index);
	darray4E getPlane(int index);
	darray3E getPlaneNormal(int index);
	ivector1D getFaceConn(int index);
	ivector1D getVertConn(int index);

	bool checkBelongPolygon(darray3E & point, int planeIndex);

	virtual void tessellate(int );
	virtual void triangulate(int );
	virtual void tetrahedralize(int );

private:
	void resetVertices();
	void setFaceConnectivity();
	void setVertConnectivity();
	void setPlanes();
};  

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Ratzinger Joseph
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Generic Convex Hull
 *
 *	Class for handling the 3D convex HULL of a shape (derived from BASE_ElementalShapes)
 *  	STILL WORKING 
 */
//class ConvexHull: public BASE_ElementalShapes{};	


#endif 



