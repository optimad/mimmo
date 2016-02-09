#ifndef CAMILOCOREPATCHES_HH
#define CAMILOCOREPATCHES_HH

// local libs
#include "libRA_BasicEleShapes.hpp"
#include "libRA_BasicMeshes.hpp"
#include "lib_STL_IO.hpp"

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Virtual class for Basic Patches Representation 
 *
 *	Base class for Volumetric Core Element, suitable for interaction with Triangulation data Structure.
 */
 
class BASE_PATCH {
  
protected:
  std::string patchType;
  
public:
    BASE_PATCH();
    ~BASE_PATCH();
    
    BASE_PATCH(const BASE_PATCH &);
    BASE_PATCH & operator=(const BASE_PATCH &);
    
    std::string getPatchType();
    
    ivector1D includeTriangulation(Class_SurfTri * );
    ivector1D excludeTriangulation(Class_SurfTri * );
    
    ivector1D includeCloudPoints(dvecarr3E &);
    ivector1D excludeCloudPoints(dvecarr3E &);
    
    ivector1D includePIDTriangulation(SHAPE *, int );
    ivector1D excludePIDTriangulation(SHAPE *, int );

    ivector1D includePIDTriangulation(SHAPE *, ivector1D & );
    ivector1D excludePIDTriangulation(SHAPE *, ivector1D & );
    
    bool isTriangleIncluded(dvecarr3E &);
    bool isTriangleIncluded(Class_SurfTri *, int indexT);
    
    virtual bool isPointIncluded(darray3E)=0;
    virtual bool isPointIncluded(Class_SurfTri *, int indexV)=0;
};

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Volumetric Cubic Patch class  
 *
 *	Volumetric Core Element, of cubic shape, suitable for interaction with Triangulation data Structure.
 */
 
class PATCH_CUBE : public BASE_PATCH, public ElementalCube{
  
 
public:
    PATCH_CUBE();
    PATCH_CUBE(darray3E & origin_, double w_, double h_, double d_); 
    ~PATCH_CUBE();
    
    PATCH_CUBE(const PATCH_CUBE &);
    PATCH_CUBE & operator=(const PATCH_CUBE &);
    
    bool isPointIncluded(darray3E);
    bool isPointIncluded(Class_SurfTri *, int indexV);

    void plotVTUpatch(std::string folder, std::string name, int number);
};

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Volumetric Cylindric Patch class  
 *
 *	Volumetric Core Element, of cylindric shape, suitable for interaction with Triangulation data Structure.
 */
 
class PATCH_CYLINDER : public BASE_PATCH, public ElementalCylinder{
 
public:
    PATCH_CYLINDER();
    PATCH_CYLINDER(darray3E & origin_, double r_, double h_); 
    ~PATCH_CYLINDER();
    
    PATCH_CYLINDER(const PATCH_CYLINDER &);
    PATCH_CYLINDER & operator=(const PATCH_CYLINDER &);
   
    bool isPointIncluded(darray3E);
    bool isPointIncluded(Class_SurfTri *, int indexV);
   
};

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Volumetric Spherical Patch class  
 *
 *	Volumetric Core Element, of spherical shape, suitable for interaction with Triangulation data Structure.
 */
 
class PATCH_SPHERE : public BASE_PATCH, public ElementalSphere{
 
public:
    PATCH_SPHERE();
    PATCH_SPHERE(darray3E & origin_, double r_); 
    ~PATCH_SPHERE();
    
    PATCH_SPHERE(const PATCH_SPHERE &);
    PATCH_SPHERE & operator=(const PATCH_SPHERE &);
   
    bool isPointIncluded(darray3E);
    bool isPointIncluded(Class_SurfTri *, int indexV);
};

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Volumetric Hexahedrical Patch class  
 *
 *	Volumetric Core Element, of generic hexahedrical shape, suitable for interaction with Triangulation data Structure.
 */
 
class PATCH_HEXA : public BASE_PATCH, public ElementalHexaHedron{
 
public:
    PATCH_HEXA();
    PATCH_HEXA(dvecarr3E & vertList_); 
    ~PATCH_HEXA();
    
    PATCH_HEXA(const PATCH_HEXA &);
    PATCH_HEXA & operator=(const PATCH_HEXA &);
   
    bool isPointIncluded(darray3E);
    bool isPointIncluded(Class_SurfTri *, int indexV);

    void plotVTUpatch(std::string folder, std::string name, int number);
};

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Volumetric Cubic HUll class  
 *
 *	Cartesian Structured Mesh, of cubic shape, suitable for interaction with Triangulation data Structure.
 */
 
class HULL_CUBE : public BASE_PATCH, public UCubicMesh{
 
public:
    HULL_CUBE();
    HULL_CUBE(darray3E origin_, double spanX_, double spanY_, double spanZ_, int nx_, int ny_, int nz_); 
    ~HULL_CUBE();
    
    HULL_CUBE(const HULL_CUBE &);
    HULL_CUBE & operator=(const HULL_CUBE &);
   
    bool isPointIncluded(darray3E);
    bool isPointIncluded(Class_SurfTri *, int indexV);
};

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Volumetric Cylindric HUll class  
 *
 *	Structured Mesh, of cylindrical shape and reference frame, suitable for interaction with Triangulation data Structure.
 */
 
class HULL_CYLINDER : public BASE_PATCH, public UCylindricalMesh{
 
public:
    HULL_CYLINDER();
    HULL_CYLINDER(darray3E origin_, double spanR_, double spanZ_, dvector1D & thetalim_, int nr_, int nt_, int nz_); 
    ~HULL_CYLINDER();
    
    HULL_CYLINDER(const HULL_CYLINDER &);
    HULL_CYLINDER & operator=(const HULL_CYLINDER &);
   
    bool isPointIncluded(darray3E);
    bool isPointIncluded(Class_SurfTri *, int indexV);
};

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Volumetric Cylindric HUll class  
 *
 *	Structured Mesh, of cylindrical shape and reference frame, suitable for interaction with Triangulation data Structure.
 */
 
class HULL_SPHERE : public BASE_PATCH, public USphericalMesh{
 
public:
    HULL_SPHERE();
    HULL_SPHERE(darray3E origin_, double spanR_, dvector1D thetalim_, dvector1D philim_, int nr_, int nt_, int np_); 
    ~HULL_SPHERE();
    
    HULL_SPHERE(const HULL_SPHERE &);
    HULL_SPHERE & operator=(const HULL_SPHERE &);
   
    bool isPointIncluded(darray3E);
    bool isPointIncluded(Class_SurfTri *, int indexV);
};

#endif

