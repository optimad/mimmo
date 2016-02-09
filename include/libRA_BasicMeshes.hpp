#ifndef CAMILOBASICMESHES_HH
#define CAMILOBASICMESHES_HH

// local libs
#include "Operators.hpp"
#include "customOperators.hpp"
#include "libRA_VTKInterfaces.hpp"

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Amos Tori
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Virtual class for 3D uniform structured mesh --> based on ucartmesh paradigm of Bitpit 
 *
 *	Base class for uniform structured grids, suitable for derivation of meshes on pure cartesian, cylindrical or spherical
 *	coordinates system.
 */

   class BASE_UStructMesh{
 
    protected:
	  darray3E origin; /**< Mesh origin*/
	  darray3E span;  /**< Mesh span */
	  double dx, dy,dz;	      /**< Mesh spacing in each direction */	
	  int nx, ny,nz;	      /**< Mesh number of cells in each direction */
	  dvector1D xnode, ynode, znode; /**<Lists holding the center cells coordinates of the mesh, in local reference sistem */
	  dvector1D xedge, yedge, zedge; /**<Lists holding the point coordinates of the mesh, in local reference system */
	  
    public:	     
  //Building stuffs	    
    BASE_UStructMesh();
    ~BASE_UStructMesh();
  
  // Copy constructor & assignment operators
    BASE_UStructMesh(const BASE_UStructMesh & other);
    BASE_UStructMesh & operator=(const BASE_UStructMesh & other);
  
  //generic mantenaince - get/set stuffs  
    void clearMesh();  
    void resizeMesh();
    void destroyNodalStructure();
    void reshapeNodalStructure();

    darray3E getOrigin();
    darray3E getSpan();
    darray3E getSpacing();
    ivector1D getDimension();

    darray3E getGridCCell(int);
    darray3E getGridCCell(int, int, int);
    darray3E getGridPoint(int);
    darray3E getGridPoint(int, int, int);

    darray3E getGlobalCCell(int);
    darray3E getGlobalCCell(int, int, int);
    darray3E getGlobalPoint(int);
    darray3E getGlobalPoint(int, int, int);
    
    ivector1D getCellNeighs(int);
    ivector1D getCellNeighs(int, int, int);
    
  //access mesh methods
    void returnCellID(darray3E & point, int &i, int &j, int &k);
    void returnCellID(dvector1D & point, int &i, int &j, int &k);
    int  accessCellData(int i, int j, int k);
    void accessCellData(int N_, int &i, int &j, int &k);
    int  accessPointData(int i, int j, int k);
    void accessPointData(int N_, int &i, int &j, int &k);
  
  // Mesh transforming  
    void translateMesh(double, double, double);
    virtual darray3E transfToGlobal( darray3E & point) = 0;
    virtual dvector1D transfToGlobal( dvector1D & point) = 0;
    virtual dvecarr3E transfToGlobal( dvecarr3E & list_points) = 0;    
    virtual darray3E transfToLocal( darray3E & point) = 0;
    virtual dvector1D transfToLocal( dvector1D & point) = 0;
    virtual dvecarr3E transfToLocal( dvecarr3E & list_points) = 0;    
    
  // various  
    virtual double getCellVolume(int index) = 0;
    virtual darray3E getLocalScaling(int index) = 0;
    virtual darray3E getLocalScaling(int i_, int j_, int k_) = 0;
    
  // interpolators  
    double interpolateCellData(darray3E & point, dvector1D & celldata);
    int interpolateCellData(darray3E & point, ivector1D & celldata);
    darray3E interpolateCellData(darray3E & point, dvecarr3E & celldata);

    double interpolatePointData(darray3E & point, dvector1D & pointdata);
    int interpolatePointData(darray3E & point, ivector1D & pointdata);
    darray3E interpolatePointData(darray3E & point, dvecarr3E & pointdata);

    void plotCloud( std::string & , std::string, int , bool,  dvecarr3E * extPoints=NULL);
    void plotCloud( std::string & , std::string, int , bool,  ivector1D & vertexList, dvecarr3E * extPoints=NULL);
    void plotGrid(std::string &, std::string , int, bool, dvecarr3E * extPoints=NULL);
    void plotGrid(std::string &, std::string , int, bool, ivector1D & cellList, dvecarr3E * extPoints=NULL);
  };
  
/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Amos Tori
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Class for 3D uniform Cartesian Mesh 
 *
 *	Derived class from UStructMesh for uniform cartesian grids.
 */
  
 class UCubicMesh : public BASE_UStructMesh {
   
 protected:
   
   std::string classType;
 
 public:
   
   UCubicMesh();
   UCubicMesh(darray3E origin_, double spanX_, double spanY_, double spanZ_, int nx_, int ny_, int nz_);
   ~UCubicMesh();
  
  // Copy constructor & assignment operators
   UCubicMesh(const UCubicMesh & other);
   UCubicMesh & operator=(const UCubicMesh & other);
   
   virtual std::string getClassType();
  
   void setMesh(darray3E origin_, double spanX_, double spanY_, double spanZ_, int nx_, int ny_, int nz_);
  // Mesh transforming  
    void scaleMesh(double, double, double);
    darray3E transfToGlobal( darray3E & point);
    dvector1D transfToGlobal( dvector1D & point);
    dvecarr3E transfToGlobal( dvecarr3E & list_points);    
    darray3E transfToLocal( darray3E & point);
    dvector1D transfToLocal( dvector1D & point);
    dvecarr3E transfToLocal( dvecarr3E & list_points);    
    
  // various  
    double getCellVolume(int index);
    darray3E getLocalScaling(int index);
    darray3E getLocalScaling(int i_, int j_, int k_);
    
}; 
  
/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Amos Tori
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Class for 3D uniform Cylindrical Mesh 
 *
 *	Derived class from UStructMesh for uniform cartesian grids, in cylindrical reference frame.
 */
  
 class UCylindricalMesh : public BASE_UStructMesh {
   
 protected:
   
   std::string classType;
   double thetaOrigin;
 public:
   
   UCylindricalMesh();
   UCylindricalMesh(darray3E origin_, double spanR_, double spanZ_, dvector1D & thetalim_, int nr_, int nt_, int nz_);
   ~UCylindricalMesh();
  
  // Copy constructor & assignment operators
   UCylindricalMesh(const UCylindricalMesh & other);
   UCylindricalMesh & operator=(const UCylindricalMesh & other);
   
   virtual std::string getClassType();
   void setMesh(darray3E origin_, double spanR_, double spanZ_, dvector1D & thetalim_, int nr_, int nt_, int nz_);
   
  // Mesh transforming  
    void scaleMesh(double, double);
    darray3E transfToGlobal( darray3E & point);
    dvector1D transfToGlobal( dvector1D & point);
    dvecarr3E transfToGlobal( dvecarr3E & list_points);    
    darray3E transfToLocal( darray3E & point);
    dvector1D transfToLocal( dvector1D & point);
    dvecarr3E transfToLocal( dvecarr3E & list_points);    
    
    double getAzimuthalOrigin();
    void  shiftAzimuthalOrigin(double);
    
  // various  
    double getCellVolume(int index);
    darray3E getLocalScaling(int index);
    darray3E getLocalScaling(int i_, int j_, int k_);

    
};  

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Amos Tori
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Class for 3D uniform Spherical Mesh 
 *
 *	Derived class from UStructMesh for uniform cartesian grids, in spherical reference frame.
 */
  
 class USphericalMesh : public BASE_UStructMesh {
   
 protected:
   
   std::string classType;
   double thetaOrigin;
   double phiOrigin;
 public:
   
   USphericalMesh();
   USphericalMesh(darray3E origin_, double spanR_, dvector1D thetalim_, dvector1D philim_, int nr_, int nt_, int np_);
   ~USphericalMesh();
  
  // Copy constructor & assignment operators
   USphericalMesh(const USphericalMesh & other);
   USphericalMesh & operator=(const USphericalMesh & other);
   
   virtual std::string getClassType();
   void setMesh(darray3E origin_, double spanR_, dvector1D thetalim_, dvector1D philim_, int nr_, int nt_, int np_);

  // Mesh transforming  
    void scaleMesh(double);
    darray3E transfToGlobal( darray3E & point);
    dvector1D transfToGlobal( dvector1D & point);
    dvecarr3E transfToGlobal( dvecarr3E & list_points);    
    darray3E transfToLocal( darray3E & point);
    dvector1D transfToLocal( dvector1D & point);
    dvecarr3E transfToLocal( dvecarr3E & list_points);    
    
    double getAzimuthalOrigin();
    void  shiftAzimuthalOrigin(double);

    double getPolarOrigin();
    void  shiftPolarOrigin(double);
    
  // various  
    double getCellVolume(int index);
    darray3E getLocalScaling(int index);
    darray3E getLocalScaling(int i_, int j_, int k_);
    
};

#endif