#ifndef CAMILOBASICMESHES2D_HH
#define CAMILOBASICMESHES2D_HH

// local libs
#include "Operators.hpp"
#include "customOperators.hpp"
#include "libRA_BasicMeshes.hpp"

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Lachance Steve
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Virtual class for 2D uniform structured mesh --> based on ucartmesh paradigm of Bitpit 
 *
 *	Base class for uniform structured grids, suitable for derivation of 2D meshes on pure cartesian or polar
 *	coordinates system.
 */

   class BASE_UStructMesh2D{
 
    protected:
	  darray2E origin; /**< Mesh origin*/
	  darray2E span;  /**< Mesh span */
	  double dx, dy;	      /**< Mesh spacing in each direction */	
	  int nx, ny;	      /**< Mesh number of cells in each direction */
	  dvector1D xnode, ynode; /**<Lists holding the center cells coordinates of the mesh, in local reference sistem */
	  dvector1D xedge, yedge;/**<Lists holding the point coordinates of the mesh, in local reference system */
	  
    public:	     
  //Building stuffs	    
    BASE_UStructMesh2D();
    ~BASE_UStructMesh2D();
  
  // Copy constructor & assignment operators
    BASE_UStructMesh2D(const BASE_UStructMesh2D & other);
    BASE_UStructMesh2D & operator=(const BASE_UStructMesh2D & other);
  
  //generic mantenaince - get/set stuffs  
    void clearMesh();  
    void resizeMesh();
    void destroyNodalStructure();
    void reshapeNodalStructure();

    darray2E getOrigin();
    darray2E getSpan();
    darray2E getSpacing();
    ivector1D getDimension();

    darray2E getGridCCell(int);
    darray2E getGridCCell(int, int);
    darray2E getGridPoint(int);
    darray2E getGridPoint(int, int);

    darray2E getGlobalCCell(int);
    darray2E getGlobalCCell(int, int);
    darray2E getGlobalPoint(int);
    darray2E getGlobalPoint(int, int);
    
    ivector1D getCellNeighs(int);
    ivector1D getCellNeighs(int, int);
    
  //access mesh methods
    void returnCellID(darray2E &point, int &i, int &j);
    void returnCellID(dvector1D &point, int &i, int &j);
    int  accessCellData(int i, int j);
    void accessCellData(int N_, int &i, int &j);
    int  accessPointData(int i, int j);
    void accessPointData(int N_, int &i, int &j);
  
  // Mesh transforming  
    void translateMesh(double, double);
    
    virtual darray2E transfToGlobal( darray2E & point) = 0;
    virtual dvector1D transfToGlobal( dvector1D & point) = 0;
    virtual dvecarr2E transfToGlobal( dvecarr2E & list_points) = 0;    
    virtual darray2E transfToLocal( darray2E & point) = 0;
    virtual dvector1D transfToLocal( dvector1D & point) = 0;
    virtual dvecarr2E transfToLocal( dvecarr2E & list_points) = 0;    
    
  // various  
    virtual double getCellVolume(int index) = 0;
    virtual darray2E getLocalScaling(int index) = 0;
    virtual darray2E getLocalScaling(int i_, int j_) = 0;
    
  // interpolators  
    double interpolateCellData(darray2E & point, dvector1D & celldata);
    int interpolateCellData(darray2E & point, ivector1D & celldata);
    darray3E interpolateCellData(darray2E & point, dvecarr3E & celldata);

    double interpolatePointData(darray2E & point, dvector1D & pointdata);
    int interpolatePointData(darray2E & point, ivector1D & pointdata);
    darray3E interpolatePointData(darray2E & point, dvecarr3E & pointdata);

    void plotCloud( std::string & , std::string, int , bool,  dvecarr3E * extPoints=NULL);
    void plotCloud( std::string & , std::string, int , bool,  ivector1D & vertexList, dvecarr3E * extPoints=NULL);
    void plotGrid(std::string &, std::string , int, bool, dvecarr3E * extPoints=NULL);
    void plotGrid(std::string &, std::string , int, bool, ivector1D & cellList, dvecarr3E * extPoints=NULL);

    darray3E vecEnlarge(darray2E);
    darray2E vecShrink(darray3E);
  };
  
/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Lachance Steve
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Class for 2D uniform Cartesian Mesh 
 *
 *	Derived class from BASE_UStructMesh2D for uniform cartesian grids.
 */
  
 class UQuadMesh : public BASE_UStructMesh2D {
   
protected:
	std::string classType;
 
public:
	
	UQuadMesh();
	UQuadMesh(darray2E origin_, double spanX_, double spanY_, int nx_, int ny_);
	~UQuadMesh();
	
	// Copy constructor & assignment operators
	UQuadMesh(const UQuadMesh & other);
	UQuadMesh & operator=(const UQuadMesh & other);
	
	virtual std::string getClassType();
	
	void setMesh(darray2E origin_, double spanX_, double spanY_,  int nx_, int ny_);
	// Mesh transforming  
	void scaleMesh(double, double);
	
	// various  
	double getCellVolume(int index);
	darray2E getLocalScaling(int index);
	darray2E getLocalScaling(int i_, int j_);
	
	darray2E transfToGlobal( darray2E & point);
	dvector1D transfToGlobal( dvector1D & point);
	dvecarr2E transfToGlobal( dvecarr2E & list_points);    
	darray2E transfToLocal( darray2E & point);
	dvector1D transfToLocal( dvector1D & point);
	dvecarr2E transfToLocal( dvecarr2E & list_points);    

}; 

#endif