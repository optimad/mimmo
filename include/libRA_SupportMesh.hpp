#ifndef __SUPPORTMESH__HH__
#define __SUPPORTMESH__HH__

// Standard C++ library
# include <vector>
# include <sstream>
# include <iostream>
# include <algorithm>
# include <map>

// local libs
# include "CG.hpp"
# include "lib_STL_IO.hpp"



/*!
 *	\date			31/1/2016
 *	\authors		Federico Gallizio
 * 	\authors 		Scafroglia Ugo
 *	\authors		Rocco Arpa
 *	\copyright		Copyright 2016 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Service Data Structure for portion of Unstructured Triangulated Meshes 
 *
 *	BASE_Support is a Service Data Structure class for Unstructured Meshes, directly derived from Class_SurfTri BitPit object of libMesh.
 *	It has no purposes of Data Structure Container, but is used to store informations of portion/patch of a preexistent Mother 
 * 	Unstructured Meshes. Retains tools for adapting and managing locally the Mother Mesh, without altering its original layout,
 *	as well as tools for extracting/handling the boundary contour of the patch.
 * 
 */
class BASE_Support: protected Class_SurfTri 
{
     
protected:
	//members
	std::string     name;			/*!< proper name of the support selected */
	std::string     classType;		/*!< string identifying type of class */
	Class_SurfTri  *mother;			/*!< pointer to mother structure ->input of the class */
	ivector1D	currentSelection;	/*!< Collection of simplicies selected from Mother Mesh -> input for the class */
	ivector1D 	mapVertMother;		/*!< list of Vertices of Mother Tesselation composing Support*/
	ivector1D       minInvConn;		/*!< Minimal inverse connectivity map, for private fast access purpose only */
public:
	//constructors
	BASE_Support();
	BASE_Support(SHAPE *mother_, int flag);
	BASE_Support(Class_SurfTri *mother_, ivector1D & list);
	BASE_Support(const BASE_Support &other);
	~BASE_Support();

	BASE_Support & operator=(const BASE_Support &other);
	
	//cleaning and updating
	virtual void cleanSupport();
	
	//set methods
	virtual void setSupport(SHAPE * mother, int flag);
	virtual void setSupport(Class_SurfTri *mother_, ivector1D &list);
	void setName(std::string );
	
	//get methods
	std::string getName();
	std::string getClassType();
	Class_SurfTri * getMother();
	ivector1D getVNeighs( int iVert);
	ivector1D getSNeighs( int iVert);
	
	ivector1D getSimplexMap();
	ivector1D getVertexMap();
	
	dvecarr3E getNormals();
	dvecarr3E getVNormals();
	
	darray3E getNormals(int sIndex);
	darray3E getVNormals(int iVert);
	
	ivector1D getWholeBoundaryIndex();
	ivector1D getConstrBoundaryIndex();
	ivector1D getFreeBoundaryIndex();
	
	dvecarr3E getWholeBoundary();
	dvecarr3E getConstrBoundary();
	dvecarr3E getFreeBoundary();
	
	//generic utilities
	dvecarr3E indexToPoints(ivector1D & map);
	int simp_LocToMother(int);
	int simp_MotherToLoc(int);
	
	ivector1D simp_LocToMother(ivector1D &);
	ivector1D simp_MotherToLoc(ivector1D &);
	
	ivector1D cleanVertexList(ivector1D & list);
	ivector1D cleanSimplexList(ivector1D & list);
	
	
	
	//plot utilities
	void plotVTU(std::string dir, std::string name, bool flag, dvecarr3E * vertices=NULL);
	void plotScalarVTU(std::string dir, std::string name, bool flag, int counter, dvector1D & field, dvecarr3E * vertices=NULL);
	void plotVectorVTU(std::string dir, std::string name, bool flag, int counter, dvecarr3E & field, dvecarr3E * vertices=NULL);
	void plotSTL(string filename,  dvecarr3E & vertices, bool flag);
	void plotDGF(string filename,  dvecarr3E & vertices);
	
	void plotWholeBoundary(std::string, std::string, bool, std::string);
	void plotConstrBoundary(std::string, std::string, bool, std::string);
	void plotFreeBoundary(std::string, std::string, bool,std::string);
	void plotBoundary(std::string, std::string, bool, std::string, ivector1D & list);

protected: 
	void compMinInvConn();
	void setClassType(std::string);
private:
	 ivector1D getMotherBoundaryVIndex();
	 Class_SurfTri extractWholeTess(dvecarr3E * vertices = NULL);
	 Class_SurfTri createBoundary3DCurve(ivector1D & list); 
//  	/*! TODO da rivedere. Secondo me non dovrebbe servire ad un cazzo di niente */
// 	  void cleanCornerNodes();
}; //end of class SupportMesh	


/*!
 *	\date			31/1/2016
 *	\authors		Federico Gallizio
 * 	\authors 		Scafroglia Ugo
 *	\authors		Rocco Arpa
 *	\copyright		Copyright 2016 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Service Data Structure for portion of Unstructured Triangulated Meshes, with boundary splitting methods 
 *
 *	SupportSplitB is a class derived from BASE_Support. It adds further features to the base class, for splitting the 
 * 	boundary contour in four sub-branches, whenever is possible.
 * 
 */
class SupportSplitB: public BASE_Support{

protected:
	//members
	ivector1D     boundSegmentedMap;	/*!< map of segmented boundaries of the patch */
	Class_SurfTri boundCurve	;	/*!<boundary curve grid manager */
	ivector1D mapOrderBoundary;		/*!<ordering map of the curve */ 
	int typeSegmentation;			/*!< interface on type of segmentation tech you want ot implement.0 is autosplitting, 1 adjacent geometries. 2 is manual. No other type supplied */
	darray3E seed;				/*!< seed point for autosplitting of boundaries */
	svector1D list_geo;			/*!< list of adjacent geometry file names, for identification of boundary branches */
	ivector1D splitPoints;			/*!< list of ordered points for boundary splitting purposes */
	bool switchBoundary;                     /*!< switch on/off boundary manipulation methods */
	bool closedLoop;			/*!< boolean. if true inhibits boundary manipulation, if support surface is a closed loop*/
public:
	SupportSplitB();
	SupportSplitB(SHAPE *mother_, int flag);
	SupportSplitB(Class_SurfTri *mother_, ivector1D & list);
	
	SupportSplitB(const SupportSplitB &other);
	~SupportSplitB();
	     
	SupportSplitB & operator=(const SupportSplitB &other);
	
	virtual void setSupport(SHAPE * mother, int flag);
	virtual void setSupport(Class_SurfTri * mother, ivector1D & list);
	
	void cleanSplitBoundaries();
	void cleanSupport();
	
	bool splitBoundary();
	virtual bool splitBoundary_geoList(svector1D & geoList);
	virtual bool splitBoundary_auto(darray3E seed, int N_, int n_eff);
	
	bool addSplitPoint(darray3E);
	bool addSplitPoint(int);
	bool removeSplitPoint(int);
	bool removeSplitPoint(darray3E);
	
	bool addSplitPoint(darray3E, bool );
	virtual bool addSplitPoint(int, bool );
	virtual bool removeSplitPoint(int, bool);
	bool removeSplitPoint(darray3E, bool);
	
	bool addSplitPoints(dvecarr3E & pointlist);
	bool addSplitPoints(ivector1D & pointlist);
	bool removeSplitPoints(ivector1D & pointlist);
	bool removeSplitPoints(dvecarr3E & pointlist);
	
	int whichBranch(darray3E );
	int whichBranch(int );
	
	int orientSegmentation(int branchIndex);
		
	int sortBoundary(darray3E &startVert );
	int sortBoundary(int startVert );
	
	int getNearestBPointIndex(darray3E);
	
	int getSplittingType();
	dvecarr3E getSplitPoints();
	ivector1D getSplitPointsIndex();
	darray3E getSeed();
	svector1D getGeoList();
	
	int getBoundaryBranchSize();
	ivector1D getBoundaryBranchIndex(int flag);
	dvecarr3E getBoundaryBranch(int flag);
	ivector2D getBoundaryBranchMap();

	Class_SurfTri * getBoundaryCurve();
	
	void switchBoundaryManipulation(bool);
	bool getBoundaryManipulationStatus();
	void plotBoundaryBranches(std::string, std::string, bool,std::string);
	
private:
	ivector1D regularizeBoundary(dvector1D & angles);
	ivector1D findSPoints(ivector1D & map, dvector1D & angles, int guessN);
	Class_SurfTri extractBoundaryFromGeoFile(std::string & filename);
	ivector1D compareCurves(Class_SurfTri & mother_,Class_SurfTri & daughter_);
	void adjust3DCurveOrientation(int);
	void adjustSeedOrientation(int);
	void createBoundary(); 
}; //end of class SupportSplitB

#endif
