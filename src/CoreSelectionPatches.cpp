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
#include "CoreSelectionPatches.hpp"

using namespace bitpit;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//BaseSelPatch IMPLEMENTATION 
/*
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
 
/*! Basic Constructor */
BaseSelPatch::BaseSelPatch(){};
/*! Basic Destructor */
BaseSelPatch::~BaseSelPatch(){};

/*! Copy Constructor 
 * \param[in] other BaseSelPatch object where copy from
 */
BaseSelPatch::BaseSelPatch(const BaseSelPatch & other){
  m_patchType = other.m_patchType;
};

/*! Copy Operator 
 * \param[in] other BaseSelPatch object where copy from
 */
BaseSelPatch & BaseSelPatch::operator=(const BaseSelPatch & other){
   m_patchType = other.m_patchType;
};

/*! Get type of patch actually instantiated */
std::string BaseSelPatch::getPatchType(){return(m_patchType);};


/*! Given a bitpit class bitpit::Patch tessellation, return cell identifiers of those simplex inside the volume of
 * BaseSelPatch 
 * \param[in] tri target tessellation
 * \param[out] result list-by-ids of simplicies included in the volumetric patch
 */
ivector1D BaseSelPatch::includeGeometry(bitpit::Patch * tri ){
  
  int nCells = tri->getCellCount();	
  ivector1D result(nCells); 
  int counter=0;
  
  for(auto &cell : tri->cells()){
   
    if(isSimplexIncluded(tri, cell.get_id())){
      result[counter] = cell.get_id();
      ++counter;
    }
  }
  result.resize(counter);
  return(result);
};

/*! Given a bitpit class bitpit::Patch tessellation, return cell identifiers of those simplex outside the volume of
 * BaseSelPatch 
 * \param[in] tri target tesselation
 * \param[out] result list-by-ids of simplicies outside the volumetric patch
 */
ivector1D BaseSelPatch::excludeGeometry(bitpit::Patch * tri){
  
	int nCells = tri->getCellCount();	
	ivector1D result(nCells); 
	int counter=0;
	
	for(auto &cell : tri->cells()){
		
		if(!isSimplexIncluded(tri, cell.get_id())){
			result[counter] = cell.get_id();
			++counter;
		}
	}
	result.resize(counter);
	return(result);
};

/*! Given a list of vertices of a point cloud, return indices of those vertices included into 
 * the volume of BaseSelPatch object 
 * \param[in] list list of cloud points
 * \param[out] result list-by-indices of vertices included in the volumetric patch
 */
ivector1D BaseSelPatch::includeCloudPoints(dvecarr3E & list){

  int size = list.size();
  ivector1D result(size); 
  int counter=0;
  
  for(int i=0; i<size; ++i){
   
    if(isPointIncluded(list[i])){
      result[counter] = i;
      ++counter;
    }
  }
  result.resize(counter);
  return(result);
};

/*! Given a list of vertices of a point cloud, return indices of those vertices outside 
 * the volume of BaseSelPatch object 
 * \param[in] list list of cloud points
 * \param[out] result list-by-indices of vertices outside the volumetric patch
 */
ivector1D BaseSelPatch::excludeCloudPoints(dvecarr3E & list){
  
  int size = list.size();
  ivector1D result(size); 
  int counter=0;
  
  for(int i=0; i<size; ++i){
   
    if(!isPointIncluded(list[i])){
      result[counter] = i;
      ++counter;
    }
  }
  result.resize(counter);
  return(result);
  
};

/*! Given a bitpit class bitpit::Patch point cloud, return identifiers of those points inside the volume of
 * BaseSelPatch  
 * \param[in] list list of cloud points
 * \param[out] result list-by-indices of vertices included in the volumetric patch
 */
ivector1D BaseSelPatch::includeCloudPoints(bitpit::Patch * tri){
	
	int nVert = tri->getVertexCount();	
	ivector1D result(nVert); 
	int counter=0;
	
	for(auto &vertex : tri->vertices()){
		
		if(isPointIncluded(tri, vertex.get_id())){
			result[counter] = vertex.get_id();
			++counter;
		}
	}
	result.resize(counter);
	return(result);
};

/*! Given a bitpit class bitpit::Patch point cloud, return identifiers of those points outside the volume of
 * BaseSelPatch  
 * \param[in] list list of cloud points
 * \param[out] result list-by-indices of vertices outside the volumetric patch
 */
ivector1D BaseSelPatch::excludeCloudPoints(bitpit::Patch * tri){
	
	int nVert = tri->getVertexCount();	
	ivector1D result(nVert); 
	int counter=0;
	
	for(auto &vertex : tri->vertices()){
		
		if(!isPointIncluded(tri, vertex.get_id())){
			result[counter] = vertex.get_id();
			++counter;
		}
	}
	result.resize(counter);
	return(result);	
};

/* Given a pidded part of a triangulated tesselation, return those pidded triangles inside the patch 
 * \param[in] sh target triangulated tesselation
 * \param[in] pidT target PID number
 * \param[out] result list-by-indices of pidded triangles inside the volumetric patch
 */
// ivector1D BaseSelPatch::includePIDTriangulation(SHAPE * sh, int pidT ){
//   
//   ivector1D temp = sh->extractPidded(pidT);
//   int size = temp.size();
//   ivector1D result(size); 
//   int counter=0;
//   
//   for(int i=0; i<size; ++i){
//     if(isTriangleIncluded(sh, temp[i]) ){
//       result[counter] = temp[i];
//       ++counter;
//     }
//   }
//   result.resize(counter);
//   return(result);
// };

/* Given a pidded part of a triangulated tesselation, return those pidded triangles outside the patch 
 * \param[in] sh target triangulated tesselation
 * \param[in] pidT target PID number
 * \param[out] result list-by-indices of pidded triangles outside the volumetric patch
 */
// ivector1D BaseSelPatch::excludePIDTriangulation(SHAPE * sh, int pidT){
// 
//   ivector1D temp = sh->extractPidded(pidT);
//   int size = temp.size();
//   ivector1D result(size); 
//   int counter=0;
//   
//   for(int i=0; i<size; ++i){
//    
//     if(!isTriangleIncluded(sh, temp[i]) ){
//       result[counter] = temp[i];
//       ++counter;
//     }
//   }
//   result.resize(counter);
//   return(result);
// };

/* Given one or more pidded part of a triangulated tesselation, return those pidded triangles inside the patch 
 * \param[in] sh target triangulated tesselation
 * \param[in] pidT target PID number
 * \param[out] result list-by-indices of pidded triangles inside the volumetric patch
 */
// ivector1D BaseSelPatch::includePIDTriangulation(SHAPE * sh, ivector1D & pidList ){
//   
//   int size = pidList.size();
//   ivector1D result(sh->nSimplex);
//   
//   int counterGlobal = 0;
//   for(int i=0; i<size; i++){
//    
//     ivector1D temp = includePIDTriangulation(sh, pidList[i]);
//     
//     int tempSize = temp.size();
//     for(int j=0; j<tempSize; ++j){
// 	result[counterGlobal] = temp[j];
// 	++counterGlobal;
//     } 
//   }
//     result.resize(counterGlobal);
//     return(result);
// };

/* Given one or more pidded part of a triangulated tesselation, return those pidded triangles outside the patch 
 * \param[in] sh target triangulated tesselation
 * \param[in] pidT target PID number
 * \param[out] result list-by-indices of pidded triangles outside the volumetric patch
 */
// ivector1D BaseSelPatch::excludePIDTriangulation(SHAPE * sh, ivector1D & pidList ){
// 
//   int size = pidList.size();
//   ivector1D result(sh->nSimplex);
//   
//   int counterGlobal = 0;
//   for(int i=0; i<size; i++){
//    
//     ivector1D temp = excludePIDTriangulation(sh, pidList[i]);
//     
//     int tempSize = temp.size();
//     for(int j=0; j<tempSize; ++j){
// 	result[counterGlobal] = temp[j];
// 	++counterGlobal;
//     } 
//   }
//     result.resize(counterGlobal);
//     return(result);
// };


/*! Return True if at least one vertex of a given triangle is included in the volumetric patch
 * \param[in] simplexVert 3 vertices of the given Triangle
 * \param[out] result boolean
 */   
bool BaseSelPatch::isSimplexIncluded(dvecarr3E & simplexVert){
  
  bool check = false;
  for(int i=0; i<simplexVert.size(); ++i){
   check = check || isPointIncluded(simplexVert[i]); 
  }
  return(check);
};

/*! Return True if at least one vertex of a given triangle is included in the volumetric patch
 * \param[in] tri pointer to a Class_SurfTri tesselation 
 * \param[in] indexT triangle index of tri.
 * \param[out] result boolean
 */ 
bool BaseSelPatch::isSimplexIncluded(bitpit::Patch * tri, int indexT){

  Cell cell = tri->getCell(indexT);
//  CellIterator it = tri->getCell(indexT);
  long * conn = cell.getConnect();
  int nVertices = cell.getVertexCount();
  bool check = false;
  for(int i=0; i<nVertices; ++i){ //recover vertex index
   check = check || isPointIncluded(tri, conn[i]); 
  }
  return(check);
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//HullCube IMPLEMENTATION 
/*
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
   
/*! Basic Constructor */
HullCube::HullCube(){
	m_patchType="hullCube";
};

 /*! Custom Constructor. Set mesh origin, span size, dimension and nodal data structure
   * \param[in] origin_ point origin in global reference system
   * \param[in] spanX_ width span -> must be > 0;
   * \param[in] spanY_ height span -> must be > 0;
   * \param[in] spanZ_ depth span -> must be > 0;
   * \param[in] nx_   number of cells in x-direction
   * \param[in] ny_   number of cells in y-direction
   * \param[in] nz_   number of cells in z-direction
   */
HullCube::HullCube(darray3E origin_, double spanX_, double spanY_, double spanZ_, int nx_, int ny_, int nz_):
		      UCubicMesh(origin_,spanX_,spanY_, spanZ_, nx_, ny_, nz_){
			      m_patchType="hullCube";
		}; 

/*! Basic Destructor */
HullCube::~HullCube(){};

/*! Copy Constructor 
 * \param[in] other HullCube object where copy from
 */
HullCube::HullCube(const HullCube & other){
  *(static_cast<BaseSelPatch*>(this))= *(static_cast<const BaseSelPatch*> (&other));
  *(static_cast<UCubicMesh*>(this))= *(static_cast<const UCubicMesh*> (&other));
};

/*! Copy Operator 
 * \param[in] other HullCube object where copy from
 */
HullCube & HullCube::operator=(const HullCube & other){
  *(static_cast<BaseSelPatch*>(this))= *(static_cast<const BaseSelPatch*> (&other));
  *(static_cast<UCubicMesh*>(this))= *(static_cast<const UCubicMesh*> (&other));
  return(*this);
};


/*! Return True if the given point is included in the volumetric patch
 * \param[in] point given vertex
 * \param[out] result boolean
 */
bool HullCube::isPointIncluded(darray3E point){

  bool check = true;
  darray3E temp = transfToLocal(point);
  darray3E pSp = getSpan();
  
  for(int i=0; i< pSp.size(); ++i){   
    temp[i] = temp[i]/pSp[i];
    
    check = check && ((temp[i] >= 0.0) && (temp[i]<=1.0));
  }
  
  return(check);  
};

/*! Return True if the given point is included in the volumetric patch
 * \param[in] tri pointer to a bitpit::Patch tesselation / point cloud
 * \param[in] indexV id of a vertex belonging to tri;
 * \param[out] result boolean
 */
bool HullCube::isPointIncluded(bitpit::Patch * tri, int indexV){

  bool check = true;
  darray3E coords = tri->getVertex(indexV).getCoords();
  darray3E temp = transfToLocal(coords);
  darray3E pSp = getSpan();
  
  for(int i=0; i< pSp.size(); ++i){   
    temp[i] = temp[i]/pSp[i];
    
    check = check && ((temp[i] >= 0.0) && (temp[i]<=1.0));
  }
  
  return(check);  
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//HullCylinder IMPLEMENTATION 
/*
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

/*! Basic Constructor */
HullCylinder::HullCylinder(){
	m_patchType="hullCylinder";
};

/*!Custom Constructor. Set mesh origin, span size, dimension and nodal data structure
 * \param[in] origin_ point origin in global reference system
 * \param[in] spanR_ max base radius span: must be >0; 
 * \param[in] spanZ_ max cylinder height span: must be >0; 
   \param[in] thetalim lower and upper limits on the angular coordinate 
 * \param[in] nr_   number of cells in r-direction
 * \param[in] nt_   number of cells in theta-direction
 * \param[in] nz_   number of cells in z-direction
 */
HullCylinder::HullCylinder(darray3E origin_, double spanR_, double spanZ_, dvector1D & thetalim_, int nr_, int nt_, int nz_):
		      UCylindricalMesh(origin_, spanR_, spanZ_, thetalim_, nr_, nt_, nz_){
			      m_patchType="hullCylinder";
		}; 

/*! Basic Destructor */
HullCylinder::~HullCylinder(){};

/*! Copy Constructor 
 * \param[in] other HullCylinder object where copy from
 */
HullCylinder::HullCylinder(const HullCylinder & other){
  *(static_cast<BaseSelPatch*>(this))= *(static_cast<const BaseSelPatch*> (&other));
  *(static_cast<UCylindricalMesh*>(this))= *(static_cast<const UCylindricalMesh*> (&other));
};

/*! Copy Operator 
 * \param[in] other PATCH_CYLINDER object where copy from
 */
HullCylinder & HullCylinder::operator=(const HullCylinder & other){
  *(static_cast<BaseSelPatch*>(this))= *(static_cast<const BaseSelPatch*> (&other));
  *(static_cast<UCylindricalMesh*>(this))= *(static_cast<const UCylindricalMesh*> (&other));
  return(*this);
};


/*! Return True if the given point is included in the volumetric patch
 * \param[in] point given vertex
 * \param[out] result boolean
 */
bool HullCylinder::isPointIncluded(darray3E point){

  bool check = true;
  darray3E pOr = getOrigin();
  darray3E temp2 = transfToLocal(point);
  darray3E pSp = getSpan();
  
    double radius_norm = temp2[0]/pSp[0];
    double height_norm = temp2[2]/pSp[2];
    
    check = check && (radius_norm <=1.0);
    check = check && ((height_norm >= 0.0) && (height_norm<=1.0));
  
  return(check);  
};

/*! Return True if the given point is included in the volumetric patch
 * \param[in] tri pointer to a bitpit::Patch tessellation / point cloud
 * \param[in] indexV id of a vertex belonging to tri.
 * \param[out] result boolean
 */
bool HullCylinder::isPointIncluded(bitpit::Patch * tri, int indexV){

 bool check = true;
  darray3E pOr = getOrigin();
  darray3E coords = tri->getVertex(indexV).getCoords();
  darray3E temp2 = transfToLocal(coords);
  darray3E pSp = getSpan();
  
    double radius_norm = temp2[0]/pSp[0];
    double height_norm = temp2[2]/pSp[2];
    
    check = check && (radius_norm <=1.0);
    check = check && (temp2[1]>=0 && temp2[1]<=pSp[1]);
    check = check && ((height_norm >= 0.0) && (height_norm<=1.0));
  
  return(check);  
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//HullSphere IMPLEMENTATION 
/*
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

/*! Basic Constructor */
HullSphere::HullSphere(){
	m_patchType="hullSphere";
};

/*! Custom Constructor. Set mesh origin, span size, dimension and nodal data structure
 * \param[in] origin_ point origin in global reference system
 * \param[in] spanR_  max radius span: must be >0;
 * \param[in] thetalim_ lower and upper limits on the polar coordinate
 * \param[in] philim_ lower and upper limits on the azimuthal coordinate
 * \param[in] nr_   number of cells in r-direction
 * \param[in] nt_   number of cells in theta-direction
 * \param[in] np_   number of cells in phi-direction
 */
HullSphere::HullSphere(darray3E origin_, double spanR_, dvector1D thetalim_, dvector1D philim_, int nr_, int nt_, int np_):
			USphericalMesh(origin_, spanR_, thetalim_, philim_, nr_, nt_, np_){
				m_patchType="hullSphere";
			}; 

/*! Basic Destructor */			
HullSphere::~HullSphere(){};

/*! Copy Constructor 
 * \param[in] other HullSphere object where copy from
 */
HullSphere::HullSphere(const HullSphere & other){
  *(static_cast<BaseSelPatch*>(this))= *(static_cast< const BaseSelPatch*> (&other));
  *(static_cast<USphericalMesh*>(this))= *(static_cast<const USphericalMesh*> (&other));  
};

/*! Copy Operator 
 * \param[in] other HullSphere object where copy from
 */
HullSphere & HullSphere::operator=(const HullSphere & other){
  *(static_cast<BaseSelPatch*>(this))= *(static_cast<const BaseSelPatch*> (&other));
  *(static_cast<USphericalMesh*>(this))= *(static_cast<const USphericalMesh*> (&other));
   return(*this);  
};   

/*! Return True if the given point is included in the volumetric patch
 * \param[in] point given vertex
 * \param[out] result boolean
 */
bool HullSphere::isPointIncluded(darray3E point){
   bool check = false;
   darray3E temp2 = transfToLocal(point);
   darray3E pSp = getSpan();
  
    double radius_norm = temp2[0]/pSp[0];
    
    check = (radius_norm <=1.0);
  
  return(check);  
};

/*! Return True if the given point is included in the volumetric patch
 * \param[in] tri pointer to a bitpit::Patch tessellation / point cloud
 * \param[in] indexV id of a vertex belonging to tri.
 * \param[out] result boolean
 */
bool HullSphere::isPointIncluded(bitpit::Patch * tri, int indexV){
   
  bool check = false;
   darray3E coords = tri->getVertex(indexV).getCoords();
   darray3E temp2 = transfToLocal(coords);
   darray3E pSp = getSpan();
  
    double radius_norm = temp2[0]/pSp[0];
    check = (radius_norm <=1.0);
    check = check && (temp2[1]>=0 && temp2[1]<=pSp[1]);
    check = check && (temp2[2]>=0 && temp2[2]<=pSp[2]);
    
  return(check);
  
};
