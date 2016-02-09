#ifndef __RA_VTKINTERFACES__
#define __RA_VTKINTERFACES__

//std libs
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>
#include <map>

// local libs
#include "Operators.hpp"
#include "customOperators.hpp"
#include "Class_VTK.hpp"
#include "Class_SurfTri.hpp"
#include "Class_UCartMesh.hpp"

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Allen Woody
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief VTK Class Interface for I/O managing of Unstructured Triangulated Meshes  
 *
 *	VTK_SHAPE is an I/O manager class for reading and writing Unstructured Meshes and data attached to them, 
 *	in VTK native .vtu format.
 */
class VTK_SHAPE: public VTK_UnstructuredGrid <VTK_SHAPE>
{
private:
	ivector1D simplTypes; /**< support list for simplicial types of tassellation. Useful while reading */
protected :
	Class_SurfTri * tri;        /**< Link class  to a surfTri external structure */
	dvecarr3E	  * ext_vertex; /**< Link class to an  external vertex list */
	dvector1D	  * scalarField;/**< Link data to an external scalar field */
	dvecarr3E	  * vectField;  /**< Link data to an external vector field */
public:
	VTK_SHAPE();
	VTK_SHAPE(Class_SurfTri *, std::string filename_);
	VTK_SHAPE(Class_SurfTri *, std::string dir_, std::string name_, std::string cod_);
	~VTK_SHAPE();
	
	void      Flush(  fstream &str, string codex_, string name  ) ; //CRTP
	void      Absorb( fstream &str, string codex_, string name  ) ; //CRTP
	void      linkTriMesh(Class_SurfTri *);
	void      linkExternalVertexList(dvecarr3E &);
	void 	  linkScalarField(dvector1D &);
	void 	  linkVectorField(dvecarr3E &);
	void      unlinkTriMesh();
	void      unlinkExternalVertexList();
	void 	  unlinkScalarField();
	void 	  unlinkVectorField();
};

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Allen Woody
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief VTK Class Interface for output managing of LevelSet on 3D Cartesian Structured Meshes  
 *
 *	VTK_LSet is an Output manager class for writing Level Set scalar field and Level Set gradient field on  
 *	Cartesian Background Mesh in VTK native .vtr format.
 */
class VTK_LSet: public VTK_RectilinearGrid <VTK_LSet>
{
	//protected:
public:
	Class_UCartMesh3D * mesh3D; /**< Link class to UcartMesh external structure */
	dvector1D * scalarField; /**< Link class to an external double scalar field */
	dvecarr3E * vectorField; /**< Link class to an external double vector field */
	
public:
	VTK_LSet();
	VTK_LSet(std::string dir, std::string name);
	VTK_LSet(std::string dir, std::string name, std::string codex);
	~VTK_LSet();
	
	void linkData(Class_UCartMesh3D &);
	void linkScalarField(dvector1D &);
	void linkVectorField(dvecarr3E &);
	void unlinkData();
	void unlinkScalarField();
	void unlinkVectorField();
	
	void Flush(  fstream &str, string codex_, string name  ) ; //CRTP
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
 *	\brief Class for VTK Output purpose 
 *
 *	Class for managing output on file of deformable structured grids as unstructured ones, based on VTK's formats (*.vtu)
 */

class VTK_BASICMESH: public VTK_UnstructuredGrid<VTK_BASICMESH>
{
protected:
	dvecarr3E * points; /*!< pointer to list of vertices */
	ivector2D * connectivity; /*!< pointer to cell-vertex connectivity */
public:
	VTK_BASICMESH();
	VTK_BASICMESH(std::string dir_, std::string name_, std::string cod_, int nP_, int nC_, int nConn_);
	~VTK_BASICMESH();
	
	void linkData(dvecarr3E &, ivector2D & );
	void unlinkData();
	void Flush(  fstream &str, string codex_, string name  ) ; //CRTP
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
 *	\brief Class for VTK Output purpose 
 *
 *	Class for managing output on file of point clouds, based on VTK's formats (*.vtu)
 */
class VTK_BASICCLOUD: public VTK_UnstructuredGrid<VTK_BASICCLOUD>
{
protected:
	dvecarr3E * points; /*!< pointer lo cloud points list */
public:
	VTK_BASICCLOUD();
	VTK_BASICCLOUD(std::string dir_, std::string name_, std::string cod_, int nP_);
	~VTK_BASICCLOUD();
	
	void linkData(dvecarr3E &);
	void unlinkData();
	void Flush(  fstream &str, string codex_, string name  ) ; //CRTP
};

/*!
 *	\date			30/1/2016
 *	\authors		Gallizio Federico
 *	\authors 		Suzanne Vega
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief VTK Class Interface for Output managing of Boundary Contours of Unstructured Surface Meshes  
 *
 *	VTK_3DCURVE is an Output manager class for writing Cloud Points related to Unstructured Surface Meshes
 *	boundary contours, handled by a Class_SurfTri object,in VTK native .vtu format.
 */
class VTK_3DCURVE: public VTK_UnstructuredGrid<VTK_3DCURVE>
{
	
protected:
	Class_SurfTri *curve;
	
public:
	VTK_3DCURVE();
	VTK_3DCURVE(Class_SurfTri & triP, std::string dir_, std::string name_, std::string cod_);
	VTK_3DCURVE(Class_SurfTri & triP, std::string filename_);
	~VTK_3DCURVE();
	
	void linkData(Class_SurfTri &);
	void unlinkData();
	void Flush(  fstream &str, string codex_, string name  ) ; //CRTP
};

#endif
