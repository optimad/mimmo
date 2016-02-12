#ifndef __MiMMO_VTKINTERFACES__
#define __MiMMO_VTKINTERFACES__

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
#include "bitpit_IO.hpp" //?? what is its new interface in bitpit? please check it out


/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Amos Tori
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Class for VTK Output of structured grids  
 *
 *	Class for managing output on file of deformable structured grids as unstructured ones, based on VTK's formats (*.vtu)
 */

class VTK_BASICMESH: public VTKUnstructuredGrid<VTK_BASICMESH>
{
protected:
	dvecarr3E * m_points; /*!< pointer to list of vertices */
	ivector2D * m_connectivity; /*!< pointer to cell-vertex connectivity */
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
 *	\brief Class for VTK Output of Cloud points 
 *
 *	Class for managing output on file of point clouds, based on VTK's formats (*.vtu)
 */
class VTK_BASICCLOUD: public VTK_UnstructuredGrid<VTK_BASICCLOUD>
{
protected:
	dvecarr3E * m_points; /*!< pointer lo cloud points list */
public:
	VTK_BASICCLOUD();
	VTK_BASICCLOUD(std::string dir_, std::string name_, std::string cod_, int nP_);
	~VTK_BASICCLOUD();
	
	void linkData(dvecarr3E &);
	void unlinkData();
	void Flush(  fstream &str, string codex_, string name  ) ; //CRTP
};

#endif // __MiMMO_VTKINTERFACES__
