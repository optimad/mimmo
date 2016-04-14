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
#ifndef __MIMMOOBJECT_HPP__
#define __MIMMOOBJECT_HPP__

#include "bitpit.hpp"
#include "MiMMO_TypeDef.hpp"

namespace mimmo{

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief MimmoObject is the base container of geometry for MiMMo library
 *
 *	MiMMO container has the information about the geometry mesh:
 *	- the type of the linked mesh of class bitpit::Patch (generic patch, surface patch or volume patch);
 *	- the pointer to the mesh of class bitpit::Patch;
 *	- the boolean indicator of the nature of the patch (true if the patch is created within the MiMMo container).
 *
 */
class MimmoObject{
public:

	//members
	int						m_type;				/**<Type of geometry (0 = generic patch, 1 = surface mesh, 2 = volume mesh). */
	bitpit::PatchKernel*	m_geometry;			/**<Reference geometry. */
	bool					m_internalPatch;	/**<If the geometry is internally created. */
	livector1D				m_mapData;			/**<Map of vertex ids actually set, for aligning external vertex data to bitpit::Patch ordering */
	liimap					m_mapDataInv;		/**<Inverse of Map of vertex ids actually set, for aligning external vertex data to bitpit::Patch ordering */
	liimap					m_mapCellInv;		/**<Inverse of Map of cell ids actually set, for aligning external vertex data to bitpit::Patch ordering */

public:
	MimmoObject(int type = 1);
	MimmoObject(int type, dvecarr3E & vertex, ivector2D * connectivity = NULL);
	MimmoObject(int type, bitpit::PatchKernel* geometry);
	~MimmoObject();

	MimmoObject(const MimmoObject & other);
	MimmoObject & operator=(const MimmoObject & other);

	void 		clear();
	bool		isEmpty();
	int			getType();
	long		getNVertex();
	long		getNCells();
	dvecarr3E	getVertex();
	darray3E	getVertex(long i);
	ivector1D	getConnectivity(long i);
	ivector2D	getConnectivity();
	bitpit::PatchKernel*		getGeometry();
	livector1D&	getMapData();
	long		getMapData(int i);
	liimap&		getMapDataInv();
	int			getMapDataInv(long id);
	liimap&		getMapCellInv();
	
	bool		setVertex(dvecarr3E & vertex);
	bool		setVertex(darray3E & vertex);
	bool		modifyVertex(darray3E & vertex, long id);
	bool		setConnectivity(ivector2D * connectivity);
	bool		setGeometry(int type, bitpit::PatchKernel* geometry);
	bool		setMapData();
	bool		setMapCell();
	bool		cleanGeometry();

	livector1D 	getVertexFromCellList(livector1D cellList);
	ivector1D	convertVertexIDtoLocal(livector1D vertList);
	
	void		write(std::string filename);
};

}

#endif /* __MIMMOOBJECT_HPP__ */
