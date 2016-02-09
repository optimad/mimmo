/*---------------------------------------------------------------------------*\
 *
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief MimmoObject is the base container of geometry for MiMMo library
 *
 *	Bla Bla
 *
 */
class MimmoObject{
private:

	//members
	bool			m_type			/**<Type of geometry (0 = surface mesh, 1 = volume mesh). */
	Patch*			m_geometry		/**<Reference geometry. */

//TODO CAPIRE SE SURFTRI/VOLTRI/PATCH HANNO METODI SET (E GET MANCANTI) PER MODIFICARE E CREARE VERTICI E CONNETTIVITA'
public:
	MimmoObject(bool type = 0);
	MimmoObject(bool type, dvecarr3E & vertex, ivector2D * connectivity = NULL);
	MimmoObject(bool type, Patch* geometry);
	~MimmoObject();

	void 		clear();
	bool		isEmpty();
	bool		getType();
	long		getNVertex();
	long		getNCells();
	dvecarr3E*	getVertex();
	ivector2D*	getConnectivity();
	Patch*		getGeometry();

	bool		setVertex(dvecarr3E & vertex);
	bool		setVertex(int index, darray3E & vertex);
	bool		setConnectivity(ivector2D & connectivity);
	bool		setGeometry(bool type, Patch* geometry);

};

#endif /* __MIMMOOBJECT_HPP__ */
