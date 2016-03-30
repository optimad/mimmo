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
#ifndef __MIMMOGEOMETRY_HPP__
#define __MIMMOGEOMETRY_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *	\date			30/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief MimmoGeometry is the class that wrap a geometry Mimmo Object .
 *
 *	The parameter of linked geometry are given by wrapper functions of this BaseManipulation
 *	derived class.
 *
 */
class MimmoGeometry: public BaseManipulation{
public:
	bool		m_write; 		/*< If true it writes the geometry on file during the execution.*/
	std::string	m_filename;		/*< Name of file to write the geometry.*/

	MimmoGeometry();
	~MimmoGeometry();

	MimmoGeometry(const MimmoGeometry & other);
	MimmoGeometry & operator=(const MimmoGeometry & other);

	int			getType();
	long		getNVertex();
	long		getNCells();
	dvecarr3E	getVertex();
	darray3E	getVertex(long i);
	ivector1D	getConnectivity(long i);

	bool		setVertex(dvecarr3E & vertex);
	bool		setVertex(darray3E & vertex);
	bool		modifyVertex(darray3E & vertex, long id);
	bool		setConnectivity(ivector2D * connectivity);
	void		setWrite(bool write);
	void		setFilename(std::string filename);

	void		write();

	void 		execute();

};

}

#endif /* __MIMMOGEOMETRY_HPP__ */
