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

enum FileType{STL};		/**Extension of file to read the geometry.*/

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

	FileType	m_type;			/**<Extension of file to read the geometry.*/
	bool		m_read; 		/**< If true it reads the geometry from file during the execution.*/
	std::string	m_rfilename;	/**< Name of file to read the geometry.*/

	bool		m_write; 		/**< If true it writes the geometry on file during the execution.*/
	std::string	m_wfilename;	/**< Name of file to write the geometry.*/

private:
	bool		m_local;		/**<Is the geometry locally instantiated?.*/

public:
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
	void		setFileType(FileType type);
	void		setFileType(int type);
	void		setRead(bool read);
	void		setReadFilename(std::string filename);
	void		setWrite(bool write);
	void		setWriteFilename(std::string filename);

	bool		write();
	void		read();

	void 		execute();

};

}

#endif /* __MIMMOGEOMETRY_HPP__ */
