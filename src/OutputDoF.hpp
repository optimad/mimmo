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
#ifndef __OUTPUTDOF_HPP__
#define __OUTPUTDOF_HPP__

#include "BaseManipulation.hpp"
#include <string>

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief OutputDoF is the class that write the degrees of freedom.
 *
 *	OutputDoF is derived from BaseManipulation class.
 *	It uses and it write the base members m_ndeg and m_displacements.
 *	After the execution of an object OutputDoF, the number of degrees of freedom and their initial
 *	displacements are write on a text file (ascii).
 */
class OutputDoF: public BaseManipulation{
private:
	//members
	std::string		m_filename;		/**<Name of the output file. The file will be an ascii text file.*/

public:
	OutputDoF(BaseManipulation* parent = NULL);
	OutputDoF(std::string filename = "output.txt", BaseManipulation* parent = NULL);
	~OutputDoF();

	OutputDoF(const OutputDoF & other);
	OutputDoF & operator=(const OutputDoF & other);

	void setFilename(std::string filename);

	//relationship methods
protected:
	void	recoverDisplacements();   //called in exec
public:
	void 	exec();

};

#endif /* __OUTPUTDOF_HPP__ */
