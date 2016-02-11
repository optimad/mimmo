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
#ifndef __INPUTDOF_HPP__
#define __INPUTDOF_HPP__

#include <string>

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief InputDoF is the class that import the initialization of the degrees of freedom.
 *
 *	InputDoF is derived from BaseManipulation class.
 *	It uses and it sets the base members m_ndeg and m_displacements.
 *	After the execution of an object InputDoF, the number of degrees of freedom and their initial
 *	displacements are set. The values can be set by the user in construction or by set methods of base class.
 *	Moreover the initial values can be read from an input file of text (ascii)
 *	by setting the related flag.
 *
 */
class InputDoF: public BaseManipulation{
private:
	//members
	bool			m_readFromFile;	/**<True if the object reads the values from file.*/
	std::string		m_filename;		/**<Name of the input file. The file has to be an ascii text file.*/

public:
	InputDoF(bool readFromFile = false);
	InputDoF(std::string filename);
	InputDoF(uint32_t ndeg, dvecarr3E & displacements);
	~InputDoF();

	void setReadFromFile(std::string filename);
	void setFilename(std::string filename);


	//relationship methods
protected:
	void	recoverDisplacements();   //called in exec
public:
	void 	exec();

};

#endif /* __INPUTDOF_HPP__ */
