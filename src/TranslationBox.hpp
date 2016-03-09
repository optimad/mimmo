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
#ifndef __TRANSLATIONBOX_HPP__
#define __TRANSLATIONBOX_HPP__

#include "BaseManipulation.hpp"
#include "FFDLattice.hpp"

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief TranslationBox is the class that applies the a translation to the box of a latticeBox object.
 *
 *	The degrees of freedom is the translation value (ndeg = 1) and the direction of the translation
 *	is a parameter of the object.
 *	The displacements has to be one term. As the displacements are array3E only the first term is used,
 *	i.e if alpha is the value of the translation in the chosen direction the only array of displacement
 *	is set m_displ[0]={aplha, 0, 0};
 *
 */
class TranslationBox: public BaseManipulation{
private:
	//members
	darray3E			m_direction;	/**<Components of the translation axis.*/
	darray3E			m_origin;		/**<Origin of box to be deformed (recovered in recoverInfo and used in useInfo).*/


public:
	TranslationBox(darray3E direction = { {0, 0, 0} });
	~TranslationBox();

	TranslationBox(const TranslationBox & other);
	TranslationBox & operator=(const TranslationBox & other);

	void setDirection(darray3E direction);
	void setTranslation(double alpha);

	//relationship methods
protected:

public:
	void 	useInfo();
	void 	execute();

};

#endif /* __TRANSLATIONBOX_HPP__ */
