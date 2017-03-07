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

namespace mimmo{

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief TranslationBox is the class that applies the a translation to a point.
 *
 *	Here the point is supposed to be an origin of a reference system.
 *	The point is translated over a direction and for a quantity set by the user or
 *	an external input.
 *	Result of the translation are saved in result of base class and
 *	in the modified member m_origin.
 * 
 *
 *	=========================================================
 * ~~~
 *	|---------------------------------------------------------------|
 *	|                 Port Input                                   	|
 *	|-------|----------|-------------------|------------------------|
 *	|PortID | PortType | variable/function | DataTypes		       	|
 *	|-------|----------|-------------------|------------------------|
 *	| 20    | M_POINT  | m_origin          | (ARRAY3, FLOAT)	   	|
 *	| 21    | M_AXIS   | m_direction       | (ARRAY3, FLOAT)	   	|
 *	| 30    | M_VALUED | m_alpha           | (SCALAR, FLOAT)	  	|
 *	|-------|----------|-------------------|------------------------|
 *
 *
 *	|---------------------------------------|-----------------------|
 *	|            Port Output                |						|
 *	|-------|-----------|-------------------|-----------------------|
 *	|PortID | PortType  | variable/function | DataTypes		       	|
 *	|-------|-----------|-------------------|-----------------------|
 *	| 20    | M_POINT   | getOrigin         | (ARRAY3, FLOAT)	   	|
 *	| 22    | M_AXIS    | getDirection      | (ARRAY3, FLOAT)	   	|
 *	| 30    | M_VALUED  | getTranslation    | (SCALAR, FLOAT)	   	|
 *	|-------|-----------|-------------------|-----------------------|
 * ~~~
 *	=========================================================
 *
 */
class TranslationBox: public BaseManipulation{
private:
	//members
	darray3E	m_direction;	/**<Components of the translation axis.*/
	darray3E	m_origin;		/**<Origin of box to be deformed.*/
	double		m_alpha;		/**<Value of translation.*/

public:
	TranslationBox(darray3E direction = { {0, 0, 0} });
	TranslationBox(const bitpit::Config::Section & rootXML);
	~TranslationBox();

	TranslationBox(const TranslationBox & other);
	TranslationBox & operator=(const TranslationBox & other);

	void		buildPorts();

	darray3E 	getDirection();
	darray3E 	getOrigin();
	double 		getTranslation();
	void 		setDirection(darray3E direction);
	void 		setOrigin(darray3E origin);
	void 		setTranslation(double alpha);

	void 	execute();

	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	
};

REGISTER(BaseManipulation, TranslationBox, "MiMMO.TranslationBox")
};

#endif /* __TRANSLATIONBOX_HPP__ */
