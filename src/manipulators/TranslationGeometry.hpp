/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
 *
 *  mimmo is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmo is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#ifndef __TRANSLATIONGEOMETRY_HPP__
#define __TRANSLATIONGEOMETRY_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *	\class TranslationGeometry
 *	\brief TranslationGeometry is the class that applies a translation to a given geometry patch.
 *
 *	The used parameters are the translation value and the direction of the translation axis.
 *
 *	=========================================================
 * ~~~
 *	|--------------------------------------------------------------|
 *	|                 Port Input                                   |
 *	|-------|----------|-------------------|-----------------------|
 *	|PortID | PortType | variable/function | DataType              |
 *	|-------|----------|-------------------|-----------------------|
 *	| 21    | M_AXIS   | setDirection      | (ARRAY3, FLOAT)       |
 *	| 30    | M_VALUED | setTranslation    | (SCALAR, FLOAT)       |
 *	| 12    | M_FILTER | setFilter         | (VECTOR, FLOAT)       |
 *	| 99    | M_GEOM   | setGeometry       | (SCALAR, MIMMO_)      |
 *	|-------|----------|-------------------|-----------------------|
 *
 *
 *	|---------------------------------------|-----------------------|
 *	|            Port Output                |                       |
 *	|-------|-----------|-------------------|-----------------------|
 *	|PortID | PortType  | variable/function | DataType              |
 *	|-------|-----------|-------------------|-----------------------|
 *	| 11    | M_GDISPLS | getDisplacements  | (VECARR3, FLOAT)      |
 *	|-------|-----------|-------------------|-----------------------|
 * ~~~
 *	=========================================================
 *
 */
class TranslationGeometry: public BaseManipulation{
private:
	//members
	darray3E	m_direction;	/**<Components of the translation axis.*/
	double		m_alpha;        /**<Angle of translation in radiant. */
	dvector1D   m_filter;       /**<Filter field for displacements modulation. */
    dvecarr3E   m_displ;        /**<Resulting displacements of geometry vertex.*/

public:
	TranslationGeometry(darray3E direction = { {0, 0, 0} });
	TranslationGeometry(const bitpit::Config::Section & rootXML);
	~TranslationGeometry();

	TranslationGeometry(const TranslationGeometry & other);
	TranslationGeometry & operator=(const TranslationGeometry & other);

	void        buildPorts();

	void        setAxis(darray3E direction);
	void        setDirection(darray3E direction);
    void        setTranslation(double alpha);
    void        setFilter(dvector1D filter);

    dvecarr3E   getDisplacements();

	void 	    execute();
	
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");
};

REGISTER(BaseManipulation, TranslationGeometry, "mimmo.TranslationGeometry")

};

#endif /* __ROTATIONGEOMETRY_HPP__ */
