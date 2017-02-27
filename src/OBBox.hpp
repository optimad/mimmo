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
#ifndef __OBBox_HPP__
#define __OBBox_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *	\date			24/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief Oriented Bounding Box calculator.
 *
 *	Builds the oriented bounding box of a 3D object (Point Clouds or superficial tessellations), passed as MimmoObject. 
 *
 *	=========================================================
 * ~~~
 *	|-------------------------------------------------------------------------------------|
 *	|                    Port Input                                                       |
 *	|-------|-------------|---------------------------------------|-----------------------|
 *	|PortID | PortType    | variable/function                     | DataType		      |
 *	|-------|-------------|---------------------------------------|-----------------------|
 *	| 99    | M_GEOM      | m_geometry                            |(SCALAR, MIMMO_)       |
 *	|-------|-------------|---------------------------------------|-----------------------|
 * 
 *
 *  |---------------------------------------------------------------|
 *	|               Port Output               						|
 *	|-------|-------------|-------------------|---------------------|
 *	|PortID | PortType    | variable/function | DataType		  	|
 *	|-------|-------------|-------------------|---------------------|
 *	| 20    | M_POINT     | getOrigin         |	(ARRAY3, FLOAT)		|
 *	| 22    | M_AXES      | getAxes 	      |	(ARR3ARR3, FLOAT)	|
 *	| 23    | M_SPAN      | getSpan           |	(ARRAY3, FLOAT)		|
 *	|-------|-------------|-------------------|---------------------|
 * ~~~
 *	=========================================================
 *
 */
class OBBox: public BaseManipulation {

protected:
	darray3E	m_origin;		/**< Origin of the OBB.*/
	darray3E    m_span;         /**< Span of the OBB. */
	dmatrix33E	m_axes;			/**< reference system of the bbox, ordered aaccording maximum shape variance */
	
public:
	OBBox();
	virtual ~OBBox();

	//copy operators/constructors
	OBBox(const OBBox & other);
	OBBox & operator=(const OBBox & other);

	void buildPorts();

	//clean structure;
	void 		clearOBBox();

	//internal methods
	darray3E	getOrigin();
	darray3E	getSpan();
	dmatrix33E	getAxes();
	void		setGeometry(MimmoObject*);
	//plotting wrappers
	void		plot(std::string directory, std::string filename, int counter, bool binary);

	//building method
	void execute();

	virtual void absorbSectionXML(bitpit::Config::Section & slotXML, std::string name = "");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");
	
	
protected:
	virtual void plotOptionalResults();
	
private:
	dmatrix33E 		eigenVectors( dmatrix33E &, darray3E & eigenValues);
	void			evaluateCovarianceMatrix(dmatrix33E &, darray3E &);
	dmatrix33E		evalCovTriangle(dvecarr3E & vv);
	dmatrix33E 		createMatrix(darray3E v1, darray3E v2);
	void 			adjustBasis( dmatrix33E &, darray3E & eigenValues);
};

}

#endif /* __LATTICE_HPP__ */
