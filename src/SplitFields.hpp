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
#ifndef __SPLITFIELDS_HPP__
#define __SPLITFIELDS_HPP__

#include "MimmoObject.hpp"
#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *	\date			30/nov/2016
 *	\authors		Rocco Arpa
 *	\authors		
 *
 *	\brief SplitField is an abstract executable block class capable of  
 *         splitting fields associated to a stitched MimmoObject, once the
 *		   original geometries of the stitched object and its division maps are set	  
 *
 *	SplitField is an abstract class. To use its features take a look to its specializations,
 *  here presented as derived class, SplitScalarField and SplitVectorField.
 * 
 *	=========================================================
 * ~~~
 *	|-------------------------------------------------------------------------------|
 *	|                 Port Input                           					        |
 *	|-------|----------|------------------------------------|-----------------------|
 *	|PortID | PortType | variable/function 					| DataType		        |
 *	|-------|----------|------------------------------------|-----------------------|
 *  | 99  	| M_GEOM   | setGeometry      					| (SCALAR, MIMMO_)	    |
 *	| 100   | M_VECGEOM| setSplittedGeometries				| (VECTOR, MIMMO_)	    |
 *  | 104   | M_MAPDCELL| setCellDivisionMap   				| (UN_MAP, LONGPAIRINTLONG) |
 * 	| 105	| M_MAPDVERT| setVertDivisionMap				| (UN_MAP, LONGPAIRINTLONG) |	
 *	|-------|----------|------------------------------------|-----------------------|
 *
 *
 *	|--------------------------------------------------------|-----------------------|
 *	|            Port Output                				 |                       |
 *	|-------|-----------|------------------------------------|-----------------------|
 *	|PortID | PortType  | variable/function 				 | DataType		      	 | 
 *	|-------|-----------|------------------------------------|-----------------------|
 *	|-------|-----------|------------------------------------|-----------------------|
 * ~~~
 *	=========================================================
 *
 */

class SplitField: public BaseManipulation{
	
private:
	int m_topo; 								/**< flag for geometrical topology control*/
protected:	
	std::vector<MimmoObject*>	m_originals;	/**< pointers to external geometries*/
	
	std::unordered_map<long, std::pair<int, long> > m_mapCellDivision; /**< division map of actual ID-cell, part Id, original ID-cell*/
	std::unordered_map<long, std::pair<int, long> > m_mapVertDivision; /**< division map of actual ID-vertex, part Id, original ID-vertex*/
	
public:
	SplitField(int topo);
	virtual ~SplitField();

	SplitField(const SplitField & other);
	SplitField & operator=(const SplitField & other);

	void buildPorts();
	
	void		setGeometry(MimmoObject * geo);
	void		setSplittedGeometries( std::vector<MimmoObject *> list);
	void		setCellDivisionMap(std::unordered_map<long, std::pair<int, long> > map);
	void		setVertDivisionMap(std::unordered_map<long, std::pair<int, long> > map);
	
	bool 		isEmpty();
	void 		clear();
	int 		getTopo();
	void 		execute();
	
	virtual void absorbSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

private:
	virtual bool split()=0;
protected:	
	bitpit::VTKElementType desumeElement(ivector2D & conn );
};

/*!
 *	\date			30/nov/2016
 *	\authors		Rocco Arpa
 *	\authors		
 *
 *	\brief SplitScalarField is specialized derived class of SplitField to split a   
 *         scalar field of doubles.
 * 
 *	=========================================================
 * ~~~
 *	|-------------------------------------------------------------------------------|
 *	|                 Port Input                           					        |
 *	|-------|--------------|------------------------------------|-----------------------|
 *	|PortID | PortType     | variable/function 					| DataType		        |
 *	|-------|--------------|------------------------------------|-----------------------|
 *  | 19  	| M_SCALARFIELD| setField      					    | (VECTOR, FLOAT)	    |
 *	|-------|--------------|------------------------------------|-----------------------|
 *
 *
 *	|--------------------------------------------------------|-----------------------|
 *	|            Port Output                				 |                       |
 *	|-------|-----------|------------------------------------|-----------------------|
 *	|PortID | PortType  | variable/function 				 | DataType		      	 | 
 *	|-------|-----------|------------------------------------|-----------------------|
 *	| 106   | M_UMGEOSFD| getSplittedData					 | (UN_MAP, MIMMO_VECFLOAT_)|
 *  |-------|-----------|------------------------------------|-----------------------|
 * ~~~
 *	=========================================================
 *
 */

class SplitScalarField: public SplitField{
private:	
	dvector1D m_field;
	std::vector< dvector1D > m_result;

public:	
	SplitScalarField(int topo);
	virtual ~SplitScalarField();
	
	void buildPorts();
	std::unordered_map<MimmoObject*,dvector1D* >	 getSplittedData(); 
	void 	setField(dvector1D field);
	
	void clear();
	
	void 	plotOptionalResults();	
private:
	bool split();
};

/*!
 *	\date			30/nov/2016
 *	\authors		Rocco Arpa
 *	\authors		
 *
 *	\brief SplitVectorField is specialized derived class of SplitField to split a   
 *         scalar field of array<double,3>.
 * 
 *	=========================================================
 * ~~~
 *	|-------------------------------------------------------------------------------|
 *	|                 Port Input                           					        |
 *	|-------|--------------|------------------------------------|-----------------------|
 *	|PortID | PortType     | variable/function 					| DataType		        |
 *	|-------|--------------|------------------------------------|-----------------------|
 *  | 11  	| M_GDISPLS    | setField      					    | (VECARR3, FLOAT)	    |
 *	|-------|--------------|------------------------------------|-----------------------|
 *
 *
 *	|--------------------------------------------------------|-----------------------|
 *	|            Port Output                				 |                       |
 *	|-------|-----------|------------------------------------|-----------------------|
 *	|PortID | PortType  | variable/function 				 | DataType		      	 | 
 *	|-------|-----------|------------------------------------|-----------------------|
 *	| 107   | M_UMGEOVFD| getSplittedData					 | (UN_MAP, MIMMO_VECARR3FLOAT_)|
 *  |-------|-----------|------------------------------------|-----------------------|
 * ~~~
 *	=========================================================
 *
 */
class SplitVectorField: public SplitField{
private:	
	dvecarr3E m_field;
	std::vector< dvecarr3E > m_result;
public:	
	
	SplitVectorField(int topo);
	virtual ~SplitVectorField();
	
	void buildPorts();
	std::unordered_map< MimmoObject*,dvecarr3E* >	 getSplittedData(); 
	void 	setField(dvecarr3E field);
	
	void clear();
	
	void 	plotOptionalResults();	
private:
	bool split();
};

}

#endif /* __SPLITFIELDS_HPP__ */
