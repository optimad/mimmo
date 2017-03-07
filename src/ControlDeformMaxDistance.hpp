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
#ifndef __CONTROLDEFORMMAXDISTANCE_HPP__
#define __CONTROLDEFORMMAXDISTANCE_HPP__

#include "BaseManipulation.hpp"
#include "MimmoObject.hpp"

namespace mimmo{

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief ControlDeformationMaxDistance is a class that check a deformation field associated to a MimmoObject geometry,
 * 		once a maximum limit distance of deformation is fixed, w.r.t. the undeformed state.
 *
 *	ControlDeformationMaxDistance is derived from BaseManipulation class.
 *	It needs a maximum, isotropic limit distance d w.r.t. geometry undeformed state, which is used to evaluate the isolevel d
 *  of the target geometry. 
 *	Returns a double value V, namely the maximum signed distance from constraint iso-level amongst all field points, 
 *  reporting how much the current deformation field violate the constraint itself.
 *  if V >0 a violation occurs. if V=0, a contact occurs, otherwise if V<0 no violation occurs. 
 *  No optional result are plot. 
 *  Class absorbs/flushes its parameters from/to xml dictionaries
 *
 *	=========================================================
 * ~~~
 *	|----------------------------------------------------------------|
 *	|                 Port Input                       		 		 |
 *	|-------|----------|-------------------|-------------------------|
 *	|PortID | PortType | variable/function | DataType 		 		 |
 *	|-------|----------|-------------------|-------------------------|
 *  | 11    | M_GDISPLS| setDefField	   | (VECARR3E, FLOAT)		 |
 *  | 30    | M_VALUED | setLimitDistance  | (SCALAR, FLOAT)		 |
 *	| 99    | M_GEOM   | setGeometry       | (SCALAR, MIMMO_)		 |
 *	|-------|----------|-------------------|-------------------------|
 *
 *
 *	|-------------------------------------------------------------------------|
 *	|            Port Output                                                  |
 *	|-------|---------------|-------------------|-----------------------------|
 *	|PortID | PortType      | variable/function | DataType 			          |
 *	|-------|---------------|-------------------|-----------------------------|
 *	| 19    | M_SCALARFIELD | getViolationField | (VECTOR, FLOAT)			  | 
 *	| 30    | M_VALUED 		| getViolation      | (SCALAR, FLOAT)			  |
 *  | 82    | M_VIOLATION   | getViolationPair  | (PAIR,PAIRMIMMO_OBJFLOAT_)  |
 *	|-------|---------------|-------------------|-----------------------------|
 * ~~~
 *	=========================================================
 *
 */
class ControlDeformMaxDistance: public BaseManipulation{
private:
	double						m_maxDist;		/**<Limit Distance*/
	dvector1D					m_violationField;	/**<Violation Distance Field */
	dvecarr3E					m_defField; 	/**<Deformation field*/
	
public:
	ControlDeformMaxDistance();
	ControlDeformMaxDistance(const bitpit::Config::Section & rootXML);
	~ControlDeformMaxDistance();

	ControlDeformMaxDistance(const ControlDeformMaxDistance & other);
	ControlDeformMaxDistance & operator=(const ControlDeformMaxDistance & other);

	void	buildPorts();

	double 									getViolation();
	dvector1D								getViolationField();
	std::pair<BaseManipulation*, double>	getViolationPair();
	
	void	setDefField(dvecarr3E field);
	void	setLimitDistance(double dist);

	void 	execute();
	
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    void plotOptionalResults();	

	
};

REGISTER(BaseManipulation, ControlDeformMaxDistance, "MiMMO.ControlDeformMaxDistance")
}

#endif /* __CONTROLDEFORMMAXDISTANCE_HPP__ */
