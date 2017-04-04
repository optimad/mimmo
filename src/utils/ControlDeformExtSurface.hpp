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
#ifndef __CONTROLDEFORMEXTSURFACE_HPP__
#define __CONTROLDEFORMEXTSURFACE_HPP__

#include "BaseManipulation.hpp"
#include "MimmoGeometry.hpp"

namespace mimmo{

/*!
 * \class ControlDeformExtSurface
 *
 * \brief ControlDeformExtSurface is a class that check a deformation field associated to a MimmoObject geometry,
 *   for eventual penetrations,  w.r.t. one or more external constraint surface meshes.
 *
 * ControlDeformExtSurface is derived from BaseManipulation class.
 * Needs one or more external surface meshes, representing the constraint of your deformed object. 
 * Returns a double value V, namely the maximum signed distance from constraint surfaces amongst all field points, 
 * reporting how much the current deformation field violate the constraint itself.
 * if V >0 a violation occurs. if V=0, a contact occurs, otherwise if V<0 no violation occurs. 
 * No optional result are plot. 
 * Class absorbs/flushes its parameters from/to xml dictionaries
 *
 *	=========================================================
 * ~~~
 *	|----------------------------------------------------------------|
 *	|                 Port Input                                     |
 *	|-------|----------|-------------------|-------------------------|
 *	|PortID | PortType | variable/function | DataType                |
 *	|-------|----------|-------------------|-------------------------|
 *	| 11    | M_GDISPLS| setDefField       | (VECARR3E, FLOAT)       |
 *	| 99    | M_GEOM   | setGeometry       | (SCALAR, MIMMO_)        |
 *	|-------|----------|-------------------|-------------------------|
 *
 *
 *	|-------------------------------------------------------------------------|
 *	|            Port Output                                                  |
 *	|-------|---------------|-------------------|-----------------------------|
 *	|PortID | PortType      | variable/function | DataType                    |
 *	|-------|---------------|-------------------|-----------------------------|
 *	| 19    | M_SCALARFIELD | getViolationField | (VECTOR, FLOAT)             | 
 *	| 30    | M_VALUED      | getViolation      | (SCALAR, FLOAT)             |
 *	| 82    | M_VIOLATION   | getViolationPair  | (PAIR,PAIRMIMMO_OBJFLOAT_)  |
 *	|-------|---------------|-------------------|-----------------------------|
 * ~~~
 *	=========================================================
 *
 */
class ControlDeformExtSurface: public BaseManipulation{
private:
	std::unordered_map<std::string, std::pair<double, int> > m_geolist; /**< list of file for geometrical proximity check*/
	dvector1D					m_violationField;	/**<Violation Field as distance from constraint */
	dvecarr3E					m_defField; 	/**<Deformation field*/
	int 						m_cellBackground; /**< Number of cells N to determine background grid spacing */ 						
	std::unordered_set<int>		m_allowed; /**< list of currently file format supported by the class*/
public:
	ControlDeformExtSurface();
	ControlDeformExtSurface(const bitpit::Config::Section & rootXML);
	~ControlDeformExtSurface();

	ControlDeformExtSurface(const ControlDeformExtSurface & other);
	ControlDeformExtSurface & operator=(const ControlDeformExtSurface & other);

	void	buildPorts();

	double 									getViolation();
	std::pair<BaseManipulation*, double>	getViolationPair();
	dvector1D								getViolationField();
	double 									getToleranceWithinViolation(std::string);
	int 									getBackgroundDetails();
	
	void	setDefField(dvecarr3E field);
	void 	setGeometry(MimmoObject * geo);
	void 	setBackgroundDetails(int nCell=50);
	const 	std::unordered_map<std::string, std::pair<double, int> > & 	getFiles() const;
	void	setFiles(std::unordered_map<std::string,std::pair<double, int> > list );
	void 	addFile(std::string file, double tol, int format);
	void 	addFile(std::string file, double tol, FileType format);
	void 	removeFile(std::string);
	void 	removeFiles();
	
	void clear();
	
	void 	execute();
	
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:	
	void plotOptionalResults();
	
private:
	void readGeometries(std::vector<std::unique_ptr<MimmoGeometry> > & extGeo, std::vector<double> & tols);
	svector1D extractInfo(std::string file);
	double evaluateSignedDistance(darray3E &point, mimmo::MimmoObject * geo, long & id, darray3E & normal, double &initRadius);
};

REGISTER(BaseManipulation, ControlDeformExtSurface,"mimmo.ControlDeformExtSurface")

}

#endif /* __CONTROLDEFORMEXTSURFACE_HPP__ */
