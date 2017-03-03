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
 \ *---------------------------------------------------------------------------*/

#ifndef __OVERLAPPINGFIELDS_HPP__
#define __OVERLAPPINGFIELDS_HPP__

#include "ReconstructFields.hpp"

namespace mimmo{

/*!
 * Class/BaseManipulation Object overlapping one or more double scalar fields associated to a MimmoObject geometry.
 * The class handles more possible overlapping referring to different geoemtries at the same time. Returns a list of
 * overlapped fields associated to their geometry. 
 * 
 * PORTS AVAILABLE IN OverlapScalarFields Class 
 * 
 *	=========================================================
 * ~~~
 *	|---------------------------------------------------------------------------------------------------|
 *	|                   Port Input                                                        				|
 *	|-------|----------------|---------------------------------------|----------------------------------|
 *	|PortID | PortType       | variable/function                     | DataTypes						|
 *	|-------|----------------|---------------------------------------|----------------------------------|
 *	| 81    | M_PAIRSCAFIELD | setAddDataField                       | (PAIR, MIMMO_VECFLOAT_)			|
 *	| 106   | M_UMGEOSFD     | setDataFieldMap                       | (UMAP, MIMMO_VECFLOAT_)			|
 *  | 200	| M_VECPAIRSF    | setDataFieldList						 | (VECTOR, PAIRMIMMO_VECFLOAT_)	|
 *	|-------|----------------|---------------------------------------|----------------------------------|
 * 
 *
 *
 *	|----------------------------------------------------------------------------|
 *	|             Port Output                                         	         |
 *	|-------|----------------|--------------------|------------------------------|
 *	|PortID | PortType       | variable/function  | DataTypes			         |
 *	|-------|----------------|--------------------|------------------------------|
 *	| 106   | M_UMGEOSFD     | getDataFieldMap    | (UMAP, MIMMO_VECFLOAT_)      |
 *  | 200   | M_VECPAIRSF    | getDataFieldList   | (VECTOR, PAIRMIMMO_VECFLOAT_)|
 *	|-------|----------------|--------------------|------------------------------|
 *
 * ~~~
 *	=========================================================
 */

class OverlapScalarFields: public mimmo::BaseManipulation {
	
private:
	
	OverlapMethod m_overlapCriterium;
	std::unordered_map < mimmo::MimmoObject*, std::vector<dvector1D *> > m_originals;
	std::unordered_map < mimmo::MimmoObject*, dvector1D > m_results;
	
public:
	OverlapScalarFields();
	OverlapScalarFields(const bitpit::Config::Section & rootXML);
	virtual ~OverlapScalarFields();
	OverlapScalarFields(const OverlapScalarFields & other);
	OverlapScalarFields & operator=(const OverlapScalarFields & other);
	
	//get-set methods
	OverlapMethod			getOverlapCriteriumENUM();
	int 					getOverlapCriterium();
	dvector1D 				getResultData(mimmo::MimmoObject * patch );
	int 					getNEffectiveFields();
	int						getNLinkedFields();
	
	std::vector<MimmoObject*>							whichGeometriesLinked();
	std::unordered_map<MimmoObject*, dvector1D* >		getDataFieldMap();
	std::vector<std::pair<MimmoObject*, dvector1D* > >	getDataFieldList();
	
	void 		setOverlapCriteriumENUM( OverlapMethod);
	void 		setOverlapCriterium( int);
	void		setAddDataField( std::pair<MimmoObject*, dvector1D*> field );
	void 		setDataFieldMap(std::unordered_map<MimmoObject*, dvector1D*> fieldMap );
	void 		setDataFieldList(std::vector<std::pair<MimmoObject*, dvector1D*> > fieldList );

	void		removeData(mimmo::MimmoObject* );
	void		removeAllData();

	void		buildPorts();
	//cleaners
	void clear();
	
	//plotting
	void	plotData(std::string dir, std::string name, bool flag, int counter, MimmoObject * geo);
	void	plotAllData(std::string dir, std::string name, bool flag);
	
	//execute
	void		execute();
	
	//XML utilities from reading writing settings to file
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:

	virtual void plotOptionalResults();
	
private:
	double 	overlapFields(dvector1D & locField);
	bitpit::VTKElementType	desumeElement(int typeGeom, ivector2D & conn); //TODO need to be moved in MimmoObject.
};	
	
/*!
 * Class/BaseManipulation Object overlapping one or more 3 double elements vector fields associated to a MimmoObject geometry.
 * The class handles more possible overlapping referring to different geoemtries at the same time. Returns a list of
 * overlapped fields associated to theri geometry. 
 * 
 * PORTS AVAILABLE IN OverlapVectorFields Class 
 * 
 *	=========================================================
 * ~~~
 *	|----------------------------------------------------------------------------------------------------|
 *	|                   Port Input                                                        				 |
 *	|-------|----------------|---------------------------------------|-----------------------------------|
 *	|PortID | PortType       | variable/function                     | DataTypes						 |
 *	|-------|----------------|---------------------------------------|-----------------------------------|
 *	| 80    | M_PAIRVECFIELD | setAddDataField                       | (PAIR, MIMMO_VECARR3FLOAT_)		 |
 *	| 107   | M_UMGEOSVD     | setDataFieldMap                       | (UMAP, MIMMO_VECARR3FLOAT_)		 |
 *  | 201	| M_VECPAIRVF    | setDataFieldList						 | (VECTOR, PAIRMIMMO_VECARR3EFLOAT_)|
 *	|-------|----------------|---------------------------------------|-----------------------------------|
 * 
 *
 *
 *	|---------------------------------------------------------------------------------|
 *	|             Port Output                                         	              |
 *	|-------|----------------|--------------------|-----------------------------------|
 *	|PortID | PortType       | variable/function  | DataTypes			              |
 *	|-------|----------------|--------------------|-----------------------------------|
 *	| 107   | M_UMGEOVFD     | getDataFieldMap    | (UMAP, MIMMO_VECARR3EFLOAT_)      |
 *  | 201   | M_VECPAIRVF    | getDataFieldList   | (VECTOR, PAIRMIMMO_VECARR3EFLOAT_)|
 *	|-------|----------------|--------------------|-----------------------------------|
 *
 * ~~~
 *	=========================================================
 */

class OverlapVectorFields: public mimmo::BaseManipulation {
	
private:
	
	OverlapMethod m_overlapCriterium;
	std::unordered_map < mimmo::MimmoObject*, std::vector<dvecarr3E *> > m_originals;
	std::unordered_map < mimmo::MimmoObject*, dvecarr3E > m_results;
	
public:
	OverlapVectorFields();
	OverlapVectorFields(const bitpit::Config::Section & rootXML);
	virtual ~OverlapVectorFields();
	OverlapVectorFields(const OverlapVectorFields & other);
	OverlapVectorFields & operator=(const OverlapVectorFields & other);
	
	//get-set methods
	OverlapMethod			getOverlapCriteriumENUM();
	int 					getOverlapCriterium();
	dvecarr3E 				getResultData(mimmo::MimmoObject * patch );
	int 					getNEffectiveFields();
	int						getNLinkedFields();
	
	std::vector<MimmoObject*>							whichGeometriesLinked();
	std::unordered_map<MimmoObject*, dvecarr3E* >		getDataFieldMap();
	std::vector<std::pair<MimmoObject*, dvecarr3E* > >	getDataFieldList();
	
	void 		setOverlapCriteriumENUM( OverlapMethod);
	void 		setOverlapCriterium( int);
	void		setAddDataField( std::pair<MimmoObject*, dvecarr3E*> field );
	void 		setDataFieldMap(std::unordered_map<MimmoObject*, dvecarr3E*> fieldMap );
	void 		setDataFieldList(std::vector<std::pair<MimmoObject*, dvecarr3E*> > fieldList );
	
	void		removeData(mimmo::MimmoObject* );
	void		removeAllData();
	
	void		buildPorts();
	//cleaners
	void clear();
	
	//plotting
	void	plotData(std::string dir, std::string name, bool flag, int counter, MimmoObject * geo);
	void	plotAllData(std::string dir, std::string name, bool flag);
	
	//execute
	void		execute();
	
	//XML utilities from reading writing settings to file
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	
protected:
	
	virtual void plotOptionalResults();
	
private:
	darray3E 	overlapFields(dvecarr3E & locField);
	std::unordered_map<long, dvecarr3E>		checkOverlapping();
	bitpit::VTKElementType	desumeElement(int typeGeom, ivector2D & conn); //TODO need to be moved in MimmoObject.
};	

};

#endif /* __OVERLAPPINGFIELDS_HPP__ */
