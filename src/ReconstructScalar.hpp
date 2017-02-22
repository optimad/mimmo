/*---------------------------------------------------------------------------*\
 * 
 *  CAMILO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License Commercial (//TODO Temporary header of license)
 *  This file is part of CAMILO.
 *
 *  CAMILO is a commercial software: you do not own rights to redistribute it 
 * 	and/or modify it both in source or pre-build formats
 *  Please contact Optimad offices for any further informations				
 *
 *  You should have received a copy of the Camilo Commercial License
 *  along with CAMILO, as well as the key to unlock the software.
 *
 \*----------------*-----------------------------------------------------------*/

#ifndef __RECONSTRUCTSCALAR_HPP__
#define __RECONSTRUCTSCALAR_HPP__

#include "BaseManipulation.hpp"
#include "CAMILO_Ports.hpp"
#include <unordered_map>
#include <utility>

namespace mimmino{

/*!
 * @enum OverlapMethod
 * @brief class for setting overlapping criterium for two different scalar fields:
 * 			1) MAX = get max(with sign) between available field values 
 * 			2) MIN = get min(with sign) between available field values 
 * 			3) AVERAGE = get simple average between available values
 * 			4) SUM = take sum of both values between overlapped fields
 */	
enum class OverlapMethod{
	MAX = 1 /**< take max values between overlapped fields*/,
	MIN = 2 /**< take min values between overlapped fields*/,
	AVERAGE = 3 /**< take averaged values between overlapped fields*/,
	SUM = 4 /**< take sum of both values between overlapped fields*/
};

/*!
 * Class/BaseManipulation Object reconstructing a scalar field on a mimmo::MimmoObjectmesh, from several
 * scalar fields defined on sub-patch of the target mesh. Field values are defined on nodes.
 * Reconstructed field is provided in m_result member of the class.
 * 
 * PORTS AVAILABLE IN ReconstructScalar Class (legend M_<> MiMMO native ports, C_<> CAMiLO own ports)
 * 
 *	=========================================================
 * ~~~
 *	|---------------------------------------------------------------------------------------------------|
 *	|                   Port Input                                                        				|
 *	|-------|----------------|---------------------------------------|----------------------------------|
 *	|PortID | PortType       | variable/function                     | DataTypes						|
 *	|-------|----------------|---------------------------------------|----------------------------------|
 *	| 81    | M_PAIRSCAFIELD | setData                               | (PAIR, MIMMO_VECFLOAT_)			|
 *	| 99    | M_GEOM         | m_geometry                            | (SCALAR, MIMMO_)					|
 *	| 1001  | C_OVERLAPMTH   | setOverlapCriterium                   | (SCALAR, INT)					|
 *  | 1002	| C_VECPAIRSF    | setData	    						 | (VECTOR, PAIRMIMMO_VECFLOAT_)	|
 *	|-------|----------------|---------------------------------------|----------------------------------|
 * 
 *
 *
 *	|-------------------------------------------------------------------|
 *	|             Port Output                                         	|
 *	|-------|----------------|--------------------|---------------------|
 *	|PortID | PortType       | variable/function  | DataTypes			|
 *	|-------|----------------|--------------------|---------------------|
 *	| 19    | M_SCALARFIELD  | getResultField     | (VECTOR, FLOAT)		|
 *  | 99    | M_GEOM         | getGeometry        | (SCALAR, MIMMO_)	|
 *	| 1001  | C_OVERLAPMTH   | getOverlapCriterium| (SCALAR, INT)		|
 *	|-------|----------------|--------------------|---------------------|
 *
 * ~~~
 *	=========================================================
 */

class ReconstructScalar: public mimmo::BaseManipulation {
	
private:
	
	OverlapMethod m_overlapCriterium;
	std::unordered_map < mimmo::MimmoObject*, dvector1D * > m_subpatch;
	dvector1D	m_result;
	
public:
	typedef std::pair<mimmo::MimmoObject*, dvector1D *>  pField; /**< Internal definition of the class */
	ReconstructScalar();
	virtual ~ReconstructScalar();
	ReconstructScalar(const ReconstructScalar & other);
	ReconstructScalar & operator=(const ReconstructScalar & other);
	
	//get-set methods
	OverlapMethod			getOverlapCriteriumENUM();
	int 					getOverlapCriterium();
	dvector1D 				getData(mimmo::MimmoObject * patch );
	int 					getNData();
	dvector1D				getResultField();
	
	std::vector< mimmo::MimmoObject	* >	whichSubMeshes();
	
	void 		setOverlapCriteriumENUM( OverlapMethod);
	void 		setOverlapCriterium( int);
	void		setData( pField );
	void 		setData(std::vector<pField> );
	void		removeData(mimmo::MimmoObject* );
	void		removeAllData();
	void		buildPorts();
	//cleaners
	void clear();
	
	//plotting

	void	plotData(std::string dir, std::string name, bool flag, int counter);
	
	//execute
	void		execute();
	
	//XML utilities from reading writing settings to file
	virtual void absorbSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
protected:
		virtual void plotOptionalResults();
	
private:
	double 	overlapFields(dvector1D & locField);
	std::unordered_map<long, dvector1D>		checkOverlapping();
	bitpit::VTKElementType	desumeElement(int typeGeom, ivector2D & conn); //TODO need to be moved in MimmoObject.
};	
	
/*!
 * Class/BaseManipulation Object reconstructing a vector field on a mimmo::MimmoObject mesh, from several
 * vector fields defined on sub-patch of the target mesh. Field values are defined on nodes.
 * Reconstructed field is provided in m_result member of the class.
 * 
 * PORTS AVAILABLE IN ReconstructScalar Class (legend M_<> MiMMO native ports, C_<> CAMiLO own ports)
 * 
 *	=========================================================
 * ~~~
 *	|-------------------------------------------------------------------------------------------------------|
 *	|                   Port Input                                                        					|
 *	|-------|----------------|---------------------------------------|--------------------------------------|
 *	|PortID | PortType       | variable/function                     | DataTypes							|
 *	|-------|----------------|---------------------------------------|--------------------------------------|
 *	| 80    | M_PAIRVECFIELD | setData                               | (PAIR, MIMMO_VECARR3FLOAT_)			|
 *	| 99    | M_GEOM         | m_geometry                            | (SCALAR, MIMMO_)						|
 *	| 1001  | C_OVERLAPMTH   | setOverlapCriterium                   | (SCALAR, INT)						|
 *  | 1003	| C_VECPAIRVF    | setData	    						 | (VECTOR, PAIRMIMMO_VECARR3FLOAT_)	|
 *	|-------|----------------|---------------------------------------|--------------------------------------|
 * 
 *
 *
 *	|-------------------------------------------------------------------|
 *	|             Port Output                                         	|
 *	|-------|----------------|--------------------|---------------------|
 *	|PortID | PortType       | variable/function  | DataTypes			|
 *	|-------|----------------|--------------------|---------------------|
 *	| 11    | M_GDISPL       | getResultField     | (VECARR3, FLOAT)	|
 *  | 99    | M_GEOM         | getGeometry        | (SCALAR, MIMMO_)	|
 *	| 1001  | C_OVERLAPMTH   | getOverlapCriterium| (SCALAR, INT)		|
 *	|-------|----------------|--------------------|---------------------|
 * 
 * ~~~
 *	=========================================================
 */

class ReconstructVector: public mimmo::BaseManipulation {
	
private:
	
	OverlapMethod m_overlapCriterium;
	std::unordered_map < mimmo::MimmoObject*, dvecarr3E * > m_subpatch;
	dvecarr3E	m_result;

public:
	typedef std::pair<mimmo::MimmoObject*, dvecarr3E *>  pVector; /**< Internal definition of the class */
	ReconstructVector();
	virtual ~ReconstructVector();
	ReconstructVector(const ReconstructVector & other);
	ReconstructVector & operator=(const ReconstructVector & other);
	
	//get-set methods
	OverlapMethod			getOverlapCriteriumENUM();
	int						getOverlapCriterium();
	dvecarr3E 				getData(mimmo::MimmoObject * patch);
	int 					getNData();
	dvecarr3E				getResultField();
	
	std::vector< mimmo::MimmoObject	* >	whichSubMeshes();
	
	void 		setOverlapCriteriumENUM( OverlapMethod);
	void		setOverlapCriterium(int );
	void		setData( pVector );
	void 		setData(std::vector<pVector> );
	void		removeData(mimmo::MimmoObject* );
	void		removeAllData();
	void		buildPorts();
	//cleaners
	void clear();
	

	//plotting
	void	plotData(std::string dir, std::string name, bool flag, int counter);
	
	//execute
	void		execute();

	//XML utilities from reading writing settings to file
	virtual void absorbSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	
protected:
	virtual void plotOptionalResults();
	
private:
	darray3E	overlapFields(dvecarr3E & locField);
	std::unordered_map<long, dvecarr3E>		checkOverlapping();
	bitpit::VTKElementType	desumeElement(int typeGeom, ivector2D & conn); //TODO need to be moved in MimmoObject.
};	

};

#endif /* __RECONSTRUCTSCALAR_HPP__ */
