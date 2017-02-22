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
 \ *----------------*-----------------------------------------------------------*/
#ifndef __CREATESEEDSONSURFACE_HPP__
#define __CREATESEEDSONSURFACE_HPP__

#include "BaseManipulation.hpp"
#include "OBBox.hpp"

using namespace mimmo;
namespace mimmino{

/*!
 * Enum class for engine  choiche to set up initial points on a 3D surface.
 */	
enum class CSeedSurf{
	RANDOM = 0 /**< Engine type, sows randomly points on surface */,
	LEVELSET = 1 /**< Engine type, sows points around a seed on surface,using geodesic distance between points */,
	CARTESIANGRID=2 /**< Engine type, sows points projecting a 3D Cartesian grid on surface */
	
};

/*!
 * Class/BaseManipulation Object to position an initial set of points on a 3D surface.
 * Three type of engines to compute point position are available:
 * 0)	CSeedSurf::RANDOM : sows points randomly on your surface, trying to displace them
 * 		at maximum euclidean distance possible on the surface.
 * 1)	CSeedSurf::LEVELSET : starting from an initial seed, sows points around it, trying 
 * 		to displace them at maximum geodesic distance possible on the surface.
 * 2)	CSeedSurf::CARTESIANGRID	evaluate points by projection of a volumetric 
 * 		cartesian grid of surface and decimating the list up the desired value of points, 
 * 		trying to displace them at maximum euclidean distance possible on the surface
 * 
 * Default engine is CARTESIANGRID
 * PORTS AVAILABLE IN CreateSeedsOnSurface Class (legend M_<> MiMMO native ports, C_<> CAMiLO own ports)
 *
 *	=========================================================
 * ~~~
 *	|----------------------------------------------------------------------------------------|
 *	|                   Port Input                                                           |
 *	|-------|----------------|---------------------------------------|-----------------------|
 *	|PortID | PortType       | variable/function                     | DataTypes	         |
 *	|-------|----------------|---------------------------------------|-----------------------|
 *	| 20    | M_POINT        | setSeed		                         | (ARRAY3, FLOAT)       |
 *	| 31    | M_VALUEI       | setNPoints	                         | (SCALAR, INT)	     | 
 *	| 32    | M_VALUEB       | setMassCenterAsSeed                   | (SCALAR, BOOL)	     | 
 *	| 99    | M_GEOM         | setGeometry                           | (SCALAR, MIMMO_)      | 
 *	| 1004  | C_SEEDSURFENG  | setEngine	                         | (SCALAR, INT)	     | 
 *	| 1900  | C_VALUEI       | setRandomFixed                        | (SCALAR, INT)	     |
 *	|-------|----------------|---------------------------------------|-----------------------|
 *
 *
 *	|-----------------------------------------------------------------------|
 *	|             Port Output                     							|
 *	|-------|----------------|--------------------|-------------------------|
 *	|PortID | PortType       | variable/function  | DataTypes	         	|
 *	|-------|----------------|--------------------|-------------------------|
 *	| 0     | M_COORDS  	 | getPoints	      | (VECARR3E, FLOAT)     	|
 *	|-------|----------------|--------------------|-------------------------|
 *
 * ~~~
 *	=========================================================
 *
 */

class CreateSeedsOnSurface: public mimmo::BaseManipulation {
	
private:
	
	dvecarr3E 	m_points;		/**< resulting points of class computation */
	int 	 	m_nPoints;		/**< total number of desired points */
	double 		m_minDist;		/**< minimum distance of tolerance */
	darray3E	m_seed;			/**< inital seed point */
	CSeedSurf   m_engine;		/**< choose kernel type for points positioning computation */ 
	bool		m_seedbaricenter; /**< bool activate mass center as starting seed */
	int 		m_randomFixed;    /**< signature for freezing random engine result*/
	//utility members
	std::unique_ptr<mimmo::OBBox> bbox;		/**<pointer to an oriented Bounding box */
	ivector1D m_deads; /**! service structure */
	
public:
	
	CreateSeedsOnSurface();
	virtual ~CreateSeedsOnSurface();
	CreateSeedsOnSurface(const CreateSeedsOnSurface & other);
	CreateSeedsOnSurface & operator=(const CreateSeedsOnSurface & other);
	
	void	buildPorts();

	//get methods
	int			getNPoints();
	dvecarr3E	getPoints();
	CSeedSurf	getEngineENUM();
	int			getEngine();
	darray3E	getSeed();
	bool		isSeedMassCenter();
	double 		getMinDistance();
	int	    	getRandomSignature();
	
	//set methods
	void		setNPoints( int);
	void		setEngineENUM( CSeedSurf);
	void		setEngine(int);
	void		setSeed(darray3E);
	void		setMassCenterAsSeed(bool );
	void		setGeometry(MimmoObject *);
	void		setRandomFixed(int signature = -1);
	
	//cleaners
	void clear();
	
	//appliers
	void 	solve(bool debug = false);
	// plotting
	void 	plotCloud(std::string dir, std::string file, int counter, bool binary);
	//execute
	void		execute();
	
	//XML utilities from reading writing settings to file
	virtual void absorbSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	
protected:
	virtual void plotOptionalResults();
	
private:

	void	solveLSet( bool debug = false);
	void	solveGrid(bool debug = false);
	void	solveRandom(bool debug = false);
	
	dvecarr3E decimatePoints(dvecarr3E &);
	
	void solveEikonal(double g, double s, std::unordered_map<long,long> & invConn, dvector1D & field);
	double updateEikonal(double g, double s, long tVert,long tCell, std::unordered_map<long int, short int> &flag, dvector1D & field);
	
	std::unordered_map<long,long>	getInverseConn();
	bool			isDeadFront(const int label);
	std::set<long>	findVertexVertexOneRing(const long &, const long & ); 
};	
	
	
};

#endif /* __CREATESEEDSONSURFACE_HPP__ */
