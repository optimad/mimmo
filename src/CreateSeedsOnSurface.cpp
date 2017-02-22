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
 \ *---------------------------------------------------------------------------*/

#include "CreateSeedsOnSurface.hpp"
#include "RCPoints.hpp"
#include "Lattice.hpp"

#include <stdlib.h>
#include <time.h>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <iostream>
#include <fstream>
#include "CAMILO_Ports.hpp"

using namespace mimmino;

//PUBLIC METHODS
/*!
 * Constructor
 */
CreateSeedsOnSurface::CreateSeedsOnSurface(){
	m_name = "MiMMiNO.CreateSeedsOnSurface";
	m_nPoints = 0;
	m_minDist = 0.0;
	m_seed = {{0.0,0.0,0.0}};
	m_engine = CSeedSurf::CARTESIANGRID;
	m_seedbaricenter = false;
	m_randomFixed = -1;
	std::unique_ptr<mimmo::OBBox> box(new mimmo::OBBox());
	bbox = std::move(box);
};

/*!
 * Destructor;
 */
CreateSeedsOnSurface::~CreateSeedsOnSurface(){
	clear();
};

/*!
 * Copy constructor
 */
CreateSeedsOnSurface::CreateSeedsOnSurface(const CreateSeedsOnSurface & other){
	*this=other;
};

/*!
 * Copy operator of the class
 */
CreateSeedsOnSurface & CreateSeedsOnSurface::operator=(const CreateSeedsOnSurface & other){
	
	*(static_cast<BaseManipulation *> (this)) = *(static_cast<const BaseManipulation *> (&other));
	m_points = other.m_points;
	m_nPoints = other.m_nPoints;
	m_minDist = other.m_minDist;
	m_seed = other.m_seed;
	m_engine = other.m_engine;
	m_seedbaricenter = other.m_seedbaricenter;
	m_randomFixed = other.m_randomFixed;
	m_deads = other.m_deads;
	return(*this);
};


/*!
 * It builds the input/output ports of the object
 */
void CreateSeedsOnSurface::buildPorts(){

	bool built = true;

	//input
	built = (built && createPortIn<MimmoObject *, CreateSeedsOnSurface>(this, &mimmino::CreateSeedsOnSurface::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	built = (built && createPortIn<darray3E, CreateSeedsOnSurface>(this, &mimmino::CreateSeedsOnSurface::setSeed, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<int, CreateSeedsOnSurface>(this, &mimmino::CreateSeedsOnSurface::setNPoints, PortType::M_VALUEI, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::INT));
	built = (built && createPortIn<int, CreateSeedsOnSurface>(this, &mimmino::CreateSeedsOnSurface::setRandomFixed, CAMILOPortType::C_VALUEI, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::INT));
	built = (built && createPortIn<bool, CreateSeedsOnSurface>(this, &mimmino::CreateSeedsOnSurface::setMassCenterAsSeed, PortType::M_VALUEB, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::BOOL));
	built = (built && createPortIn<int, CreateSeedsOnSurface>(this, &mimmino::CreateSeedsOnSurface::setEngine, CAMILOPortType::C_SEEDSURFENG, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::INT));

	//output
	built = (built && createPortOut<dvecarr3E, CreateSeedsOnSurface>(this, &mimmino::CreateSeedsOnSurface::getPoints, PortType::M_COORDS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	
	m_arePortsBuilt = built;
};

/*!
 * Return number of points to be distributed over 3D surface.
 * See setNPoints documentation. 
 */
int	CreateSeedsOnSurface::getNPoints(){
	return m_nPoints;
};

/*!
 * Return list of points distributed over 3D surface, after class execution.
 */
dvecarr3E	CreateSeedsOnSurface::getPoints(){
	return m_points;
};

/*!
 * Return engine type for point distribution calculation.
 * ENUM overloading.
 */
CSeedSurf	CreateSeedsOnSurface::getEngineENUM(){
	return m_engine;
};


/*!
 * Return engine type for point distribution calculation.
 * 
 */
int	CreateSeedsOnSurface::getEngine(){
	return static_cast<int>(m_engine);
};

/*!
 * Return seed point used for calculation. If geometry baricenter is
 * used as initial seed (see setMassCenterAsSeed), it can be 
 * correctly visualized after the execution of the object.
 */
darray3E	CreateSeedsOnSurface::getSeed(){
	return m_seed;
}

/*!
 * Return true, if the option to use geometry baricenter as initial seed is 
 * activated. (see setMassCenterAsSeed method documentation).
 */
bool	CreateSeedsOnSurface::isSeedMassCenter(){
	return m_seedbaricenter;
};


/*!
 * Return minimum absolute distance during the point distribution calculation.
 */
double	CreateSeedsOnSurface::getMinDistance(){
	return m_minDist;
}

/*!
 * Return the signature of the current ramdom distribution of points on target surface,
 * whenever is fixed or not, for result replication. See setRandomFixed method.  This option will
 * only make sense if a CSeedSurf::RANDOM engine is employed. Otherwise it is ignored. 
 * \return signature. 
 */
int	CreateSeedsOnSurface::getRandomSignature(){
	return m_randomFixed;
}

/*!
 * Set the number of points to be distributed.
 * \param[in]	val	number of total points
 */
void	CreateSeedsOnSurface::setNPoints(int val){
	m_nPoints = std::max(val,0);
}

/*!
 * Set the engine type for point distribution evaluation
 * \param[in]	eng engine type passed as enum CSeedSurf
 */
void	CreateSeedsOnSurface::setEngineENUM(CSeedSurf eng){
	m_engine = eng;
}

/*!
 * Set the engine type for point distribution evaluation
 * \param[in]	eng engine type passed as integer. See CSeedSurf enum.
 */
void	CreateSeedsOnSurface::setEngine(int eng){
	if(eng <0 && eng >2)	eng = 2;
	setEngineENUM(static_cast<CSeedSurf>(eng));
}

/*!
 * Set the starting seed point for points distribution on surface evaluation.
 * The assigned seed is ignored if the geometry baricenter is set as initial seed
 * (see setMassCenterAsSeed).
 * \param[in]	seed starting 3D point in the neighborhood of the surface
 */
void	CreateSeedsOnSurface::setSeed(darray3E seed){
	m_seed = seed;
}

/*!
 * Set the starting seed point for points distribution on surface evaluation as the mass center
 * of the surface itself. Any seed passed with setSeed method will be ignored if this option is active.
 * \param[in]	flag true activate the option, false deactivate it.
 */
void	CreateSeedsOnSurface::setMassCenterAsSeed(bool flag){
	m_seedbaricenter = flag;
}

/*!
 * Set Geometry, check if its search bvtree is built, if not, build it. 
 * Point Cloud geometries, or pure Volume meshes are currently not supported by this block.
 * Reimplemented from BaseManipulation::setGeometry()
 */
void	CreateSeedsOnSurface::setGeometry(MimmoObject * geo){
	if(geo == NULL)	return;
	if(geo->getType() != 1)	return;
	
	BaseManipulation::setGeometry(geo);
	if(!geo->isBvTreeBuilt())			getGeometry()->buildBvTree();
	if(!geo->areAdjacenciesBuilt() )	getGeometry()->buildAdjacencies();
	bbox->setGeometry(geo);
}

/*!
 * Set the signature (each integer >= 0) of your random distribution. Same signature will be able to reproduce 
 * the exact random distribution in multiple runs. If signature is < 0 (default), point distribution will randomly 
 * vary run by run.It is possible to get the current signature after each random execution using the getRandomSignature method.
 * This option will only make sense if a CSeedSurf::RANDOM engine is employed. Otherwise it is ignored.
 *\param[in] signature integer number   
 */
void	CreateSeedsOnSurface::setRandomFixed( int signature){
	 m_randomFixed = signature;
}

/*!
 * Clear contents of the class
 */
void 	CreateSeedsOnSurface::clear(){
	m_points.clear();
	m_nPoints = 0;
	m_minDist = 0.0;
	m_seed = {{0.0,0.0,0.0}};
	m_engine = CSeedSurf::CARTESIANGRID;
	m_seedbaricenter = false;
	m_randomFixed = -1;
	m_deads.clear();
	
}

/*!
 * Execution command of the class. Wrapper to solve( bool ) method
 */
void 	CreateSeedsOnSurface::execute(){
	solve(false);
}

/*!
 * Plot Optional results of the class, that is the point cloud distribution on surface
 * as a list of xyz data points w/ identifier mark and as a point cloud in *.vtu format
 */
void 	CreateSeedsOnSurface::plotOptionalResults(){
	
	std::string dir = m_outputPlot;
	std::string nameGrid  = m_name + "CLOUD";
	plotCloud(dir, nameGrid, getClassCounter(), true );
}

/*!
 * Apply command of the class;
 */
void  CreateSeedsOnSurface::solve(bool debug){
	
	if(getGeometry() == NULL || m_nPoints< 1){
		if(debug)	std::cout<<"No geometry linked, or not enough total point set in "<<m_name<<" object. Doing Nothing"<<std::endl;
		return;
	}
	m_points.clear();
	bbox->execute();
	if(m_seedbaricenter)	m_seed = bbox->getOrigin();
	
	switch(m_engine){
		case CSeedSurf::RANDOM :
			solveRandom(debug);
			break;
		
		case CSeedSurf::LEVELSET :
			solveLSet(debug);
			break;

		case CSeedSurf::CARTESIANGRID :
			solveGrid(debug);
			break;
		
		default: //never been reached
		break;
	}
};

/*!
 * Find your optimal distribution starting from a seed point and calculating geodesic distance from point of each 
 * triangulated surface node. Add the most distant point from seed to list of candidates, thus update geodesic distance field from the two points, 
 * and repeat the process up the desired number of candidates.
 *\param[in]	debug	flag to activate logs of solver execution   
 */
void 	CreateSeedsOnSurface::solveLSet(bool debug){
	
	if(debug)	std::cout<<m_name<<" : started LevelSet engine"<<std::endl;
	dvecarr3E initList; 
	m_deads.reserve(m_nPoints);
	
	getGeometry()->getPatch()->buildAdjacencies();
	//find the nearest point of triagulation to the seed
	if(!(getGeometry()->isKdTreeBuilt())) getGeometry()->buildKdTree();
	
	bitpit::SurfUnstructured * tri = static_cast<bitpit::SurfUnstructured * >(getGeometry()->getPatch()); 
	auto map    = getGeometry()->getMapData();
	
	double distance = 0.0;
	for(auto &cell : tri->getCells()){
		distance += tri->evalCellArea(cell.getId());
	}
	distance /= double(tri->getCellCount());
	distance = std::pow(distance, 0.5);
	
	livector1D neighs, excl;
	darray3E projSeed = bvTreeUtils::projectPoint(&m_seed, getGeometry()->getBvTree());
	bitpit::Vertex vertSeed(0, projSeed);
	int nSize = 0;
	while( nSize < 1){
		getGeometry()->getKdTree()->hNeighbors(&vertSeed, distance, &neighs, & excl );   
		nSize = neighs.size();
		distance *= 1.1;
	}	
	
	long candidate = neighs[0];
	double minDist = norm2(projSeed - tri->getVertexCoords(candidate));
	for(int i=1; i<nSize; ++i){
		double val = norm2(projSeed - tri->getVertexCoords(neighs[i]));
		if(val < minDist){
			candidate = neighs[i];
			minDist = val;
		}
	}
	
	m_deads.push_back(getGeometry()->getMapDataInv(candidate));
	int deadSize = m_deads.size();
	if(debug)	std::cout<<m_name<<" : projected seed point"<<std::endl;

	std::unordered_map<long,long> invConn = getInverseConn();
	if(debug)	std::cout<<m_name<<" : created geometry inverse connectivity"<<std::endl;
	   
	while(deadSize < m_nPoints){
		
		dvector1D field(tri->getVertexCount(), 1.0E18);
		for(auto & dd : m_deads)	field[dd] = 0.0;
		
		solveEikonal(1.0,1.0, invConn, field);
		
		double maxField;
		maxval(field, maxField);
		dvector1D::iterator itF = std::find(field.begin(), field.end(), maxField);
		
		m_deads.push_back(std::distance(field.begin(), itF));

		deadSize = m_deads.size();
		if(debug)	std::cout<<m_name<<" : geodesic distance field for point "<<deadSize-1<<" found"<<std::endl;
	}
	
	//store result in m_points.
	m_points.reserve(deadSize);
	for(auto val: m_deads){
		m_points.push_back(tri->getVertexCoords(map[val]));
	}
	
	m_minDist = 1.E18;
	
	for(int i=0; i<deadSize; ++i){
		for(int j=i+1; j<deadSize; ++j){
			m_minDist = std::fmin(m_minDist,norm2(m_points[i] - m_points[j]));
		}
	}

	m_deads.clear();
	if(debug)	std::cout<<m_name<<" : distribution of point successfully found w/ LevelSet engine "<<std::endl;
};

/*!
 * Find your optimal distribution, projecting an initial 3D cartesian grid on the 3D object surface.
 * The final number of points is reached decimating projected points up to desired value m_points, trying
 * to maximaze euclidean distance between them
 */
void 	CreateSeedsOnSurface::solveGrid(bool debug){
	
	if(debug)	std::cout<<m_name<<" : started CartesianGrid engine"<<std::endl;
	iarray3E dim;
	
	//get the seed and project it on surface
	darray3E projSeed = bvTreeUtils::projectPoint(&m_seed, getGeometry()->getBvTree());
	if(debug)	std::cout<<m_name<<" : projected seed point"<<std::endl;
	if (m_nPoints == 1)	{
		m_points.clear();
		m_points.push_back(projSeed);
		return;
	}
	
	m_minDist = norm2(bbox->getSpan());
	double dx = m_minDist/std::pow(double(m_nPoints),0.5);
	{
		//check dimension
		std::vector<std::pair<double,int> > mapDimension(3);
		std::vector<std::pair<double,int> >::iterator itm;
		std::vector<std::pair<double,int> >::reverse_iterator ritm;
		
		mapDimension[0] = std::make_pair(bbox->getSpan()[0]/dx, 0);
		mapDimension[1] = std::make_pair(bbox->getSpan()[1]/dx, 1);
		mapDimension[2] = std::make_pair(bbox->getSpan()[2]/dx, 2);

		std::sort(mapDimension.begin(), mapDimension.end());
		
		for(itm = mapDimension.begin(); itm !=mapDimension.end(); ++itm){
			dim[itm->second] = std::max(int(itm->first + 0.5),1);
		}
		
		int cumCells = dim[0]*dim[1]*dim[2];
		ritm=mapDimension.rbegin();
		while (cumCells < m_nPoints){
		
			dim[ritm->second] += 1;
			ritm++;
			if(ritm == mapDimension.rend())	ritm = mapDimension.rbegin();
			
			cumCells = dim[0]*dim[1]*dim[2];
		}
	}
	
	//raise dimension to number of effective nodes for each direction.
	dim[0] +=1;
	dim[1] +=1;
	dim[2] +=1;
	
	//create a lattice on it
	mimmo::Lattice * grid = new Lattice();
	grid->setShape(mimmo::ShapeType::CUBE);
	grid->setOrigin(bbox->getOrigin());
	grid->setSpan(bbox->getSpan());
	grid->setRefSystem(bbox->getAxes());
	grid->setDimension(dim);
	grid->execute();
//	grid->plotGrid("./", "lattice",0,false);
	
	if(debug)	std::cout<<m_name<<" : build volume cartesian grid wrapping 3D surface"<<std::endl;
	//find narrow band cells and extracting their centroids
	dvecarr3E centroids = grid->getGlobalCellCentroids();
	
	dvecarr3E initList;
	initList.reserve(centroids.size());
	
	double distR = 0.5*norm2(grid->getSpacing());
	double dist, dummy;
	long id;
	darray3E normal;
	for(auto & p : centroids){
		dummy=distR;
		dist =  mimmo::bvTreeUtils::signedDistance(&p,getGeometry()->getBvTree(),id, normal, dummy);
		if(std::abs(dist) < distR){
			initList.push_back(p-dist*normal);
		}	
		normal.fill(0.0);	
	}
	if(debug)	std::cout<<m_name<<" : found grid cell centers in the narrow band of 3D surface and projected them on it "<<std::endl;
	
	
	//rearrange the list, putting the most nearest point to the seed on top
	// and decimate points up to desired value.
	int initListSize = initList.size();
	if( initListSize > m_nPoints){
		
		dvecarr3E secondList;
		secondList.reserve(initList.size());
		{
			int counter =0;
			dvecarr3E::iterator it, itMaxNorm = initList.begin();
			double normPoint, minValue = 1.E18;

			for(it= initList.begin(); it !=initList.end(); ++it){
				normPoint= norm2(*it - m_seed);
				if(normPoint < minValue){
					minValue = normPoint;
					itMaxNorm= it;
				}
			}
		
			darray3E temp = *itMaxNorm;
			initList.erase(itMaxNorm);
			secondList.push_back(temp);
			secondList.insert(secondList.end(),initList.begin(), initList.end());
			initList.clear();
		}
		
		initList = decimatePoints(secondList);
		if(debug)	std::cout<<m_name<<" : candidates decimated "<<std::endl;
	}
	//store result in m_points.
	m_points = initList;
	
	if(debug)	std::cout<<m_name<<" : distribution of point successfully found w/ CartesianGrid engine "<<std::endl;
	delete grid; grid = NULL;
};


/*!
 * Find distribution randomly. Regularize the distribution so that each node fullfills 
 * the maximum distance possible between points requirement.
 *\param[in]	debug	flag to activate logs of solver execution 
 */
void 	CreateSeedsOnSurface::solveRandom(bool debug){
	
	if(debug)	std::cout<<m_name<<" : started Random engine"<<std::endl;
											 dvecarr3E initList(getNPoints()); 
	
	//get the seed and project it on surface
	initList[0] = bvTreeUtils::projectPoint(&m_seed, getGeometry()->getBvTree());
	if(debug)	std::cout<<m_name<<" : projected seed point"<<std::endl;
	if (m_nPoints == 1)	{
		 m_points.clear();
		 m_points = initList;
		 return;
	}
	
	m_minDist = norm2(bbox->getSpan()) / 2.0;
	
    //fill the oriented bounding box of the figure randomly with max of 100 and 5*m_nPoints
	dvecarr3E tentative;
	{
		darray3E span = bbox->getSpan();
		dmatrix33E axes = bbox->getAxes();
		darray3E minP = bbox->getOrigin();
		for(int i=0; i<3; ++i) {minP += - 0.5*span[i]*axes[i];}

		
		if (m_randomFixed <0 ){
			m_randomFixed = (unsigned int)time(NULL);
		}
		srand(static_cast<unsigned int>(m_randomFixed));
		
		int nTent = std::max(100,5*m_nPoints);
		tentative.resize(nTent+1);
		tentative[0] = initList[0];
		
		for(int i = 0; i<nTent; ++i){
			tentative[i+1] = minP;
			for(int j=0; j<3; ++j){
				
				double valrand = (std::rand()%100)/99.0;
				tentative[i+1] +=  valrand * span[j]*axes[j]; 
			}
		}
		
		//project tentative points on surface.
		for(int i = 0; i<nTent; ++i){
			minP = bvTreeUtils::projectPoint(&tentative[i+1], getGeometry()->getBvTree()); 
			tentative[i+1] = minP;
		}
	}
	

	if(debug)	std::cout<<m_name<<" : found random points"<<std::endl;
	//decimate points up to desired value
	initList = decimatePoints(tentative);
	if(debug)	std::cout<<m_name<<" : decimated random points"<<std::endl;
	//store result in m_points.
	m_points = initList;
	if(debug)	std::cout<<m_name<<" : distribution of point successfully found w/ Random engine "<<std::endl;
};


/*!
 * Plot point distribution as *.vtu Cloud point.
 * \param[in] dir folder path
 * \param[in] file file name without tag
 * \param[in] counter identifier integer for the file
 * \param[in] binary flag to write ASCII-false, binary-true
 */
void CreateSeedsOnSurface::plotCloud(std::string dir, std::string file, int counter, bool binary){
	
	bitpit::VTKFormat codex = bitpit::VTKFormat::ASCII;
	if(binary){codex=bitpit::VTKFormat::APPENDED;}
	
	ivector1D conn(m_nPoints);
	for(int i=0; i<m_nPoints; i++){
		conn[i] = i;
	}
	
	bitpit::VTKUnstructuredGrid vtk(dir, file, bitpit::VTKElementType::VERTEX);
	vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, m_points) ;
	vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, conn) ;
	vtk.setDimensions( m_nPoints, m_nPoints);
	vtk.setCodex(codex);
	if(counter>=0){vtk.setCounter(counter);}
	
	vtk.write();
}

/*!
 * Decimate a cloud of points with number greater than m_nPoints, to desired value.
 * Regularize distribution of nodes to meet minimum distance requirements of the class.
 * First point of the list is meant as the starting seed of decimation.
 *\param[in] list of 3D points in space, with size > m_nPoints 
 */
dvecarr3E CreateSeedsOnSurface::decimatePoints(dvecarr3E & list){
	
	dvecarr3E result(m_nPoints);
	int listS = list.size();
	
	bitpit::KdTree<3,darray3E,long>	kdT;
	//build a kdtree of points
	kdT.nodes.resize(list.size() + kdT.MAXSTK);
	long label=0;
	for(auto & val : list ){
		kdT.insert(&val, label);
		++label;
	}
		
	long candidate = 0; //starting seed;
	
	//use kdTree to search points within a radius of m_minDist from the seed.
	// if no matches are found, choose randomly the next candidate
	livector1D neighs, excl,effective;
	std::set<long>	visited;
	darray3E temp;
	double valdum;
	std::map<double, long> chooseCand; 
	livector1D finalCandidates;
	
	finalCandidates.reserve(m_nPoints);
	int candSize = finalCandidates.size();
	while( candSize < m_nPoints){
		
		finalCandidates.push_back(candidate);
		candSize++;
		visited.insert(candidate);
		excl.insert(excl.end(), visited.begin(), visited.end());
		kdT.hNeighbors(&list[candidate], m_minDist, &neighs, & excl );

		visited.insert(neighs.begin(), neighs.end());
		
		int sizeN = visited.size();
		neighs.clear();
		excl.clear();
		
		
		while(sizeN == listS){
			m_minDist *= 0.9;
			visited.clear();
			visited.insert(finalCandidates.begin(), finalCandidates.end());
			for(auto ind : finalCandidates){
				excl.insert(excl.end(), visited.begin(), visited.end());
				kdT.hNeighbors(&list[ind], m_minDist, &neighs, & excl );
				visited.insert(neighs.begin(), neighs.end());
			}
			sizeN = visited.size();
			excl.clear();
			neighs.clear();
		}
			
		effective.reserve(list.size() - visited.size());
		std::set<long>::iterator it1 = visited.begin();
		for(int i=0; i< listS; ++i){
			if(it1 != visited.end() || i != *it1 ){
				effective.push_back(i);
				
			}else{ ++it1;}
		}
			
		for(auto & ind : effective){
			dvector1D norms(finalCandidates.size());
			int countNorms=0;
			for(auto indEx: finalCandidates){
				norms[countNorms] = norm2(list[ind] - list[indEx]);
				++countNorms;
			}
				
			double norm = 0.0, variance=0.0;
			for(auto &val : norms)	norm += val;
			norm /= double(norms.size());
			for(auto &val : norms)	val  = std::abs(val/norm - 1.0);
			maxval(norms, variance);
			norm /= std::pow((variance+1.0),2);
				
			chooseCand[norm] = ind;
		}
			
		candidate = (chooseCand.rbegin())->second;
		
		effective.clear();
		chooseCand.clear();
    }
	
	int counter = 0;
	for(auto index : finalCandidates){
		result[counter] = list[index];
		++counter;
	}
	
	return result;
};


//PRIVATE METHODS
/*!
 * Update m_sdf distance field on a target node of a superficial tessellation solving 
 * the Eikonal equation |grad(u)| = g, using  a fast marching method. Tessellation must be
 * mandatorily a triangular one. 
 * \param[in] g			propagation speed of the Eikonal equation      
 * \param[in] s			flag for inwards/outwards propagation (s = -+1)
 * \param[in] Itarget   id of the target node in bitpit::PatchKernel indexing
 * \param[in] Ttarget   id of the cell which the Itarget belongs to in bitpit::PatchKernel indexing 
 * \param[in] flag	 	flag vector reporting eikonal front advancing status on nodes. Using dead(= 0), alive(= 1),and far away(= 2) identifiers.    
 * @param[in] field     reference distance field 
 * \return	result, updated value of the m_sdf distance field on the target node.
 */
double CreateSeedsOnSurface::updateEikonal(double g, double s, long tVert,long tCell, std::unordered_map<long int, short int> &flag, dvector1D & field){
	//get the pointer to reference geometry
	liimap & vmap = getGeometry()->getMapDataInv();
	bitpit::PatchKernel * tri = getGeometry()->getPatch();
	
	//todo get dimension? of Simplex. Get a check for non triangular cell or generalize it
	
	livector1D				oneRing;
	long					I, U, V, W;
	
	int						k, m;
	
	int						select = -1;
	double					value(1.0e+18);
	double					dVU, dVW, dWU, dVP;
	std::array<double,3>	eVU, eVW, eWU, eVP, P;
	
	
	double					a, b, c, A, B, C, K, discr ;
	double					phi_U, phi_W, phi_P;
	int 					discrType;
	
	double					xi1, xi2;
	double					tempVal1, tempVal2;
	
	//get current field value of the node;
	value = std::abs(field[vmap[tVert]]);
	
	V = tVert;
	
	{
		// find 1-Ring cells around target node
		bitpit::Cell & targetCell = tri->getCell(tCell);
		int locVert = targetCell.findVertex(tVert);
		if(locVert == bitpit::Vertex::NULL_ID) return value;
		oneRing = tri->findCellVertexOneRing(tCell, locVert);
	}
	
	// Loop over cells in the 1-Ring --------------------------------------------------- //
	for (auto && oneIndex : oneRing) {
		
		// Cell data get id of vertex composing triangular cell
		I = oneIndex;
		bitpit::Cell & cellI = tri->getCell(I);
		
		k = cellI.findVertex(V);
		k = (k + 1) % cellI.getVertexCount();
		U = cellI.getVertex(k);
		m = (k + 1) % cellI.getVertexCount();
		W = cellI.getVertex(m);
		
		// discriminate case, according to flag vector of deads, alives and far-aways
		if ((flag[U] == 0) && (flag[W] == 0)) {
			select = 2;
		}
		else {
			if ((flag[U] == 0) || (flag[W] == 0)) {
				select = 1;
				if (flag[W] == 0) {
					U = W;
				}
			}
			else {
				select = 0;
			}
		}
		
		//Compute solution to the 2D Eikonal equation
		switch (select){
			
			case 1 : //  with 1 dead node
				eVU = tri->getVertexCoords(V) - tri->getVertexCoords(U);
				dVU = norm2(eVU);
				value = std::min(value, std::abs(field[vmap[U]]) + g*dVU); ///??????????????????????????????????????????
				break;
				
			case 2 : // with 2 dead nodes
				
				eVW = tri->getVertexCoords(V) - tri->getVertexCoords(W);
				dVW = norm2(eVW);
				eVW = eVW/dVW;
				eVU = tri->getVertexCoords(V) - tri->getVertexCoords(U);
				dVU = norm2(eVU);
				eVU = eVU/dVU;
				eWU = tri->getVertexCoords(W) - tri->getVertexCoords(U);
				dWU = norm2(eWU);
				eWU = eWU/dWU;
				
				// Coeffs -------------------------------------------------------------------- //
				phi_U = std::abs(field[vmap[U]]);
				phi_W = std::abs(field[vmap[W]]);
				K = phi_W - phi_U;
				a = pow(dWU, 2);
				b = -dWU*dVU*dotProduct(eWU, eVU);
				c = pow(dVU, 2);
				A = a*(pow(K, 2) - a);
				B = b*(pow(K, 2) - a);
				C = (pow(K, 2)*c - pow(b, 2));
				discr = pow(B, 2) - A*C;
				
				// Find optimal solution ----------------------------------------------------- //
				discrType = (discr < -1.0e-12) + 2*(std::abs(A) > 1.0e-12) +3*((std::abs(A) < 1.0e-12) && (std::abs(A) >= 0.0)) -1;
				
				switch(discrType){
					case 1: //2 distinct solutions
						discr = std::abs(discr);
						
						xi1 = (-B - sqrt(discr))/A;
						xi2 = (-B + sqrt(discr))/A;
						
						// Restriction of solutions onto [0, 1]
						xi1 = std::min(1.0, std::max(0.0, xi1));
						xi2 = std::min(1.0, std::max(0.0, xi2));
						
						// Solution #1
						P = (1.0 - xi1) * tri->getVertexCoords(U)  +  xi1 * tri->getVertexCoords(W);
						eVP = tri->getVertexCoords(V) - P;
						dVP = norm2(eVP);
						eVP = eVP/dVP;
						phi_P = (1.0 - xi1) * phi_U + xi1 * phi_W;
						tempVal1 = phi_P + g * dVP;
						
						// Solution #2
						P = (1.0 - xi2) * tri->getVertexCoords(U)  +  xi2 * tri->getVertexCoords(W);
						eVP = tri->getVertexCoords(V) - P;
						dVP = norm2(eVP);
						eVP = eVP/dVP;
						phi_P = (1.0 - xi2) * phi_U + xi2 * phi_W;
						tempVal2 = phi_P + g * dVP;
						
						break;
						
					case 2: // coincident solutions	
						discr = std::abs(discr);
						
						// Solution #1
						P = tri->getVertexCoords(U);
						eVP = tri->getVertexCoords(V) - P;
						dVP = norm2(eVP);
						eVP = eVP/dVP;
						phi_P = phi_U;
						tempVal1 = phi_P + g*dVP;
						
						// Solution #2
						P = tri->getVertexCoords(W);
						eVP = tri->getVertexCoords(V) - P;
						dVP = norm2(eVP);
						eVP = eVP/dVP;
						phi_P = phi_W;
						tempVal2 = phi_P + g*dVP;
						
						break;
						
					default: //no real solutions indeed
						tempVal1 = value;
						tempVal2 = value;
						break;
				}//end on switch discrType	
				
				// Update solution for case 2:
				value = std::min(value, std::min(tempVal1, tempVal2));
				break;
				
					default: //doing really nothing. No dead nodes to hang out.
						break;
		}//end switch select
		
	} //loop on oneRing
	
	tri = NULL;
	
	return(value);
}; 

/*!
 * Solve the 3D Eikonal equation |grad(u)| = g, using  a fast marching method, on a target triangulation 
 * associated unstructured superficial grid. Tessellation must be mandatorily a triangular one.
 * (SurfaceConstraints::m_sdf values in the unknown region must be set to the value 1.0e+18.
 *  Dead vertices of front to be propagated must be set to zero, initially)
 * @param[in] g Propagation speed.
 * @param[in] s Velocity sign (+1 --> propagate outwards, -1 --> propagate inwards).
 * @param[in] invConn inverse connectivity of yout current geometry
 * @param[out] field field to be computed, already allocated.
 */
void CreateSeedsOnSurface::solveEikonal(double g, double s, std::unordered_map<long,long> & invConn, dvector1D & field ){	
	//recover bitpit::PatchKernel
	liimap & vmap = getGeometry()->getMapDataInv();
	bitpit::PatchKernel * tri = getGeometry()->getPatch();
	
		
	// declare total size and support structure
	long 	N(tri->getVertexCount());
	std::unordered_map<long int, short int> active;
	
	{ //FLAG DEAD/ALIVE/FAR-AWAY VERTICES                                                   
		
		long myId;
		bool check;
		std::set<long>				neighs;
		std::set<long>::iterator	it, itend;
		
		//set active vector size and mark its position  with original geometry ids
		active.reserve(N);
		for ( const auto &vertex : tri->getVertices() ){ 
			myId           = vertex.getId() ;
			active[myId] = 2 ;
		}
		
		//fill active vector
		for ( const auto &vertex : tri->getVertices() ){ 
			myId 	=	vertex.getId();
			
			// Dead vertices
			if( isDeadFront(vmap[myId]) ){ 
				active[myId] = 0;
				
			}else{
				
				//loop over neighbors
				check = false;
				neighs = findVertexVertexOneRing(invConn[myId], myId); 
				it = neighs.begin();
				itend = neighs.end();
				while(!check && it !=itend){
					
					check = s*field[vmap[*it]] >= 0.0 && field[vmap[*it]] < 1.0E+18;
					++it;
				};
				
				active[myId] = 2 - (int) check;
			}
		}
	}	
	
	
	{ // Construct min heap data structure 
		long                            m(0), I(0), myId, J ;
		double                          value ;
		
		std::set<long>				neighs;
		std::set<long>::iterator	it,itbeg, itend;
		
		std::vector<std::array<int,2>>  map(N), *mapPtr = &map;
		
		bitpit::MinPQueue<double, long> heap(N, true, mapPtr);
		
		// Inserting alive vertices in  the heap
		
		for(const auto & vertex : tri->getVertices()){
			
			myId = vertex.getId();
			if(active[myId] == 1) {
				//assign a value to your actual vertex
				value = updateEikonal(s, g, myId, invConn[myId], active, field);
				
				//store it into heap
				map[m][0] = vmap[myId];
				map[vmap[myId]][1] = m;
				
				heap.keys[m] = value;
				heap.labels[m] = myId;
				
				//update counter
				++m;
			}
			++I;
		}//next vertex	
		
		// Build min-heap
		heap.heap_size = m;		
		heap.buildHeap();
		
		
		// FAST MARCHING                                                                       //
		while (heap.heap_size > 0) {
			
			// Extract root from min-heap
			heap.extract(value, myId);
			
			// Update level set value
			//value =  s*updateEikonal(s, g, myId, invConn[myId], active);
			field[vmap[myId]] = value; 
			
			// Update flag to dead;
			active[myId] = 0;
			
			//update neighbours
			neighs = findVertexVertexOneRing(invConn[myId], myId); 
			itbeg = neighs.begin();
			itend = neighs.end();
			
			for(it=itbeg; it != itend; ++it) {
				J = *it;
				
				if(active[J] == 1){
					
					//update local value;
					value = updateEikonal(s,g,J,invConn[J],active, field);
					
					//update its value in the min-heap
					I = vmap[J];
					heap.modify( map[I][1],value,J );
					
				}else if( active[J] == 2){
					
					//update local value;
					value = updateEikonal(s,g,J,invConn[J],active, field);
					
					//reflag vertex as alive vertex
					active[J] = 1;
					I = vmap[J];
					
					// Insert neighbor into the min heap
					map[heap.heap_size][0] = I ;
					map[I][1] = heap.heap_size;
					
					heap.insert(value, J);
				}
			}
		}//end while	
	};
	
	tri = NULL;
	return;
};	

/*!
 * Get a minimal inverse connectivity of the target geometry mesh associated to the class.
 * Each vertex (passed by Id) is associated at list to one of the possible simplex 
 * (passed by Id) which it belongs. This is returned in an unordered_map having as key the 
 * vertex Id and as value the Cell id. Id is meant as the unique label identifier associated
 * to bitpit::PatchKernel original geometry
 *\return	unordered_map of vertex ids (key) vs cell-belonging-to ids(value) 
 */
std::unordered_map<long,long> CreateSeedsOnSurface::getInverseConn(){
	
	bitpit::PatchKernel * tri = getGeometry()->getPatch();
	std::unordered_map<long,long> invConn ;
	
	long cellId, vertID;
	const long * locConn;
	int sizeConn, counter;
	
	for(const auto &cell : tri->getCells()){
		cellId = cell.getId();
		sizeConn = cell.getVertexCount();
		locConn = cell.getConnect();
		for(counter=0; counter<sizeConn; ++counter) invConn[locConn[counter]] = cellId;
	}
	
	locConn  = NULL;
	tri = NULL;
	
	return(invConn);
};
/*!
 * Return true if a given vertex belongs to the current constrained boundary front of your patch
 * \param[in]	label index of vertex, in sequential mimmo::MimmoObject notation
 * \return boolean, true if vertex belongs to constrained set, false if not 
 */
bool CreateSeedsOnSurface::isDeadFront(const int label){
	
	ivector1D::iterator got = std::find(m_deads.begin(), m_deads.end(), label);
	if(got == m_deads.end()) return false;
	return true;
}
/*!
 * Return VertexVertex One Ring of a specified target vertex
 * \param[in]	cellId 	bitpit::PatchKernel Id of a cell which target belongs to  
 * \param[in]	vertId	bitpit::PatchKernel Id of the target vertex
 * \return		list of all vertex in the One Ring of the target, by their bitpit::PatchKernel Ids 
 */
std::set<long> CreateSeedsOnSurface::findVertexVertexOneRing(const long & cellId, const long & vertexId){
	
	std::set<long> result;
	bitpit::PatchKernel * tri = getGeometry()->getPatch();
	bitpit::Cell &cell =  tri->getCell(cellId);
	
	int loc_target = cell.findVertex(vertexId);
	if(loc_target ==bitpit::Vertex::NULL_ID) return result;
	
	livector1D list = tri->findCellVertexOneRing(cellId, loc_target);
	
	long connSize;
	for(auto && index : list){
		bitpit::Cell & cell = tri->getCell(index);
		connSize = cell.getVertexCount();
		
		for(int i=0; i<connSize; ++i){
			result.insert(cell.getVertex(i));
		}
	}
	
	result.erase(vertexId);
	tri = NULL;
	return result;
}

/*!
 * Get infos from a XML bitpit::Config::section. The parameters that can absorb are
 * 
 * 1) NPoints - total points to distribute 
 * 2) Engine  - type of distribution engine 0-Random,2-CartesianGrid,1-Levelset;
 * 3) Seed    - initial seed point;
 * 4) MassCenterAsSeed - boolean, if true use geometry mass center sa seed
 * 5) RandomFixed- get signature to fix distribution pattern when 0-RANDOM engine is selected
 * 6) PlotInExecution - plot optional result during object execution
 * 7) OutputPlot  - set path to store optional result at 6) 
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void CreateSeedsOnSurface::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name ){
	
	if(slotXML.hasOption("NPoints")){
		std::string input = slotXML.get("NPoints");
		input = bitpit::utils::trim(input);
		int value = 0;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
			value = std::max(value,0);
		}
		setNPoints(value);
	}
	
	if(slotXML.hasOption("Engine")){
		std::string input = slotXML.get("Engine");
		input = bitpit::utils::trim(input);
		int value = 2;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
			value = std::min(std::max(value,0),2);
		}
		setEngine(value);
	}
	
	
	if(slotXML.hasOption("Seed")){
		std::string input = slotXML.get("Seed");
		input = bitpit::utils::trim(input);
		darray3E temp;
		temp.fill(0.0);
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> temp[0]>>temp[1]>>temp[2];
		}
		setSeed(temp);
	}

	if(slotXML.hasOption("MassCenterAsSeed")){
		std::string input = slotXML.get("MassCenterAsSeed");
		input = bitpit::utils::trim(input);
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setMassCenterAsSeed(value);
	}
	
	if(slotXML.hasOption("RandomFixed")){
		std::string input = slotXML.get("RandomFixed");
		input = bitpit::utils::trim(input);
		int value = -1;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setRandomFixed(value);
	}
	
	if(slotXML.hasOption("PlotInExecution")){
		std::string input = slotXML.get("PlotInExecution");
		input = bitpit::utils::trim(input);
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setPlotInExecution(value);
	}
	
	if(slotXML.hasOption("OutputPlot")){
		std::string input = slotXML.get("OutputPlot");
		input = bitpit::utils::trim(input);
		std::string temp = ".";
		if(!input.empty())	setOutputPlot(input);
		else			  	setOutputPlot(temp);
	}
	
return;	
};

/*!
 * Plot class infos to a XML bitpit::Config::section. The parameters that can be flushed are
 * 
 * 1) NPoints - total points to distribute 
 * 2) Engine  - type of distribution engine 0-Random,1-CartesianGrid,2-Levelset;
 * 3) Seed    - initial seed point;
 * 4) MassCenterAsSeed - boolean, if true use geometry mass center sa seed
 * 5) RandomFixedSeed - get signature to fix distribution pattern when 0-RANDOM engine is selected
 * 6) PlotInExecution - plot optional result during object execution
 * 7) OutputPlot  - set path to store optional result at 6) 
 *
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void CreateSeedsOnSurface::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
	if(m_nPoints != 0 ){
		slotXML.set("NPoints", std::to_string(m_nPoints));
	}
	
	int engint = getEngine();
	if(engint != 2 ){
		slotXML.set("Engine", std::to_string(engint));
	}
	
	darray3E seed = getSeed();
	if(norm2(seed) != 0.0 ){
		
		std::stringstream ss;
		ss<<std::scientific<<seed[0]<<'\t'<<seed[1]<<'\t'<<seed[2];
		slotXML.set("Seed", ss.str());
	}
	
    if(isSeedMassCenter() ){
		slotXML.set("MassCenterAsSeed", std::to_string(1));
	}	
	
	int signat = getRandomSignature();
	if( signat != -1 && m_engine == mimmino::CSeedSurf::RANDOM){
		slotXML.set("RandomFixed", std::to_string(signat));
	}
	
	if(isPlotInExecution()){
		slotXML.set("PlotInExecution", std::to_string(1));
	}
	
	if(m_outputPlot != "."){
		slotXML.set("OutputPlot", m_outputPlot);
	}
	
};


