/*---------------------------------------------------------------------------*\
 * 
 *  mimmo
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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
 \ *---------------------------------------------------------------------------*/

#include "CreateSeedsOnSurface.hpp"
#include "Lattice.hpp"
#include "SkdTreeUtils.hpp"

#include <stdlib.h>
#include <time.h>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <iostream>
#include <fstream>
#include <cmath>

namespace mimmo{

/*!
 * Constructor
 */
CreateSeedsOnSurface::CreateSeedsOnSurface(){
    m_name = "mimmo.CreateSeedsOnSurface";
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
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
CreateSeedsOnSurface::CreateSeedsOnSurface(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.CreateSeedsOnSurface";
    m_nPoints = 0;
    m_minDist = 0.0;
    m_seed = {{0.0,0.0,0.0}};
    m_engine = CSeedSurf::CARTESIANGRID;
    m_seedbaricenter = false;
    m_randomFixed = -1;
    std::unique_ptr<mimmo::OBBox> box(new mimmo::OBBox());
    bbox = std::move(box);

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.CreateSeedsOnSurface"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Destructor;
 */
CreateSeedsOnSurface::~CreateSeedsOnSurface(){
    clear();
};

/*!
 * Copy constructor
 */
CreateSeedsOnSurface::CreateSeedsOnSurface(const CreateSeedsOnSurface & other):BaseManipulation(other){
    m_points = other.m_points;
    m_nPoints = other.m_nPoints;
    m_minDist = other.m_minDist;
    m_seed = other.m_seed;
    m_engine = other.m_engine;
    m_seedbaricenter = other.m_seedbaricenter;
    m_randomFixed = other.m_randomFixed;
    m_deads = other.m_deads;
    m_sensitivity = other.m_sensitivity;
    bbox = std::move(std::unique_ptr<mimmo::OBBox>(new mimmo::OBBox(*(other.bbox.get()))));
};

/*!
 * Assignment operator
 */
CreateSeedsOnSurface & CreateSeedsOnSurface::operator=(CreateSeedsOnSurface other){
    swap(other);
    return *this;
}

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void CreateSeedsOnSurface::swap(CreateSeedsOnSurface & x) noexcept
{
    std::swap(m_points, x.m_points);
    std::swap(m_nPoints, x.m_nPoints);
    std::swap(m_minDist, x.m_minDist);
    std::swap(m_seed, x.m_seed);
    std::swap(m_engine, x.m_engine);
    std::swap(m_seedbaricenter, x.m_seedbaricenter);
    std::swap(m_randomFixed, x.m_randomFixed);
    std::swap(m_deads, x.m_deads);
//     std::swap(m_sensitivity, x.m_sensitivity);
    m_sensitivity.swap(x.m_sensitivity);
    std::swap(bbox, x.bbox);
    BaseManipulation::swap(x);
}
/*!
 * It builds the input/output ports of the object
 */
void
CreateSeedsOnSurface::buildPorts(){

    bool built = true;

    //input
    built = (built && createPortIn<MimmoObject *, CreateSeedsOnSurface>(this, &mimmo::CreateSeedsOnSurface::setGeometry,M_GEOM, true));
    built = (built && createPortIn<darray3E, CreateSeedsOnSurface>(this, &mimmo::CreateSeedsOnSurface::setSeed, M_POINT));
    built = (built && createPortIn<int, CreateSeedsOnSurface>(this, &mimmo::CreateSeedsOnSurface::setNPoints, M_VALUEI));
    built = (built && createPortIn<int, CreateSeedsOnSurface>(this, &mimmo::CreateSeedsOnSurface::setRandomFixed, M_VALUEI2 ));
    built = (built && createPortIn<bool, CreateSeedsOnSurface>(this, &mimmo::CreateSeedsOnSurface::setMassCenterAsSeed, M_VALUEB));
    built = (built && createPortIn<dmpvector1D, CreateSeedsOnSurface>(this, &mimmo::CreateSeedsOnSurface::setSensitivityMap, M_FILTER));
    
    //output
    built = (built && createPortOut<dvecarr3E, CreateSeedsOnSurface>(this, &mimmo::CreateSeedsOnSurface::getPoints, M_COORDS));

    m_arePortsBuilt = built;
};

/*!
 * Return number of points to be distributed over 3D surface.
 * \return number of points
 */
int
CreateSeedsOnSurface::getNPoints(){
    return m_nPoints;
};

/*!
 * Return list of points distributed over 3D surface, after class execution.
 * \return list of points
 */
dvecarr3E
CreateSeedsOnSurface::getPoints(){
    return m_points;
};

/*!
 * Return CSeedSurf engine type for point distribution calculation.
 * \return engine type
 */
CSeedSurf
CreateSeedsOnSurface::getEngineENUM(){
    return m_engine;
};


/*!
 * Return engine type for point distribution calculation.
 * \return engine type
 */
int
CreateSeedsOnSurface::getEngine(){
    return static_cast<int>(m_engine);
};

/*!
 * Return seed point used for calculation. If geometry baricenter is
 * used as initial seed (see setMassCenterAsSeed), it can be 
 * correctly visualized after the execution of the object.
 * \return seed point
 */
darray3E
CreateSeedsOnSurface::getSeed(){
    return m_seed;
}

/*!
 * Return true, if the option to use geometry mass center as initial seed is
 * activated. (see setMassCenterAsSeed method documentation).
 * \return mass center flag
 */
bool
CreateSeedsOnSurface::isSeedMassCenter(){
    return m_seedbaricenter;
};


/*!
 * Return minimum absolute distance during the point distribution calculation.
 * \return minimum distance between points
 */
double
CreateSeedsOnSurface::getMinDistance(){
    return m_minDist;
}

/*!
 * Return the signature of the current random distribution of points on target surface,
 * whenever is fixed or not, for result replication. See setRandomFixed method.  This option will
 * only make sense if a CSeedSurf::RANDOM engine is employed. Otherwise it is ignored. 
 * \return signature. 
 */
int
CreateSeedsOnSurface::getRandomSignature(){
    return m_randomFixed;
}

/*!
 * Set the number of points to be distributed.
 * \param[in]    val    number of total points
 */
void
CreateSeedsOnSurface::setNPoints(int val){
    m_nPoints = std::max(val,0);
}

/*!
 * Set the engine type for point distribution evaluation
 * \param[in]    eng engine type passed as enum CSeedSurf
 */
void
CreateSeedsOnSurface::setEngineENUM(CSeedSurf eng){
    m_engine = eng;
}

/*!
 * Set the engine type for point distribution evaluation
 * \param[in]    eng engine type passed as integer. See CSeedSurf enum.
 */
void
CreateSeedsOnSurface::setEngine(int eng){
    if(eng <0 || eng >2)    eng = 2;
    setEngineENUM(static_cast<CSeedSurf>(eng));
}

/*!
 * Set the starting seed point for points distribution on surface evaluation.
 * The assigned seed is ignored if the geometry baricenter is set as initial seed
 * (see setMassCenterAsSeed).
 * \param[in]    seed starting 3D point in the neighborhood of the surface
 */
void
CreateSeedsOnSurface::setSeed(darray3E seed){
    m_seed = seed;
}

/*!
 * Set the starting seed point for points distribution on surface evaluation as the mass center
 * of the surface itself. Any seed passed with setSeed method will be ignored if this option is active.
 * \param[in]    flag true activate the option, false deactivate it.
 */
void
CreateSeedsOnSurface::setMassCenterAsSeed(bool flag){
    m_seedbaricenter = flag;
}

/*!
 * Set Geometry, check if its search bvtree is built, if not, build it. 
 * Point Cloud geometries, or pure Volume meshes are currently not supported by this block.
 * Reimplemented from BaseManipulation::setGeometry()
 * \param[in] geo pointer to target geometry
 */
void
CreateSeedsOnSurface::setGeometry(MimmoObject * geo){
    if(geo == NULL)    return;
    if(geo->isEmpty())  return;
    if(geo->getType() != 1)    return;

    BaseManipulation::setGeometry(geo);
    if(!geo->isSkdTreeSync())            getGeometry()->buildSkdTree();
    if(!geo->areAdjacenciesBuilt() )    getGeometry()->buildAdjacencies();
    bbox->setGeometry(geo);
}

/*!
 * Set the signature (each integer >= 0) of your random distribution. Same signature will be able to reproduce 
 * the exact random distribution in multiple runs. If signature is < 0 (default), point distribution will randomly 
 * vary run by run.It is possible to get the current signature after each random execution using the getRandomSignature method.
 * This option will only make sense if a CSeedSurf::RANDOM engine is employed. Otherwise it is ignored.
 *\param[in] signature integer number   
 */
void
CreateSeedsOnSurface::setRandomFixed( int signature){
    m_randomFixed = signature;
}

/*!
 * Set a sensitivity scalar field, referred to the target geometry linked, to drive placement of the seeds points
 * on the most sensitive part of the geometry. The sensitivity field MUST be defined on geometry vertices.  
 *\param[in] field sensitivity  
 */
void
CreateSeedsOnSurface::setSensitivityMap( dmpvector1D field){
    m_sensitivity = field;
}

/*!
 * Clear contents of the class
 */
void
CreateSeedsOnSurface::clear(){
    m_points.clear();
    m_nPoints = 0;
    m_minDist = 0.0;
    m_seed = {{0.0,0.0,0.0}};
    m_engine = CSeedSurf::CARTESIANGRID;
    m_seedbaricenter = false;
    m_randomFixed = -1;
    m_deads.clear();
    m_sensitivity.clear();

}

/*!
 * Execution command of the class. Wrapper to solve( bool ) method
 */
void
CreateSeedsOnSurface::execute(){
    solve(false);
}

/*!
 * Plot Optional results of the class, that is the point cloud distribution on surface
 * as a list of xyz data points w/ identifier mark and as a point cloud in *.vtu format
 */
void
CreateSeedsOnSurface::plotOptionalResults(){

    std::string dir = m_outputPlot;
    std::string nameGrid  = m_name + "CLOUD";
    plotCloud(dir, nameGrid, getId(), true );
}

/*!
 * Apply command of the class;
 * \param[in] debug flag to activate logs of solver execution
 */
void
CreateSeedsOnSurface::solve(bool debug){

    if(getGeometry() == NULL){
        throw std::runtime_error(m_name + " : NULL pointer to linked geometry.");
    }
    
    if(getGeometry()->isEmpty()|| m_nPoints< 1){
        throw std::runtime_error(m_name + " : empty geometry or not seeding intial point defined.");
    }
    m_points.clear();
    bbox->execute();
    if(m_seedbaricenter)    m_seed = bbox->getOrigin();
    
    checkField();
    normalizeField();
    
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
 *\param[in]    debug    flag to activate logs of solver execution
 */
void
CreateSeedsOnSurface::solveLSet(bool debug){

    if(debug)    (*m_log)<<m_name<<" : started LevelSet engine"<<std::endl;
    dvecarr3E initList;
    m_deads.reserve(m_nPoints);

    //understand if the class is a pure triangulation or not
    std::unique_ptr<MimmoObject> objTriangulated;
    MimmoObject * workgeo;
    dmpvector1D worksensitivity;
    if(!checkTriangulation()){
        objTriangulated = std::move(triangulate());
        workgeo = objTriangulated.get();
        worksensitivity = m_sensitivity;
    }else{
        workgeo   = getGeometry();
        worksensitivity = m_sensitivity_triangulated;
    }

    if(!(workgeo->areAdjacenciesBuilt()) ) workgeo->buildAdjacencies();
    if(!(workgeo->isSkdTreeSync())) workgeo->buildSkdTree();
    if(!(workgeo->isKdTreeSync())) workgeo->buildKdTree();
    bitpit::SurfUnstructured * tri = static_cast<bitpit::SurfUnstructured * >(workgeo->getPatch());

    double distance = 0.0;
    for(const auto &cell : tri->getCells()){
        distance += tri->evalCellArea(cell.getId());
    }
    distance /= double(tri->getCellCount());
    distance = std::pow(distance, 0.5);

    livector1D neighs, excl;
    darray3E projSeed = skdTreeUtils::projectPoint(&m_seed, workgeo->getSkdTree());
    bitpit::Vertex vertSeed(0, projSeed);
    //find the vertex of the mesh nearest to projSeed.
    int nSize = 0;
    while( nSize < 1){
        workgeo->getKdTree()->hNeighbors(&vertSeed, distance, &neighs, & excl );
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

    m_deads.push_back(candidate);
    int deadSize = m_deads.size();
    if(debug)    (*m_log)<<m_name<<" : projected seed point"<<std::endl;

    std::unordered_map<long,long> invConn = getInverseConn(*(workgeo->getPatch()));
    if(debug)    (*m_log)<<m_name<<" : created geometry inverse connectivity"<<std::endl;

    while(deadSize < m_nPoints){

        dmpvector1D field;
        for (const auto & v : tri->getVertices()){
            field.insert(v.getId(), 1.0E+18);
        }
        for(const auto & dd : m_deads)    field[dd] = 0.0;

        solveEikonal(1.0,1.0, *(workgeo->getPatch()), invConn, field);

        //modulate field with current working sensitivity field
        auto itSE=worksensitivity.end();
        for(auto itSX =worksensitivity.begin(); itSX != itSE; ++itSX){
            field[itSX.getId()] *= *itSX;
        }
        
        double maxField= 0.0;
        long candMax =0;
        auto itE=field.end();
        for(auto itX =field.begin(); itX != itE; ++itX){
            if(*itX > maxField){
              maxField = *itX;
              candMax = itX.getId();
            }
        }
        
        m_deads.push_back(candMax);

        deadSize = m_deads.size();
        if(debug)    (*m_log)<<m_name<<" : geodesic distance field for point "<<deadSize-1<<" found"<<std::endl;
    }

    //store result in m_points.
    m_points.reserve(deadSize);
    for(const auto & val: m_deads){
        m_points.push_back(tri->getVertexCoords(val));
    }

    m_minDist = 1.E18;

    for(int i=0; i<deadSize; ++i){
        for(int j=i+1; j<deadSize; ++j){
            m_minDist = std::fmin(m_minDist,norm2(m_points[i] - m_points[j]));
        }
    }

    m_deads.clear();
    if(debug)    (*m_log)<<m_name<<" : distribution of point successfully found w/ LevelSet engine "<<std::endl;
};

/*!
 * Find your optimal distribution, projecting an initial 3D cartesian grid on the 3D object surface.
 * The final number of points is reached decimating projected points up to desired value m_points, trying
 * to maximaze euclidean distance between them
 * \param[in] debug flag to activate logs of solver execution
 */
void
CreateSeedsOnSurface::solveGrid(bool debug){

    if(debug)    (*m_log)<<m_name<<" : started CartesianGrid engine"<<std::endl;
    iarray3E dim;
    
    if(!(getGeometry()->isSkdTreeSync())) getGeometry()->buildSkdTree();
    
    //get the seed and project it on surface
    darray3E projSeed = skdTreeUtils::projectPoint(&m_seed, getGeometry()->getSkdTree());
    if(debug)    (*m_log)<<m_name<<" : projected seed point"<<std::endl;
    if (m_nPoints == 1)    {
        m_points.clear();
        m_points.push_back(projSeed);
        return;
    }

    m_minDist = norm2(bbox->getSpan());
    double dx = m_minDist/int(std::pow(double(m_nPoints),0.5)+ 0.5);
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
            if(ritm == mapDimension.rend())    ritm = mapDimension.rbegin();

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
    //grid->plotGrid("./", "lattice",0,false);

    if(debug)    (*m_log)<<m_name<<" : build volume cartesian grid wrapping 3D surface"<<std::endl;
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
        dist =  mimmo::skdTreeUtils::signedDistance(&p,getGeometry()->getSkdTree(),id, normal, dummy);
        if(std::abs(dist) < distR){
            initList.push_back(p-dist*normal);
        }
        normal.fill(0.0);
    }
    if(debug)    (*m_log)<<m_name<<" : found grid cell centers in the narrow band of 3D surface and projected them on it "<<std::endl;


    //rearrange the list, putting the most nearest point to the seed on top
    // and decimate points up to desired value.
    int initListSize = initList.size();
    if( initListSize > m_nPoints){

        dvecarr3E secondList;
        secondList.reserve(initList.size());
        {
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
        if(debug)    (*m_log)<<m_name<<" : candidates decimated "<<std::endl;
    }
    //store result in m_points.
    m_points = initList;
    m_nPoints = (int)m_points.size();
    if(debug)    (*m_log)<<m_name<<" : distribution of point successfully found w/ CartesianGrid engine "<<std::endl;
    delete grid; grid = NULL;
};


/*!
 * Find distribution randomly. Regularize the distribution so that each node fullfills 
 * the maximum distance possible between points requirement.
 * \param[in]    debug    flag to activate logs of solver execution
 */
void
CreateSeedsOnSurface::solveRandom(bool debug){

    if(debug)    (*m_log)<<m_name<<" : started Random engine"<<std::endl;
    dvecarr3E initList(getNPoints());

    if(!(getGeometry()->isSkdTreeSync())) getGeometry()->buildSkdTree();
    
    //get the seed and project it on surface
    initList[0] = skdTreeUtils::projectPoint(&m_seed, getGeometry()->getSkdTree());
    if(debug)    (*m_log)<<m_name<<" : projected seed point"<<std::endl;
    if (m_nPoints == 1)    {
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
            minP = skdTreeUtils::projectPoint(&tentative[i+1], getGeometry()->getSkdTree());
            tentative[i+1] = minP;
        }
    }


    if(debug)    (*m_log)<<m_name<<" : found random points"<<std::endl;
    //decimate points up to desired value
    initList = decimatePoints(tentative);
    if(debug)    (*m_log)<<m_name<<" : decimated random points"<<std::endl;
    //store result in m_points.
    m_points = initList;
    m_nPoints = (int)m_points.size();
    if(debug)    (*m_log)<<m_name<<" : distribution of point successfully found w/ Random engine "<<std::endl;
};


/*!
 * Plot point distribution as *.vtu Cloud point.
 * \param[in] dir folder path
 * \param[in] file file name without tag
 * \param[in] counter identifier integer for the file
 * \param[in] binary flag to write ASCII-false, binary-true
 */
void
CreateSeedsOnSurface::plotCloud(std::string dir, std::string file, int counter, bool binary){

    if(m_points.empty())    return;
    
    bitpit::VTKFormat codex = bitpit::VTKFormat::ASCII;
    if(binary){codex=bitpit::VTKFormat::APPENDED;}

    ivector1D conn(m_nPoints);
    for(int i=0; i<m_nPoints; i++){
        conn[i] = i;
    }

    dvector1D sens(m_nPoints, 1.0);
    for(int i=0; i<m_nPoints; ++i){
        sens[i] = interpolateSensitivity(m_points[i]);
    }
    
    bitpit::VTKUnstructuredGrid vtk(dir, file, bitpit::VTKElementType::VERTEX);
    vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, m_points) ;
    vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, conn) ;
    vtk.setDimensions( m_nPoints, m_nPoints);
    vtk.addData("sensitivity", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT,sens);
    vtk.setCodex(codex);
    if(counter>=0){vtk.setCounter(counter);}

    vtk.write();
}

/*!
 * Decimate a cloud of points with number greater than m_nPoints, to desired value.
 * Regularize distribution of nodes to meet minimum distance & max sensitivity possible requirements of the class.
 * First point of the list is meant as the starting seed of decimation.
 * \param[in] list of 3D points in space, with size > m_nPoints
 * \return decimated list of points
 */
dvecarr3E
CreateSeedsOnSurface::decimatePoints(dvecarr3E & list){

    dvecarr3E result(m_nPoints);
    int listS = list.size();

    //reference all coordinate list to the seed candidate 0;
    // modulate the coordinate of each point with its respective sensitivity.
    dvecarr3E listRefer(listS, {{0.0,0.0,0.0}});
    for(int j=0; j<listS; ++j){
        listRefer[j] = (list[j] - list[0])*interpolateSensitivity(list[j]);
    }

    bitpit::KdTree<3,darray3E,long>    kdT;
    //build a kdtree of points
    kdT.nodes.resize(listS + kdT.MAXSTK);
    long label=0;
    for(auto & val : listRefer ){
        kdT.insert(&val, label);
        ++label;
    }

    long candidate = 0; //starting seed;

    //use kdTree to search points within a radius of m_minDist from the seed.
    // if no matches are found, choose randomly the next candidate
    livector1D neighs, excl,effective;
    std::set<long>    visited;
    std::map<double, long> chooseCand;
    livector1D finalCandidates;

    finalCandidates.reserve(m_nPoints);
    int candSize = finalCandidates.size();
    while( candSize < m_nPoints){

        finalCandidates.push_back(candidate);
        candSize++;
        visited.insert(candidate);
        excl.insert(excl.end(), visited.begin(), visited.end());
        kdT.hNeighbors(&listRefer[candidate], m_minDist, &neighs, & excl );

        visited.insert(neighs.begin(), neighs.end());

        int sizeN = visited.size();
        neighs.clear();
        excl.clear();


        while(sizeN == listS){
            m_minDist *= 0.9;
            visited.clear();
            visited.insert(finalCandidates.begin(), finalCandidates.end());
            for(const auto & ind : finalCandidates){
                excl.insert(excl.end(), visited.begin(), visited.end());
                kdT.hNeighbors(&listRefer[ind], m_minDist, &neighs, & excl );
                visited.insert(neighs.begin(), neighs.end());
            }
            sizeN = visited.size();
            excl.clear();
            neighs.clear();
        }

        effective.reserve(list.size() - visited.size());
        std::set<long>::iterator it1 = visited.begin();
        for(int i=0; i< listS; ++i){
            if( !visited.count(i) )
                effective.push_back(i);
        }

        for(const auto & ind : effective){
            dvector1D norms(finalCandidates.size());
            int countNorms=0;
            for(const auto & indEx: finalCandidates){
                norms[countNorms] = norm2(list[ind] - list[indEx]);
                ++countNorms;
            }

            double norm = 0.0, variance=0.0;
            for(auto &val : norms)    norm += val;
            norm /= double(norms.size());
            for(auto &val : norms)    val  = std::abs(val/norm - 1.0);
            maxval(norms, variance);
            norm /= std::pow((variance+1.0),2);

            chooseCand[norm] = ind;
        }

        candidate = (chooseCand.rbegin())->second;

        effective.clear();
        chooseCand.clear();
    }

    int counter = 0;
    for(const auto & index : finalCandidates){
        result[counter] = list[index];
        ++counter;
    }

    return result;
};


/*!
 * Update m_sdf distance field on a target node of a superficial tessellation solving 
 * the Eikonal equation |grad(u)| = g, using  a fast marching method. Tessellation must be
 * mandatorily a triangular one. 
 * \param[in] g       propagation speed of the Eikonal equation
 * \param[in] s       flag for inwards/outwards propagation (s = -+1)
 * \param[in] tVert   id of the target node in bitpit::PatchKernel indexing
 * \param[in] tCell   id of the cell which the Itarget belongs to in bitpit::PatchKernel indexing
 * \param[in] tri     reference to target triangulated surface.
 * \param[in] flag    flag vector reporting eikonal front advancing status on nodes. Using dead(= 0), alive(= 1),and far away(= 2) identifiers.
 * @param[in] field   reference distance field 
 * \return    updated value of the m_sdf distance field on the target node.
 */
double
CreateSeedsOnSurface::updateEikonal(double g, double s, long tVert,long tCell, bitpit::PatchKernel &tri, std::unordered_map<long int, short int> &flag, dmpvector1D & field){

    BITPIT_UNUSED(s);

    livector1D                oneRing;
    long                    I, U, V, W;

    int                        k, m;

    int                        select = -1;
    double                    value(1.0e+18);
    double                    dVU, dVW, dWU, dVP;
    std::array<double,3>    eVU, eVW, eWU, eVP, P;


    double                    a, b, c, A, B, C, K, discr ;
    double                    phi_U, phi_W, phi_P;
    int                     discrType;

    double                    xi1, xi2;
    double                    tempVal1, tempVal2;

    //get current field value of the node;
    value = std::abs(field[tVert]);

    V = tVert;

    {
        // find 1-Ring cells around target node
        bitpit::Cell & targetCell = tri.getCell(tCell);
        int locVert = targetCell.findVertex(tVert);
        if(locVert == bitpit::Vertex::NULL_ID) return value;
        oneRing = tri.findCellVertexOneRing(tCell, locVert);
    }

    // Loop over cells in the 1-Ring --------------------------------------------------- //
    for (auto && oneIndex : oneRing) {

        // Cell data get id of vertex composing triangular cell
        I = oneIndex;
        bitpit::Cell & cellI = tri.getCell(I);
        long * connCellI = cellI.getConnect();
        k = cellI.findVertex(V);
        k = (k + 1) % cellI.getVertexCount();
        U = connCellI[k];
        m = (k + 1) % cellI.getVertexCount();
        W = connCellI[m];

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
            eVU = tri.getVertexCoords(V) - tri.getVertexCoords(U);
            dVU = norm2(eVU);
            value = std::min(value, std::abs(field[U]) + g*dVU); ///??????????????????????????????????????????
            break;

        case 2 : // with 2 dead nodes

            eVW = tri.getVertexCoords(V) - tri.getVertexCoords(W);
            dVW = norm2(eVW);
            eVW = eVW/dVW;
            eVU = tri.getVertexCoords(V) - tri.getVertexCoords(U);
            dVU = norm2(eVU);
            eVU = eVU/dVU;
            eWU = tri.getVertexCoords(W) - tri.getVertexCoords(U);
            dWU = norm2(eWU);
            eWU = eWU/dWU;

            // Coeffs -------------------------------------------------------------------- //
            phi_U = std::abs(field[U]);
            phi_W = std::abs(field[W]);
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
                P = (1.0 - xi1) * tri.getVertexCoords(U)  +  xi1 * tri.getVertexCoords(W);
                eVP = tri.getVertexCoords(V) - P;
                dVP = norm2(eVP);
                eVP = eVP/dVP;
                phi_P = (1.0 - xi1) * phi_U + xi1 * phi_W;
                tempVal1 = phi_P + g * dVP;

                // Solution #2
                P = (1.0 - xi2) * tri.getVertexCoords(U)  +  xi2 * tri.getVertexCoords(W);
                eVP = tri.getVertexCoords(V) - P;
                dVP = norm2(eVP);
                eVP = eVP/dVP;
                phi_P = (1.0 - xi2) * phi_U + xi2 * phi_W;
                tempVal2 = phi_P + g * dVP;

                break;

            case 2: // coincident solutions
                discr = std::abs(discr);

                // Solution #1
                P = tri.getVertexCoords(U);
                eVP = tri.getVertexCoords(V) - P;
                dVP = norm2(eVP);
                eVP = eVP/dVP;
                phi_P = phi_U;
                tempVal1 = phi_P + g*dVP;

                // Solution #2
                P = tri.getVertexCoords(W);
                eVP = tri.getVertexCoords(V) - P;
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

    return(value);
}; 

/*!
 * Solve the 3D Eikonal equation |grad(u)| = g, using  a fast marching method, on a target triangulation 
 * associated unstructured superficial grid. Tessellation must be mandatorily a triangular one.
 * (SurfaceConstraints::m_sdf values in the unknown region must be set to the value 1.0e+18.
 *  Dead vertices of front to be propagated must be set to zero, initially)
 * \param[in] g Propagation speed.
 * \param[in] s Velocity sign (+1 --> propagate outwards, -1 --> propagate inwards).
 * \param[in] tri reference to target triangulated surface.
 * \param[in] invConn inverse connectivity of target triangulated surface
 * \param[in,out] field field to be computed, already allocated.
 */
void
CreateSeedsOnSurface::solveEikonal(double g, double s,bitpit::PatchKernel &tri, std::unordered_map<long,long> & invConn, dmpvector1D & field ){
    
    // declare total size and support structure
    long     N(tri.getVertexCount());
    std::unordered_map<long int, short int> active;

    std::unordered_map<long, int> vmap;
    int countV = 0;
    for(const auto & vert: tri.getVertices()){
        vmap[vert.getId()] = countV;
        ++countV;
    }
    
    { //FLAG DEAD/ALIVE/FAR-AWAY VERTICES

        long myId;
        bool check;
        std::set<long>                neighs;
        std::set<long>::iterator    it, itend;

        //set active vector size and mark its position  with original geometry ids
        active.reserve(N);
        for ( const auto &vertex : tri.getVertices() ){
            myId           = vertex.getId() ;
            active[myId] = 2 ;
        }

        //fill active vector
        for ( const auto &vertex : tri.getVertices() ){
            myId     =    vertex.getId();

            // Dead vertices
            if( isDeadFront(myId) ){
                active[myId] = 0;

            }else{

                //loop over neighbors
                check = false;
                neighs = findVertexVertexOneRing(tri,invConn[myId], myId);
                it = neighs.begin();
                itend = neighs.end();
                while(!check && it !=itend){

                    check = s*field[*it] >= 0.0 && field[*it] < 1.0E+18;
                    ++it;
                };

                active[myId] = 2 - (int) check;
            }
        }
    }


    { // Construct min heap data structure
        long                            m(0), I(0), myId, J ;
        double                          value ;

        std::set<long>                neighs;
        std::set<long>::iterator    it,itbeg, itend;

        std::vector<std::array<int,2>>  map(N), *mapPtr = &map;

        bitpit::MinPQueue<double, long> heap(N, true, mapPtr);

        // Inserting alive vertices in  the heap

        for(const auto & vertex : tri.getVertices()){

            myId = vertex.getId();
            if(active[myId] == 1) {
                //assign a value to your actual vertex
                value = updateEikonal(s, g, myId, invConn[myId], tri, active, field);

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
            field[myId] = value;

            // Update flag to dead;
            active[myId] = 0;

            //update neighbours
            neighs = findVertexVertexOneRing(tri, invConn[myId], myId);
            itbeg = neighs.begin();
            itend = neighs.end();

            for(it=itbeg; it != itend; ++it) {
                J = *it;

                if(active[J] == 1){

                    //update local value;
                    value = updateEikonal(s,g,J,invConn[J], tri, active, field);

                    //update its value in the min-heap
                    I = vmap[J];
                    heap.modify( map[I][1],value,J );

                }else if( active[J] == 2){

                    //update local value;
                    value = updateEikonal(s,g,J, invConn[J],tri, active, field);

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
};

/*!
 * Get a minimal inverse connectivity of the target geometry mesh associated to the class.
 * Each vertex (passed by Id) is associated at list to one of the possible simplex 
 * (passed by Id) which it belongs. This is returned in an unordered_map having as key the 
 * vertex Id and as value the Cell id. Id is meant as the unique label identifier associated
 * to bitpit::PatchKernel original geometry
 * \param[in] geo reference to target surface geometry
 *\return    unordered_map of vertex ids (key) vs cell-belonging-to ids(value)
 */
std::unordered_map<long,long>
CreateSeedsOnSurface::getInverseConn(bitpit::PatchKernel & geo){

    std::unordered_map<long,long> invConn ;

    long cellId;
    for(const auto &cell : geo.getCells()){
        cellId = cell.getId();
        auto vList = cell.getVertexIds();
        for(const auto & idV : vList) invConn[idV] = cellId;
    }

    return(invConn);
};

/*!
 * Return true if a given vertex belongs to the current constrained boundary front of your patch
 * \param[in]    label index of vertex, in sequential mimmo::MimmoObject notation
 * \return boolean, true if vertex belongs to constrained set, false if not 
 */
bool CreateSeedsOnSurface::isDeadFront(const long int label){

    livector1D::iterator got = std::find(m_deads.begin(), m_deads.end(), label);
    if(got == m_deads.end()) return false;
    return true;
}

/*!
 * Return VertexVertex One Ring of a specified target vertex
 * \param[in]    geo        target surface geometry
 * \param[in]    cellId     bitpit::PatchKernel Id of a cell which target belongs to
 * \param[in]    vertexId    bitpit::PatchKernel Id of the target vertex
 * \return        list of all vertex in the One Ring of the target, by their bitpit::PatchKernel Ids
 */
std::set<long>
CreateSeedsOnSurface::findVertexVertexOneRing(bitpit::PatchKernel &geo, const long & cellId, const long & vertexId){

    std::set<long> result;
    bitpit::Cell &cell =  geo.getCell(cellId);

    int loc_target = cell.findVertex(vertexId);
    if(loc_target == bitpit::Vertex::NULL_ID) return result;

    livector1D list = geo.findCellVertexOneRing(cellId, loc_target);

    long connSize;
    for(const auto & index : list){
        bitpit::Cell & cell = geo.getCell(index);
        auto vList = cell.getVertexIds();
        for(const auto & idV : vList){
            result.insert(idV);
        }
    }

    result.erase(vertexId);
    return result;
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
CreateSeedsOnSurface::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name ){

    BITPIT_UNUSED(name);
    
    BaseManipulation::absorbSectionXML(slotXML, name);
    
    if(slotXML.hasOption("NPoints")){
        std::string input = slotXML.get("NPoints");
        input = bitpit::utils::string::trim(input);
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
        input = bitpit::utils::string::trim(input);
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
        input = bitpit::utils::string::trim(input);
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
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setMassCenterAsSeed(value);
    }

    if(slotXML.hasOption("RandomFixed")){
        std::string input = slotXML.get("RandomFixed");
        input = bitpit::utils::string::trim(input);
        int value = -1;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setRandomFixed(value);
    }


};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
CreateSeedsOnSurface::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);
    
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
    if( signat != -1 && m_engine == mimmo::CSeedSurf::RANDOM){
        slotXML.set("RandomFixed", std::to_string(signat));
    }


};

/*!
 * Interpolate sensitivity field on a 3D point
 * \param[in] point belonging to the target geometry surface.
 */
double
CreateSeedsOnSurface::interpolateSensitivity(darray3E & point){

    MimmoObject* geo = getGeometry();
    long supportCell = skdTreeUtils::locatePointOnPatch(point, *(geo->getSkdTree()));
    if(supportCell == bitpit::Cell::NULL_ID)    return 0.0;

    bitpit::Cell & cell = geo->getPatch()->getCell(supportCell);
    bitpit::ConstProxyVector<long> nVList = cell.getVertexIds();
    dvector1D weights(nVList.size(), 0), val(nVList.size(),0);
    double wtot = 0.0;
    int countV = 0;
    for(const auto & idLoc : nVList){
        if(m_sensitivity.exists(idLoc)){
            val[countV] = m_sensitivity[idLoc];
        }else{
            val[countV] = 0.0;
        }
        double valdist = norm2(geo->getVertexCoords(idLoc) - point);
        if ( valdist< 1.E-18){ 
            return val[countV];
        }
        
        weights[countV] = 1.0/valdist;
        wtot += weights[countV];
        ++countV;
    }
    
    weights /= wtot;
    
    double result = 0.0;
    for(int i=0; i<countV; ++i){
        result +=  weights[i]*val[i];
    }

    return result;
}

/*!
 * Check your current point data sensitivity field associated to linked geometry. 
 * Do nothing if class linked geometry is a null pointer or empty.
 * If field geometry is not linked or empty or if field geometry its uncoherent with class linked 
 * geometry or not referred to MPVLocation point or if the field itself does not carry any value, create a default unitary field.
 * If the field is still coherent but miss values on some geometry nodes, complete it assigning zero value on missing ids.
 * If no exception occurs, leave the field as it is.
 */
void CreateSeedsOnSurface::checkField(){
    
    if(getGeometry() == NULL)   return;
    if(getGeometry()->isEmpty()) return;
    dmpvector1D defaultField;
    defaultField.setGeometry(getGeometry());
    defaultField.setDataLocation(MPVLocation::POINT);
    //create unity field;
    for(const auto vert: getGeometry()->getVertices()){
        defaultField.insert(vert.getId(), 1.0);
    }
    
    m_log->setPriority(bitpit::log::Verbosity::DEBUG);
    if(getGeometry() != m_sensitivity.getGeometry()){
        m_sensitivity = defaultField;
        (*m_log)<<"warning in "<<m_name<<" : Not suitable data field connected. Reference geometry linked by the class and by MimmoPiercedvector field differs. Using default field"<<std::endl;
        return;
    }else if(m_sensitivity.getDataLocation() != MPVLocation::POINT){
        m_sensitivity = defaultField;
        (*m_log)<<"warning in "<<m_name<<" : Wrong data field connected. Data Location of MimmoPiercedvector field is not referred geometry Vertex/Point. Using default field"<<std::endl;
        return;
    }else{    
        if(!m_sensitivity.completeMissingData(0.0)){
            m_sensitivity = defaultField;
            (*m_log)<<"warning in "<<m_name<<" : Not coherent data field connected. Data Ids of MimmoPiercedvector field are not aligned with geometry Vertex/Point field. Using default field"<<std::endl;
            return;
        }
    }
    m_log->setPriority(bitpit::log::Verbosity::NORMAL);
    
}


/*!
 * Normalize your current point data sensitivity field associated to linked geometry. 
 * It is assumed at this point that the field is fully coeherent, as result of checkField method check.
 */
void CreateSeedsOnSurface::normalizeField(){
    
    double minSense=0.0,maxSense=0.0;
    minval(m_sensitivity.getRawDataAsVector(), minSense);
    //operate translation.
    if (!std::isnan(minSense)){
        for(auto &val: m_sensitivity){
            val += -1.0*minSense;
        }
    }
    maxval(m_sensitivity.getRawDataAsVector(), maxSense);
    //operate normalization.
    if (!std::isnan(maxSense) || maxSense != 0.0){
        for(auto &val: m_sensitivity){
            val /= maxSense;
        }
    }

    m_log->setPriority(bitpit::log::Verbosity::DEBUG);
    if(std::isnan(minSense) || std::isnan(maxSense) || maxSense == 0.0){
        (*m_log)<<"warning in "<<m_name<<" : Not valid data of sensitivity field detected. Using default unitary field"<<std::endl;
        dmpvector1D defaultField;
        defaultField.setGeometry(getGeometry());
        defaultField.setDataLocation(MPVLocation::POINT);
        //create unity field;
        for(const auto vert: getGeometry()->getVertices()){
            defaultField.insert(vert.getId(), 1.0);
        }
        m_sensitivity = defaultField;
    }
    m_log->setPriority(bitpit::log::Verbosity::NORMAL);
}


/*!
 * \return true if the linked geometry is a homogeneous triangular 3D surface, false otherwise
 */
bool
CreateSeedsOnSurface::checkTriangulation(){
    if(getGeometry()->getType() != 1) return false;
    bool check = true;
    bitpit::PatchKernel::CellIterator it = getGeometry()->getPatch()->cellBegin();
    bitpit::PatchKernel::CellIterator itEnd = getGeometry()->getPatch()->cellEnd();
    while(check && it != itEnd){
        check = ( (*it).getType() == bitpit::ElementType::TRIANGLE );
        ++it;
    }
    return check;
}

/*!
 * \return a homogeneous triangulated and indipendent clone of the current target geometry linked to the class.
 */
std::unique_ptr<MimmoObject>
CreateSeedsOnSurface::triangulate(){

    std::unique_ptr<MimmoObject> temp = getGeometry()->clone();
    m_sensitivity_triangulated = m_sensitivity;
    m_sensitivity_triangulated.setGeometry(temp.get());

    long maxID, newID, newVertID;
    
    const auto orderedCellID = temp->getCells().getIds(true);
    maxID = orderedCellID[(int)orderedCellID.size()-1];
    newID = maxID+1;
    {
        const auto orderedVertID = temp->getVertices().getIds(true);
        newVertID = orderedVertID[(int)orderedCellID.size()-1] +1;
    }
    
    bitpit::ElementType eletype;
    bitpit::ElementType eletri = bitpit::ElementType::TRIANGLE;
    livector1D connTriangle(3);
    for(const auto &idcell : orderedCellID){
        
        livector1D conn = temp->getCellConnectivity(idcell);
        eletype = temp->getPatch()->getCell(idcell).getType();
        short pid = temp->getPatch()->getCell(idcell).getPID();
        
        switch (eletype){
            case bitpit::ElementType::QUAD:
            case bitpit::ElementType::PIXEL:
            {
                temp->getPatch()->deleteCell(idcell);
                for(std::size_t i=0; i<2; ++i){
                    connTriangle[0] = conn[0];
                    connTriangle[1] = conn[i+1];
                    connTriangle[2] = conn[i+2];
                    temp->addConnectedCell(connTriangle, eletri, pid, newID);
                    ++newID;
                }
            }
                break;
            case bitpit::ElementType::POLYGON:
            {
                std::size_t startIndex = 1;
                std::size_t nnewTri = conn.size() - startIndex;
                //calculate barycenter and add it as new vertex
                darray3E barycenter = temp->getPatch()->evalCellCentroid(idcell);
                temp->addVertex(barycenter, newVertID);
                // adding new vertex, adding also a sensitivity exstimation on new point.
                double sens_new = interpolateSensitivity(barycenter);
                m_sensitivity_triangulated.insert(newVertID, sens_new);
                //delete current polygon
                temp->getPatch()->deleteCell(idcell);
                //insert new triangles from polygon subdivision
                for(std::size_t i=0; i<nnewTri; ++i){
                    connTriangle[0] = newVertID;
                    connTriangle[1] = conn[ startIndex + std::size_t( i % nnewTri) ];
                    connTriangle[2] = conn[ startIndex + std::size_t( (i+1) % nnewTri ) ];
                    temp->addConnectedCell(connTriangle, eletri, pid, newID);
                    ++newID;
                }
                //increment label of vertices
                ++newVertID;
                
            }
                break;
            case bitpit::ElementType::TRIANGLE:
                //do nothing
                break;    
            default:
                throw std::runtime_error("unrecognized cell type in 3D surface mesh of CreateSeedsOnSurface");
                break;
        }
    }
    return std::move(temp);
}

}