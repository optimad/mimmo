/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
#include "ControlDeformExtSurface.hpp"
#include "SkdTreeUtils.hpp"
#include <volcartesian.hpp>
#include <CG.hpp>

namespace mimmo{

/*!
 * Default constructor of ControlDeformExtSurface
 */
ControlDeformExtSurface::ControlDeformExtSurface(){
    m_name = "mimmo.ControlDeformExtSurface";
    m_tolerance = 0.0;

    m_allowed.insert((FileType::_from_string("STL"))._to_integral());
    m_allowed.insert((FileType::_from_string("SURFVTU"))._to_integral());
    m_allowed.insert((FileType::_from_string("NAS"))._to_integral());

};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ControlDeformExtSurface::ControlDeformExtSurface(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.ControlDeformExtSurface";
    m_tolerance = 0.0;
    m_allowed.insert((FileType::_from_string("STL"))._to_integral());
    m_allowed.insert((FileType::_from_string("SURFVTU"))._to_integral());
    m_allowed.insert((FileType::_from_string("NAS"))._to_integral());

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.ControlDeformExtSurface"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!Default destructor of ControlDeformExtSurface
 */
ControlDeformExtSurface::~ControlDeformExtSurface(){};

/*!Copy constructor of ControlDeformExtSurface. Deformation field referred to geometry
 * and result violation field are not copied.
 */
ControlDeformExtSurface::ControlDeformExtSurface(const ControlDeformExtSurface & other):BaseManipulation(other){
    m_allowed = other.m_allowed;
    m_geoList = other.m_geoList;
    m_geoFileList = other.m_geoFileList;
    m_tolerance = other.m_tolerance;
};

/*!
 * Assignment operator of ControlDeformExtSurface. Deformation field referred to geometry
 * and result violation field are not copied.
 */
ControlDeformExtSurface & ControlDeformExtSurface::operator=(ControlDeformExtSurface other){
    swap(other);
    return(*this);
};

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void ControlDeformExtSurface::swap(ControlDeformExtSurface & x) noexcept
{
    std::swap(m_allowed, x.m_allowed);
    std::swap(m_geoList, x.m_geoList);
    std::swap(m_geoFileList, x.m_geoFileList);
    std::swap(m_tolerance, x.m_tolerance);
    m_violationField.swap(x.m_violationField);
    m_defField.swap(x.m_defField);
    BaseManipulation::swap(x);
}
/*! It builds the input/output ports of the object
 */
void
ControlDeformExtSurface::buildPorts(){
    bool built = true;

    built = (built && createPortIn<dmpvecarr3E*, ControlDeformExtSurface>(this, &mimmo::ControlDeformExtSurface::setDefField, M_GDISPLS));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, ControlDeformExtSurface>(this, &mimmo::ControlDeformExtSurface::setGeometry, M_GEOM, true));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, ControlDeformExtSurface>(this, &mimmo::ControlDeformExtSurface::addConstraint, M_GEOM2));

    built = (built && createPortOut<double, ControlDeformExtSurface>(this, &mimmo::ControlDeformExtSurface::getViolation, M_VALUED));
    built = (built && createPortOut<dmpvector1D*, ControlDeformExtSurface>(this, &mimmo::ControlDeformExtSurface::getViolationField, M_SCALARFIELD));
    m_arePortsBuilt = built;
};

/*!
  Return the value of violation of deformed geometry, after class execution.
  If value is positive or at least zero, constraint violation occurs, penetrating or touching at
  least in one point the constraint sourface. Otherwise, returning negative values means
  that no violation occurs.
  \return violation value
 */
double
ControlDeformExtSurface::getViolation(){

    double result = -1.0*std::numeric_limits<double>::max();
    for(const auto & val : m_violationField){
        result = std::fmax(result, val);
    }
#if MIMMO_ENABLE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_MAX, m_communicator);
#endif
    return    result;
};

/*!
  Return the violation field relative to deformed target geometry, onto the geometry itself.
  The info is available after class execution.
  If value is positive or at least zero, constraint violation occurs, penetrating or
  touching at least in one point of the constraint sourface. Otherwise, returning negative values
  means that no violation occurs.
  \return violation field values
 */
dmpvector1D *
ControlDeformExtSurface::getViolationField(){
    return  &m_violationField;
};


/*!
 * Get proximity tolerance for constraints fixed into the class.
  Tolerance acts as an offset treshold, so that a deformed body can be
  considered as violating constraints when it's dangerously close,
  even if it's not really touching/penetrating them.
  \return     fixed proximity tolerance
 */
double
ControlDeformExtSurface::getTolerance(){
    return m_tolerance;
}

/*!
 * Return the actual list of external geometry files selected as constraint to check your deformation.
   Only constraints specified by files are returned.
 * \return list of external geometries files with their type
 */
const ControlDeformExtSurface::fileListWithType &
ControlDeformExtSurface::getConstraintFiles() const{
    return    m_geoFileList;
}

/*!
 * Set the deformative field associated to each point of the target geometry.
 * Field resize occurs in execution, if point dimension between field and geoemetry does not match.
 * \param[in]    field of deformation
 */
void
ControlDeformExtSurface::setDefField(dmpvecarr3E *field){
    if(!field)  return;
    m_defField.clear();
    m_defField = *field;
};

/*!
 * Set link to target geometry for your selection.
 * \param[in] target pointer to target geometry
 */
void
ControlDeformExtSurface::setGeometry( mimmo::MimmoSharedPointer<mimmo::MimmoObject> target){
    if(target == nullptr)       return;
    if(target->getType() != 1){
        (*m_log)<<"warning: ControlDeformExtSurface cannot support current geometry. It works only w/ 3D surface."<<std::endl;
        return;
    }
    m_geometry = target;

};

/*!
 * Set proximity tolerance for constraints (double >= 0).
  Tolerance acts as an offset treshold, so that a deformed body can be
  considered as violating constraints when it's dangerously close,
  even if it's not really touching/penetrating them.
  \param[in]   tol proximity tolerance
 */
void
ControlDeformExtSurface::setTolerance(double tol){
    m_tolerance = std::max(0.0, tol);
}

/*!
 * Add a surface geometry as constraint for violation control.
  For MPI version, this is the only method to provide partitioned external constraint
  geometries
 * \param[in] constraint pointer to constraint geometry
 */
void
ControlDeformExtSurface::addConstraint( mimmo::MimmoSharedPointer<mimmo::MimmoObject> constraint){
    if(constraint == nullptr)       return;
    if(constraint->getType() != 1){
        (*m_log)<<"warning: ControlDeformExtSurface cannot support current constraint geometry. It works only w/ 3D surface."<<std::endl;
        return;
    }
    m_geoList.insert(constraint);
};



/*!
 * Set a list of external geometry files to check your deformation violation.
   The int format of the file (see FileType enum) must be present as second argument.
   This is a method to pass constraint geometries to the class, alternative to addConstraint.
   In case of MPI version, please be aware that since these geometries are read from
   file, they will be retained on the rank 0 processor only and not properly distributed
   among ranks. For a full control of distributed geometries, partition them externally and
   add them to the class using addConstraint method.
 * \param[in] files list external geometries to be read
 */
void
ControlDeformExtSurface::setConstraintFiles(ControlDeformExtSurface::fileListWithType files){
    for(auto & val : files){
        addConstraintFile(val.first, val.second);
    }
};

/*!
 * As setConstraintFiles, but adding a new constraint geometry file once at the time.
 *\param[in] file of external geometry to be read
 *\param[in] format  type of file as integer (see Filetype enum). only nas, stl, surfvtu and nas are supported.
 */
void
ControlDeformExtSurface::addConstraintFile(std::string file, int format){
    if(m_allowed.count(format)>0){
        m_geoFileList[file] = format;
    }
};

/*!
* As setConstraintFiles, but adding a new constraint geometry file once at the time.
 *\param[in] file of external geometry to be read
 *\param[in] format  type of file as  Filetype enum. only nas, stl, stvtu and sqvtu are supported.
 */
void     ControlDeformExtSurface::addConstraintFile(std::string file, FileType format){
    addConstraintFile(file,format._to_integral());
};

/*!
 * Remove an existent file from the list of external geometry files. If not in the list, do nothing
 * \param[in] file to be removed from the list
 */
void
ControlDeformExtSurface::removeConstraintFile(std::string file){
    if(m_geoFileList.count(file) >0)    m_geoFileList.erase(file);
};

/*!
 * Empty your list of external constraint geometry files
 */
void
ControlDeformExtSurface::removeConstraintFiles(){
    m_geoFileList.clear();
};

/*!
 * Reset class to default parameters.
 */
void
ControlDeformExtSurface::clear(){
    removeConstraintFiles();
    m_geoList.clear();
    m_defField.clear();
    m_violationField.clear();
    m_tolerance = 0.0;
    BaseManipulation::clear();
};

/*!Execution command.
  Calculate violation value and store it in the class member m_violationField
 */
void
ControlDeformExtSurface::execute(){

    MimmoSharedPointer<MimmoObject> geo = getGeometry();
    if(geo == nullptr){
        (*m_log)<<"Error "+ m_name + " : null pointer to linked geometry found"<<std::endl;
        throw std::runtime_error("Error "+m_name + " : null pointer to linked geometry found");

    }
    if(geo->isEmpty()){
        (*m_log)<<"Warning in " + m_name +" : empty linked geometry found"<<std::endl;
    }

    bool check = m_defField.getGeometry() == geo;
    check = check && m_defField.getDataLocation() == MPVLocation::POINT;
    if(!check){
        (*m_log)<<"Warning in "<<m_name<<": Unsuitable deformation field linked"<<std::endl;
    }
    m_defField.completeMissingData({{0.0,0.0,0.0}});

    m_violationField.clear();
    m_violationField.initialize(geo, MPVLocation::POINT, -1.0E25);

    std::vector<long> pointIds = geo->getVerticesIds();
    dvecarr3E points(pointIds.size());

    //adding deformation to points and calculate their local
    //barycenter
    int count = 0;
    darray3E dbMin, dbMax;
    dbMin.fill(std::numeric_limits<double>::max());
    dbMax.fill(-1.0*std::numeric_limits<double>::max());
    for (long id : pointIds){
        pointIds[count] = id;
        points[count] = geo->getVertexCoords(id) + m_defField[id];
        for (int i=0; i<3; ++i){
            dbMin[i] = std::min(dbMin[i], points[count][i]);
            dbMax[i] = std::max(dbMax[i], points[count][i]);
        }
        ++count;
    }
    darray3E loc_dcenter = 0.5*(dbMin + dbMax);


    //evaluate global AABB of geometry and save its center
    darray3E pbMin, pbMax;
    getGlobalBoundingBox(geo, pbMin, pbMax);
    darray3E originalGlobalCenter = 0.5*(pbMin + pbMax);

    //put together constrained surfaces in a unique list
    std::vector<MimmoSharedPointer<MimmoObject>> constraint_geos;
    //from files first
    readFileConstraints(constraint_geos);
    // then add those directly linked with addConstraint method.
    constraint_geos.insert(constraint_geos.end(), m_geoList.begin(), m_geoList.end());

    //check list size of constraints and rise a warning if the list is empty.
    if(constraint_geos.empty()) {
        (*m_log)<<"Warning in " + m_name +" : no valid constraint geometries are linked to the class. "<<std::endl;
    }

    // start examining one by one all external constraints
    bool checkOpen;
    darray3E bbMin, bbMax;

    for(MimmoSharedPointer<MimmoObject> localg : constraint_geos){

        //check constraint skdtree
        localg->buildSkdTree();
        //get the global bounding box of the constraint surface
        getGlobalBoundingBox(localg, bbMin, bbMax);
        //evaluate a guess of the search radius as distance using AABBs of the global
        //target and the global constraint.
        darray3E insMin, insMax;
        double suppval(1.0E-08);
        if(bitpit::CGElem::intersectBoxBox(pbMin,pbMax, bbMin, bbMax, insMin, insMax)){
           suppval =  std::max(1.0E-08, 0.5*norm2(insMax - insMin) );
        }
        double searchRadius = std::max( suppval, norm2(originalGlobalCenter - 0.5*(bbMin + bbMax)) );

        //calculate where the center of the global target is located w.r.t to constraint, to keep up the
        //right sign for the distance.
        double refsign = 1.0;
        std::vector<double> distances;
        distances.resize(1);
        evaluateSignedDistance(std::vector<darray3E>(1,originalGlobalCenter), localg, searchRadius, distances);
        if(distances[0] > 0.0) refsign *= -1.0;

        //Calculate now the real distances of the deformed points set.
        distances.resize(points.size());
        //evaluate a first guess of the search radius as distance using AABBs of the local
        //point set and the global constraint.
        suppval = 1.0E-08;
        if(bitpit::CGElem::intersectBoxBox(dbMin,dbMax, bbMin, bbMax, insMin, insMax)){
           suppval =  std::max(1.0E-08, 0.5*norm2(insMax - insMin) );
        }
        searchRadius = std::max( suppval, norm2(loc_dcenter - 0.5*(bbMin + bbMax)) );

        evaluateSignedDistance(points, localg, searchRadius, distances);

        distances *= refsign;

        //push values directly in m_violationField.
        count = 0;
        for(long id : pointIds){
            m_violationField[id] = std::max( m_violationField[id], (distances[count] + m_tolerance) );
            ++count;
        }

    }//end looping on constraint geometries.

    //write resume file.
    writeLog();

};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
ControlDeformExtSurface::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::absorbSectionXML(slotXML, name);

    std::unordered_map<std::string, int > mapp;
    if(slotXML.hasSection("Files")){

        const bitpit::Config::Section & filesXML = slotXML.getSection("Files");

        for(auto & subfile : filesXML.getSections()){
            std::string path;
            std::string tag, tolstring;
            double value = 1.E-8;

            if(subfile.second->hasOption("fullpath"))    {
                path = subfile.second->get("fullpath");
                path = bitpit::utils::string::trim(path);
            }
            if(subfile.second->hasOption("tag")){
                tag = subfile.second->get("tag");
                tag = bitpit::utils::string::trim(tag);
                //check tag;
                auto maybe_tag = FileType::_from_string_nothrow(tag.c_str());
                if(!maybe_tag)    tag.clear();
                else    tag = maybe_tag->_to_string();
            }
            if(!path.empty() && !tag.empty()){
                mapp[path] = int(FileType::_from_string(tag.c_str()));
            }
        }

        setConstraintFiles(mapp);

    }

    if(slotXML.hasOption("Tolerance")){
        std::string input = slotXML.get("Tolerance");
        input = bitpit::utils::string::trim(input);
        double value = 0.0;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setTolerance(value);
    }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
ControlDeformExtSurface::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);
    bitpit::Config::Section & filesXML = slotXML.addSection("Files");

    int counter = 0;
    for(auto & file : m_geoFileList){
        std::string name = "file"+std::to_string(counter);
        bitpit::Config::Section & local = filesXML.addSection(name);
        local.set("fullpath", file.first);
        std::string typetag = (FileType::_from_integral(file.second))._to_string();
        local.set("tag", typetag);
        ++counter;
    }

    slotXML.set("Tolerance", std::to_string(m_tolerance));

};

/*!
 * Read external constraints provided by files (m_geoFileList) and return it
 * in a list of shared pointers to MimmoObjects.
   For MPI versions, every rank will have its geometry container, mesh info
   will be retained on rank 0 only
 * \param[in,out] extFileGeos list of geometry containers, empty in input, filled in output
 */
void
ControlDeformExtSurface::readFileConstraints(std::vector<MimmoSharedPointer<MimmoObject> > & extFileGeos){

    extFileGeos.resize(m_geoFileList.size());

    int counter = 0;
    std::unique_ptr<MimmoGeometry> geo (new MimmoGeometry(MimmoGeometry::IOMode::READ));

    for(auto & geoinfo : m_geoFileList){

        svector1D info = extractInfo(geoinfo.first);
        geo->setDir(info[0]);
        geo->setFilename(info[1]);
        geo->setFileType(geoinfo.second);
        geo->execute();

        MimmoSharedPointer<MimmoObject> locMesh = geo->getGeometry();
        if (locMesh == nullptr) continue;

        extFileGeos[counter] = locMesh;

        extFileGeos[counter]->updateAdjacencies();
        extFileGeos[counter]->buildSkdTree();

        ++counter;
    }
    extFileGeos.resize(counter);
};

/*!
 * Extract root dir/filename/tag from an absolute file pattern
 * \return dir/filename/tag
 */
svector1D
ControlDeformExtSurface::extractInfo(std::string file){

    std::string root, name, tag,temp;
    std::string key1=".", key2="/\\";

    std::size_t found = file.find_last_of(key2);
    root = file.substr(0, found);
    temp = file.substr(found+1);

    found = temp.find_last_of(key1);
    name = temp.substr(0,found);
    tag = temp.substr(found+1);

    svector1D result(3);
    result[0] = root;
    result[1] = name;
    result[2] = tag;

    return     result;
}

/*!
 * Internal custom Wrapper to skdTreeUtils::signedDistance. Calculate signed distance of a point cloud set
   w.r.t to a constraint geometry.
   An initial homogeneous search radius is provided, that can be enlarged step by step until all the points
   has valid distance in output.
   MPI/serial version handling is done internally.

 * \param[in] points pointer to 3D targets points list
 * \param[in] geo    target geometry w/ skdTree in it
 * \param[in] initRadius guess initial search radius.
 * \param[out] distances with sign of each point from target surface.
 */
void
ControlDeformExtSurface::evaluateSignedDistance(const std::vector<darray3E> &points, MimmoSharedPointer<MimmoObject> &geo,
                                                double initRadius, std::vector<double> & distances)
{

    geo->buildSkdTree();

    double dist = std::numeric_limits<double>::max();
    double rate = 0.05;
    int kmax = 1000;
    int kiter = 0;
    double sRadius(initRadius);
    distances.resize(points.size());

    std::vector<darray3E> work = points;
    std::unordered_map<std::size_t,std::size_t> mapPosIndex;
    for(std::size_t i=0; i<points.size(); ++i){
        mapPosIndex[i] = i;
    }

    bool checkEmptyWork = work.empty();
#if MIMMO_ENABLE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &checkEmptyWork, 1, MPI_C_BOOL, MPI_LAND, m_communicator);
#endif

    while(!checkEmptyWork && kiter < kmax){

        std::vector<long> suppCellIds(work.size());
        std::vector<double> distanceWork(work.size());
        std::vector<darray3E> normals(work.size());

#if MIMMO_ENABLE_MPI
        std::vector<int> suppCellRanks(work.size());
        skdTreeUtils::signedGlobalDistance(work.size(), work.data(), geo->getSkdTree(), suppCellIds.data(), suppCellRanks.data(), normals.data(), distanceWork.data(), sRadius, false);
#else
        skdTreeUtils::signedDistance(work.size(), work.data(), geo->getSkdTree(), suppCellIds.data(), distanceWork.data(), normals.data(), sRadius);
#endif

        //get all points with distances not calculated.
        std::vector<std::array<double,3>> failedPoints;
        std::unordered_map<std::size_t, std::size_t> failedPosIndex;
        failedPoints.reserve(work.size());
        int count(0);
        for(long & val : suppCellIds){
            if(val == bitpit::Cell::NULL_ID){
                failedPoints.push_back(work[count]);
                failedPosIndex[failedPoints.size()-1] = mapPosIndex[count];
            }else{
                distances[mapPosIndex[count]] = distanceWork[count];
            }
            ++count;
        }

        std::swap(failedPoints, work);
        std::swap(failedPosIndex, mapPosIndex);

        //check work if empty
        checkEmptyWork = work.empty();
#if MIMMO_ENABLE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &checkEmptyWork, 1, MPI_C_BOOL, MPI_LAND, m_communicator);
#endif
        //increase search radius
        sRadius *= (1.0 + rate);
        //increase iteration
        kiter++;
    }

    //check for residual unchecked distances, and put them to -1.0E+25
    // just to be manageable by visualizers like paraview.
    for(auto & tuple : mapPosIndex){
        distances[tuple.second] = -1.0E+25;
    }

}

/*!
 * Plot optional results in execution, that is the violation distance field on current deformed geometry.
 * Reimeplemented from BaseManipulation::plotOptionalResults;
 */
void
ControlDeformExtSurface::plotOptionalResults(){

    m_defField.setName("Deformation Field");
    m_violationField.setName("Violation Distance Field");
    write(getGeometry(), m_defField, m_violationField);

}

/*!
 * Write log file
 */
void
ControlDeformExtSurface::writeLog(){

    double violationMax = getViolation();
    std::size_t geoListSize = m_geoList.size();

#if MIMMO_ENABLE_MPI
    if(m_rank == 0){
#endif
        std::string logname = m_name+std::to_string(getId())+"_violation.log";
        std::ofstream log;
        log.open(logname);
        log<<"mimmo "<<m_name<<" resume file"<<std::endl;
        log<<std::endl;
        log<<std::endl;
        log<<" violation value : " << violationMax << std::endl;
        log<<std::endl;
        log<<" directly linked constraints (number) : " << geoListSize << std::endl;
        log<<std::endl;
        log<<" additional constraints from files : " ;
        for (auto ig : m_geoFileList){
            log<< ig.first <<"  "<<std::endl;
        }
        log.close();

#if MIMMO_ENABLE_MPI
    }
    MPI_Barrier(m_communicator); // other ranks stalled here, waiting 0 to finish writing.
 #endif

}


/*!
 * Wrapper to MimmoObject bounding box geometry calculation. No matter the version, return always the
   global bounding box of the object
   \param[in] geo target geometry
   \param[out] bMin min point of AABB
   \param[out] bMax max point of AABB
 */
void
ControlDeformExtSurface::getGlobalBoundingBox(MimmoSharedPointer<MimmoObject> & geo, darray3E & bMin, darray3E & bMax){
    if(geo->getBoundingBoxSyncStatus() != mimmo::SyncStatus::SYNC){
        geo->update();
    }
    geo->getBoundingBox(bMin, bMax, true);
}

}//end mimmo namespace
