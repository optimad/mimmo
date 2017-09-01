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
\*---------------------------------------------------------------------------*/
#include "ControlDeformExtSurface.hpp"
#include "SkdTreeUtils.hpp"
#include <cmath>

namespace mimmo{

/*!
 * Default constructor of ControlDeformExtSurface
 */
ControlDeformExtSurface::ControlDeformExtSurface(){
    m_name = "mimmo.ControlDeformExtSurface";
    m_cellBackground = 50;
    m_allowed.insert((FileType::_from_string("STL"))._to_integral());
    m_allowed.insert((FileType::_from_string("STVTU"))._to_integral());
    m_allowed.insert((FileType::_from_string("SQVTU"))._to_integral());
    m_allowed.insert((FileType::_from_string("NAS"))._to_integral());

};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ControlDeformExtSurface::ControlDeformExtSurface(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.ControlDeformExtSurface";
    m_cellBackground = 50;
    m_allowed.insert((FileType::_from_string("STL"))._to_integral());
    m_allowed.insert((FileType::_from_string("STVTU"))._to_integral());
    m_allowed.insert((FileType::_from_string("SQVTU"))._to_integral());
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

/*!Copy constructor of ControlDeformExtSurface.
 */
ControlDeformExtSurface::ControlDeformExtSurface(const ControlDeformExtSurface & other):BaseManipulation(){
    m_allowed = other.m_allowed;
    m_geolist = other.m_geolist;
    m_cellBackground = other.m_cellBackground;
};

/*!
 * Assignement operator of ControlDeformExtSurface. Create an exact copy of the class,
 * except for the deformation field referred to the target geometry.
 */
ControlDeformExtSurface & ControlDeformExtSurface::operator=(const ControlDeformExtSurface & other){
    *(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
    m_allowed = other.m_allowed;
    m_geolist = other.m_geolist;
    m_cellBackground = other.m_cellBackground;
    //deformation and violation field are not copied
    return(*this);
};

/*! It builds the input/output ports of the object
 */
void
ControlDeformExtSurface::buildPorts(){
    bool built = true;

    built = (built && createPortIn<dvecarr3E, ControlDeformExtSurface>(this, &mimmo::ControlDeformExtSurface::setDefField, PortType::M_GDISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT, true));
    built = (built && createPortIn<MimmoObject*, ControlDeformExtSurface>(this, &mimmo::ControlDeformExtSurface::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_, true));

    built = (built && createPortOut<double, ControlDeformExtSurface>(this, &mimmo::ControlDeformExtSurface::getViolation, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<dvector1D, ControlDeformExtSurface>(this, &mimmo::ControlDeformExtSurface::getViolationField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<std::pair<BaseManipulation*, double>, ControlDeformExtSurface>(this, &mimmo::ControlDeformExtSurface::getViolationPair, PortType::M_VIOLATION, mimmo::pin::containerTAG::PAIR, mimmo::pin::dataTAG::PAIRMIMMO_OBJFLOAT_));
    m_arePortsBuilt = built;
};

/*! 
 * Return the value of violation of deformed geometry, after class execution. 
 *  If value is positive or at least zero, constraint violation occurs, penetrating or touching at least in one point the 
 *  constraint sourface. Otherwise, returning negative values means that no violation occurs 
 *  delimited by the plane (where plane normal pointing), true the exact opposite.
 * \return violation value
 */
double 
ControlDeformExtSurface::getViolation(){

    double result = -1.0E+18;
    for(auto & val : m_violationField){
        result = std::fmax(result, val);
    }

    return    result;
};


/*! 
 *  Return the value of violation of deformed geometry, after class execution. 
 *  A BaseManipulation object pointer, representing the sender of geometry to which violation is referred, 
 *  is attached to violation value; 
 *  See getViolation method doc for further details. 
 * \return std::pair structure reporting geometry sender pointer and violation value.
 */
std::pair<BaseManipulation*, double> 
ControlDeformExtSurface::getViolationPair(){

    //get map of Input ports of the class.
    std::map<short int, mimmo::PortIn*> mapPorts = getPortsIn();

    //get class who send geometry here - portID = 99 -> M_GEOM

    std::vector<BaseManipulation*> senders = mapPorts[99]->getLink();

    std::string etiq;
    if(senders.size() == 0){
        return    std::make_pair(this, getViolation());
    }else{
        return    std::make_pair(senders[0], getViolation());
    }

};

/*! 
 * Return the violation distances of each point of deformed geometry, on the geometry itself. The info is available after class execution. 
 *  If value is positive or at least zero, constraint violation occurs, penetrating or touching at least in one point the 
 *  constraint sourface. Otherwise, returning negative values means that no violation occurs 
 *  delimited by the plane (where plane normal pointing), true the exact opposite.
 * \return violation field values
 */
dvector1D 
ControlDeformExtSurface::getViolationField(){
    return(m_violationField);
};


/*!
 * Get Tolerance fixed for each external files listed into the class. Tolerance is needed
 * to consider a deformed body close to target constraints not touching/penetrating them.
 * \param[in] file target external file listed into the class
 * \return     fixed tolerance (if file is not listed return 0.0)
 */
double
ControlDeformExtSurface::getToleranceWithinViolation(std::string file){
    if(m_geolist.count(file) > 0)    return m_geolist[file].first;
    else                            return 0.0;
}

/*!
 * Get number of cells used for determining background grid spacing.
 * See ControlDeformExtSurface::setBackgroundDetails() method documentation;
 * \return     number of cells (default 20)
 */
int
ControlDeformExtSurface::getBackgroundDetails(){
    return m_cellBackground;
}

/*!
 * Set the deformative field associated to each point of the target geometry. 
 * Field resize occurs in execution, if point dimension between field and geoemetry does not match.
 * \param[in]    field of deformation
 */
void
ControlDeformExtSurface::setDefField(dvecarr3E field){
    m_defField.clear();
    m_violationField.clear();
    m_defField = field;
    m_violationField.resize(field.size(),-1.E+18);
};

/*!
 * Set link to target geometry for your selection. Reimplementation of 
 * GenericSelection::setGeometry();
 * \param[in] target pointer to target geometry
 */
void
ControlDeformExtSurface::setGeometry( MimmoObject * target){
    if(target->isEmpty())    return;

    if(target->getType() != 1){
        (*m_log)<<"warning: ControlDeformExtSurface cannot support current geometry. It works only w/ 3D surface."<<std::endl;
        return;
    }
    m_geometry = target;

};


/*!
 * Set the number of Cells NC necessary to determine the spacing of a background grids.
 * A background axis-aligned cartesian volume grid wrapping each constraint + deformed body is created to evaluate 
 * distances of each deformed point from constraint surface. Spacing of such grid is
 * calculated as dh = D/NC where dh is the spacing, D is the diagonal of grid bounding-box.  
 * \param[in] nCell number of cells for cartesian background grid
 */
void
ControlDeformExtSurface::setBackgroundDetails(int nCell){
    m_cellBackground = std::fmax(2, nCell);
}


/*!
 * Return the actual list of external geometry files selected as constraint to check your deformation.
 * Fixed tolerance for each file are reported as first argument of the map, integer identifying format 
 * type of the file (as FileType enum) as second argument of the map.
 * \return list of external geometries files with their tolerance and type
 */
const std::unordered_map<std::string, std::pair<double, int> >
& ControlDeformExtSurface::getFiles() const{
    return    m_geolist;
}

/*!
 * Set a list of external geometry files w/ their relative offset tolerance as constraint 
 * to check your deformation violation. The int format of the file (see FileType enum) must be 
 * present as second argument
 * \param[in] files list external geometries to be read
 */
void
ControlDeformExtSurface::setFiles(std::unordered_map<std::string, std::pair<double, int> >  files){
    for(auto && val : files){
        addFile(val.first, val.second.first, val.second.second);
    }
};

/*!
 * Add a new file with its offset tolerance to a list of external constraint geometry files 
 *\param[in] file of external geometry to be read
 *\param[in] tol  offset tolerance, negative or positive, to consider a deformed body close to target constraints already not touching or penetrating them
 *\param[in] format  type of file as integer (see Filetype enum). only nas, stl, stvtu and sqvtu are supported. 
 */
void
ControlDeformExtSurface::addFile(std::string file, double tol, int format){
    if(m_allowed.count(format)>0){
        m_geolist[file] = std::make_pair(tol, format);
    }
};

/*!
 * Add a new file with its offset tolerance to a list of external constraint geometry files 
 *\param[in] file of external geometry to be read
 *\param[in] tol  offset tolerance, negative or positive, to consider a deformed body close to target constraints already not touching or penetrating them
 *\param[in] format  type of file as  Filetype enum. only nas, stl, stvtu and sqvtu are supported. 
 */
void     ControlDeformExtSurface::addFile(std::string file, double tol, FileType format){
    addFile(file,tol,format._to_integral());
};

/*!
 * Remove an existent file to a list of external geometry files. If not in the list, do nothing 
 * \param[in] file to be removed from the list 
 */
void
ControlDeformExtSurface::removeFile(std::string file){
    if(m_geolist.count(file) >0)    m_geolist.erase(file);
};

/*!
 * Empty your list of external constraint geometry files 
 */
void
ControlDeformExtSurface::removeFiles(){
    m_geolist.clear();
};

/*!
 * Clear your class
 */
void
ControlDeformExtSurface::clear(){
    removeFiles();
    m_defField.clear();
    m_violationField.clear();
    BaseManipulation::clear();
    m_cellBackground = 50;
};

/*!Execution command. Calculate violation value and store it in the class member m_violation
 */
void
ControlDeformExtSurface::execute(){

    MimmoObject * geo = getGeometry();
    if(geo->isEmpty()) return;

    int nDFS = m_defField.size();
    m_defField.resize(geo->getNVertex(),darray3E{{0.0,0.0,0.0}});
    m_violationField.resize(nDFS,-1.E+18);

    dvector1D violationField;
    violationField.resize(nDFS);

    dvecarr3E points = geo->getVertexCoords();
    dvecarr3E pointsOR = points;

    //evaluate deformable geometry barycenter************************
    darray3E geoBary = {{0.0,0.0,0.0}};
    double distBary = 0.0;

    for(auto & p: points)    geoBary += p;

    geoBary /= (double)geo->getNVertex();

    for(auto & p: points){
        distBary= std::fmax(distBary,norm2(p-geoBary));
    }
    //***************************************************************

    //adding deformation to points **********************************
    points+= m_defField;
    //***************************************************************

    //read external surfaces*****************************************
    std::vector<std::unique_ptr<MimmoGeometry> > extgeo;
    dvector1D tols;
    readGeometries(extgeo, tols);
    //***************************************************************

    if(extgeo.size() < 1)    return;

    //evaluating bounding box of deformation*************************
    darray3E bbMinDef = {{0,0,0}}, bbMaxDef={{0,0,0}};
    darray3E bbMin, bbMax;
    {
        int count2 = 0;
        for(auto & p : points){
            if(count2 == 0)    {
                bbMinDef = p;
                bbMaxDef = p;
            }else{
                for(int i=0; i<3; ++i ){
                    bbMinDef[i] = std::fmin(bbMinDef[i], p[i]);
                    bbMaxDef[i] = std::fmax(bbMaxDef[i], p[i]);
                }
            }
            count2++;
        }
        count2=0;
        for(auto & p : pointsOR){
            for(int i=0; i<3; ++i ){
                bbMinDef[i] = std::fmin(bbMinDef[i], p[i]);
                bbMaxDef[i] = std::fmax(bbMaxDef[i], p[i]);
            }
            count2++;
        }
    }
    //***************************************************************
    int counterExtGeo =0;
    // start examining all external geometries***********************
    for(auto &gg : extgeo){

        //check constraints properties ******************************
        MimmoObject * local = gg->getGeometry();
        if(!(local->isBvTreeBuilt()))    local->buildBvTree();
        bool checkOpen = local->isClosedLoop();

        double dist;
        double radius, radius_old;
        long id;
        darray3E normal;
        int count;
        //************************************************************

        //add local constraints to the bounding box compute **********
        bbMin = bbMinDef;
        bbMax = bbMaxDef;

        for(auto & p : local->getVertexCoords()){
            for(int i=0; i<3; ++i ){
                bbMin[i] = std::fmin(bbMin[i], p[i]);
                bbMax[i] = std::fmax(bbMax[i], p[i]);
            }
        }
        //*************************************************************

        //Evaluate background grid dimension for Level set ************
        //calculate dh as a tot part of bb diagonal    *******************

        darray3E span = bbMax - bbMin;
        bbMin += -1.0*(0.05*span);
        span *= 1.1;
        double dh = norm2(span)/(double)m_cellBackground;
        iarray3E dim;
        for(int i=0; i<3; ++i){
            dim[i] = (int)(span[i]/dh + 0.5);
            dim[i] = std::max(2, dim[i]);
        }
        //*************************************************************

        //check if its more convenient evaluate distances with a direct BvTree
        //evaluation of distance on target deformation points or using a background grid +
        //interpolation.
        double nReq = double(dim[0]+1)*double(dim[1]+1)*double(dim[2]+1);
        double nAva = 0.8*nDFS;

        if( nReq > nAva){

            //going to use direct evaluation.

            radius = distBary;
            radius_old = radius;
            dvector1D refsigns(nDFS, 1.0);

            //get the actual sign of distance of the undeformed cloud w.r.t constraints
            if(checkOpen){
                count = 0;
                for(auto &p : pointsOR){
                    dist = evaluateSignedDistance(p, local, id, normal, radius);
                    if(radius < 1.0E-12)    radius = radius_old;
                    else                     radius_old = radius;
                    if(dist < 0.0)    refsigns[count] = -1.0;
                    count++;
                }
            }else{
                dist = evaluateSignedDistance(geoBary, local, id, normal, radius);

                if(dist < 0.0) refsigns *= -1.0;
                if(radius <1.0E-12)    radius = radius_old;
            }

            //evaluate distance of the deformation cloud w.r.t. constraints
            count = 0;
            for(auto &p : points){
                dist = evaluateSignedDistance(p, local, id, normal, radius);
                if(radius < 1.0E-12)    radius = radius_old;
                else                     radius_old = radius;
                violationField[count] = -1.0*refsigns[count]*dist;
                count++;
            }

        }else{

            //going to use background grid to evaluate distances

            //instantiate a VolCartesian;
            bitpit::VolCartesian * mesh = new bitpit::VolCartesian(0,3,bbMin,span, dim);
            mesh->update();

            //calculate distance on point of background grid.
            dvecarr3E background_points(mesh->getVertexCount());
            int count = 0;
            for(auto & v: mesh->getVertices()){
                background_points[count] = v.getCoords();
                count++;
            }
            dvector1D background_LS(background_points.size());

            radius = 0.5*norm2(span);
            radius_old = radius;
            count = 0;
            for(auto &p : background_points){
                dist = evaluateSignedDistance(p, local, id, normal, radius);
                if(radius < 1.0E-12)    radius = radius_old;
                else                     radius_old = radius;
                background_LS[count] = dist;
                count++;
            }

            //interpolate on background to obtain distances

            //evaluate sign of distances of the undeformed cloud w.r.t to constraints.
            dvector1D refsigns(nDFS, 1.0);
            if(checkOpen){
                count=  0;
                for(auto & p : pointsOR){
                    ivector1D stencil;
                    dvector1D weights;
                    int res = mesh->linearVertexInterpolation(p, stencil,weights);

                    dist = 0.0;
                    for(int j=0; j<res; ++j){
                        dist += background_LS[stencil[j]]*weights[j];
                    }

                    if(dist < 0.0)    refsigns[count] = -1.0;
                    ++count;
                }
            }else{
                dist = evaluateSignedDistance(geoBary, local, id, normal, radius);
                if(dist< 0)    refsigns *= -1.0;
                if(radius < 1.0E-12)    radius = radius_old;
            }

            //evaluate violation field
            count = 0;
            for(auto & p : points){
                ivector1D stencil;
                dvector1D weights;

                int res = mesh->linearVertexInterpolation(p, stencil,weights);

                violationField[count] = 0.0;
                for(int j=0; j<res; ++j){
                    violationField[count] += background_LS[stencil[j]]*weights[j];
                }

                violationField[count] *= -1.0*refsigns[count];
                ++count;

            }
            delete mesh;
            mesh = NULL;
        }//end if else

        for(int i=0; i< nDFS; ++i){
            m_violationField[i] = std::fmax(m_violationField[i], (violationField[i] + tols[counterExtGeo]));
        }
        ++counterExtGeo;
    }
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
    
    std::unordered_map<std::string, std::pair<double, int> > mapp;
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

            if(subfile.second->hasOption("tolerance"))    {
                tolstring = subfile.second->get("tolerance");
                tolstring = bitpit::utils::string::trim(tolstring);
                if(!tolstring.empty()){
                    std::stringstream ss(tolstring);
                    ss >> value;
                }
            }

            if(!path.empty() && !tag.empty()){
                mapp[path] = std::make_pair(value,(int) FileType::_from_string(tag.c_str()));
            }
        }

        setFiles(mapp);

    }

    if(slotXML.hasOption("BGDetails")){
        std::string input = slotXML.get("BGDetails");
        input = bitpit::utils::string::trim(input);
        int value = 50;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setBackgroundDetails(value);
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
    for(auto & file : m_geolist){
        std::string name = "file"+std::to_string(counter);
        bitpit::Config::Section & local = filesXML.addSection(name);
        local.set("fullpath", file.first);
        std::string typetag = (FileType::_from_integral(file.second.second))._to_string();
        local.set("tag", typetag);
        std::string tolerance;
        {
            std::stringstream ss;
            ss<<std::scientific<<file.second.first;
            tolerance = ss.str();
        }
        local.set("tolerance", tolerance);

        ++counter;
    }

    slotXML.set("BGDetails", std::to_string(m_cellBackground));

};

/*!
 * Read all external geoemetries from files (whose name is stored in m_geolist) and return it 
 * in a list of unique pointers pointing to MimmoGeometry objects.
 * \param[in,out] extGeo list of read external constraint geoemetries.
 * \param[in,out] tols   tolerance for each effective geometry read
 */
void
ControlDeformExtSurface::readGeometries(std::vector<std::unique_ptr<MimmoGeometry> > & extGeo, std::vector<double>& tols){

    extGeo.resize(m_geolist.size());
    tols.resize(m_geolist.size());

    int counter = 0;
    for(auto & geoinfo : m_geolist){

        svector1D info = extractInfo(geoinfo.first);
        std::unique_ptr<MimmoGeometry> geo (new MimmoGeometry());
        geo->setIOMode(IOMode::READ);
        geo->setDir(info[0]);
        geo->setFilename(info[1]);
        geo->setFileType(geoinfo.second.second);
        geo->setBuildBvTree(true);
        geo->execute();

        if(geo->getGeometry()->getNVertex() == 0 || geo->getGeometry()->getNCells() == 0 || !geo->getGeometry()->isBvTreeSupported()){
            (*m_log)<<"warning: failed to read geometry in ControlDeformExtSurface::readGeometries. Skipping file..."<<std::endl;
        }else{
            if (!geo->getGeometry()->areAdjacenciesBuilt()) geo->getGeometry()->getPatch()->buildAdjacencies();
            extGeo[counter] = std::move(geo);
            tols[counter] = geoinfo.second.first;
            ++counter;
        }
    }

    extGeo.resize(counter);
    tols.resize(counter);
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
 * Evaluate Signed Distance for a point from given BvTree of a open/closed geometry 3D surface. 
 * Return distance from target geometry with sign. Positive distance is returned, 
 * if the point is placed on the same side of simplex support normal. Negative otherwise.
 * Wrapper to bvTreeUtils::signedDistance
 *
 * \param[in] point 3D target point
 * \param[in] geo    target geometry w/ bvTree in it
 * \param[in,out] id returns id of support geometry simplex from which distance is evaluated 
 * \param[in,out] normal returns normal of support geometry simplex from which distance is evaluated 
 * \param[in,out] initRadius guess initial distance.
 * 
 * \return signed distance from target surface.  
 */
double
ControlDeformExtSurface::evaluateSignedDistance(darray3E &point, mimmo::MimmoObject * geo, long & id, darray3E & normal, double &initRadius){

    double dist = 1.0E+18;
    double rate = 0.02;
    int kmax = 1000;
    int kiter = 0;
    bool flag = true;

    if(!geo->isBvTreeBuilt())    geo->buildBvTree();

    while(flag && kiter < kmax){
        dist = skdTreeUtils::signedDistance(&point, geo->getBvTree(), id, normal, initRadius);
        flag = (dist == 1.0E+18);
        if(flag)    initRadius *= (1.0+ rate*((double)flag));
        kiter++;
    }

    return dist;
}

/*!
 * Plot optional results in execution, that is the violation distance field on current deformed geometry.
 * Reimeplemented from BaseManipulation::plotOptionalResults;
 */
void
ControlDeformExtSurface::plotOptionalResults(){
    if(getGeometry()->isEmpty())    return;

    dvecarr3E  points = getGeometry()->getVertexCoords();
    m_defField.resize(points.size());
    points+=m_defField;
    ivector2D connectivity = getGeometry()->getCompactConnectivity();

    bitpit::VTKElementType  elDM = bitpit::VTKElementType::TRIANGLE;

    std::string name = m_name +std::to_string(getClassCounter())+ "_ViolationField";
    bitpit::VTKUnstructuredGrid output(m_outputPlot, name, elDM);
    output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points) ;
    output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity) ;
    output.setDimensions(connectivity.size(), points.size());


    m_violationField.resize(points.size(), -1.e+18);
    std::string sdfstr = "Violation Distance Field";
    output.addData(sdfstr, bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, m_violationField);
    output.write();
}

}
