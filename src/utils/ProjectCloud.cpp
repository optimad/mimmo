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
#include "ProjectCloud.hpp"
#include "SkdTreeUtils.hpp"

namespace mimmo{

/*!
 * Default constructor of ProjectCloud
 */
ProjectCloud::ProjectCloud(){
    m_name = "mimmo.ProjectCloud";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ProjectCloud::ProjectCloud(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.ProjectCloud";

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.ProjectCloud"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}


/*!Default destructor of ProjectCloud
 */
ProjectCloud::~ProjectCloud(){};

/*!
 * Copy constructor of ProjectCloud. Result Projected points are not copied
 */
ProjectCloud::ProjectCloud(const ProjectCloud & other):BaseManipulation(other){
    m_points = other.m_points;
};


/*!
 * Assignment operator of ProjectCloud. Result Projected points are not copied
 */
ProjectCloud & ProjectCloud::operator=(ProjectCloud other){
    swap(other);
    return *this;
};

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void ProjectCloud::swap(ProjectCloud & x) noexcept
{
    std::swap(m_points, x.m_points);
    std::swap(m_proj, x.m_proj);
    std::swap(m_labels, x.m_labels);
    std::swap(m_internalPC, x.m_internalPC);

    BaseManipulation::swap(x);
}
/*! It builds the input/output ports of the object
 */
void
ProjectCloud::buildPorts(){
    bool built = true;

    built = (built && createPortIn<dvecarr3E, ProjectCloud>(this, &mimmo::ProjectCloud::setCoords, M_COORDS, true,1));
    built = (built && createPortIn<MimmoObject * , ProjectCloud>(this, &mimmo::ProjectCloud::setGeometry, M_GEOM, true));
    built = (built && createPortIn<MimmoObject * , ProjectCloud>(this, &mimmo::ProjectCloud::setCoords, M_GEOM2, true,1));

    built = (built && createPortOut<dvecarr3E, ProjectCloud>(this, &mimmo::ProjectCloud::getProjectedCoords, M_COORDS));
    built = (built && createPortOut<MimmoObject*, ProjectCloud>(this, &mimmo::ProjectCloud::getProjectedCloud, M_GEOM));

    m_arePortsBuilt = built;
};

/*!
    It gets the coordinates of points stored in the object.
    \param[out] labels return also the labels attached to points if structure is not null.
 *  \return coordinates of points stored in the object.
 */
dvecarr3E
ProjectCloud::getOriginalCoords(livector1D * labels){

    if(labels ){
        *labels = m_labels;
    }
    return(m_points);
};

/*!
    It gets the resulting points after class functionality execution
 * \return projected points stored in the object.
 */
dvecarr3E
ProjectCloud::getProjectedCoords(){
    return getProjectedCoords(nullptr);
};

/*!
    It gets the resulting points after class functionality execution
    \param[out] labels return also the labels attached to points if structure is not null.
 * \return projected points stored in the object.
 */
dvecarr3E
ProjectCloud::getProjectedCoords(livector1D * labels){

    if(labels){
        *labels = m_labels;
         labels->resize(m_proj.size());
    }
    return  m_proj;
};

/*!
    It gets the resulting points after class functionality execution
    in the form of a MimmoObject point cloud.
 * \return point cloud of projected point
 */
MimmoObject *
ProjectCloud::getProjectedCloud(){

    m_internalPC = std::unique_ptr<MimmoObject>(new MimmoObject(3)); //new point cloud class
    m_internalPC->getPatch()->reserveVertices(m_proj.size());
    int counter = 0;
    for(darray3E & val : m_proj){
        m_internalPC->addVertex(val, m_labels[counter]);
        ++counter;
    }
    return  m_internalPC.get();
};


/*!
   It sets the coordinates of original points to be processed.
 * \param[in] coords coordinates of points to be used .
 */
void
ProjectCloud::setCoords(dvecarr3E coords){
    m_points = coords;
    //update labels
    m_labels.resize(m_points.size());
    int counter = 0;
    for(long & label : m_labels){
        label = counter;
        ++counter;
    }
};

/*!
   Read the coordinates of original points to be processed from any valid MimmoObject
   and treat them as a point cloud.
 * \param[in] coordsPC generic MimmoObject geoemtry
 */
void
ProjectCloud::setCoords(MimmoObject * coordsPC){

    m_points.clear();
    m_labels.clear();

    if(coordsPC){
        m_points.reserve(coordsPC->getPatch()->getVertexCount());
        m_labels.reserve(m_points.size());
        for(bitpit::Vertex & vert : coordsPC->getVertices()){
            m_points.push_back(vert.getCoords());
            m_labels.push_back(vert.getId());
        }
    }
};

/*!
   Set target reference surface geometry where the point cloud needs to be projected.
   Only surface geometry of MimmoObject type 1 are allowed.
 * \param[in] refSurface reference geometry
 */
void
ProjectCloud::setGeometry(MimmoObject * refSurface){
    if(!refSurface) return;
    if(refSurface->getType() != 1) return;

    BaseManipulation::setGeometry(refSurface);
};

/*!
 * clear contents of the class
 */
void ProjectCloud::clear(){
    BaseManipulation::clear();
    m_points.clear();
    m_proj.clear();
    m_labels.clear();
    m_internalPC = nullptr;
}


/*!Execution command.
 * Project list of points on the target geometry.
 */
void
ProjectCloud::execute(){

    if(getGeometry() == NULL){
        (*m_log)<<m_name + " : NULL pointer to linked geometry found"<<std::endl;
        throw std::runtime_error(m_name + "NULL pointer to linked geometry found");
    }

    if(getGeometry()->isEmpty()){
        (*m_log)<<m_name + " : empty linked geometry found"<<std::endl;
        return;
    }
    if(getGeometry()->getType() != 1){
        return;
    }

    if(!getGeometry()->isSkdTreeSync())    getGeometry()->buildSkdTree();

    //project points on surface.
    int counter = 0;
    m_proj.resize(m_points.size(), {{0.0,0.0,0.0}});
    for(darray3E &val : m_points){
        m_proj[counter]= skdTreeUtils::projectPoint(&val, getGeometry()->getSkdTree());
        ++counter;
    }
    return;
};

/*!
 * Plot optional result of the class in execution, that is the projected cloud
 * as standard vtk unstructured grid.
 */
void
ProjectCloud::plotOptionalResults(){
    bitpit::VTKFormat codex = bitpit::VTKFormat::APPENDED;

    int size = m_proj.size();
    ivector1D conn(size);
    for(int i=0; i<size; i++){
        conn[i] = i;
    }
    std::string dir = "./";
    std::string file = m_name + "_" + std::to_string(getId());

    livector1D labels = m_labels;
    labels.resize(m_proj.size());

    bitpit::VTKUnstructuredGrid vtk(dir, file, bitpit::VTKElementType::VERTEX);
    vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, m_proj) ;
    vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, conn) ;
    vtk.setDimensions(size, size);
    vtk.addData("labels", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, labels);
    vtk.setCodex(codex);

    vtk.write();
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
ProjectCloud::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name ){

    BITPIT_UNUSED(name);
    //start absorbing
    BaseManipulation::absorbSectionXML(slotXML, name);
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
ProjectCloud::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    BaseManipulation::flushSectionXML(slotXML, name);

};

}
