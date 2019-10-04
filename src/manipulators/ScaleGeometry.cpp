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
#include "ScaleGeometry.hpp"

namespace mimmo{


/*!
 * Default constructor of ScaleGeometry
 */
ScaleGeometry::ScaleGeometry(darray3E scaling){
    m_scaling = scaling;
    m_origin.fill(0.0);
    m_meanP = false;
    m_name = "mimmo.ScaleGeometry";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ScaleGeometry::ScaleGeometry(const bitpit::Config::Section & rootXML){

    m_scaling.fill(1.0);
    m_origin.fill(0.0);
    m_meanP = false;
    m_name = "mimmo.ScaleGeometry";

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.ScaleGeometry"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!Default destructor of ScaleGeometry
 */
ScaleGeometry::~ScaleGeometry(){};

/*!Copy constructor of ScaleGeometry. No result geometry displacements are copied.
 */
ScaleGeometry::ScaleGeometry(const ScaleGeometry & other):BaseManipulation(other){
    m_scaling = other.m_scaling;
    m_origin = other.m_origin;
    m_meanP = other.m_meanP;
    m_filter = other.m_filter;
};

/*!
 * Assignment operator. No result geometry displacements are copied.
 */
ScaleGeometry & ScaleGeometry::operator=(ScaleGeometry other){
    swap(other);
    return *this;
}

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void ScaleGeometry::swap(ScaleGeometry & x) noexcept
{
    std::swap(m_scaling, x.m_scaling);
    std::swap(m_origin, x.m_origin);
    std::swap(m_meanP, x.m_meanP);
//     std::swap(m_filter, x.m_filter);
//     std::swap(m_displ, x.m_displ);
    m_filter.swap(x.m_filter);
    m_displ.swap(x.m_displ);
    BaseManipulation::swap(x);
}


/*! It builds the input/output ports of the object
 */
void
ScaleGeometry::buildPorts(){
    bool built = true;
    built = (built && createPortIn<darray3E, ScaleGeometry>(&m_origin, M_POINT));
    built = (built && createPortIn<darray3E, ScaleGeometry>(&m_scaling, M_SPAN));
    built = (built && createPortIn<dmpvector1D*, ScaleGeometry>(this, &mimmo::ScaleGeometry::setFilter, M_FILTER));
    built = (built && createPortIn<MimmoObject*, ScaleGeometry>(&m_geometry, M_GEOM, true));
    built = (built && createPortOut<dmpvecarr3E*, ScaleGeometry>(this, &mimmo::ScaleGeometry::getDisplacements, M_GDISPLS));
    built = (built && createPortOut<MimmoObject*, ScaleGeometry>(this, &BaseManipulation::getGeometry, M_GEOM));
    m_arePortsBuilt = built;
};

/*!It sets the center point.
 * \param[in] origin origin for scaling transform
 */
void
ScaleGeometry::setOrigin(darray3E origin){
    m_origin = origin;
}

/*!It sets if the center point for scaling is the mean point of the geometry.
 * \param[in] meanP origin for scaling transform is the mean point?
 */
void
ScaleGeometry::setMeanPoint(bool meanP){
    m_meanP = meanP;
}

/*!It sets the scaling factors of each axis.
 * \param[in] scaling scaling factor values for x, y and z absolute axis.
 */
void
ScaleGeometry::setScaling(darray3E scaling){
    m_scaling = scaling;
}

/*!It sets the filter field to modulate the displacements of the vertices
 * of the target geometry.
 * \param[in] filter filter field defined on geometry vertices.
 */
void
ScaleGeometry::setFilter(dmpvector1D *filter){
    if(!filter) return;
    m_filter = *filter;
}

/*!
 * Return actual computed displacements field (if any) for the geometry linked.
 * \return  deformation field
 */
dmpvecarr3E*
ScaleGeometry::getDisplacements(){
    return &m_displ;
};

/*!Execution command. It perform the scaling by computing the displacements
 * of the points of the geometry. It applies a filter field eventually set as input.
 */
void
ScaleGeometry::execute(){

    if(getGeometry() == NULL){
        (*m_log)<<m_name + " : NULL pointer to linked geometry found"<<std::endl;
        throw std::runtime_error(m_name + "NULL pointer to linked geometry found");
    }

    if(getGeometry()->isEmpty()){
        (*m_log)<<m_name + " : empty linked geometry found"<<std::endl;
    }

    checkFilter();
    m_displ.clear();
    m_displ.setDataLocation(mimmo::MPVLocation::POINT);
    m_displ.reserve(getGeometry()->getNVertices());
    m_displ.setGeometry(getGeometry());


    long ID;
    darray3E value;
    //computing centroid
    int nV = m_geometry->getNVertices();
    darray3E center = m_origin;
    if (m_meanP){
        center.fill(0.0);
        for (const auto & vertex : m_geometry->getVertices()){
            center += vertex.getCoords();
        }
        center /=double(nV);
    }
    for (const auto & vertex : m_geometry->getVertices()){
        ID = vertex.getId();

        darray3E coords = vertex.getCoords();
        value = (( m_scaling*(coords - center) + center ) - coords) * m_filter[ID] ;
        m_displ.insert(ID, value);
    }
};

/*!
 * Directly apply deformation field to target geometry.
 */
void
ScaleGeometry::apply(){
    _apply(m_displ);
}

/*!
 * Check if the filter is related to the target geometry.
 * If not create a unitary filter field.
 */
void
ScaleGeometry::checkFilter(){
    bool check = m_filter.getDataLocation() == mimmo::MPVLocation::POINT;
    check = check && m_filter.completeMissingData(0.0);
    check = check && m_filter.getGeometry() == getGeometry();

    if (!check){
        m_log->setPriority(bitpit::log::Verbosity::DEBUG);
        (*m_log)<<"Not valid filter found in "<<m_name<<". Proceeding with default unitary field"<<std::endl;
        m_log->setPriority(bitpit::log::Verbosity::NORMAL);

        m_filter.clear();
        m_filter.setGeometry(m_geometry);
        m_filter.setDataLocation(mimmo::MPVLocation::POINT);
        m_filter.reserve(getGeometry()->getNVertices());
        for (const auto & vertex : getGeometry()->getVertices()){
            m_filter.insert(vertex.getId(), 1.0);
        }
    }
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
ScaleGeometry::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("MeanPoint")){
        std::string input = slotXML.get("MeanPoint");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setMeanPoint(value);
    };


    if(slotXML.hasOption("Origin")){
        std::string input = slotXML.get("Origin");
        input = bitpit::utils::string::trim(input);
        darray3E temp = {{0.0,0.0,0.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            for(auto &val : temp) ss>>val;
        }
        setOrigin(temp);
    }

    if(slotXML.hasOption("Scaling")){
        std::string input = slotXML.get("Scaling");
        input = bitpit::utils::string::trim(input);
        darray3E temp = {{0.0,0.0,0.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            for(auto &val : temp) ss>>val;
        }
        setScaling(temp);
    }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
ScaleGeometry::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);
    {
        std::stringstream ss;
        ss<<std::scientific<<m_origin[0]<<'\t'<<m_origin[1]<<'\t'<<m_origin[2];
        slotXML.set("Origin", ss.str());
    }

    {
        bool value = m_meanP;
        std::string towrite = std::to_string(value);
        slotXML.set("MeanPoint", towrite);
    }

    {
        std::stringstream ss;
        ss<<std::scientific<<m_scaling[0]<<'\t'<<m_scaling[1]<<'\t'<<m_scaling[2];
        slotXML.set("Scaling", ss.str());
    }


};

}
