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
#include "RotationGeometry.hpp"

namespace mimmo{


/*!
 * Default constructor of RotationGeometry
 * \param[in] origin origin of rotation axis
 * \param[in] direction direction of rotation axis
 */
RotationGeometry::RotationGeometry(darray3E origin, darray3E direction){
    m_origin = origin;
    m_direction = direction;
    m_alpha=  0.0;
    m_name = "mimmo.RotationGeometry";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
RotationGeometry::RotationGeometry(const bitpit::Config::Section & rootXML){

    m_origin.fill(0.0);
    m_direction.fill(0.0);
    m_alpha = 0.0;
    m_name = "mimmo.RotationGeometry";

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.RotationGeometry"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!Default destructor of RotationGeometry
 */
RotationGeometry::~RotationGeometry(){};

/*!
 * Copy constructor of RotationGeometry. No result geometry displacements are copied.
 */
RotationGeometry::RotationGeometry(const RotationGeometry & other):BaseManipulation(other){
    m_origin = other.m_origin;
    m_direction = other.m_direction;
    m_alpha = other.m_alpha;
    m_filter = other.m_filter;
};

/*!
 * Assignment operator. No result geometry displacements are copied.
 */
RotationGeometry & RotationGeometry::operator=(RotationGeometry other){
    swap(other);
    return *this;
}

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void RotationGeometry::swap(RotationGeometry & x) noexcept
{
    std::swap(m_origin, x.m_origin);
    std::swap(m_direction, x.m_direction);
    std::swap(m_alpha, x.m_alpha);
    m_filter.swap(x.m_filter);
    m_displ.swap(x.m_displ);
    BaseManipulation::swap(x);
}

/*! It builds the input/output ports of the object
 */
void
RotationGeometry::buildPorts(){
    bool built = true;
    built = (built && createPortIn<darray3E, RotationGeometry>(&m_origin, M_POINT));
    built = (built && createPortIn<darray3E, RotationGeometry>(&m_direction, M_AXIS));
    built = (built && createPortIn<double, RotationGeometry>(&m_alpha, M_VALUED));
    built = (built && createPortIn<dmpvector1D*, RotationGeometry>(this, &mimmo::RotationGeometry::setFilter, M_FILTER));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, RotationGeometry>(&m_geometry, M_GEOM, true));
    built = (built && createPortOut<dmpvecarr3E*, RotationGeometry>(this, &mimmo::RotationGeometry::getDisplacements, M_GDISPLS));
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, RotationGeometry>(this, &BaseManipulation::getGeometry, M_GEOM));
    m_arePortsBuilt = built;
};

/*!It sets the origin and direction of the rotation axis.
 * \param[in] origin origin of rotation axis
 * \param[in] direction direction of rotation axis
 */
void
RotationGeometry::setAxis(darray3E origin, darray3E direction){
    m_origin = origin;
    m_direction = direction;
}

/*!It sets the origin of the rotation axis.
 * \param[in] origin origin of rotation axis
 */
void
RotationGeometry::setOrigin(darray3E origin){
    m_origin = origin;
}

/*!It sets the direction of the rotation axis.
 * \param[in] direction direction of rotation axis
 */
void
RotationGeometry::setDirection(darray3E direction){
    m_direction = direction;
    double L = norm2(m_direction);
    for (int i=0; i<3; i++)
        m_direction[i] /= L;
}

/*!It sets the value of the rotation.
 * \param[in] alpha value of rotation [rad]
 */
void
RotationGeometry::setRotation(double alpha){
    m_alpha = alpha;
}

/*!It sets the filter field to modulate the displacements of the vertices
 * of the target geometry.
 * \param[in] filter filter field defined on geometry vertices.
 */
void
RotationGeometry::setFilter(dmpvector1D *filter){
    if(!filter) return;
    m_filter = *filter;
}

/*!
 * Return actual computed displacements field (if any) for the geometry linked.
 * \return  deformation field
 */
dmpvecarr3E*
RotationGeometry::getDisplacements(){
    return &m_displ;
};

/*!Execution command. It saves in "rot"-terms the modified axes and origin, by the
 * rotation conditions. This terms can be recovered and passed by a pin to a child object
 * by the related get-methods.
 */
void
RotationGeometry::execute(){

    if(getGeometry() == nullptr){
        (*m_log)<<m_name + " : nullptr pointer to linked geometry found"<<std::endl;
        throw std::runtime_error(m_name + "nullptr pointer to linked geometry found");
    }

    if(getGeometry()->isEmpty()){
        (*m_log)<<m_name + " : empty linked geometry found"<<std::endl;
    }

    checkFilter();

    m_displ.clear();
    m_displ.setDataLocation(mimmo::MPVLocation::POINT);
    m_displ.reserve(getGeometry()->getNVertices());
    m_displ.setGeometry(getGeometry());


    //compute coefficients and constant vectors of rodriguez formula
    double a   = std::cos(m_alpha);
    darray3E b =  (1.0 - std::cos(m_alpha)) * m_direction;
    double c   = std::sin(m_alpha);

    darray3E point, rotated;
    long ID;
    darray3E value;
    for (const auto & vertex : m_geometry->getVertices()){
        point = vertex.getCoords();
        ID = vertex.getId();

        point -= m_origin;
        //rodrigues formula
        rotated = a * point +
                b * dotProduct(m_direction, point) +
                c * crossProduct(m_direction, point);

        rotated += m_origin;
        point += m_origin;

        value = (rotated-point)*m_filter[ID];
        m_displ.insert(ID, value);

    }
};

/*!
 * Directly apply deformation field to target geometry.
 */
void
RotationGeometry::apply(){
    _apply(m_displ);
}

/*!
 * Check if the filter is related to the target geometry.
 * If not create a unitary filter field.
 */
void
RotationGeometry::checkFilter(){
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
RotationGeometry::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::absorbSectionXML(slotXML, name);

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

    if(slotXML.hasOption("Direction")){
        std::string input = slotXML.get("Direction");
        input = bitpit::utils::string::trim(input);
        darray3E temp = {{0.0,0.0,0.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            for(auto &val : temp) ss>>val;
        }
        setDirection(temp);
    }

    if(slotXML.hasOption("Rotation")){
        std::string input = slotXML.get("Rotation");
        input = bitpit::utils::string::trim(input);
        double temp = 0.0;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        setRotation(temp);
    }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
RotationGeometry::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    {
        std::stringstream ss;
        ss<<std::scientific<<m_origin[0]<<'\t'<<m_origin[1]<<'\t'<<m_origin[2];
        slotXML.set("Origin", ss.str());
    }

    {
        std::stringstream ss;
        ss<<std::scientific<<m_direction[0]<<'\t'<<m_direction[1]<<'\t'<<m_direction[2];
        slotXML.set("Direction", ss.str());
    }

    slotXML.set("Rotation", std::to_string(m_alpha));

};

}
