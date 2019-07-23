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
#include "TranslationGeometry.hpp"

namespace mimmo{


/*!
 * Default constructor of TranslationGeometry
 */
TranslationGeometry::TranslationGeometry(darray3E direction){
    m_direction = direction;
    m_alpha = 0.0;
    m_name = "mimmo.TranslationGeometry";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
TranslationGeometry::TranslationGeometry(const bitpit::Config::Section & rootXML){

    m_direction.fill(0.0);
    m_alpha = 0.0;
    m_name = "mimmo.TranslationGeometry";

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.TranslationGeometry"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!Default destructor of TranslationGeometry
 */
TranslationGeometry::~TranslationGeometry(){};

/*!Copy constructor of TranslationGeometry.No result geometry displacements are copied.
 */
TranslationGeometry::TranslationGeometry(const TranslationGeometry & other):BaseManipulation(other){
    m_direction = other.m_direction;
    m_alpha = other.m_alpha;
    m_filter = other.m_filter;
};

/*!Assignment operator of TranslationGeometry.No result geometry displacements are copied.
 */
TranslationGeometry & TranslationGeometry::operator=(TranslationGeometry other){
    swap(other);
    return *this;
};

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void TranslationGeometry::swap(TranslationGeometry & x) noexcept
{
    std::swap(m_direction, x.m_direction);
    std::swap(m_alpha, x.m_alpha);
//     std::swap(m_filter, x.m_filter);
//     std::swap(m_displ, x.m_displ);
    m_filter.swap(x.m_filter);
    m_displ.swap(x.m_displ);
    BaseManipulation::swap(x);
}

/*! It builds the input/output ports of the object
 */
void
TranslationGeometry::buildPorts(){
    bool built = true;
    built = (built && createPortIn<darray3E, TranslationGeometry>(&m_direction, M_AXIS));
    built = (built && createPortIn<double, TranslationGeometry>(&m_alpha, M_VALUED));
    built = (built && createPortIn<dmpvector1D, TranslationGeometry>(this, &mimmo::TranslationGeometry::setFilter, M_FILTER));
    built = (built && createPortIn<MimmoObject*, TranslationGeometry>(&m_geometry, M_GEOM, true));
    built = (built && createPortOut<dmpvecarr3E, TranslationGeometry>(this, &mimmo::TranslationGeometry::getDisplacements, M_GDISPLS));
    built = (built && createPortOut<MimmoObject*, TranslationGeometry>(this, &BaseManipulation::getGeometry, M_GEOM));
    m_arePortsBuilt = built;
};

/*!It sets the direction of the translation axis.
 * \param[in] direction Direction of translation axis.
 */
void
TranslationGeometry::setDirection(darray3E direction){
    m_direction = direction;
    double L = norm2(m_direction);
    for (int i=0; i<3; i++)
        m_direction[i] /= L;
}

/*!It sets the value of the translation.
 * \param[in] alpha Value of translation in length unity.
 */
void
TranslationGeometry::setTranslation(double alpha){
    m_alpha = alpha;
}

/*!It sets the filter field to modulate the displacements of the vertices
 * of the target geometry.
 * \param[in] filter Filter field defined on geometry vertices.
 */
void
TranslationGeometry::setFilter(dmpvector1D filter){
    m_filter = filter;
}

/*!
 * Return actual computed displacements field (if any) for the geometry linked.
 * \return  deformation field
 */
dmpvecarr3E
TranslationGeometry::getDisplacements(){
    return m_displ;
};

/*!Execution command. It perform the translation by computing the displacements
 * of the points of the geometry. It applies a filter field eventually set as input.
 */
void
TranslationGeometry::execute(){

    if(getGeometry() == NULL){
        (*m_log)<<m_name + " : NULL pointer to linked geometry found"<<std::endl;
        throw std::runtime_error(m_name + "NULL pointer to linked geometry found");
    }


    if(getGeometry()->isEmpty()){
        (*m_log)<<m_name + " : empty linked geometry found"<<std::endl;
    }

    checkFilter();
    m_displ.clear();
    m_displ.clear();
    m_displ.setDataLocation(mimmo::MPVLocation::POINT);
    m_displ.reserve(getGeometry()->getNVertices());
    m_displ.setGeometry(getGeometry());


    long ID;
    darray3E value;
    for (const auto & vertex : m_geometry->getVertices()){
        ID = vertex.getId();
        value = m_alpha*m_direction*m_filter[ID];
        m_displ.insert(ID, value);
    }
};

/*!
 * Directly apply deformation field to target geometry.
 */
void
TranslationGeometry::apply(){

    if (getGeometry() == NULL) return;
    darray3E vertexcoords;
    long int ID;
    for (const auto & vertex : m_geometry->getVertices()){
        vertexcoords = vertex.getCoords();
        ID = vertex.getId();
        if(m_displ.exists(ID))  vertexcoords += m_displ[ID];
        getGeometry()->modifyVertex(vertexcoords, ID);
    }

#if MIMMO_ENABLE_MPI
	getGeometry()->updatePointGhostExchangeInfo();
#endif

}

/*!
 * Check if the filter is related to the target geometry.
 * If not create a unitary filter field.
 */
void
TranslationGeometry::checkFilter(){
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
TranslationGeometry::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::absorbSectionXML(slotXML, name);

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

    if(slotXML.hasOption("Translation")){
        std::string input = slotXML.get("Translation");
        input = bitpit::utils::string::trim(input);
        double temp = 0.0;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        setTranslation(temp);
    }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
TranslationGeometry::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    BaseManipulation::flushSectionXML(slotXML, name);
    {
        std::stringstream ss;
        ss<<std::scientific<<m_direction[0]<<'\t'<<m_direction[1]<<'\t'<<m_direction[2];
        slotXML.set("Direction", ss.str());
    }

    slotXML.set("Translation", std::to_string(m_alpha));
};

}
