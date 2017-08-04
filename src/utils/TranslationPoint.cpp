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
#include "TranslationPoint.hpp"

namespace mimmo{

/*!
 * Default constructor of TranslationPoint
 */
TranslationPoint::TranslationPoint(darray3E direction){
    m_direction = direction;
    m_alpha = 0.0;
    m_name = "mimmo.TranslationPoint";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
TranslationPoint::TranslationPoint(const bitpit::Config::Section & rootXML){

    m_direction.fill(0.0);
    m_alpha = 0.0;
    m_name = "mimmo.TranslationPoint";

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.TranslationPoint"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!Default destructor of TranslationPoint
 */
TranslationPoint::~TranslationPoint(){};

/*!Copy constructor of TranslationPoint.
 */
TranslationPoint::TranslationPoint(const TranslationPoint & other):BaseManipulation(other){
    m_origin = other.m_origin;
    m_alpha = other.m_alpha;
    m_direction = other.m_direction;
};

/*!
 * Swap Function
 * \param[in] x object to be swapped
 */
void TranslationPoint::swap(TranslationPoint & x) noexcept
{
    std::swap(m_origin, x.m_origin);
    std::swap(m_alpha, x.m_alpha);
    std::swap(m_direction, x.m_direction);
    BaseManipulation::swap(x);
}

/*! It builds the input/output ports of the object
 */
void
TranslationPoint::buildPorts(){
    
    PortManager::instance().addPort(M_POINT, MC_ARRAY3, MD_FLOAT);
    PortManager::instance().addPort(M_AXIS, MC_ARRAY3, MD_FLOAT);
    PortManager::instance().addPort(M_VALUED, MC_SCALAR, MD_FLOAT);
    
    bool built = true;
    built = (built && createPortIn<darray3E, TranslationPoint>(&m_origin, M_POINT));
    built = (built && createPortIn<darray3E, TranslationPoint>(&m_direction, M_AXIS));
    built = (built && createPortIn<double, TranslationPoint>(&m_alpha, M_VALUED));
    built = (built && createPortOut<darray3E, TranslationPoint>(this, &mimmo::TranslationPoint::getOrigin, M_POINT));
    built = (built && createPortOut<darray3E, TranslationPoint>(this, &mimmo::TranslationPoint::getDirection, M_AXIS));
    built = (built && createPortOut<double, TranslationPoint>(this, &mimmo::TranslationPoint::getTranslation, M_VALUED));
    m_arePortsBuilt = built;
};

/*!It gets the direction of the translation.
 * \return Direction of translation.
 */
darray3E
TranslationPoint::getDirection(){
    return(m_direction);
}

/*!It gets the value of the translation.
 * \return Value of translation.
 */
double
TranslationPoint::getTranslation(){
    return(m_alpha);
}

/*!It gets the original position of the point to be translated (before the execution)
 * or the position of the translated point (after the execution of the object).
 * \return Position of the point.
 */
darray3E
TranslationPoint::getOrigin(){
    return(m_origin);
}

/*!It sets the direction of the translation.
 * \param[in] direction Direction of translation.
 */
void
TranslationPoint::setDirection(darray3E direction){
    m_direction = direction;
}

/*!It sets the value of the translation.
 * \param[in] alpha Value of translation.
 */
void
TranslationPoint::setTranslation(double alpha){
    m_alpha = alpha;
}

/*!It sets the original coordinates of the point to be translated.
 * \param[in] origin Position of the point.
 */
void
TranslationPoint::setOrigin(darray3E origin){
    m_origin = origin;
}

/*!Execution command. It modifies the coordinates of the origin
 * with the translation conditions.
 * The result of the translation is stored in member result of base class
 * and in the member m_origin.
 * After exec() the original point coordinates will be permanently modified.
 */
void
TranslationPoint::execute(){
    for (int i=0; i<3; i++){
        m_origin[i] += m_alpha * m_direction[i];
    }
};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
TranslationPoint::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

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
TranslationPoint::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

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

    slotXML.set("Translation", std::to_string(m_alpha));

};


}

