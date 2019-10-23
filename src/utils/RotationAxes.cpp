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
#include "RotationAxes.hpp"

namespace mimmo{


/*!
 * Default constructor of RotationAxes
 */
RotationAxes::RotationAxes(darray3E origin, darray3E direction){
    setOrigin(origin);
    setDirection(direction);

    for(auto & val: m_axes)val.fill(0.0);
    m_axes[0][0] = 1.0;
    m_axes[1][1] = 1.0;
    m_axes[1][1] = 1.0;
    m_axes_origin = {{0,0,0}};

    for(auto & val: m_rotax)val.fill(0.0);
    m_rotax[0][0] = 1.0;
    m_rotax[1][1] = 1.0;
    m_rotax[1][1] = 1.0;
    m_rotax_origin = {{0,0,0}};

    m_alpha=0.0;

    m_name = "mimmo.RotationAxes";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
RotationAxes::RotationAxes(const bitpit::Config::Section & rootXML){

    m_origin.fill(0.0);
    m_direction.fill(0.0);

    for(auto & val: m_axes)
        val.fill(0.0);
    m_axes[0][0] = 1.0;
    m_axes[1][1] = 1.0;
    m_axes[1][1] = 1.0;
    m_axes_origin = {{0,0,0}};

    for(auto & val: m_rotax)
        val.fill(0.0);
    m_rotax[0][0] = 1.0;
    m_rotax[1][1] = 1.0;
    m_rotax[1][1] = 1.0;
    m_rotax_origin = {{0,0,0}};

    m_alpha=0.0;

    m_name = "mimmo.RotationAxes";

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.RotationAxes"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!Default destructor of RotationAxes
 */
RotationAxes::~RotationAxes(){};

/*!Copy constructor of RotationAxes.
 */
RotationAxes::RotationAxes(const RotationAxes & other):BaseManipulation(other){
    m_axes = other.m_axes;
    m_axes_origin = other.m_axes_origin;
    m_alpha=other.m_alpha;
    m_origin = other.m_origin;
    m_direction = other.m_direction;
    for(auto & val: m_rotax)
        val.fill(0.0);
    m_rotax[0][0] = 1.0;
    m_rotax[1][1] = 1.0;
    m_rotax[1][1] = 1.0;
    m_rotax_origin = {{0,0,0}};
};

/*!Assignement operator of RotationAxes.
 */
RotationAxes & RotationAxes::operator=(RotationAxes other){
    swap(other);
    return *this;
};


/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void RotationAxes::swap(RotationAxes & x) noexcept
{
    std::swap(m_axes, x.m_axes);
    std::swap(m_axes_origin, x.m_axes_origin);
    std::swap(m_alpha, x.m_alpha);
    std::swap(m_origin, x.m_origin);
    std::swap(m_direction, x.m_direction);
    std::swap(m_rotax, x.m_rotax);
    std::swap(m_rotax_origin, x.m_rotax_origin);
    BaseManipulation::swap(x);
}
/*! It builds the input/output ports of the object
 */
void
RotationAxes::buildPorts(){
    bool built = true;
    built = (built && createPortIn<darray3E, RotationAxes>(this, &mimmo::RotationAxes::setOrigin, M_POINT));
    built = (built && createPortIn<darray3E, RotationAxes>(this, &mimmo::RotationAxes::setDirection, M_AXIS));
    built = (built && createPortIn<double, RotationAxes>(this, &mimmo::RotationAxes::setRotation, M_VALUED));
    built = (built && createPortIn<darray3E, RotationAxes>(this, &mimmo::RotationAxes::setAxesOrigin, M_POINT2));
    built = (built && createPortIn<dmatrix33E, RotationAxes>(this, &mimmo::RotationAxes::setAxes, M_AXES));
    built = (built && createPortOut<darray3E, RotationAxes>(this, &mimmo::RotationAxes::getRotatedOrigin, M_POINT));
    built = (built && createPortOut<dmatrix33E, RotationAxes>(this, &mimmo::RotationAxes::getRotatedAxes, M_AXES));
    m_arePortsBuilt = built;
};

/*!It sets the origin and direction of the rotation axis.
 * \param[in] origin Origin of rotation axis.
 * \param[in] direction Direction of rotation axis.
 */
void
RotationAxes::setAxis(darray3E origin, darray3E direction){
    setOrigin(origin);
    setDirection(direction);
}

/*!It sets the origin of the rotation axis.
 * \param[in] origin Origin of rotation axis.
 */
void
RotationAxes::setOrigin(darray3E origin){
    m_origin = origin;
}

/*!It sets the direction of the rotation axis.
 * \param[in] direction Direction of rotation axis.
 */
void
RotationAxes::setDirection(darray3E direction){
    m_direction = direction;
    double norm = norm2(m_direction);
    if(norm > std::numeric_limits<double>::min()){
        m_direction /= norm;
    }
}

/*!It sets the value of the rotation.
 * \param[in] alpha Value of rotation axis.
 */
void
RotationAxes::setRotation(double alpha){
    m_alpha = alpha;
}

/*!It sets the reference system to be rotated.
 * \param[in] axes Original reference system.
 */
void
RotationAxes::setAxes(dmatrix33E axes){
    m_axes = axes;
}

/*!It sets the origin of the reference system to be rotated.
 * \param[in] axes_origin Origin of reference system.
 */
void
RotationAxes::setAxesOrigin(darray3E axes_origin){
    m_axes_origin = axes_origin;
}

/*!It gets the rotated reference system.
 * \return Rotated reference system.
 */
dmatrix33E
RotationAxes::getRotatedAxes(){
    return(m_rotax);
}

/*!It gets the rotated origin of the reference system.
 * \return Rotated origin of reference system.
 */
darray3E
RotationAxes::getRotatedOrigin(){
    return(m_rotax_origin);
}

/*!Execution command. It saves in "rot"-terms the modified axes and origin, by the
 * rotation conditions. This terms can be recovered and passed by a pin to a child object
 * by the related get-methods.
 */
void
RotationAxes::execute(){

    //Rotation of origin
    m_rotax_origin = {{0,0,0}};
    m_axes_origin += (-1.0)*m_origin;
    //rodrigues formula
    m_rotax_origin = m_axes_origin * std::cos(m_alpha) +
            dotProduct(m_direction, m_axes_origin) * (1.0 - std::cos(m_alpha)) * m_direction +
            crossProduct(m_direction, m_axes_origin) * std::sin(m_alpha);

    m_rotax_origin += m_origin;
    m_axes_origin += m_origin;

    //rotation of axes
    m_rotax.fill(darray3E{{0,0,0}});
    //rodrigues formula
    for (int i=0; i<3; i++){
        m_rotax[i] = m_axes[i] * std::cos(m_alpha) +
                dotProduct(m_direction, m_axes[i]) * (1.0 - std::cos(m_alpha)) * m_direction +
                crossProduct(m_direction, m_axes[i]) * std::sin(m_alpha);
    }

};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
RotationAxes::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

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

    if(slotXML.hasSection("RefSystem")){
        const bitpit::Config::Section & rfXML = slotXML.getSection("RefSystem");
        std::string rootAxis = "axis";
        std::string axis;
        dmatrix33E temp;
        temp[0].fill(0.0); temp[0][0] = 1.0;
        temp[1].fill(0.0); temp[1][1] = 1.0;
        temp[2].fill(0.0); temp[2][2] = 1.0;
        for(int i=0; i<3; ++i){
            axis = rootAxis + std::to_string(i);
            std::string input = rfXML.get(axis);
            input = bitpit::utils::string::trim(input);
            if(!input.empty()){
                std::stringstream ss(input);
                for(auto &val : temp[i]) ss>>val;
            }
        }
        setAxes(temp);
    }

    if(slotXML.hasOption("OriginRS")){
        std::string input = slotXML.get("OriginRS");
        input = bitpit::utils::string::trim(input);
        darray3E temp = {{0.0,0.0,0.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            for(auto &val : temp) ss>>val;
        }
        setAxesOrigin(temp);
    }
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
RotationAxes::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

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

    {
        bitpit::Config::Section & rsXML = slotXML.addSection("RefSystem");
        std::string rootAxis = "axis";
        std::string localAxis;
        int counter=0;
        for(auto &axis : m_axes){
            localAxis = rootAxis+std::to_string(counter);
            std::stringstream ss;
            ss<<std::scientific<<axis[0]<<'\t'<<axis[1]<<'\t'<<axis[2];
            rsXML.set(localAxis, ss.str());
            ++counter;
        }
    }

    {
        std::stringstream ss;
        ss<<std::scientific<<m_axes_origin[0]<<'\t'<<m_axes_origin[1]<<'\t'<<m_axes_origin[2];
        slotXML.set("OriginRS", ss.str());
    }

};

}
