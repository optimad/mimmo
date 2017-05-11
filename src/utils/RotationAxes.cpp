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
    m_origin = origin;
    m_direction = direction;
    m_name = "mimmo.RotationAxes";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
RotationAxes::RotationAxes(const bitpit::Config::Section & rootXML){

    m_origin.fill(0.0);
    m_direction.fill(0.0);
    m_name = "mimmo.RotationAxes";

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::trim(input);
    if(input == "mimmo.RotationAxes"){
        absorbSectionXML(rootXML);
    }else{
        std::cout<<"Warning in custom xml mimmo::RotationAxes constructor. No valid xml data found"<<std::endl;
    };
}

/*!Default destructor of RotationAxes
 */
RotationAxes::~RotationAxes(){};

/*!Copy constructor of RotationAxes.
 */
RotationAxes::RotationAxes(const RotationAxes & other):BaseManipulation(other){
    m_origin = other.m_origin;
    m_direction = other.m_direction;
};

/*!Assignement operator of RotationAxes.
 */
RotationAxes & RotationAxes::operator=(const RotationAxes & other){
    *(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
    m_origin = other.m_origin;
    m_direction = other.m_direction;
    return(*this);
};


/*! It builds the input/output ports of the object
 */
void
RotationAxes::buildPorts(){
    bool built = true;
    built = (built && createPortIn<darray3E, RotationAxes>(&m_origin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<darray3E, RotationAxes>(&m_direction, PortType::M_AXIS, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<double, RotationAxes>(&m_alpha, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<darray3E, RotationAxes>(&m_axes_origin, PortType::M_POINT2, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<dmatrix33E, RotationAxes>(&m_axes, PortType::M_AXES, mimmo::pin::containerTAG::ARR3ARR3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<darray3E, RotationAxes>(this, &mimmo::RotationAxes::getRotatedOrigin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<dmatrix33E, RotationAxes>(this, &mimmo::RotationAxes::getRotatedAxes, PortType::M_AXES, mimmo::pin::containerTAG::ARR3ARR3, mimmo::pin::dataTAG::FLOAT));
    m_arePortsBuilt = built;
};

/*!It sets the origin and direction of the rotation axis.
 * \param[in] origin Origin of rotation axis.
 * \param[in] direction Direction of rotation axis.
 */
void
RotationAxes::setAxis(darray3E origin, darray3E direction){
    m_origin = origin;
    m_direction = direction;
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
    double L = norm2(m_direction);
    for (int i=0; i<3; i++)
        m_direction[i] /= L;
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
    m_axes_origin -= m_origin;
    //rodrigues formula
    m_rotax_origin = m_axes_origin * cos(m_alpha) +
            dotProduct(m_direction, m_axes_origin) * (1 - cos(m_alpha)) * m_direction +
            crossProduct(m_direction, m_axes_origin) * sin(m_alpha);

    m_rotax_origin += m_origin;
    m_axes_origin += m_origin;

    //rotation of axes
    m_rotax.fill(darray3E{{0,0,0}});
    //rodrigues formula
    for (int i=0; i<3; i++){
        m_rotax[i] = m_axes[i] * cos(m_alpha) +
                dotProduct(m_direction, m_axes[i]) * (1 - cos(m_alpha)) * m_direction +
                crossProduct(m_direction, m_axes[i]) * sin(m_alpha);
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

    if(slotXML.hasOption("Priority")){
        std::string input = slotXML.get("Priority");
        int value =0;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss>>value;
        }
        setPriority(value);
    };

    if(slotXML.hasOption("Origin")){
        std::string input = slotXML.get("Origin");
        input = bitpit::utils::trim(input);
        darray3E temp = {{0.0,0.0,0.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            for(auto &val : temp) ss>>val;
        }
        setOrigin(temp);
    }

    if(slotXML.hasOption("Direction")){
        std::string input = slotXML.get("Direction");
        input = bitpit::utils::trim(input);
        darray3E temp = {{0.0,0.0,0.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            for(auto &val : temp) ss>>val;
        }
        setDirection(temp);
    }

    if(slotXML.hasOption("Rotation")){
        std::string input = slotXML.get("Rotation");
        input = bitpit::utils::trim(input);
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
            input = bitpit::utils::trim(input);
            if(!input.empty()){
                std::stringstream ss(input);
                for(auto &val : temp[i]) ss>>val;
            }
        }
        setAxes(temp);
    }

    if(slotXML.hasOption("OriginRS")){
        std::string input = slotXML.get("OriginRS");
        input = bitpit::utils::trim(input);
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

    slotXML.set("ClassName", m_name);
    slotXML.set("Priority", std::to_string(getPriority()));


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

