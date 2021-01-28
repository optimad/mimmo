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
\ *---------------------------------------------------------------------------*/

#include "Primitive.hpp"

namespace mimmo{

/*!
 * Basic Constructor
 */
Primitive::Primitive(){
    m_name = "mimmo.Primitive";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
Primitive::Primitive(const bitpit::Config::Section & rootXML){
    m_name = "mimmo.Primitive";
    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.Primitive"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*! Destructor */
Primitive::~Primitive(){};

/*! Copy Constructor
 *\param[in] other Primitive object
 */
Primitive::Primitive(const Primitive & other):BaseManipulation(other), UStructMesh(other){
};

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void Primitive::swap(Primitive & x) noexcept
{
    UStructMesh::swap(x);
    BaseManipulation::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void Primitive::buildPorts(){

    bool built = true;
    built = (built && createPortIn<darray3E, Primitive>(this, &mimmo::Primitive::setInfLimits, M_INFLIMITS));
    built = (built && createPortIn<dmatrix33E, Primitive>(this, &mimmo::Primitive::setRefSystem, M_AXES));
    built = (built && createPortIn<darray3E, Primitive>(this, &mimmo::Primitive::setSpan, M_SPAN));
    built = (built && createPortIn<darray3E, Primitive>(this, &mimmo::Primitive::setOrigin, M_POINT));
    built = (built && createPortIn<mimmo::ShapeType, Primitive>(this, &mimmo::Primitive::setShape, M_SHAPE));
    built = (built && createPortIn<const BasicShape *, Primitive>(this, &mimmo::Primitive::setShape, M_COPYSHAPE));
    built = (built && createPortIn<int, Primitive>(this, &mimmo::Primitive::setShape, M_SHAPEI));

    // creating output ports
    built = (built && createPortOut<BasicShape *, Primitive>(this, &mimmo::Primitive::getShape, M_COPYSHAPE));

    m_arePortsBuilt = built;
};

/*!
 * Clean every setting and data hold by the class
 */
void Primitive::clearPrimitive(){
    clear(); //base manipulation stuff clear
    clearMesh(); // structured mesh cleaned
};

/*!
 * Execute your object. Currently do nothing.
 */
void
Primitive::execute(){
};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void Primitive::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("Shape")){
        std::string input = slotXML.get("Shape");
        input = bitpit::utils::string::trim(input);

        if(input == "CYLINDER"){
            setShape(ShapeType::CYLINDER);
        }else if(input =="SPHERE"){
            setShape(ShapeType::SPHERE);
        }else if(input =="WEDGE"){
            setShape(ShapeType::WEDGE);
        }else{
            setShape(ShapeType::CUBE);
        }
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
    };

    if(slotXML.hasOption("Span")){
        std::string input = slotXML.get("Span");
        input = bitpit::utils::string::trim(input);
        darray3E temp = {{0.0,0.0,0.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            for(auto &val : temp) ss>>val;
        }
        setSpan(temp);
    };

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
        setRefSystem(temp);
    };

    if(slotXML.hasOption("InfLimits")){
        std::string input = slotXML.get("InfLimits");
        input = bitpit::utils::string::trim(input);
        darray3E temp = {{0.0,0.0,0.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            for(auto &val : temp) ss>>val;
        }
        setInfLimits(temp);
    };

}

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void Primitive::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    std::string towrite = "CUBE";

    if(getShapeType() == ShapeType::CYLINDER){
        towrite = "CYLINDER";
    } else if(getShapeType() == ShapeType::SPHERE){
        towrite = "SPHERE";
    } else if(getShapeType() == ShapeType::WEDGE){
        towrite = "WEDGE";
    }
    slotXML.set("Shape", towrite);

    {
        std::stringstream ss;
        ss<<std::scientific<<getOrigin()[0]<<'\t'<<getOrigin()[1]<<'\t'<<getOrigin()[2];
        slotXML.set("Origin", ss.str());
    }

    {
        std::stringstream ss;
        ss<<std::scientific<<getSpan()[0]<<'\t'<<getSpan()[1]<<'\t'<<getSpan()[2];
        slotXML.set("Span", ss.str());
    }

    {
        auto rs = getRefSystem();
        bitpit::Config::Section & rsXML = slotXML.addSection("RefSystem");
        std::string rootAxis = "axis";
        std::string localAxis;
        int counter=0;
        for(auto &axis : rs){
            localAxis = rootAxis+std::to_string(counter);
            std::stringstream ss;
            ss<<std::scientific<<axis[0]<<'\t'<<axis[1]<<'\t'<<axis[2];
            rsXML.set(localAxis, ss.str());
            ++counter;
        }
    }

};

}
