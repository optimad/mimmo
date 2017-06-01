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
#include "MultiApply.hpp"

namespace mimmo{

/*!
 * Default constructor of Apply
 */
MultiApply::MultiApply():BaseManipulation(){
    m_name = "mimmo.MultiApply";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
MultiApply::MultiApply(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.MultiApply";

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::trim(input);
    if(input == "mimmo.MultiApply"){
        absorbSectionXML(rootXML);
    }else{
        (*m_log)<<"Warning in custom xml mimmo::MultiApply constructor. No valid xml data found"<<std::endl;
    };
}

/*!Default destructor of Apply
 */
MultiApply::~MultiApply(){};

/*!Copy constructor of Apply.
 */
MultiApply::MultiApply(const MultiApply & other):BaseManipulation(){
    *this = other;
};

/*!Assignement operator of Apply.
 */
MultiApply & MultiApply::operator=(const MultiApply & other){
    *(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
    m_input = other.m_input;
    return(*this);
};

/*! It builds the input/output ports of the object
 */
void
MultiApply::buildPorts(){
    bool built = true;
    built = (built && createPortIn<std::unordered_map<MimmoObject*, dvecarr3E*>, MultiApply>(this, &MultiApply::setInputList, PortType::M_UMGEOVFD, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::MIMMO_VECARR3FLOAT_, true));
    m_arePortsBuilt = built;
};


/*! Set the input of the class as std::pair of target MimmoObject * and
 * relative displacement field dvecarr3E *. If target geometry already exist or it is NULL
 * no input is added to the class.
 * \param[in] input geometry and its displacements
 */
void
MultiApply::addInput(std::pair<MimmoObject*, dvecarr3E*> input){

    if (input.first == NULL || m_input.count(input.first) > 0 ) return;
    m_input[input.first] = input.second;
}


/*! Set a whole list of inputs of the class as an unorderd map of target MimmoObject * as key and
 * relative displacement field dvecarr3E * as argument. If a target geometry already exists in the member list or it is NULL
 * no input is added to the class. Multiple insertion of list of inputs are appended in the unique member list.
 * \param[in] input list of geometries and displacements
 */
void
MultiApply::setInputList(std::unordered_map<MimmoObject*, dvecarr3E*> input){

    for(auto & val : input){

        addInput(val);
    }
}

/*!
 * Clear the input list stored into the class
 */
void
MultiApply::clearList(){
    m_input.clear();
}

/*!Execution command.
 * It applies the deformations stored in the class  to their relative linked geometry.
 * After exec() the original geometries will be permanently modified.
 */
void
MultiApply::execute(){


    if (m_input.empty()) return;


    for(auto val: m_input){



        dvecarr3E vertex = val.first->getVertexCoords();
        long nv = val.first->getNVertex();
        nv = long(std::min(int(nv), int(val.second->size())));

        livector1D & idmap = val.first->getMapData();
        for (long i=0; i<nv; i++){

            vertex[i] += (*(val.second))[i];
            val.first->modifyVertex(vertex[i], idmap[i]);
        }

    }
    return;
};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void MultiApply::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    if(slotXML.hasOption("Priority")){
        input = slotXML.get("Priority");
        int value =0;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss>>value;
        }
        setPriority(value);
    };

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void MultiApply::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    slotXML.set("ClassName", m_name);
    slotXML.set("Priority", std::to_string(getPriority()));

};	


}


