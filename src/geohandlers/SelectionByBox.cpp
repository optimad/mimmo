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

#include "MeshSelection.hpp"

namespace mimmo{

/*!
 * Basic Constructor
 */
SelectionByBox::SelectionByBox(){
    m_name = "mimmo.SelectionByBox";
    m_type = SelectionType::BOX;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
SelectionByBox::SelectionByBox(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.SelectionByBox";
    m_type = SelectionType::BOX;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.SelectionByBox"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Custom Constructor
 * \param[in] origin Origin of the box->baricenter.
 * \param[in] span     Span of the box, width/height/depth.
 * \param[in] target    Pointer to MimmoObject target geometry.
 */
SelectionByBox::SelectionByBox(darray3E origin, darray3E span, mimmo::MimmoSharedPointer<MimmoObject> target){
    m_name = "mimmo.SelectionByBox";
    m_type = SelectionType::BOX;
    setGeometry(target);
    setOrigin(origin);
    setSpan(span[0],span[1],span[2]);
};

/*!
 * Destructor
 */
SelectionByBox::~SelectionByBox(){};

/*!
 * Copy Constructor
 */
SelectionByBox::SelectionByBox(const SelectionByBox & other):GenericSelection(other), Cube(other){
};

/*!
 * Copy operator - CAS idiom
 */
SelectionByBox & SelectionByBox::operator=(SelectionByBox other){
    swap(other);
    return *this;
};

/*!
 * Swap function. Assign data of this class to another of the same type and vice-versa.
 * Resulting patch of selection is not swapped, ever.
 * \param[in] x SelectionByBox object
 */
void SelectionByBox::swap(SelectionByBox & x) noexcept
{
   GenericSelection::swap(x);
   Cube::swap(x);
}

/*!
 * It builds the input/output ports of the object.
 */
void
SelectionByBox::buildPorts(){

    GenericSelection::buildPorts();
    bool built = m_arePortsBuilt;

    built = (built && createPortIn<darray3E, SelectionByBox>(this, &SelectionByBox::setOrigin, M_POINT ));
    built = (built && createPortIn<darray3E, SelectionByBox>(this, &SelectionByBox::setSpan, M_SPAN));
    built = (built && createPortIn<dmatrix33E, SelectionByBox>(this, &SelectionByBox::setRefSystem, M_AXES));

    m_arePortsBuilt = built;
};

/*!
 * Clear content of the class
 */
void
SelectionByBox::clear(){
    m_subpatch.reset();
    BaseManipulation::clear();
};



/*!
 * Extract portion of target geometry which are enough near
 * to the external geometry provided
 * \return ids of cell of target tessellation extracted
 */
livector1D
SelectionByBox::extractSelection(){
    switch(m_topo){
    case 3:
        if(m_dual)    return    excludeCloudPoints(getGeometry());
        else          return    includeCloudPoints(getGeometry());
        break;
    default:
        if(m_dual)    return    excludeGeometry(getGeometry());
        else          return    includeGeometry(getGeometry());
        break;
    }
};


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
SelectionByBox::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("Dual")){
        std::string input = slotXML.get("Dual");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setDual(value);
    }

    if(slotXML.hasOption("Origin")){
        std::string input = slotXML.get("Origin");
        input = bitpit::utils::string::trim(input);
        darray3E temp = {{0.0,0.0,0.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp[0]>>temp[1]>>temp[2];
            setOrigin(temp);
        }else{
            setOrigin(temp);
        }
    }

    if(slotXML.hasOption("Span")){
        std::string input = slotXML.get("Span");
        input = bitpit::utils::string::trim(input);
        darray3E temp = {{1.0,1.0,1.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp[0]>>temp[1]>>temp[2];
            setSpan(temp);
        }else{
            setSpan(temp);
        }
    }

    if(slotXML.hasSection("RefSystem")){

        const bitpit::Config::Section & axesXML = slotXML.getSection("RefSystem");
        dmatrix33E axes;
        for(int i=0; i<3; ++i){
            axes[i].fill(0.0);
            axes[i][i] = 1.0;
        }

        if(axesXML.hasOption("axis0")){
            std::string input = axesXML.get("axis0");
            input = bitpit::utils::string::trim(input);
            if(!input.empty()){
                std::stringstream ss(input);
                ss>>axes[0][0]>>axes[0][1]>>axes[0][2];
            }
        }

        if(axesXML.hasOption("axis1")){
            std::string input = axesXML.get("axis1");
            input = bitpit::utils::string::trim(input);
            if(!input.empty()){
                std::stringstream ss(input);
                ss>>axes[1][0]>>axes[1][1]>>axes[1][2];
            }
        }

        if(axesXML.hasOption("axis2")){
            std::string input = axesXML.get("axis2");
            input = bitpit::utils::string::trim(input);
            if(!input.empty()){
                std::stringstream ss(input);
                ss>>axes[2][0]>>axes[2][1]>>axes[2][2];
            }
        }

        setRefSystem(axes);
    }
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
SelectionByBox::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML,name);

    int value = m_dual;
    slotXML.set("Dual", std::to_string(value));

    {
        darray3E org = getOrigin();
        std::stringstream ss;
        ss<<std::scientific<<org[0]<<'\t'<<org[1]<<'\t'<<org[2];
        slotXML.set("Origin",ss.str());
    }

    {
        darray3E span = getSpan();
        std::stringstream ss;
        ss<<std::scientific<<span[0]<<'\t'<<span[1]<<'\t'<<span[2];
        slotXML.set("Span",ss.str());
    }

    {
        dmatrix33E axes = getRefSystem();
        bitpit::Config::Section & axesXML = slotXML.addSection("RefSystem");

        for(int i=0; i<3; ++i){
            std::string name = "axis"+std::to_string(i);
            std::stringstream ss;
            ss<<std::scientific<<axes[i][0]<<'\t'<<axes[i][1]<<'\t'<<axes[i][2];
            axesXML.set(name, ss.str());
        }
    }

};

}
