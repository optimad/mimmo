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
 \ *---------------------------------------------------------------------------*/

#include "FVMeshSelection.hpp"
namespace mimmo{

/*!
 * Basic Constructor.
 * Parameter topo get topology of the target MimmoFvMesh where performing extraction:
 * - 1 for volume bulk and surface boundary
 * - 2 for surface bulk and 3DCurve boundary
 * No other values are allowed
 * \param[in] topo topology of the target MimmoFvMesh.
 */
FVSelectionByBox::FVSelectionByBox(int topo): FVGenericSelection(topo), Cube(){
    m_name = "mimmo.FVSelectionByBox";
    m_type = FVSelectionType::BOX;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
FVSelectionByBox::FVSelectionByBox(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.FVSelectionByBox";
    m_type = FVSelectionType::BOX;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);

    std::string fallback_name2 = "1";
    std::string input_topo = rootXML.get("Topology", fallback_name2);
    input_topo = bitpit::utils::string::trim(input_topo);

    int topo = std::stoi(input_topo);
    topo = std::min(2,std::max(1, topo));
    m_topo = topo;

    if(input == "mimmo.FVSelectionByBox"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Custom Constructor
 * \param[in] topo topology of the target MimmoFvMesh.See Basic Constructor doxy.
 * \param[in] origin Origin of the box->baricenter.
 * \param[in] span     Span of the box, width/height/depth.
 */
FVSelectionByBox::FVSelectionByBox(int topo, darray3E origin, darray3E span): FVGenericSelection(topo), Cube(origin, span){
    m_name = "mimmo.FVSelectionByBox";
    m_type = FVSelectionType::BOX;
};

/*!
 * Destructor
 */
FVSelectionByBox::~FVSelectionByBox(){};

/*!
 * Copy Constructor
 */
FVSelectionByBox::FVSelectionByBox(const FVSelectionByBox & other):FVGenericSelection(other), Cube(other){};

/*!
 * Copy operator - CAS idiom
 */
FVSelectionByBox & FVSelectionByBox::operator=(FVSelectionByBox other){
    swap(other);
    return *this;
};

/*!
 * Swap function. Assign data of this class to another of the same type and vice-versa.
 * Resulting patch of selection is not swapped, ever.
 * \param[in] x FVSelectionByBox object
 */
void FVSelectionByBox::swap(FVSelectionByBox & x) noexcept
{
   FVGenericSelection::swap(x);
   Cube::swap(x);
}

/*!
 * It builds the input/output ports of the object.
 */
void
FVSelectionByBox::buildPorts(){

    bool built = true;

    FVGenericSelection::buildPorts();

    built = (built && createPortIn<darray3E, FVSelectionByBox>(this, &FVSelectionByBox::setOrigin, M_POINT ));
    built = (built && createPortIn<darray3E, FVSelectionByBox>(this, &FVSelectionByBox::setSpan, M_SPAN));
    built = (built && createPortIn<dmatrix33E, FVSelectionByBox>(this, &FVSelectionByBox::setRefSystem, M_AXES));

    m_arePortsBuilt = built;
};

/*!
 * Clear content of the class
 */
void
FVSelectionByBox::clear(){
    m_volpatch.reset(nullptr);
    m_bndpatch.reset(nullptr);
    m_dual = false;
    m_bndgeometry = NULL;
    BaseManipulation::clear();
};



/*!
 * Extract portion of target bulk-boundary geometry got by box
 * \param[out] bulk cell ids of target bulk extracted
 * \param[out] boundary cell ids of target boundary extracted, divided by PID if any.
 */
void
FVSelectionByBox::extractSelection(livector1D & bulk, livector1D & boundary){

    if(m_dual){
        bulk = excludeGeometry(m_geometry);
        boundary = excludeGeometry(m_bndgeometry);
    }else{
        bulk = includeGeometry(m_geometry);
        boundary = includeGeometry(m_bndgeometry);
    }
};


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
FVSelectionByBox::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("Topology")){
        std::string input = slotXML.get("Topology");
        input = bitpit::utils::string::trim(input);
        int value = 1;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        if (m_topo != value)  warningXML(m_log, m_name);
    }else{
        warningXML(m_log, m_name);
    }


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
FVSelectionByBox::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML,name);

    slotXML.set("Topology", std::to_string(m_topo));
    slotXML.set("Dual", std::to_string(int(m_dual)));
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
