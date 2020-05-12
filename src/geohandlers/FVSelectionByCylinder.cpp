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
#include <bitpit_common.hpp>

namespace mimmo{

/*!
 * Basic Constructor.
 * Parameter topo get topology of the target MimmoFvMesh where performing extraction:
 * - 1 for volume bulk and surface boundary
 * - 2 for surface bulk and 3DCurve boundary

 * No other values are allowed
 * \param[in] topo topology of the target MimmoFvMesh.
 */
FVSelectionByCylinder::FVSelectionByCylinder(int topo): FVGenericSelection(topo), Cylinder(){
    m_name = "mimmo.FVSelectionByCylinder";
    m_type = FVSelectionType::CYLINDER;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
FVSelectionByCylinder::FVSelectionByCylinder(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.FVSelectionByCylinder";
    m_type = FVSelectionType::CYLINDER;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);

    std::string fallback_name2 = "1";
    std::string input_topo = rootXML.get("Topology", fallback_name2);
    input_topo = bitpit::utils::string::trim(input_topo);

    int topo = std::stoi(input_topo);
    topo = std::min(2,std::max(1, topo));
    m_topo = topo;


    if(input == "mimmo.FVSelectionByCylinder"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Custom Constructor. Pay attention span of angular coordinate must be at most 2*pi.
 * \param[in] topo topology of the target MimmoFvMesh.See Basic Constructor doxy.
 * \param[in] origin Origin of the cylinder->baricenter
 * \param[in] span     Span of the cylinder, basis radius/ span of angular coord in radians/height
 * \param[in] infLimTheta    Starting origin of the angular coordinate. default is 0 radians.
 * \param[in] mainAxis    Orientation of the cylinder height axis
 */
FVSelectionByCylinder::FVSelectionByCylinder(int topo, darray3E origin, darray3E span,
                                            double infLimTheta, darray3E mainAxis):
                                            FVGenericSelection(topo), Cylinder(origin, span)
{
    m_name = "mimmo.FVSelectionByCylinder";
    m_type = FVSelectionType::CYLINDER;
    setInfLimits(infLimTheta,1);
    setRefSystem(2, mainAxis);
};

/*!
 * Destructor
 */
FVSelectionByCylinder::~FVSelectionByCylinder(){};

/*!
 * Copy Constructor
 */
FVSelectionByCylinder::FVSelectionByCylinder(const FVSelectionByCylinder & other):FVGenericSelection(other), Cylinder(other){};

/*!
 * Copy operator
 */
FVSelectionByCylinder & FVSelectionByCylinder::operator=(FVSelectionByCylinder other){
    swap(other);
    return *this;
};

/*!
 * Swap function. Assign data of this class to another of the same type and vice-versa.
 * Resulting patch of selection is not swapped, ever.
 * \param[in] x FVSelectionByCylinder object
 */
void FVSelectionByCylinder::swap(FVSelectionByCylinder & x) noexcept
{
    FVGenericSelection::swap(x);
    Cylinder::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
FVSelectionByCylinder::buildPorts(){

    bool built = true;

    FVGenericSelection::buildPorts();

    built = (built && createPortIn<darray3E, FVSelectionByCylinder>(this, &FVSelectionByCylinder::setOrigin,M_POINT));
    built = (built && createPortIn<darray3E, FVSelectionByCylinder>(this, &FVSelectionByCylinder::setSpan, M_SPAN));
    built = (built && createPortIn<dmatrix33E, FVSelectionByCylinder>(this, &FVSelectionByCylinder::setRefSystem, M_AXES));
    built = (built && createPortIn<darray3E, FVSelectionByCylinder>(this, &FVSelectionByCylinder::setInfLimits, M_INFLIMITS));

    m_arePortsBuilt = built;
};

/*!
 * Clear content of the class
 */
void
FVSelectionByCylinder::clear(){
    m_volpatch.reset();
    m_bndpatch.reset();
    m_dual = false;
    m_bndgeometry.reset();
    BaseManipulation::clear();
};

/*!
 * Extract portion of target bulk-boundary geometry got by cylinder
 * \param[out] bulk cell ids of target bulk extracted
 * \param[out] boundary cell ids of target boundary extracted, divided by PID if any.
 */
void
FVSelectionByCylinder::extractSelection(livector1D & bulk, livector1D & boundary){

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
FVSelectionByCylinder::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

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
        darray3E temp = {{1.0,2.0*BITPIT_PI,1.0}};
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

    if(slotXML.hasOption("InfLimits")){
        std::string input = slotXML.get("InfLimits");
        input = bitpit::utils::string::trim(input);
        darray3E temp = {{0.0,0.0,0.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp[0]>>temp[1]>>temp[2];
            setInfLimits(temp);
        }else{
            setInfLimits(temp);
        }
    }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
FVSelectionByCylinder::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

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

    {
        darray3E inflim = getInfLimits();
        std::stringstream ss;
        ss<<std::scientific<<inflim[0]<<'\t'<<inflim[1]<<'\t'<<inflim[2];
        slotXML.set("InfLimits",ss.str());
    }

};

}
