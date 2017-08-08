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

#include "MeshSelection.hpp"
#include "levelSet.hpp"
#include <cstddef>
namespace mimmo{

/*!
 * Basic Constructor
 */
SelectionBySphere::SelectionBySphere(){
    m_name = "mimmo.SelectionBySphere";
    m_type = SelectionType::SPHERE;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
SelectionBySphere::SelectionBySphere(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.SelectionBySphere";
    m_type = SelectionType::SPHERE;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.SelectionBySphere"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Custom Constructor. Pay attention span of angular and polar coords are at most 2*pi and pi respectively. 
 *  Inf limits of polar coordinate must be > 0 and < pi.
 * \param[in] origin Origin of the sphere->baricenter
 * \param[in] span     Span of the cylinder, main radius/ span of angular coord in radians/span of the polar coord in radians
 * \param[in] infLimTheta    Starting origin of the angular coordinate. default is 0 radians.
 * \param[in] infLimPhi    Starting origin of the polar coordinate. default is 0 radians.
 * \param[in] target    Pointer to a target geometry
 */
SelectionBySphere::SelectionBySphere(darray3E origin, darray3E span, double infLimTheta, double infLimPhi, MimmoObject * target){
    m_name = "mimmo.SelectionBySphere";
    m_type = SelectionType::SPHERE;
    setGeometry(target);
    setOrigin(origin);
    setSpan(span[0],span[1],span[2]);
    setInfLimits(infLimTheta,1);
    setInfLimits(infLimPhi,2);
};

/*!
 * Destructor
 */
SelectionBySphere::~SelectionBySphere(){};

/*!
 * Copy Constructor
 */
SelectionBySphere::SelectionBySphere(const SelectionBySphere & other):GenericSelection(other), Sphere(other){
};

/*!
 * Copy operator
 */
SelectionBySphere & SelectionBySphere::operator=(SelectionBySphere other){
    swap(other);
    return *this;
};

/*!
 * Swap function. Assign data of this class to another of the same type and vice-versa.
 * Resulting patch of selection is not swapped, ever.
 * \param[in] x SelectionBySphere object
 */
void SelectionBySphere::swap(SelectionBySphere & x) noexcept
{
    GenericSelection::swap(x);
    Sphere::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
SelectionBySphere::buildPorts(){

    bool built = true;

    GenericSelection::buildPorts();

    built = (built && createPortIn<darray3E, SelectionBySphere>(this, &SelectionBySphere::setOrigin,M_POINT));
    built = (built && createPortIn<darray3E, SelectionBySphere>(this, &SelectionBySphere::setSpan, M_SPAN));
    built = (built && createPortIn<dmatrix33E, SelectionBySphere>(this, &SelectionBySphere::setRefSystem, M_AXES));
    built = (built && createPortIn<darray3E, SelectionBySphere>(this, &SelectionBySphere::setInfLimits, M_INFLIMITS));

    m_arePortsBuilt = built;
};

/*!
 * Clear content of the class
 */
void SelectionBySphere::clear(){
    m_subpatch.release();
    BaseManipulation::clear();
};

/*!
 * Extract portion of target geometry which are enough near to the external geometry provided
 * \return ids of cell of target tesselation extracted 
 */
livector1D
SelectionBySphere::extractSelection(){
    switch(m_topo){
    case 3:
        if(m_dual)  return    excludeCloudPoints(getGeometry());
        else        return    includeCloudPoints(getGeometry());
        break;
    default:
        if(m_dual)  return    excludeGeometry(getGeometry());
        else        return    includeGeometry(getGeometry());
        break;
    }
};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
SelectionBySphere::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    //start absorbing
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
        darray3E temp = {{1.0,2.0*M_PI,M_PI}};
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
SelectionBySphere::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);
    
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

    {
        darray3E inflim = getInfLimits();
        std::stringstream ss;
        ss<<std::scientific<<inflim[0]<<'\t'<<inflim[1]<<'\t'<<inflim[2];
        slotXML.set("InfLimits",ss.str());
    }

};

}
