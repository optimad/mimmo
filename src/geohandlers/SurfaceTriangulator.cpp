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

#include "SurfaceTriangulator.hpp"

namespace mimmo{

/*!
 * Constructor
 */
SurfaceTriangulator::SurfaceTriangulator(){
    m_name = "mimmo.SurfaceTriangulator";
    m_workOnTarget = false;

};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
SurfaceTriangulator::SurfaceTriangulator(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.SurfaceTriangulator";
    m_workOnTarget = false;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.SurfaceTriangulator"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Destructor;
 */
SurfaceTriangulator::~SurfaceTriangulator(){};

/*!
 * Copy constructor. Internal geometry data structure is not copied.
 */
SurfaceTriangulator::SurfaceTriangulator(const SurfaceTriangulator & other):BaseManipulation(other){
    m_workOnTarget = other.m_workOnTarget;
};

/*!
 * Assignment operator of the class. Copy class data as relative copy constructor.
 */
SurfaceTriangulator & SurfaceTriangulator::operator=(SurfaceTriangulator other){
    swap(other);
    return  *this;
};

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void SurfaceTriangulator::swap(SurfaceTriangulator &x) noexcept
{
   std::swap(m_workOnTarget, x.m_workOnTarget);
   std::swap(m_intPatch, x.m_intPatch);
   BaseManipulation::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void SurfaceTriangulator::buildPorts(){

    bool built = true;
    built = (built && createPortIn<MimmoObject *, SurfaceTriangulator>(this, &mimmo::SurfaceTriangulator::setGeometry,M_GEOM, true));
    built = (built && createPortOut<MimmoObject*, SurfaceTriangulator>(this, &mimmo::SurfaceTriangulator::getGeometry, M_GEOM));
    m_arePortsBuilt = built;
};


/*!
 * \return pointer to geometry resulting by class manipulation
 */
MimmoObject*  SurfaceTriangulator::getGeometry(){
    if(m_workOnTarget)  return BaseManipulation::getGeometry();
    else                return m_intPatch.get();
};

/*!
 * \return true if the class is directly applying on target mesh, false if the class
 * is creating a newly indipendent patch.
 */
bool  SurfaceTriangulator::isWorkingOnTarget(){
    return m_workOnTarget;
};

/*!
 * Set reference geometry for your SurfaceTriangulator class. Geometry must
 * be a surface tessellation.
 * \param[in] geo pointer to target geometry
 */
void SurfaceTriangulator::setGeometry(MimmoObject * geo){
    if(!geo)    return;
    BaseManipulation::setGeometry(geo);
};

/*!
 * If true, activate option to apply triangulator directly on the target geometry object linked to the class
 * with setGeometry method. If false, create internally an indipendent clone of the target geometry and work on it.
 * False is default.
 * \param[in] flag activation flag
 */
void        SurfaceTriangulator::setWorkOnTarget(bool flag){
    m_workOnTarget = flag;
};

/*!
 * Triangulate your target non-homogeneous surface.
 */
void SurfaceTriangulator::execute(){

    MimmoObject * target = BaseManipulation::getGeometry();

    if (target == NULL){
        throw std::runtime_error (m_name + " : NULL pointer to linked geometry");
    }

    if(!m_workOnTarget){
        //create internal patch as clone of the linked geometry and set target pointing to it
        m_intPatch = std::move(target->clone());
        target = m_intPatch.get();
    }

    target->triangulate();
};

/*!
 * Plot optional results during execution, that is the newly triangulated surface
 */
void     SurfaceTriangulator::plotOptionalResults(){

    bitpit::VTKUnstructuredGrid & vtk = getGeometry()->getPatch()->getVTK();
    vtk.setDirectory(m_outputPlot + "/");
    vtk.setName(m_name +std::to_string(getId()));
    getGeometry()->getPatch()->write();
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SurfaceTriangulator::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    //start absorbing
    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("WorkOnTarget")){
        std::string input = slotXML.get("WorkOnTarget");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setWorkOnTarget(value);
    }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SurfaceTriangulator::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML,name);
    slotXML.set("WorkOnTarget", std::to_string(int(m_workOnTarget)));
};

}
