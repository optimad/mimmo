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

#include "StitchGeometry.hpp"

using namespace std;
using namespace bitpit;
namespace mimmo{

/*!
 * Default constructor of StitchGeometry.
 * Format admissible are linked to your choice of topology. See FileType enum
 * \param[in] topo    Topology of your geometries (1-surface, 2-volume, 3-points cloud)
 */
StitchGeometry::StitchGeometry(int topo){
    m_name         = "mimmo.StitchGeometry";
    m_geocount = 0;
    m_topo     = std::min(1, topo);
    if(m_topo > 3)    m_topo = 1;
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
StitchGeometry::StitchGeometry(const bitpit::Config::Section & rootXML){

    std::string fallback_name = "ClassNONE";
    std::string fallback_topo = "-1";
    std::string input_name = rootXML.get("ClassName", fallback_name);
    std::string input_topo = rootXML.get("Topology", fallback_topo);
    input_name = bitpit::utils::trim(input_name);
    input_topo = bitpit::utils::trim(input_topo);

    int topo = std::stoi(input_topo);
    m_topo = std::max(1,topo);
    if (m_topo >3) m_topo = 1;

    m_geocount = 0;
    m_name = "mimmo.StitchGeometry";


    if(input_name == "mimmo.StitchGeometry"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Default destructor of StitchGeometry.
 */
StitchGeometry::~StitchGeometry(){
    clear();
};

/*!Copy constructor of StitchGeometry.Soft Copy of MimmoObject;
 */
StitchGeometry::StitchGeometry(const StitchGeometry & other):BaseManipulation(){
    *this = other;
};

/*!
 * Assignement operator of StitchGeometry. Soft copy of MimmoObject
 */
StitchGeometry & StitchGeometry::operator=(const StitchGeometry & other){
    clear();
    *(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(&other));

    m_topo = other.m_topo;
    m_extgeo = other.m_extgeo;
    m_geocount = other.m_geocount;

    //warning the internal data structure and its division map is not copied. Relaunch the execution eventually to fill it.
    return *this;
};

/*!
 * Building ports of the class
 */
void
StitchGeometry::buildPorts(){
    bool built = true;
    built = (built && createPortIn<std::vector<MimmoObject*>, StitchGeometry>(this, &mimmo::StitchGeometry::setGeometries, PortType::M_VECGEOM, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::MIMMO_, true, 1));
    built = (built && createPortIn<MimmoObject*, StitchGeometry>(this, &mimmo::StitchGeometry::addGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_, true, 1));

    built = (built && createPortOut<MimmoObject*, StitchGeometry>(this, &mimmo::StitchGeometry::getGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    m_arePortsBuilt = built;
}

/*!
 * Return kind of topology supported by the object.
 * \return topology supported (1-surface, 2-volume, 3-points cloud)
 */
int 
StitchGeometry::getTopology(){
    return m_topo;
}

/*!
 * Get current stitched geometry pointer. Reimplementation of BaseManipulation::getGeometry,
 * \return Pointer to stitched geometry.
 */
MimmoObject *
StitchGeometry::getGeometry(){
    return m_patch.get();
}

/*!
 * Add an external geometry to be stitched.Topology of the geometry must be coeherent
 * with topology of the class;
 * \param[in] geo  Pointer to MimmoObject
 */
void
StitchGeometry::addGeometry(MimmoObject* geo){
    if(geo->isEmpty()) return;
    if(geo->getType() != m_topo)    return;
    if(m_extgeo.count(geo)    > 0)    return;

    m_extgeo.insert(std::make_pair(geo,m_geocount) );
    m_geocount++;
};

/*!
 * Set geometries to be stitched. List will be saved as is, replacing any other saved list. 
 * Reimplementation of BaseManipulation::setGeometry
 * \param[in] external Pointer to stitched geometry.
 */
void
StitchGeometry::setGeometries(std::vector<MimmoObject *> external){

    for(auto & obj: external){
        addGeometry(obj);
    }
};

/*!
 * Check if stitched geometry is present or not.
 * \return true - no geometry present, false otherwise.
 */
bool 
StitchGeometry::isEmpty(){
    return (m_patch.get() == NULL);
}

/*!
 * Clear all stuffs in your class
 */
void
StitchGeometry::clear(){
    m_extgeo.clear();
    m_geocount = 0;
    m_patch.reset(nullptr);
    BaseManipulation::clear();
};

/*!Execution command.
 * It stitches together multiple geometries in the same object.
 */
void
StitchGeometry::execute(){
    if(m_extgeo.empty())    return;

    std::unique_ptr<MimmoObject> dum(new MimmoObject(m_topo));
    long nCells = 0;
    long nVerts = 0;

    for(auto &obj:m_extgeo){
        nCells += obj.first->getNCells();
        nVerts += obj.first->getNVertex();
    }

    //reserving memory
    dum->getPatch()->reserveVertices(nVerts);
    dum->getPatch()->reserveCells(nCells);

    long cV = 0, cC = 0;

    //start filling your stitched object.
    {
        //optional vars;
        long vId, cId;
        short PID;
        bitpit::ElementInfo::Type eltype;

        for(auto &obj : m_extgeo){

            //start extracting/reversing vertices of the current obj
            for(auto & vv : obj.first->getVertices()){
                vId = vv.getId();
                dum->addVertex(obj.first->getVertexCoords(vId), cV);
                //update map;
                cV++;
            }

            //start extracting/reversing cells of the current obj
            for(auto & cc : obj.first->getCells()){
                cId = cc.getId();
                PID = cc.getPID();
                eltype = cc.getType();
                //get the local connectivity and update with new vertex numbering;
                livector1D conn = obj.first->getCellConnectivity(cId);
                dum->addConnectedCell(conn, eltype, PID, cC);
                //update map;
                cC++;
            }
        }
    }//scope for optional vars;

    dum->cleanGeometry();

    m_patch = std::move(dum);
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void StitchGeometry::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    //checking topology
    if(slotXML.hasOption("Topology")){
        std::string input = slotXML.get("Topology");
        input = bitpit::utils::trim(input);
        int temptop = -1;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temptop;
        }
        if(m_topo != temptop)    return;
    }

    BaseManipulation::absorbSectionXML(slotXML, name);
    
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void StitchGeometry::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);
    slotXML.set("Topology", m_topo);

};

/*!
 * Plot stitched geometry in *.vtu file as optional result;
 */
void 
StitchGeometry::plotOptionalResults(){
    if(isEmpty()) return;
    std::string name = m_name + "_" + std::to_string(getClassCounter()) +  "_Patch";
    m_patch->getPatch()->write(name);
}

}
