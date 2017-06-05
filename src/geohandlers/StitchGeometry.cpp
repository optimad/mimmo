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

    m_buildBvTree = other.m_buildBvTree;
    m_buildKdTree = other.m_buildKdTree;
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
    built = (built && createPortIn<std::vector<MimmoObject*>, StitchGeometry>(this, &mimmo::StitchGeometry::setGeometry, PortType::M_VECGEOM, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortIn<MimmoObject*, StitchGeometry>(this, &mimmo::StitchGeometry::setAddGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));

    built = (built && createPortOut<MimmoObject*, StitchGeometry>(this, &mimmo::StitchGeometry::getGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortOut<std::vector<MimmoObject*>, StitchGeometry>(this, &mimmo::StitchGeometry::getOriginalGeometries, PortType::M_VECGEOM, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortOut<std::unordered_map<long,std::pair<int, long> >, StitchGeometry>(this, &mimmo::StitchGeometry::getCellDivisionMap, PortType::M_MAPDCELL, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::LONGPAIRINTLONG));
    built = (built && createPortOut<std::unordered_map<long,std::pair<int, long> >, StitchGeometry>(this, &mimmo::StitchGeometry::getVertDivisionMap, PortType::M_MAPDVERT, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::LONGPAIRINTLONG));
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
 * Get list of input original geometries.
 * \return Pointers to original geometries.
 */
std::vector<MimmoObject *> 
StitchGeometry::getOriginalGeometries(){
    std::vector<MimmoObject * > res(m_extgeo.size());
    for(auto & val: m_extgeo){
        res[val.second] = val.first;
    }
    return res;
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
 * Get the complete map of the stitched object cell ids vs original geometries cell ids.
 * The map reports as key the actual id of the cell in the stitched object and a pair argument
 * with an integer id identifying the original part which belongs and the cell-id it had in the original
 * geometry structure.
 * Please note the part identifier numbering follows the order of the original geometry list m_extgeo.
 * \return map with cells IDs of the stitched geometry and pair of local and global cells IDs of input geometries
 */
std::unordered_map<long, std::pair<int,long> >
StitchGeometry::getCellDivisionMap(){
    return m_mapCellDivision;
}

/*!
 * Get the complete map of the stiched object vertex ids vs original geometries vertex ids.
 * The map reports as key the actual id of the vertex in the stitched object and a pair argument
 * with an integer id identifying the original part which belongs and the vertex-id it had in the original
 * geometry structure. Please note the part identifier numbering follows the order of the original geometry list m_extgeo.
 * \return map with vertices IDs of the stitched geometry and pair of local and global vertices IDs of input geometries
 */
std::unordered_map<long, std::pair<int,long> >
StitchGeometry::getVertDivisionMap(){
    return m_mapVertDivision;
}

/*!
 * Add an external geometry to be stitched.Topology of the geometry must be coeherent
 * with topology of the class;
 * \param[in] geo  Pointer to MimmoObject
 */
void
StitchGeometry::setAddGeometry(MimmoObject* geo){
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
StitchGeometry::setGeometry(std::vector<MimmoObject *> external){

    for(auto & obj: external){
        setAddGeometry(obj);
    }
};

/*!It sets if the BvTree of stitched geometry has to be built during execution.
 * \param[in] build If true the BvTree is built in execution and stored in
 * the related MimmoObject member.
 */
void
StitchGeometry::setBuildBvTree(bool build){
    m_buildBvTree = build;
}

/*!It sets if the KdTree of all the patch geometries has to be built during execution.
 * \param[in] build If true the KdTree is built in execution and stored in
 * the related MimmoObject member.
 */
void
StitchGeometry::setBuildKdTree(bool build){
    m_buildKdTree = build;
}

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
    m_mapCellDivision.clear();
    m_mapVertDivision.clear();
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

    //clean maps
    m_mapCellDivision.clear();
    m_mapVertDivision.clear();

    long cV = 0, cC = 0;

    //start filling your stitched object.
    //divion maps will be filled coherently

    {
        //optional vars;
        long vId, cId;
        short PID;
        bitpit::ElementInfo::Type eltype;

        for(auto &obj : m_extgeo){

            //start extracting/reversing vertices of the current obj
            std::unordered_map<long,long> mapVloc;

            for(auto & vv : obj.first->getVertices()){
                vId = vv.getId();
                dum->addVertex(obj.first->getVertexCoords(vId), cV);
                //update map;
                m_mapVertDivision[cV] = std::make_pair(obj.second, vId);
                mapVloc[vId] = cV;
                cV++;
            }

            //start extracting/reversing cells of the current obj
            for(auto & cc : obj.first->getCells()){
                cId = cc.getId();
                PID = cc.getPID();
                eltype = cc.getType();
                //get the local connectivity and update with new vertex numbering;
                livector1D conn = obj.first->getCellConnectivity(cId);
                for(auto && ee: conn)    ee = mapVloc[ee];

                dum->addConnectedCell(conn, eltype, PID, cC);
                //update map;
                m_mapCellDivision[cC] = std::make_pair(obj.second, cId);
                cC++;
            }
        }
    }//scope for optional vars;

    if(m_buildBvTree)    dum->buildBvTree();
    if(m_buildKdTree)    dum->buildKdTree();

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

    if(slotXML.hasOption("Priority")){
        input = slotXML.get("Priority");
        int value =0;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss>>value;
        }
        setPriority(value);
    };

    if(slotXML.hasOption("BvTree")){
        input = slotXML.get("BvTree");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss >> value;
        }
        setBuildBvTree(value);
    };

    if(slotXML.hasOption("KdTree")){
        input = slotXML.get("KdTree");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss >> value;
        }
        setBuildKdTree(value);
    };

    if(slotXML.hasOption("PlotInExecution")){
        std::string input = slotXML.get("PlotInExecution");
        input = bitpit::utils::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setPlotInExecution(value);
    }

    if(slotXML.hasOption("OutputPlot")){
        std::string input = slotXML.get("OutputPlot");
        input = bitpit::utils::trim(input);
        std::string temp = ".";
        if(!input.empty())    setOutputPlot(input);
        else                  setOutputPlot(temp);
    }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void StitchGeometry::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    slotXML.set("ClassName", m_name);
    slotXML.set("Priority", std::to_string(getPriority()));
    slotXML.set("Topology", m_topo);

    std::string output;

    if(m_buildBvTree){
        output = std::to_string(m_buildBvTree);
        slotXML.set("BvTree", output);
    }

    if(m_buildKdTree){
        output = std::to_string(m_buildKdTree);
        slotXML.set("KdTree", output);
    }

    if(isPlotInExecution()){
        slotXML.set("PlotInExecution", std::to_string(1));
    }

    if(m_outputPlot != "."){
        slotXML.set("OutputPlot", m_outputPlot);
    }
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
