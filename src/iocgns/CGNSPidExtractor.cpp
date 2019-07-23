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

#include "CGNSPidExtractor.hpp"
#include <unordered_map>

using namespace std;
using namespace bitpit;
namespace mimmo{

/*!
 * Default constructor of CGNSPidExtractor.
 */
CGNSPidExtractor::CGNSPidExtractor(){
    m_name         = "mimmo.CGNSPidExtractor";
    m_force = false;
    m_targetpid.insert(0);
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
CGNSPidExtractor::CGNSPidExtractor(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.CGNSPidExtractor";
    m_force = false;
    m_targetpid.insert(0);

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.CGNSPidExtractor"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Default destructor of CGNSPidExtractor.
 */
CGNSPidExtractor::~CGNSPidExtractor(){};

/*!
 * Copy constructor of CGNSPidExtractor.
 * No result patch is copied.
 */
CGNSPidExtractor::CGNSPidExtractor(const CGNSPidExtractor & other):BaseManipulation(other){
    m_force = other.m_force;
    m_targetpid = other.m_targetpid;
};

/*!
 * Assignment operator. No result patch is copied.
 */
CGNSPidExtractor & CGNSPidExtractor::operator=(CGNSPidExtractor other){
    swap(other);
    return *this;
}

/*!
 * Swap function. No result patch member is swapped.
 * \param[in] x object to be swapped.
 */
void CGNSPidExtractor::swap(CGNSPidExtractor & x) noexcept
{
    std::swap(m_force, x.m_force);
    std::swap(m_targetpid, x.m_targetpid);
    BaseManipulation::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
CGNSPidExtractor::buildPorts(){
    bool built = true;
    built = (built && createPortIn<MimmoObject*, CGNSPidExtractor>(this, &CGNSPidExtractor::setGeometry, M_GEOM, true));

    built = (built && createPortOut<MimmoObject*, CGNSPidExtractor>(this, &CGNSPidExtractor::getPatch, M_GEOM));
    m_arePortsBuilt = built;
};



/*!
 * Return the extracted surface by PID. If not valid PID is found, return an exact copy of the original target geometry.
 * \return pointer to extracted surface
 */
MimmoObject*
CGNSPidExtractor::getPatch(){
    return m_patch.get();
}

/*!
 * Return current active PIDs for selection in the class.
 * \return list of target PIDs
 */
std::set<long>
CGNSPidExtractor::whatPIDActive(){
    return m_targetpid;
}

/*!
 * Return if the class is set to reduce (in execution)
 * the resulting PID-selected patch in a homogeneous straight triangulation,
 * if eventually not.
 * \return the surface has to be triangulated
 */
bool
CGNSPidExtractor::isForcedToTriangulate(){
    return m_force;
}

/*!
 * Set current active PID for selectioning in the class. If PID is not
 * supported by the linked target geoemetry, return a stand-alone copy of the target geometry itself in execution.
 * \param[in]    val target PID to insert in the extraction list
 */
void
CGNSPidExtractor::addPID( long val){
    if(val > 0) m_targetpid.insert(val);
}

/*!
 * Set list of current active PID for selectioning in the class. If no PID in list is not
 * supported by the linked target geoemetry, return a stand-alone copy of the target geometry itself in execution.
 * \param[in]   vval list of target PIDs to use during the extraction
 */
void
CGNSPidExtractor::setPID( std::vector<long> vval){
    for(const auto & val : vval) addPID(val);
}

/*!
 * Set the class to reduce (in execution) the resulting PID-selected patch in a homogeneous straight triangulation,
 * if eventually not.
 * \param[in]   flag if true, try to triangulate homogeneously the selected patch.
 */
void
CGNSPidExtractor::setForcedToTriangulate( bool flag){
    m_force = flag;
}

/*!
 * Set current geometry to an external surface boundary mesh, from a mimmo CGNS reader.
 * \param[in]   geo pointer to input geometry.
 */
void
CGNSPidExtractor::setGeometry(MimmoObject * geo){
    if(geo == NULL)    return;
    if(geo->getType() != 1)    return;
    BaseManipulation::setGeometry(geo);
}


/*!
 * Execution command.
 * Extract patch with target PIDs.
 */
void
CGNSPidExtractor::execute(){

    MimmoObject * mother = getGeometry();

    if(!mother){
        (*m_log)<<m_name + " : NULL pointer to linked geometry found"<<std::endl;
        throw std::runtime_error("NULL pointer to linked geometry found");
        return;
    }

    if(mother->isEmpty() ){
#if MIMMO_ENABLE_MPI
        (*m_log)<<"WARNING "<<m_name << " on rank "<<m_rank<<" : empty linked geometry found"<<std::endl;
#else
        (*m_log)<<"WARNING "<<m_name<<" : empty linked geometry found"<<std::endl;
#endif
    }

    livector1D extracted;

    for(const auto & val: m_targetpid){
        livector1D temp= getGeometry()->extractPIDCells(val);
        extracted.insert(extracted.end(), temp.begin(), temp.end());
    }
    std::unique_ptr<MimmoObject> patchTemp(new MimmoObject(1));
    livector1D vertExtracted = mother->getVertexFromCellList(extracted);

    bitpit::PiercedVector<bitpit::Vertex> & motherverts = mother->getVertices();
    for(const auto & val: vertExtracted){
        patchTemp->addVertex(motherverts.at(val).getCoords(), val);
    }

    int rank;
    bitpit::PiercedVector<bitpit::Cell> & mothercells = mother->getCells();
    for(const auto &val: extracted){
        bitpit::Cell & cell = mothercells.at(val);
        rank = -1;
#if MIMMO_ENABLE_MPI
        rank = mother->getPatch()->getCellRank(val);
#endif
        patchTemp->addCell(cell, val, rank);
    }


    auto originalmap = mother->getPIDTypeListWNames();
    auto currentPIDmap = patchTemp->getPIDTypeList();


    for(const auto & val: currentPIDmap){
        patchTemp->setPIDName(val, originalmap[val]);
    }

    //check if patch is forced to be triangulated in case of polygons (nVertices > 4) more
    // nodes are added to the mesh.
    // For Now this slot works only in SERIAL mode.
    if(m_force){
#if MIMMO_ENABLE_MPI == 0

        long maxID, newID, newVertID;

        const auto orderedCellID = patchTemp->getCells().getIds(true);
        maxID = orderedCellID[(int)orderedCellID.size()-1];
        newID = maxID+1;
        {
            const auto orderedVertID = patchTemp->getVertices().getIds(true);
            newVertID = orderedVertID[(int)orderedCellID.size()-1] +1;
        }

        bitpit::ElementType eletype;
        bitpit::ElementType eletri = bitpit::ElementType::TRIANGLE;
        livector1D connTriangle(3);



        for(const auto &idcell : orderedCellID){

            livector1D conn = patchTemp->getCellConnectivity(idcell);
            eletype = patchTemp->getPatch()->getCell(idcell).getType();
            long pid = patchTemp->getPatch()->getCell(idcell).getPID();

            switch (eletype){
                case bitpit::ElementType::QUAD:
                case bitpit::ElementType::PIXEL:
                {
                    patchTemp->getPatch()->deleteCell(idcell);
                    for(std::size_t i=0; i<2; ++i){
                        connTriangle[0] = conn[0];
                        connTriangle[1] = conn[i+1];
                        connTriangle[2] = conn[i+2];
                        patchTemp->addConnectedCell(connTriangle, eletri, pid, newID);
                        ++newID;
                    }
                }
                    break;
                case bitpit::ElementType::POLYGON:
                {
                    std::size_t startIndex = 1;
                    std::size_t nnewTri = conn.size() - startIndex;
                    //calculate barycenter and add it as new vertex
                    darray3E barycenter = patchTemp->getPatch()->evalCellCentroid(idcell);
                    patchTemp->addVertex(barycenter, newVertID);
                    //delete current polygon
                    patchTemp->getPatch()->deleteCell(idcell);
                    //insert new triangles from polygon subdivision
                    for(std::size_t i=0; i<nnewTri; ++i){
                        connTriangle[0] = newVertID;
                        connTriangle[1] = conn[ startIndex + std::size_t( i % nnewTri) ];
                        connTriangle[2] = conn[ startIndex + std::size_t( (i+1) % nnewTri ) ];
                        patchTemp->addConnectedCell(connTriangle, eletri, pid, newID);
                        ++newID;
                    }
                    //increment label of vertices
                    ++newVertID;

                }
                    break;
                case bitpit::ElementType::TRIANGLE:
                        //do nothing
                    break;
                default:
                    throw std::runtime_error("unrecognized cell type in 3D surface mesh of CGNSPidExtractor");
                    break;
            }
        }
#else
    //TODO provide implementation to deal with insertion/deletion of vertices and cells in parallel
    (*m_log)<< "WARNING " <<m_name <<" : forced triangulation is not available yet in MPI compilation."<<std::endl;
#endif
    }

    m_patch = std::move(patchTemp);

#if MIMMO_ENABLE_MPI
    //delete orphan ghosts
    m_patch->buildAdjacencies();
    m_patch->deleteOrphanGhostCells();
    if(m_patch->getPatch()->countOrphanVertices() > 0){
        m_patch->getPatch()->deleteOrphanVertices();
    }
    m_patch->setPartitioned();
#endif

    if(getGeometry()->isInfoSync()) m_patch->buildPatchInfo();

#if MIMMO_ENABLE_MPI
    if(getGeometry()->arePointGhostExchangeInfoSync()) m_patch->updatePointGhostExchangeInfo();

#endif
};


/*!
 * Plot optional result of the class in execution, that is the selected patch
 * as standard vtk unstructured grid.
 */
void CGNSPidExtractor::plotOptionalResults(){
    if(!getPatch())    return;
    std::string dir = m_outputPlot+"/";
    std::string name = m_name + "_" + std::to_string(getId()) +  "_Patch";
    getPatch()->getPatch()->getVTK().setDirectory(dir);
    getPatch()->getPatch()->getVTK().setName(name);
    getPatch()->getPatch()->write();
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
CGNSPidExtractor::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    BaseManipulation::absorbSectionXML(slotXML, name);

    std::vector<long> temp;
    if(slotXML.hasOption("nPID")){
        input = slotXML.get("nPID");
        int value = 0;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
            value = std::max(0, value);
        }
        temp.resize(value, -1);
    };

    if(slotXML.hasOption("PID") && !temp.empty()){
        input = slotXML.get("PID");
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            for(auto & val: temp){
                ss >> val;
            }
        }
        setPID(temp);
    };

    if(slotXML.hasOption("ForcedToTriangulate")){
        input = slotXML.get("ForcedToTriangulate");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setForcedToTriangulate(value);
    };

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
CGNSPidExtractor::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("nPID", std::to_string(int(m_targetpid.size())));

    std::stringstream output;
    auto it = m_targetpid.begin();
    auto itE = m_targetpid.end();
    while(it != itE)    {
        output<<*it<<" ";
        ++it;
    }
    slotXML.set("PID", output.str());
    slotXML.set("ForcedToTriangulate", std::to_string(int(m_force)));


};


}
