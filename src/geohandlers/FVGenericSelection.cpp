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

namespace mimmo {

/*!
 * Basic Constructor.
 * Parameter topo get topology of the target MimmoFvMesh where performing extraction:
 * - 1 for volume bulk and surface boundary
 * - 2 for surface bulk and 3DCurve boundary

 * No other values are allowed
 * \param[in] topo topology of the target MimmoFvMesh.
 */
FVGenericSelection::FVGenericSelection(int topo){
    m_type = FVSelectionType::UNDEFINED;
    m_topo = std::min(2, std::max(1, topo)); /*default to volume bulk-surface boundary geometry*/
    m_dual = false; /*default to exact selection*/
    m_bndgeometry.reset();
};

/*!
 * Basic Destructor
 */
FVGenericSelection::~FVGenericSelection(){};

/*!
 * Copy Constructor, any already calculated selection is not copied.
 */
FVGenericSelection::FVGenericSelection(const FVGenericSelection & other):BaseManipulation(other){
    m_type = other.m_type;
    m_topo = other.m_topo;
    m_dual = other.m_dual;
    m_bndgeometry = other.m_bndgeometry;
};

/*!
 * Copy operator, any already calculated selection is not copied.
 */
FVGenericSelection & FVGenericSelection::operator=(const FVGenericSelection & other){
    *(static_cast<BaseManipulation *>(this)) = *(static_cast<const BaseManipulation * >(&other));
    m_type = other.m_type;
    m_topo = other.m_topo;
    m_dual = other.m_dual;
    m_bndgeometry = other.m_bndgeometry;
    /*m_subpatch is not copied and it is obtained in execution*/
    return *this;
};

/*!
 * Swap function. Assign data of this class to another of the same type and vice-versa.
 * Resulting patches of selection are not swapped, ever.
 * \param[in] x FVGenericSelection object
 */
void FVGenericSelection::swap(FVGenericSelection & x) noexcept
{
    std::swap(m_type, x.m_type);
    std::swap(m_topo, x.m_topo);
    std::swap(m_dual, x.m_dual);
    std::swap(m_bndgeometry, x.m_bndgeometry);
    BaseManipulation::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
FVGenericSelection::buildPorts(){

    bool built = true;

    built = (built && createPortIn<mimmo::MimmoSharedPointer<MimmoObject>, FVGenericSelection>(this, &FVGenericSelection::setGeometry, M_GEOM, true));
    built = (built && createPortIn<mimmo::MimmoSharedPointer<MimmoObject>, FVGenericSelection>(this, &FVGenericSelection::setBoundaryGeometry, M_GEOM2, true));
    built = (built && createPortIn<bool, FVGenericSelection>(this, &FVGenericSelection::setDual,M_VALUEB));

    built = (built && createPortOut<mimmo::MimmoSharedPointer<MimmoObject>, FVGenericSelection>(this, &FVGenericSelection::getVolumePatch, M_GEOM));
    built = (built && createPortOut<mimmo::MimmoSharedPointer<MimmoObject>, FVGenericSelection>(this, &FVGenericSelection::getBoundaryPatch, M_GEOM2));
    built = (built && createPortOut<mimmo::MimmoSharedPointer<MimmoObject>, FVGenericSelection>(this, &FVGenericSelection::getInternalBoundaryPatch, M_GEOM3));

    m_arePortsBuilt = built;
};


/*!
 * Return type of method for selection.
 * See FVSelectionType enum for available methods.
 * \return type of selection method
 */
FVSelectionType
FVGenericSelection::whichMethod(){
    return    m_type;
};

/*!
 * Return pointer-by-copy to bulk sub-patch extracted by the class
 * \return pointer to Bulk Mesh MimmoObject extracted sub-patch
 */
mimmo::MimmoSharedPointer<MimmoObject>
FVGenericSelection::getVolumePatch(){
    return    m_volpatch;
};

/*!
 * Return pointer-by-copy to boundary sub-patch extracted by the class
 * \return pointer to Boundary MimmoObject extracted sub-patch
 */
mimmo::MimmoSharedPointer<MimmoObject>
FVGenericSelection::getBoundaryPatch(){
    return    m_bndpatch;
};


/*!
 * Return pointer-by-copy to internal boundary sub-patch extracted by the class
 * \return pointer to Boundary MimmoObject extracted sub-patch
 */
mimmo::MimmoSharedPointer<MimmoObject>
FVGenericSelection::getInternalBoundaryPatch(){
    return    m_intbndpatch;
};
/*!
 * Return pointer-by-copy to bulk sub-patch extracted by the class
 * \return pointer to Bulk Mesh MimmoObject extracted sub-patch
 */
const mimmo::MimmoSharedPointer<MimmoObject>
FVGenericSelection::getVolumePatch() const{
    return    m_volpatch;
};

/*!
 * Return pointer-by-copy to boundary sub-patch extracted by the class
 * \return pointer to Boundary MimmoObject extracted sub-patch
 */
const mimmo::MimmoSharedPointer<MimmoObject>
FVGenericSelection::getBoundaryPatch() const{
    return    m_bndpatch;
};

/*!
 * Return pointer-by-copy to internal boundary sub-patch extracted by the class
 * \return pointer to Boundary MimmoObject extracted sub-patch
 */
const mimmo::MimmoSharedPointer<MimmoObject>
FVGenericSelection::getInternalBoundaryPatch() const{
    return    m_intbndpatch;
};


/*!
 * Set link to target bulk geometry for your selection.
 * Reimplementation of mimmo::BaseManipulation::setGeometry();
 *  \param[in] target Pointer to MimmoObject with bulk target geometry.
 */
void
FVGenericSelection::setGeometry( mimmo::MimmoSharedPointer<MimmoObject> target){
    if(target == nullptr)  return;
    int type = target->getType();
    if(m_topo == 1 && type != 2) return;
    if(m_topo == 2 && type != 1) return;
    m_geometry = target;
};

/*!
 * Set link to target boundary geometry for your selection.
 *  \param[in] target Pointer to MimmoObject with boundary target geometry.
 */
void
FVGenericSelection::setBoundaryGeometry( mimmo::MimmoSharedPointer<MimmoObject> target){
    if(target == nullptr)  return;
    int type = target->getType();
    if(m_topo == 1 && type != 1) return;
    if(m_topo == 2 && type != 4) return;
    m_bndgeometry = target;
};


/*!
 * Set your class behavior selecting a portion of a target geoemetry.
 * Given a initial set up, gets the dual result (its negative) of current selection.
 * For instance, in a class extracting geometry inside the volume of an
 * elemental shape, gets all other parts of it not included in the shape.
 * \param[in] flag true-Activate / false-Deactivate "dual" feature .
 */
void
FVGenericSelection::setDual(bool flag ){
    m_dual = flag;
}

/*!
 * Return actual status of "dual" feature of the class. See setDual method.
 * \return  true/false for "dual" feature activated or not
 */
bool
FVGenericSelection::isDual(){
    return m_dual;
};


/*!
 * Execute your object. A selection is extracted and trasferred in
 * two indipendent MimmoObject structures pointed by m_volpatch and m_bndpatch members
 */
void
FVGenericSelection::execute(){
    if(m_geometry == nullptr || m_bndgeometry == nullptr) {
        (*m_log)<<m_name + " : nullptr pointer to target bulk/boundary geometry found or both"<<std::endl;
        throw std::runtime_error (m_name + " : nullptr pointer to target bulk/boundary geometry found or both");
        return;
    }
    if(m_geometry->isEmpty() || m_bndgeometry->isEmpty() ){
        (*m_log)<<m_name + " : empty bulk/boundary geometry linked"<<std::endl;
    };

    if(!checkCoherenceBulkBoundary()){
        (*m_log)<<m_name + " : id-vertex uncoherent bulk/boundary geometry linked"<<std::endl;
    }

    m_volpatch.reset();
    m_bndpatch.reset();
    m_intbndpatch.reset();

    livector1D extractedVol;
    livector1D extractedBnd;

    extractSelection(extractedVol, extractedBnd);

    if(extractedVol.empty() || extractedBnd.empty()) {
        (*m_log)<<m_name + " : empty selection for bulk or boundary performed. check block set-up"<<std::endl;
    }

    /*Create subpatches.*/
    int topovol=2, topobnd=1;
    if(m_topo == 2){
        topovol = 1;
        topobnd = 4;
    }

    mimmo::MimmoSharedPointer<MimmoObject> tempVol(new MimmoObject(topovol));
    mimmo::MimmoSharedPointer<MimmoObject> tempBnd(new MimmoObject(topobnd));
    mimmo::MimmoSharedPointer<MimmoObject> tempInternalBnd(new MimmoObject(topobnd));

    //VOLUME PART
    {
        livector1D vertExtracted = m_geometry->getVertexFromCellList(extractedVol);
        for(const auto & idV : vertExtracted){
            tempVol->addVertex(m_geometry->getVertexCoords(idV), idV);
        }

        int rank;
        for(const auto & idCell : extractedVol){
            bitpit::Cell & cell = m_geometry->getPatch()->getCell(idCell);
            rank = -1;
#if MIMMO_ENABLE_MPI
            rank = m_geometry->getPatch()->getCellRank(idCell);
#endif
            tempVol->addCell(cell, idCell, rank);
        }
    }


    //BOUNDARY PART
    // get vertices of real boundary mesh.
    livector1D vertBndExtracted;
    long bndMaxCellId = -1;
    if(!extractedBnd.empty()){
        bndMaxCellId = *(std::max_element(extractedBnd.begin(), extractedBnd.end()));
    }
#if MIMMO_ENABLE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &bndMaxCellId, 1, MPI_LONG, MPI_MAX, m_communicator);
#endif
    ++bndMaxCellId;
    {
        vertBndExtracted = m_bndgeometry->getVertexFromCellList(extractedBnd);

        //TO AVOID INCONSISTENCY WITH VOLUME MESH, CHECK all verts extracted are
        // present into the volume mesh.
        {
            livector1D epurated_Verts;
            epurated_Verts.reserve(vertBndExtracted.size());
            bitpit::PiercedVector<bitpit::Vertex> & volv = tempVol->getVertices();
            for(long id: vertBndExtracted){
                if(volv.exists(id)) epurated_Verts.push_back(id);
            }
            std::swap(vertBndExtracted, epurated_Verts);

            // last step, recover cell strictly in the pool of this new list of epurated vertices.
            extractedBnd = m_bndgeometry->getCellFromVertexList(vertBndExtracted, true);
        }

        for(const auto & idV : vertBndExtracted){
            tempBnd->addVertex(m_bndgeometry->getVertexCoords(idV), idV);
        }

        int rank;
        for(const auto & idCell : extractedBnd){
            bitpit::Cell & cell = m_bndgeometry->getPatch()->getCell(idCell);
            rank = -1;
#if MIMMO_ENABLE_MPI
            rank = m_bndgeometry->getPatch()->getCellRank(idCell);
#endif
            tempBnd->addCell(cell, idCell, rank);
        }
    }

     // get the internal boundary, those interfaces of the volume mesh belonging to the border.
    {
        std::unordered_map<long, std::set<int>> bndFacesMap = tempVol->extractBoundaryFaceCellID(true);
        long faceCount(0);
        for(auto & tuple : bndFacesMap)     faceCount += long(tuple.second.size());

        livector1D candidateVerts = tempVol->extractBoundaryVertexID(bndFacesMap);
        tempVol->buildInterfaces();
        //fill the wrapped internal boundary mesh (can coprehend also the real boundary)
#if MIMMO_ENABLE_MPI
        {
            std::vector<long> offsetFaceCount(m_nprocs, 0);
            offsetFaceCount[m_rank] = faceCount;
            MPI_Allreduce(MPI_IN_PLACE, &offsetFaceCount, m_nprocs, MPI_LONG, MPI_MAX, m_communicator);
            for(int i=0; i<m_rank; ++i){
                bndMaxCellId += offsetFaceCount[i];
            }
        }
#endif

        for(long idV : candidateVerts){
            tempInternalBnd->addVertex(tempVol->getVertexCoords(idV), idV);
        }

        int rank;
        long PID;
        for(const auto & tuplemap : bndFacesMap){
            bitpit::Cell & cell = tempVol->getPatch()->getCell(tuplemap.first);
            const long * interfacesList = cell.getInterfaces();
            for(int locface : tuplemap.second){
                if(interfacesList[locface] < 0) continue;
                bitpit::Interface & face = tempVol->getPatch()->getInterface(interfacesList[locface]);
                rank = -1;
#if MIMMO_ENABLE_MPI
                rank = tempVol->getPatch()->getCellRank(face.getOwner());
#endif
                tempInternalBnd->addConnectedCell(std::vector<long>(face.getConnect(), face.getConnect()+face.getConnectSize()),
                                                  face.getType(), tempVol->getPatch()->getCell(face.getOwner()).getPID(), bndMaxCellId, rank);
                ++bndMaxCellId;
            }
        }

        tempVol->resetInterfaces();
        tempVol->resetAdjacencies();
        //remove cells shared with the real boundary, from this "internal" pot.
        livector1D sharedCells = tempInternalBnd->getCellFromVertexList(vertBndExtracted, true);

        //delete cells and orphans.
        tempInternalBnd->getPatch()->deleteCells(sharedCells);
        tempInternalBnd->buildAdjacencies();
        if(tempInternalBnd->getPatch()->countOrphanVertices() > 0){
            tempInternalBnd->getPatch()->deleteOrphanVertices();
        }
        tempInternalBnd->resetAdjacencies();
    }

    m_volpatch = tempVol;
    m_bndpatch = tempBnd;
    m_intbndpatch = tempInternalBnd;

    m_volpatch->cleanGeometry();
    m_bndpatch->cleanGeometry();
    m_intbndpatch->cleanGeometry();
    m_volpatch->update();
    m_bndpatch->update();
    m_intbndpatch->update();

};

/*!
 * Plot optional result of the class in execution. It plots the selected patch
 * as standard vtk unstructured grid.
 */
void
FVGenericSelection::plotOptionalResults(){

	std::string originalname = m_name;
	 if(getVolumePatch()){
		 m_name = originalname + "_Volume_Patch";
		 write(getVolumePatch());
	 }
    if(getBoundaryPatch()){
    	m_name = originalname + "_Boundary_Patch";
		 write(getBoundaryPatch());
    }
    if(getInternalBoundaryPatch()){
     m_name = originalname + "_InternalBoundary_Patch";
      write(getInternalBoundaryPatch());
    }
	m_name = originalname;

}

/*!
 * Check if boundary and bulk geometry as coherence in vertex indexing.
 */
bool
FVGenericSelection::checkCoherenceBulkBoundary(){

    livector1D vertBnd = m_bndgeometry->getVertices().getIds();
    bitpit::PiercedVector<bitpit::Vertex> & vertBulk = m_geometry->getVertices();

    for (const auto &id: vertBnd){
       if(!vertBulk.exists(id)) return false;
    }

    return true;
}

}
