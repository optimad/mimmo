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
    m_topo = std::min(2, std::max(1, topo)); /*default to volume bulk-surface boundary geometry*/
    m_dual = false; /*default to exact selection*/
    m_bndgeometry= nullptr;
    m_selectEngine = nullptr;
};

/*!
 * Basic Destructor
 */
FVGenericSelection::~FVGenericSelection(){};

/*!
 * Copy Constructor, any already calculated selection is not copied.
 */
FVGenericSelection::FVGenericSelection(const FVGenericSelection & other):BaseManipulation(other){
    m_topo = other.m_topo;
    m_dual = other.m_dual;
    m_bndgeometry = other.m_bndgeometry;
    m_selectEngine = other.m_selectEngine;
};

/*!
 * Copy operator, any already calculated selection is not copied.
 */
FVGenericSelection & FVGenericSelection::operator=(FVGenericSelection other){
    swap(other);
    return *this;
};

/*!
 * Swap function. Assign data of this class to another of the same type and vice-versa.
 * Resulting patches of selection are not swapped, ever.
 * \param[in] x FVGenericSelection object
 */
void FVGenericSelection::swap(FVGenericSelection & x) noexcept
{
    std::swap(m_topo, x.m_topo);
    std::swap(m_dual, x.m_dual);
    std::swap(m_bndgeometry, x.m_bndgeometry);
    std::swap(m_selectEngine, x.m_selectEngine);
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
 * Set GenericSelection block that performs selection coherently on bulk+boundary compound.
 * \param[in] selectBlock valid pointer to a selection block.
 */
void
FVGenericSelection::setSelection(MimmoSharedPointer<GenericSelection> selectBlock){
    if(selectBlock == nullptr) return;
    m_selectEngine = selectBlock;
}


/*!
 * Return pointer to bulk sub-patch extracted by the class
 * \return pointer to Bulk Mesh MimmoObject extracted sub-patch
 */
mimmo::MimmoSharedPointer<MimmoObject>
FVGenericSelection::getVolumePatch(){
    return    m_volpatch;
};

/*!
 * Return pointer to boundary sub-patch extracted by the class
 * \return pointer to Boundary MimmoObject extracted sub-patch
 */
mimmo::MimmoSharedPointer<MimmoObject>
FVGenericSelection::getBoundaryPatch(){
    return    m_bndpatch;
};


/*!
 * Return pointer to internal boundary sub-patch extracted by the class
 * \return pointer to Boundary MimmoObject extracted sub-patch
 */
mimmo::MimmoSharedPointer<MimmoObject>
FVGenericSelection::getInternalBoundaryPatch(){
    return    m_intbndpatch;
};

/*!
 * Return pointer to bulk sub-patch extracted by the class
 * \return pointer to Bulk Mesh MimmoObject extracted sub-patch
 */
const mimmo::MimmoSharedPointer<MimmoObject>
FVGenericSelection::getVolumePatch() const{
    return    m_volpatch;
};

/*!
 * Return pointer to boundary sub-patch extracted by the class
 * \return pointer to Boundary MimmoObject extracted sub-patch
 */
const mimmo::MimmoSharedPointer<MimmoObject>
FVGenericSelection::getBoundaryPatch() const{
    return    m_bndpatch;
};

/*!
 * Return pointer to internal boundary sub-patch extracted by the class
 * \return pointer to Boundary MimmoObject extracted sub-patch
 */
const mimmo::MimmoSharedPointer<MimmoObject>
FVGenericSelection::getInternalBoundaryPatch() const{
    return    m_intbndpatch;
};



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
        throw std::runtime_error (m_name + " : nullptr pointer to one or both bulk/boundary targets found");
        return;
    }

    if(m_selectEngine == nullptr) {
        throw std::runtime_error (m_name + " : no valid selection engine found");
        return;
    }

    if(!checkCoherenceBulkBoundary()){
        (*m_log)<<"Warning in "+m_name + " : id-vertex uncoherent bulk/boundary geometry linked"<<std::endl;
    }

    //selecting bulk
    m_selectEngine->setGeometry(m_geometry);
    m_selectEngine->execute();
    m_volpatch = m_selectEngine->getPatch();

    //selecting boundary
    m_selectEngine->setGeometry(m_bndgeometry);
    m_selectEngine->execute();
    m_bndpatch = m_selectEngine->getPatch();

    //clean up the boundary. THis is necessary since the selection on bulk and boundary
    // may lead to a situation in which some vertex nodes in the boundary patch selected are not
    // in the pot of the bulk patch selected
    cleanUpBoundaryPatch();

    //create the internal patch;
    m_intbndpatch = createInternalBoundaryPatch();



};

/*! Once bulk and boundary are filled with selections, ensure that boundary patch
    has only vertices shared with the bulk borders. Clean up the mesh eventually.*/
void
FVGenericSelection::cleanUpBoundaryPatch(){
    //check boundary subpatch vertices and retain only those inside bulk subpatch.
    //local cleaning on rank;
    livector1D boundaryClearedVerts;
    bitpit::PiercedVector<bitpit::Vertex> & bulkV = m_volpatch->getVertices();
    boundaryClearedVerts.reserve(m_bndpatch->getNVertices());
    long id;
    for(const bitpit::Vertex & vert : m_bndpatch->getVertices()){
        id = vert.getId();
        if(bulkV.exists(id))    boundaryClearedVerts.push_back(id);
    }

    std::unordered_set<long> survivedCells;
    {
        livector1D temp = m_bndpatch->getCellFromVertexList(boundaryClearedVerts, true); // true strictly cell defined by this set.
        survivedCells.insert(temp.begin(), temp.end());
    }
    livector1D toDeleteCells;
    livector1D bndC = m_bndpatch->getCells().getIds();

    toDeleteCells.reserve(bndC.size() - survivedCells.size());
    for(long id: bndC){
        if(survivedCells.count(id) > 0) continue;
        toDeleteCells.push_back(id);
    }
    toDeleteCells.shrink_to_fit();

    bool checkAllClean = toDeleteCells.empty();
#if MIMMO_ENABLE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &checkAllClean,1, MPI_C_BOOL, MPI_LAND, m_communicator);
#endif

    if(!checkAllClean){
        //unsync all mimmoobject tracking structures
        m_bndpatch->setUnsyncAll();
        //delete the selected cells
        m_bndpatch->getPatch()->deleteCells(toDeleteCells);
#if MIMMO_ENABLE_MPI
        //delete possible isolated ghosts
        m_bndpatch->deleteOrphanGhostCells();
#endif
        //check for orphans vertices and delete it
        if(m_bndpatch->getPatch()->countOrphanVertices() > 0){
            m_bndpatch->getPatch()->deleteOrphanVertices();
        }
        m_bndpatch->getPatch()->squeezeVertices();
        m_bndpatch->getPatch()->squeezeCells();

        m_bndpatch->update();
    }
}


/*! Once bulk and boundary patch selections are filled and clean, retrieve the internal
    boundary patch.
    \return the internal boundary patch
*/
MimmoSharedPointer<MimmoObject>
FVGenericSelection::createInternalBoundaryPatch(){

    //get the vertices retained by the real boundary
    livector1D realBoundaryVertices = m_bndpatch->getVertices().getIds();

    //check interfaces of volpatch
    bool resetVolInterfaces = (m_volpatch->getInterfacesSyncStatus() == SyncStatus::NONE);
    m_volpatch->updateInterfaces();

    MimmoSharedPointer<MimmoObject> m_patch = m_volpatch->extractBoundaryMesh();

    if(resetVolInterfaces){
        m_volpatch->destroyInterfaces();
    }

    //remove cells shared with the real boundary, from this "internal" pot.
    livector1D sharedCells = m_patch->getCellFromVertexList(realBoundaryVertices, true); //strict=true means cells defined by this vertices list.

    //unsync all mimmoobject tracking structures
    m_patch->setUnsyncAll();
    //delete cells and orphans.
    m_patch->getPatch()->deleteCells(sharedCells);
#if MIMMO_ENABLE_MPI
    m_patch->deleteOrphanGhostCells();
#endif
    if(m_patch->getPatch()->countOrphanVertices() > 0){
        m_patch->getPatch()->deleteOrphanVertices();
    }
    m_patch->getPatch()->squeezeVertices();
    m_patch->getPatch()->squeezeCells();

    m_patch->update();

    return m_patch;
}
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

    m_geometry->update();
    m_bndgeometry->update();

    //check if boundary and bulk are globally empty or not.
    std::array<bool,2> emptyCompound = {{m_geometry->isEmpty(), m_bndgeometry->isEmpty()}};
#if MIMMO_ENABLE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &emptyCompound,2, MPI_C_BOOL, MPI_LAND, m_communicator);
#endif
    //if they differ return false.
    if(emptyCompound[0] != emptyCompound[1])    return false;

    bool check = true;
    //analyze vertices
    livector1D boundaryIds = m_bndgeometry->getVerticesIds(true); //only internals nodes.
    std::unordered_set<long> bulkset;
    {
        livector1D bulkIds = m_geometry->extractBoundaryVertexID(false); //only rank internals
        bulkset.insert(bulkIds.begin(), bulkIds.end());
    }
    if(!boundaryIds.empty()){
        for(long id: boundaryIds){
            check = check && (bulkset.count(id) > 0);
        }
    }

#if MIMMO_ENABLE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &check, 1, MPI_C_BOOL, MPI_LAND, m_communicator);
#endif


    return check;
}

}
