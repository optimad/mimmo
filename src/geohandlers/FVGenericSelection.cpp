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
#include <cstddef>
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
    m_bndgeometry = NULL;
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

    built = (built && createPortIn<MimmoObject *, FVGenericSelection>(this, &FVGenericSelection::setGeometry, M_GEOM, true));
    built = (built && createPortIn<MimmoObject *, FVGenericSelection>(this, &FVGenericSelection::setBoundaryGeometry, M_GEOM2, true));
    built = (built && createPortIn<bool, FVGenericSelection>(this, &FVGenericSelection::setDual,M_VALUEB));

    built = (built && createPortOut<MimmoObject *, FVGenericSelection>(this, &FVGenericSelection::getVolumePatch, M_GEOM));
    built = (built && createPortOut<MimmoObject *, FVGenericSelection>(this, &FVGenericSelection::getBoundaryPatch, M_GEOM2));
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
MimmoObject*
FVGenericSelection::getVolumePatch(){
    return    m_volpatch.get();
};

/*!
 * Return pointer-by-copy to boundary sub-patch extracted by the class
 * \return pointer to Boundary MimmoObject extracted sub-patch
 */
MimmoObject*
FVGenericSelection::getBoundaryPatch(){
    return    m_bndpatch.get();
};

/*!
 * Return pointer-by-copy to bulk sub-patch extracted by the class
 * \return pointer to Bulk Mesh MimmoObject extracted sub-patch
 */
const MimmoObject*
FVGenericSelection::getVolumePatch() const{
    return    m_volpatch.get();
};

/*!
 * Return pointer-by-copy to boundary sub-patch extracted by the class
 * \return pointer to Boundary MimmoObject extracted sub-patch
 */
const MimmoObject*
FVGenericSelection::getBoundaryPatch() const{
    return    m_bndpatch.get();
};

/*!
 * Set link to target bulk geometry for your selection.
 * Reimplementation of mimmo::BaseManipulation::setGeometry();
 *  \param[in] target Pointer to MimmoObject with bulk target geometry.
 */
void
FVGenericSelection::setGeometry( MimmoObject * target){
    if(target == NULL)  return;
    if(target->isEmpty()) return;
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
FVGenericSelection::setBoundaryGeometry( MimmoObject * target){
    if(target == NULL)  return;
    if(target->isEmpty()) return;
    int type = target->getType();
    if(m_topo == 1 && type != 1) return;
    if(m_topo == 2 && type != 4) return;
    m_bndgeometry = target;
};


/*!
 * Set your class behavior selecting a portion of a target geoemetry.
 * Given a initial set up, gets the dual
 * result (its negative) of current selection.
 * For instance, in a class extracting geometry inside the volume of an
 * elemental shape, gets all other parts of it not included in the shape.
 * For a class extracting portions of geometry
 * matching an external tessellation, gets all the other parts not matching it.
 * For a class extracting portion of
 * geometry identified by a PID list, get all other parts not marked by such list.
 * \param[in] flag Active/Inactive "dual" feature true/false.
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
    if(m_geometry == NULL || m_bndgeometry == NULL) {
        throw std::runtime_error (m_name + " : NULL pointer to target bulk/boundary geometry found or both");
    }
    if(m_geometry->isEmpty() || m_bndgeometry->isEmpty() ){
        throw std::runtime_error (m_name + " : empty bulk/boundary geometry linked");
    };

    if(!checkCoherenceBulkBoundary()){
        throw std::runtime_error (m_name + " : connected bulk and boundary geometries are not coherent");
    }
    
    m_volpatch.reset(nullptr);
    m_bndpatch.reset(nullptr);
    
    livector1D extractedVol;
    short maxPID = 0;
    std::unordered_map<short,livector1D> extractedBndVert;
    {
        std::unordered_map<short,livector1D> extractedBnd;
        
        extractSelection(extractedVol, extractedBnd);
        
        for (const auto & touple: extractedBnd){
            extractedBndVert[touple.first] = m_bndgeometry->getVertexFromCellList(touple.second);
            maxPID = std::max(touple.first, maxPID);
        }
    }

    if(extractedVol.empty()) {
        throw std::runtime_error (m_name + " : empty selection performed. check block set-up");
    }

    /*Create subpatches.*/
    int topovol=2, topobnd=1;
    if(m_topo == 2){
        topovol = 1;
        topobnd = 4;
    }
    
    std::unique_ptr<MimmoObject> tempVol(new MimmoObject(topovol));
    std::unique_ptr<MimmoObject> tempBnd(new MimmoObject(topobnd));

    {
        livector1D TT;
        bitpit::ElementType eltype;
        short PID;

        livector1D vertExtracted = m_geometry->getVertexFromCellList(extractedVol);
        for(const auto & idV : vertExtracted){
            tempVol->addVertex(m_geometry->getVertexCoords(idV), idV);
        }

        for(const auto & idCell : extractedVol){
            bitpit::Cell & cell = m_geometry->getPatch()->getCell(idCell);
            eltype = cell.getType();
            PID = (short)cell.getPID();
            TT = m_geometry->getCellConnectivity(idCell);
            tempVol->addConnectedCell(TT,eltype, PID, idCell);
            TT.clear();
        }
        
        tempVol->buildInterfaces();
    }
    
    //now create boundary mesh of selection.
    {
        std::set<long> boundaryInterfaces;
        std::set<long> boundaryVertices;
        bitpit::ConstProxyVector<long> vcount;
        for(const auto & interf : tempVol->getInterfaces()){
            if(interf.isBorder()){
                boundaryInterfaces.insert(interf.getId());
                vcount = interf.getVertexIds();
                boundaryVertices.insert(vcount.begin(), vcount.end());
            }
        }
        
        //fill new boundary
        bitpit::PiercedVector<bitpit::Interface> & bulkInterf = tempVol->getInterfaces();
        
        for(const auto & idV: boundaryVertices){
            tempBnd->addVertex(tempVol->getVertexCoords(idV),idV);
        }
        
        short PID = maxPID+1;
        bitpit::ElementType eltype;
        livector1D conn;
        for(const auto & idI: boundaryInterfaces){
            
            int size = bulkInterf[idI].getConnectSize();
            conn.resize(size);
            long * cc = bulkInterf[idI].getConnect();
            for(int i=0; i<size; ++i){
                conn[i] = cc[i];
            }
            eltype = bulkInterf[idI].getType();
            tempBnd->addConnectedCell(conn, eltype, PID, idI);
        }
    }
    
    //post-processing boundary (assign existent PID if any).
    for(auto & touple : extractedBndVert){
        livector1D cellids = tempBnd->getCellFromVertexList(touple.second, true);
        for(const long & cid: cellids){
            tempBnd->setPIDCell(cid, touple.first);
        }
    }
    
    m_volpatch = std::move(tempVol);
    m_bndpatch = std::move(tempBnd);
};

/*!
 * Plot optional result of the class in execution. It plots the selected patch
 * as standard vtk unstructured grid.
 */
void
FVGenericSelection::plotOptionalResults(){
    if(getVolumePatch() == NULL || getBoundaryPatch() == NULL) return;
    if(getVolumePatch()->isEmpty() || getBoundaryPatch()->isEmpty()) return;

    std::string dir = m_outputPlot;
    std::string namevol = m_name + "_Volume_Patch."+ std::to_string(getId());
    std::string namebnd = m_name + "_Boundary_Patch."+ std::to_string(getId());
    
    getVolumePatch()->getPatch()->write(namevol);
    getBoundaryPatch()->getPatch()->write(namebnd);
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
