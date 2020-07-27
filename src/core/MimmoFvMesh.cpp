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

#include "MimmoFvMesh.hpp"

namespace mimmo{

/*!
 * Default constructor
 */
InfoBoundaryPatch::InfoBoundaryPatch(){
    name = "";
    type = "NONE";
}

/*!
 * Default destructor
 */
InfoBoundaryPatch::~InfoBoundaryPatch(){}


/*!
 * Default Constructor
 */
MimmoFvMesh::MimmoFvMesh(){
    m_name = "";
    m_bulk.reset();
    m_boundary.reset();
}

/*!
 * Default Destructor
 */
MimmoFvMesh::~MimmoFvMesh(){}

/*!
 * Custom constructor. Pass bulk and boundary of the mesh.
 * The ownership of the argument will be transferred to the internal
 * members of the class.
 * Please note, The class admits only the following combinations as bulk/boundary pair:
 *
 *  - bulk MimmoObject of type 2-Volume and boundary MimmoObject of type 1-Surface
 *  - bulk MimmoObject of type 1-Surface and boundary MimmoObject of type 4-3DCurve.
 *
 *\param[in] bulk unique pointer to the bulk MimmoObject
 *\param[in] boundary unique pointer to the boundary MimmoObject
 */
MimmoFvMesh::MimmoFvMesh(MimmoSharedPointer<MimmoObject> bulk, MimmoSharedPointer<MimmoObject> boundary){
    //check couples first
    if(!bulk && !boundary){
        throw std::runtime_error("Error in MimmoFvMesh Constructor: nullptr linked as bulk mesh or boundary mesh");
    }
    std::string key = std::to_string(bulk->getType()) + std::to_string(boundary->getType());

    if(key !="21" && key !="14"){
        throw std::runtime_error("Error in MimmoFvMesh Constructor: not valid types in pair bulk-boundary meshes");
    }

    m_bulk = bulk;
    m_boundary = boundary;
}

/*!
 * Copy Constructor. Softly copy all internal meshes of argument class
 * as external meshes for the current one.
   \param[in] other object ot be copied
 */
MimmoFvMesh::MimmoFvMesh(const MimmoFvMesh & other):BaseManipulation(other){
    m_bulk = other.m_bulk;
    m_boundary = other.m_boundary;
}


/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void MimmoFvMesh::swap(MimmoFvMesh & x) noexcept{
    std::swap(m_bulk, x.m_bulk);
    std::swap(m_boundary, x.m_boundary);
    std::swap(m_infoBoundary, x.m_infoBoundary);
    BaseManipulation::swap(x);
}

/*!
 * Build ports of the class
 */
void MimmoFvMesh::buildPorts(){

    bool built = true;
    //Don't push mandatory ports here, this class is thought more as an intermediate interface,
    //than a proper mimmo executable block. See IOOFoam families, for example.
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, MimmoFvMesh>(this, &MimmoFvMesh::setGeometry, M_GEOM ));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, MimmoFvMesh>(this, &MimmoFvMesh::setBoundaryGeometry, M_GEOM2));
    // creating output ports
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, MimmoFvMesh>(this, &MimmoFvMesh::getGeometry, M_GEOM));
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, MimmoFvMesh>(this, &MimmoFvMesh::getBoundaryGeometry, M_GEOM2));

    m_arePortsBuilt = built;
};

/*!
 * Set bulk geometry. The class softly set it, i.e. simply copy pointer
   to the candidate bulk geometry.
 * MimmoObject type admitted are 2-Volume or 1-Surface.
 * \param[in] bulk geometry to be linked
 */
void MimmoFvMesh::setGeometry(MimmoSharedPointer<MimmoObject> bulk){
    if (bulk == nullptr) return;
    if (bulk->getType() != 2 && bulk->getType() != 1) return;
    m_bulk = bulk;
}

/*!
 * Set boundary geometry. The class softly set it, i.e. simply copy pointer
   to the candidate boundary geometry.
 * MimmoObject type admitted are 1-Surface or 4-3dCurve.
 * \param[in] boundary geometry to be linked
 */
void MimmoFvMesh::setBoundaryGeometry(MimmoSharedPointer<MimmoObject> boundary){
    if (boundary == nullptr) return;
    if (boundary->getType() != 1 && boundary->getType() != 4) return;
    m_boundary = boundary;
}

/*!
 * \return current bulk geometry
 */
MimmoSharedPointer<MimmoObject>
MimmoFvMesh::getGeometry(){
    return m_bulk;
}

/*!
 * \return current boundary geometry
 */
MimmoSharedPointer<MimmoObject>
MimmoFvMesh::getBoundaryGeometry(){
    return m_boundary;
}

/*!
 * Add optional information on a boundary patch identified by a Part IDentifier
 * /TODO: method does nothing at the moment.
 * \param[in] PID patch identifier
 * \param[in] info struct holding option info on the patch.
 */
void
MimmoFvMesh::addInfoBoundaryPatch(const long & PID, const InfoBoundaryPatch & info){
    BITPIT_UNUSED(PID);
    BITPIT_UNUSED(info);
}


/*!
 * \return optional information on on a boundary patch identified by a Part IDentifier
 * /TODO: method does nothing at the moment.
 * \param[in] PID patch identifier
 */
InfoBoundaryPatch
MimmoFvMesh::getInfoBoundaryPatch(const long & PID){
    BITPIT_UNUSED(PID);
    return InfoBoundaryPatch();
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void MimmoFvMesh::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
    BITPIT_UNUSED(name);
    BaseManipulation::absorbSectionXML(slotXML, name);
}

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void MimmoFvMesh::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    BaseManipulation::flushSectionXML(slotXML, name);
}


/*!
 * Create boundary mesh starting from a valid bulk present in your
   current bulk slot(linked or owned).
 */
void MimmoFvMesh::createBoundaryMesh(){

    MimmoSharedPointer<MimmoObject> bulk = getGeometry();
    if (!bulk->areInterfacesBuilt()) bulk->buildInterfaces();

    std::set<long> boundaryInterfaces;
    std::set<long> boundaryVertices;
    bitpit::ConstProxyVector<long> vcount;
    for(const auto & interf : bulk->getInterfaces()){
        if(interf.isBorder()){
            boundaryInterfaces.insert(interf.getId());
            vcount = interf.getVertexIds();
            boundaryVertices.insert(vcount.begin(), vcount.end());
        }
    }

    //fill new boundary
    int type = int(bulk->getType() == 2) + 4*int(bulk->getType() == 1);
    MimmoSharedPointer<MimmoObject> temp(new MimmoObject(type));

    bitpit::PiercedVector<bitpit::Interface> & bulkInterf = bulk->getInterfaces();

    for(const auto & idV: boundaryVertices){
        temp->addVertex(bulk->getVertexCoords(idV),idV);
    }


    long PID = 0;
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
        PID = long(bulkInterf[idI].getPID());

#if MIMMO_ENABLE_MPI
        int rank = bulk->getPatch()->getCellRank(bulkInterf[idI].getOwner());
        temp->addConnectedCell(conn, eltype, PID, idI, rank);
#else
        temp->addConnectedCell(conn, eltype, PID, idI);
#endif

    }

    m_boundary = temp;
}


/*!
 * Check coeherence between bulk and boundary mesh, i.e.:
 *
 * - bulk non-null and non-empty
 * - if boundary non null and non empty, it checks out:
 *   - bulk/boundary pair must be of MimmoObject types 2-Volume bulk and 1-Surface boundary or 1-Surface bulk and 4-3Dcurve boundary
 *   - ids tight matching between border interfaces on the bulk and cells of the boundary
 * - otherwise:
 *   - create internally a boundary patch, without pid segmentation informations, and return true.
 *
 * \return false if one of the check failed.
 */
bool  MimmoFvMesh::checkMeshCoherence(){

    MimmoSharedPointer<MimmoObject> bulk = getGeometry();
    if(bulk == nullptr)  return false;

    MimmoSharedPointer<MimmoObject> boundary = getBoundaryGeometry();

    bool checkBoundaries = true;
    if(boundary == nullptr)  checkBoundaries = false;

    if(checkBoundaries){
        if(!bulk->areInterfacesBuilt())  return false;
        std::string key = std::to_string(bulk->getType()) + std::to_string(boundary->getType());
        if(key !="21" && key !="14")    return false;

        std::set<long> idBorderInterf;
        for(const auto & interf : bulk->getInterfaces()){
            if(interf.isBorder())   idBorderInterf.insert(interf.getId());
        };
        std::vector<long> idInterf;
        idInterf.insert(idInterf.end(), idBorderInterf.begin(), idBorderInterf.end());
        std::vector<long> idCell   = boundary->getCells().getIds(true);

        return idInterf == idCell;
    }else{
        createBoundaryMesh();
        return true;
    }
}

}
