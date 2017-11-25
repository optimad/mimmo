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

#include "MimmoObject.hpp"
#include "Operators.hpp"
#include <set>

using namespace std;
using namespace bitpit;

namespace mimmo{

/*!
 * MimmoSurfUnstructured default constructor
 */
MimmoSurfUnstructured::MimmoSurfUnstructured():bitpit::SurfUnstructured(){}

/*!
 * MimmoSurfUnstructured custom constructor
 * \param[in] patch_dim dimensionality of elements (2D-triangles/quads/polygons, 1D-lines) 
 */
MimmoSurfUnstructured::MimmoSurfUnstructured(int patch_dim):
                       bitpit::SurfUnstructured(int(patch_dim),int(3)){}

/*!
 * MimmoSurfUnstructured custom constructor
 * \param[in] id custom identification label of the mesh
 * \param[in] patch_dim dimensionality of elements (2D-triangles/quads/polygons, 1D-lines) 
 */
MimmoSurfUnstructured::MimmoSurfUnstructured(const int & id, int patch_dim):
                       bitpit::SurfUnstructured(id, int(patch_dim), int(3)){}

/*!
 * MimmoSurfUnstructured custom constructor
 * \param[in] stream input stream where reading from
 */
MimmoSurfUnstructured::MimmoSurfUnstructured(std::istream & stream):
                       bitpit::SurfUnstructured(stream){}


/*!
 * Basic Destructor
 */
MimmoSurfUnstructured::~MimmoSurfUnstructured(){}

/*!
 * Cloning a MimmoSurfUnstructured in an independent object of base type 
 * bitpit::PatchKernel
 */
std::unique_ptr<bitpit::PatchKernel>
MimmoSurfUnstructured::clone() const{
    return std::unique_ptr<MimmoSurfUnstructured>(new MimmoSurfUnstructured(*this));
}

/*!
 * MimmoVolUnstructured default constructor
 * \param[in] dimension dimensionality of elements (3D - tetrahedra/hexahedra ..., 2D-triangles/quads/polygons)
 */
MimmoVolUnstructured::MimmoVolUnstructured(const int& dimension):
bitpit::VolUnstructured(dimension){}

/*!
 * MimmoVolUnstructured custom constructor
 * \param[in] id custom identification label of the mesh
 * \param[in] dimension dimensionality of elements (3D - tetrahedra/hexahedra ..., 2D-triangles/quads/polygons)
 */
MimmoVolUnstructured::MimmoVolUnstructured(const int & id, const int & dimension):
bitpit::VolUnstructured(id, dimension){}

/*!
 * Basic Destructor
 */
MimmoVolUnstructured::~MimmoVolUnstructured(){}

/*!
 * Cloning a MimmoVolUnstructured in an independent object of base type 
 * bitpit::PatchKernel
 */
std::unique_ptr<bitpit::PatchKernel>
MimmoVolUnstructured::clone() const{
    return std::unique_ptr<MimmoVolUnstructured>(new MimmoVolUnstructured(*this));
}

/*!
 * MimmoPointCloud basic constructor
 */
MimmoPointCloud::MimmoPointCloud():
bitpit::SurfUnstructured(int(2),int(3)){}

/*!
 * MimmoPointCloud custom constructor
 * \param[in] id custom identification label of the mesh
 */
MimmoPointCloud::MimmoPointCloud(const int & id):
bitpit::SurfUnstructured(id, int(2), int(3)){}

/*!
 * Basic Destructor
 */
MimmoPointCloud::~MimmoPointCloud(){}

/*!
 * Cloning a MimmoPointCloud in an independent object of base type 
 * bitpit::PatchKernel
 */
std::unique_ptr<bitpit::PatchKernel>
MimmoPointCloud::clone() const{
    return std::unique_ptr<MimmoPointCloud>(new MimmoPointCloud(*this));
}



/*!
 * Default constructor of MimmoObject.
 * It requires a int flag identifying the type of mesh meant to be created:
 *  - surface unstructured mesh = 1
 *  - volume unstructured mesh  = 2
 *  - 3D Cloud Point            = 3
 *  - 3D tessellated Curve      = 4
 *
 * \param[in] type type of mesh
 */
MimmoObject::MimmoObject(int type){

    m_log = &bitpit::log::cout(MIMMO_LOG_FILE);
    m_type = max(type,1);
    if (m_type > 4){
        (*m_log)<<"Error MimmoObject: unrecognized data structure type in class construction. Switch to DEFAULT 1-Surface"<<std::endl;
        throw std::runtime_error ("MimmoObject : unrecognized mesh type in class construction");
    }
    switch(m_type){
        case 1:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoSurfUnstructured(2)));
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new bitpit::SurfaceSkdTree(dynamic_cast<SurfaceKernel*>(m_patch.get()))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 2:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoVolUnstructured(3)));
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new VolumeSkdTree(dynamic_cast<VolumeKernel*>(m_patch.get()))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 3:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoPointCloud()));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 4:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoSurfUnstructured(1)));
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new bitpit::SurfaceSkdTree(dynamic_cast<SurfaceKernel*>(m_patch.get()))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        default:
            //never been reached
            break;
    }

    m_internalPatch = true;
    m_extpatch = NULL;

    m_skdTreeSupported = (m_type != 3);
    m_skdTreeSync = false;
    m_kdTreeSync = false;
    m_AdjBuilt = false;
    m_IntBuilt = false;
}

/*!
 * Custom constructor of MimmoObject.
 * This constructor builds a generic mesh from given vertex list and its related
 * local connectivity. Mesh-element details are desumed by local cell-vertices connectivity size.
 * Connectivities supported can be of any kind. Mixed elements are allowed.
 * Cell element types supported are all those of bitpit::ElementType list. In particular:
 *  - any 2D cell for surface meshes 
 *  - any 3D cell for volume meshes
 *  - VERTEX only for 3D point clouds
 *  - LINE cells only for 3D curves
 *
 * Cloud points are always supported (no connectivity field must be provided)
 * Type of meshes supported are described in the default MimmoObject constructor documentation.
 * Please note, despite type argument, if null connectivity structure is provided, MimmoObject will be built 
 * as a standard cloud point of vertices provided by vertex.
 * Connectivity for each standard cell is simply defined as a list of vertex indices which compose it
 * (indices are meant as positions in the provided vertex list argument).
 * Polygonal and Polyhedral cells requires a special convention to define their connectivity:
 * - polygons: require on top of the list the total number of vertices which compose it,
 *             followed by the indices in the vertex list.(ex. polygon with 5 vertices: 5 i1 i2 i3 i4 i5)
 * - polyhedra: require on top the total number of faces, followed by the number of vertices which composes the local face + the vertex 
 *   indices which compose the faces, and so on for all the faces which compose the polyhedron (ex. polyhedron with 3 faces, 
 *   first face is a triangle, second is a quad etc...  3 3 i1 i2 i3 4 i2 i3 i4 i5 ...)
 * \param[in] type type of meshes.
 * \param[in] vertex Coordinates of geometry vertices.
 * \param[in] connectivity pointer to mesh connectivity list (optional).
 */
MimmoObject::MimmoObject(int type, dvecarr3E & vertex, livector2D * connectivity){

    m_log = &bitpit::log::cout(MIMMO_LOG_FILE);
    
    if(connectivity == NULL){
        m_type = 3;
    }else{
        m_type = max(type,1);
        if (m_type > 4){
            (*m_log)<<"Error MimmoObject: unrecognized data structure type in class construction. Switch to DEFAULT 1-Surface"<<std::endl;
            throw std::runtime_error ("MimmoObject : unrecognized mesh type in class construction");
        }
    }

    switch(m_type){
        case 1:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoSurfUnstructured(2)));
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new bitpit::SurfaceSkdTree(dynamic_cast<SurfaceKernel*>(m_patch.get()))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 2:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoVolUnstructured(3)));
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new VolumeSkdTree(dynamic_cast<VolumeKernel*>(m_patch.get()))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 3:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoPointCloud()));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 4:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoSurfUnstructured(1)));
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new bitpit::SurfaceSkdTree(dynamic_cast<SurfaceKernel*>(m_patch.get()))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        default:
            //never been reached
            break;
    }
    
    m_internalPatch = true;
    m_extpatch = NULL;

    m_skdTreeSupported = (m_type != 3);
    m_skdTreeSync = false;
    m_kdTreeSync = false;
    m_AdjBuilt = false;
    m_IntBuilt = false;

    bitpit::ElementType eltype;
    std::size_t sizeVert = vertex.size();

    m_patch->reserveVertices(sizeVert);



    for(const auto & vv : vertex)   addVertex(vv);

    m_log->setPriority(bitpit::log::Priority::DEBUG);

    if(m_type != 3){
        livector1D temp;
        std::size_t sizeCell = connectivity->size();
        m_patch->reserveCells(sizeCell);
        for(auto const & cc : *connectivity){
            eltype = desumeElement(cc);
            if(eltype != bitpit::ElementType::UNDEFINED){
                addConnectedCell(cc, eltype);
            }else{
                (*m_log)<<"warning: in MimmoObject custom constructor. Undefined cell type detected and skipped."<<std::endl; 
            }
        }

        m_pidsType.insert(0);

    }else{
        (*m_log)<<"Not supported connectivity found for MimmoObject"<<std::endl;
        (*m_log)<<"Proceeding as Point Cloud geometry"<<std::endl;
    }	
    m_log->setPriority(bitpit::log::Priority::NORMAL);
};

/*!
 * Custom constructor of MimmoObject.
 * This constructor builds a geometry data structure soft linking to an external bitpit::PatchKernel object, that is
 * it does not own the geometry data structure, but simple access it, while it is instantiated elsewhere.
 * Search Trees will be referred to this linked geometry.
 * If a null geometry patch is linked, a standard MimmoObject is built.
 * The mesh type needs to be specified (see default constructor MimmoObject(int type) doc).
 * \param[in] type type of mesh
 * \param[in] geometry pointer to a geometry of class PatchKernel to be linked.
 */
MimmoObject::MimmoObject(int type, bitpit::PatchKernel* geometry){

    m_log = &bitpit::log::cout(MIMMO_LOG_FILE);
    m_type = max(type,1);
    if (m_type > 4 || geometry == NULL){
        (*m_log)<<"Error MimmoObject: unrecognized data structure type or NULL argument in class construction."<<std::endl;
        throw std::runtime_error ("MimmoObject : unrecognized mesh type or NULL argument in class construction");
    }

    m_internalPatch = false;
    m_extpatch = geometry;
    switch(m_type){
        case 1:
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new bitpit::SurfaceSkdTree(dynamic_cast<SurfaceKernel*>(geometry))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 2:
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new VolumeSkdTree(dynamic_cast<VolumeKernel*>(geometry))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 3:
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 4:
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new bitpit::SurfaceSkdTree(dynamic_cast<SurfaceKernel*>(geometry))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        default:
            //never been reached
            break;
    }

    m_skdTreeSupported = (m_type != 3);
    m_skdTreeSync = false;
    m_kdTreeSync = false;
    m_AdjBuilt = false;
    m_IntBuilt = false;
    
    for(const auto &cell : getPatch()->getCells()){
        auto PID = cell.getPID();
        m_pidsType.insert((short)PID);
    }
}

/*!
 * Default destructor
 */
MimmoObject::~MimmoObject(){

};

/*!
 * Copy constructor of MimmoObject. 
 * Internal allocated PatchKernel is copied from the argument object as a soft link 
 * (copied by its pointer only, geometry data structure is treated as external in the new copied class). 
 * Search trees are instantiated, but their containts (if any) are not copied by default.
 */
MimmoObject::MimmoObject(const MimmoObject & other){

    m_log = &bitpit::log::cout(MIMMO_LOG_FILE);
    
    m_extpatch = other.m_extpatch;
    if(other.m_internalPatch){
        m_extpatch      = other.m_patch.get();
    }
    m_internalPatch = false;
    
    m_type              = other.m_type;
    m_pidsType          = other.m_pidsType;
    m_skdTreeSupported  = other.m_skdTreeSupported;
    m_AdjBuilt          = other.m_AdjBuilt;
    m_IntBuilt          = other.m_IntBuilt;
    
    m_skdTreeSync    = false;
    m_kdTreeSync    = false;
    
    //instantiate empty trees:
    switch(m_type){
        case 1:
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new bitpit::SurfaceSkdTree(dynamic_cast<SurfaceKernel*>(m_extpatch))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 2:
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new VolumeSkdTree(dynamic_cast<VolumeKernel*>(m_extpatch))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 3:
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 4:
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new bitpit::SurfaceSkdTree(dynamic_cast<SurfaceKernel*>(m_extpatch))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        default:
            //never been reached
            break;
    }
};

/*!
* Assignement operator of MimmoObject.
* Please be careful when using it, because internal bitpit::PatchKernel structure
* is only copied by their pointer (soft copy). If you destroy the original MimmoObject
* the copied class will have an internal patch pointing to nullptr.
* Search trees are instantiated, but their containts (if any) are not copied by default.
* \param[in] other reference to another MimmoObject
*/
MimmoObject & MimmoObject::operator=(MimmoObject other){
    swap(other);
    return *this;
};

/*!
 * Swap function. Swap data between current object and one of the same type. 
 * \param[in] x object to be swapped.
 */
void MimmoObject::swap(MimmoObject & x) noexcept
{
    std::swap(m_patch, x.m_patch);
    std::swap(m_extpatch, x.m_extpatch);
    std::swap(m_internalPatch, x.m_internalPatch);
    std::swap(m_type, x.m_type);
    std::swap(m_pidsType, x.m_pidsType);
    std::swap(m_skdTreeSupported, x.m_skdTreeSupported);
    std::swap(m_AdjBuilt, x.m_AdjBuilt);
    std::swap(m_IntBuilt, x.m_IntBuilt);
    std::swap(m_skdTree, x.m_skdTree);
    std::swap(m_kdTree, x.m_kdTree);
    std::swap(m_skdTreeSync, x.m_skdTreeSync);
    std::swap(m_kdTreeSync, x.m_kdTreeSync);
}


/*!
 * Is the object empty?
 * \return True/False if geometry data structure is empty.
 */
bool
MimmoObject::isEmpty(){
    bool check = getNVertex() > 0;
    if(m_type != 3) check = check && (getNCells() > 0);
    return !check;
};

/*!
 * Is the object empty? const overloading.
 * \return True/False if geometry data structure is empty.
 */
bool
MimmoObject::isEmpty() const{
    bool check = getNVertex() > 0;
    if(m_type != 3) check = check && (getNCells() > 0);
    return !check;
};


/*! 
 * Is the skdTree (former bvTree) ordering supported w/ the current geometry?
 * \return True for connectivity-based meshes, false por point clouds
 */
bool
MimmoObject::isBvTreeSupported(){
    return isSkdTreeSupported();
};

/*! 
 * Is the skdTree (former bvTree) ordering supported w/ the current geometry?
 * \return True for connectivity-based meshes, false por point clouds
 */
bool
MimmoObject::isSkdTreeSupported(){
    return m_skdTreeSupported;
};

/*!
 * Return the type of mesh currently hold by the class 
 * (for type of mesh allowed see MimmoObject(int type) documentation).
 * \return integer flag for mesh type
 */
int
MimmoObject::getType(){
    return m_type;
};

/*!
 * Return the total number of vertices within the data structure.
 * \return number of mesh vertices 
 */
long
MimmoObject::getNVertex() const {
    const auto p = getPatch();
    return p->getVertexCount();
};

/*!
 * Return the total number of cells within the data structure.
 * \return number of mesh cells.
 */
long
MimmoObject::getNCells() const {
    const auto p = getPatch();
    return p->getCellCount();
};

/*!
 * Return the compact list of vertices hold by the class
 * \return coordinates of mesh vertices 
 * \param[out] mapDataInv pointer to inverse of Map of vertex ids, for aligning external vertex data to bitpit::Patch ordering
 */
dvecarr3E
MimmoObject::getVertexCoords(liimap* mapDataInv){
    dvecarr3E result(getNVertex());
    int  i = 0;

    auto pvert = getVertices();

    if (mapDataInv != NULL){
        for (auto const & vertex : pvert){
            result[i] = vertex.getCoords();
            (*mapDataInv)[vertex.getId()] = i;
            ++i;
        }
    }
    else{
        for (auto const & vertex : pvert){
            result[i] = vertex.getCoords();
            ++i;
        }
    }
    return result;
};

/*!
 * Gets the vertex coordinates marked with a unique label i (bitpit::PatchKernel id);
 * If i is not found, return a [1.e18,1.e18,1.e18] default value.
 * \param[in] i bitpit::PatchKernel id of the vertex in the mesh.
 * \return Coordinates of the i-th vertex of geometry mesh.
 */
darray3E
MimmoObject::getVertexCoords(long i){
    if(!(getVertices().exists(i)))	return darray3E({{1.e18,1.e18,1.e18}});
    return 	getPatch()->getVertexCoords(i);
};

/*!
 * Return reference to the PiercedVector structure of vertex hold 
 * by bitpit::PatchKernel class member
 * \return  Pierced Vector structure of vertices
 */
bitpit::PiercedVector<bitpit::Vertex> &
MimmoObject::getVertices(){
    return getPatch()->getVertices();
}

/*!
 * Return const reference to the PiercedVector structure of vertex hold 
 * by bitpit::PatchKernel class member
 * \return  const Pierced Vector structure of vertices
 */
const bitpit::PiercedVector<bitpit::Vertex> &
MimmoObject::getVertices() const {
    const auto p = getPatch();
    return p->getVertices();
}

/*!
 * Get the geometry compact connectivity.
 * "Compact" means that in the connectivity matrix the vertex indexing is referred to 
 * the local, compact and sequential numbering of vertices, as it 
 * gets them from the internal method getVertexCoords().
 * Special connectivity is returned for polygons and polyhedra support. See getCellConnectivity
 * doxy for further info.
 * \return local connectivity list 
 * \param[in] mapDataInv inverse of Map of vertex ids actually set, for aligning external vertex data to bitpit::Patch ordering
 */
livector2D
MimmoObject::getCompactConnectivity(liimap & mapDataInv){

    livector2D connecti(getNCells());
    int np, counter =0;

    for(auto const & cell : getCells()){
        np = cell.getConnectSize();
        const long * conn_ = cell.getConnect();
        connecti[counter].resize(np);
        bitpit::ElementType eltype = cell.getType();
        
        if(eltype == bitpit::ElementType::POLYGON){
            connecti[counter][0] = conn_[0];
            for (int i=1; i<np; ++i){
                connecti[counter][i] = mapDataInv[conn_[i]];
            }
        }else if(eltype == bitpit::ElementType::POLYHEDRON){
            connecti[counter][0] = conn_[0];
            for(int nF = 0; nF < conn_[0]-1; ++nF){
                int facePos = cell.getFaceStreamPosition(nF);
                int beginVertexPos = facePos + 1;
                int endVertexPos   = facePos + 1 + conn_[facePos];
                connecti[counter][facePos] = conn_[facePos]; 
                for (int i=beginVertexPos; i<endVertexPos; ++i){
                    connecti[counter][i] = mapDataInv[conn_[i]];
                }
            }
        }else{
            for (int i=0; i<np; ++i){
                connecti[counter][i] = mapDataInv[conn_[i]];
            }
        }
        ++counter;
    }
    return connecti;
}

/*!
 * It gets the connectivity of the cells of the linked geometry. Index of vertices in
 * connectivity matrix are returned according to bitpit::PatchKernel unique indexing.
 * Special connectivity is returned for polygons and polyhedra support. See getCellConnectivity
 * doxy for further info.
 * \return cell-vertices connectivity 
 */
livector2D
MimmoObject::getConnectivity(){

    livector2D connecti(getNCells());
    int np, counter =0;

    for(auto const & cell : getCells()){
        np = cell.getConnectSize();
        const long * conn_ = cell.getConnect();
        connecti[counter].resize(np);
        for (int i=0; i<np; ++i){
            connecti[counter][i] = conn_[i];
        }
        ++counter;
    }

    return connecti;
};

/*!
 * It gets the connectivity of a cell, with vertex id's in bitpit::PatchKernel unique indexing.
 * Connectivity of polygons is returned as (nV,V1,V2,V3,V4,...) where nV is the number of vertices
 * defining the polygon and V1,V2,V3... are indices of vertices.
 * Connectivity of polyhedrons  is returned as (nF,nF1V, V1, V2, V3,..., nF2V, V2,V4,V5,V6,...,....)
 * where nF is the total number of faces of polyhedros, nF1V is the number of vertices composing the face 1,
 * followed by the indices of vertices which draws it. Face 2,3,..,n are defined in the same way.
 * Any other standard cell element is uniquely defined by its list of vertex indices.
 * if i cell does not exists, return an empty list.
 * \param[in] i bitpit::PatchKernel ID of the cell.
 * \return i-th cell-vertices connectivity, in bitpit::PatchKernel vertex indexing.
 */
livector1D
MimmoObject::getCellConnectivity(long i){
    if (!(getCells().exists(i)))    return livector1D(0);

    bitpit::Cell & cell = getPatch()->getCell(i); 
    int np = cell.getConnectSize();
    const long * conn_ = cell.getConnect();
    livector1D connecti(np);

    for (int j=0; j<np; j++){
        connecti[j] = conn_[j];
    }
    return connecti;
};

/*!
 * Return reference to the PiercedVector structure of cells hold 
 * by bitpit::PatchKernel class member
 * \return  Pierced Vector structure of cells
 */
bitpit::PiercedVector<bitpit::Cell> &
MimmoObject::getCells(){
    return getPatch()->getCells();
}

/*!
 * Return const reference to the PiercedVector structure of cells hold 
 * by bitpit::PatchKernel class member
 * \return  const Pierced Vector structure of cells
 */
const bitpit::PiercedVector<bitpit::Cell> &
MimmoObject::getCells()  const{
    const auto p = getPatch();
    return p->getCells();
}

/*!
 * Return reference to the PiercedVector structure of interfaces hold 
 * by bitpit::PatchKernel class member
 * \return  Pierced Vector structure of interfacess
 */
bitpit::PiercedVector<bitpit::Interface> &
MimmoObject::getInterfaces(){
    return getPatch()->getInterfaces();
}

/*!
 * Return const reference to the PiercedVector structure of interfaces hold 
 * by bitpit::PatchKernel class member
 * \return  const Pierced Vector structure of interfaces
 */
const bitpit::PiercedVector<bitpit::Interface> &
MimmoObject::getInterfaces()  const{
    const auto p = getPatch();
    return p->getInterfaces();
}

/*!
 * Return a compact vector with the Ids of cells given by
 * bitpit::PatchKernel unique-labeled indexing
 * \return id of cells
 */
livector1D
MimmoObject::getCellsIds(){
    return getPatch()->getCells().getIds();
};

/*!
 * \return pointer to bitpit::PatchKernel structure hold by the class.
 */
PatchKernel*
MimmoObject::getPatch(){
    if(!m_internalPatch) return m_extpatch;
    else return m_patch.get();
};

/*!
 * \return pointer to bitpit::PatchKernel structure hold by the class.
 */
const PatchKernel*
MimmoObject::getPatch() const{
    if(!m_internalPatch) return m_extpatch;
    else return m_patch.get();
};

/*!
 * Return the indexing vertex map, to pass from local, compact indexing
 * to bitpit::PatchKernel unique-labeled indexing
 * \return local/unique-id map
 */
livector1D
MimmoObject::getMapData(){
    livector1D mapData(getNVertex());
    int i = 0;
    for (auto const & vertex : getVertices()){
        mapData[i] = vertex.getId();
        ++i;
    }
    return mapData;
};

/*!
 * Return the indexing vertex map, to pass from bitpit::PatchKernel unique-labeled indexing
 * to local, compact indexing;
 * \return unique-id/local map
 */
liimap
MimmoObject::getMapDataInv(){
    liimap mapDataInv;
    int i = 0;
    for (auto const & vertex : getVertices()){
        mapDataInv[vertex.getId()] = i;
        ++i;
    }
    return mapDataInv;
};

/*!
 * Return the indexing cell map, to pass from local, compact indexing
 * to bitpit::PatchKernel unique-labeled indexing
 * \return local/unique-id map
 */
livector1D
MimmoObject::getMapCell(){
    livector1D mapCell(getNCells());
    int i = 0;
    for (auto const & cell : getCells()){
        mapCell[i] = cell.getId();
        ++i;
    }
    return mapCell;
};

/*!
 * Return the indexing cell map, to pass from bitpit::PatchKernel unique-labeled indexing
 * to local, compact indexing;
 * \return unique-id/local map
 */
liimap
MimmoObject::getMapCellInv(){
    liimap mapCellInv;
    int i = 0;
    for (auto const & cell : getCells()){
        mapCellInv[cell.getId()] = i;
        ++i;
    }
    return mapCellInv;
};

/*!
 * \return the list of PID types actually present in your geometry.
 * If empty list is returned, pidding is actually not supported for this geometry
 */
std::unordered_set<short> 	&
MimmoObject::getPIDTypeList(){
    return m_pidsType;
};

/*!
 * \return the list of PID associated to each cell of tessellation in compact 
 * sequential ordering
 * If empty list is returned, pidding is not supported for this geometry
 */
shivector1D 	
MimmoObject::getCompactPID() {
    if(!m_skdTreeSupported || m_pidsType.empty())	return shivector1D(0);
    shivector1D result(getNCells());
    int counter=0;
    for(auto const & cell : getCells()){
        result[counter] = (short)cell.getPID();
        ++counter;
    }
    return(result);
};

/*!
 * \return the map of the PID (argument) associated to each cell unique-id (key) in tessellation.
 * If empty map is returned, pidding is not supported for this geometry.
 */
std::unordered_map<long,short> 	
MimmoObject::getPID() {
    if(!m_skdTreeSupported || m_pidsType.empty())	return std::unordered_map<long,short>();
    std::unordered_map<long,short> 	result;
    for(auto const & cell : getCells()){
        result[cell.getId()] = (short) cell.getPID();
    }
    return(result);
};

/*!
 * \return true if the skdTree ordering structure for cells is built/synchronized 
 * with your current geometry
 */
bool
MimmoObject::isBvTreeBuilt(){
    return isSkdTreeSync();
}

/*!
 * \return true if the skdTree ordering structure for cells is built/synchronized 
 * with your current geometry
 */
bool
MimmoObject::isBvTreeSync(){
    return isSkdTreeSync();
}

/*!
 * \return true if the skdTree ordering structure for cells is built/synchronized 
 * with your current geometry
 */
bool
MimmoObject::isSkdTreeSync(){
    return m_skdTreeSync;
}

/*!
 * \return pointer to geometry skdTree search structure
 */
bitpit::PatchSkdTree*
MimmoObject::getBvTree(){
    return getSkdTree();
}

/*!
 * \return pointer to geometry skdTree search structure
 */
bitpit::PatchSkdTree*
MimmoObject::getSkdTree(){
    if (!m_skdTreeSupported) return NULL;
    if (!isSkdTreeSync()) buildSkdTree();
    return m_skdTree.get();
}


/*!
 * \return true if the kdTree vertices ordering structure is 
 * built/synchronized with your current geometry
 */
bool
MimmoObject::isKdTreeBuilt(){
    return isKdTreeSync();
}

/*!
 * \return true if the kdTree vertices ordering structure is 
 * built/synchronized with your current geometry
 */
bool
MimmoObject::isKdTreeSync(){
    return m_kdTreeSync;
}

/*!
 * \return pointer to geometry KdTree internal structure
 */
bitpit::KdTree<3, bitpit::Vertex, long >*
MimmoObject::getKdTree(){
    if (!m_kdTreeSync) buildKdTree();
    return m_kdTree.get();
}


/*!
 * Set the vertices structure of the class, clearing any previous vertex list stored.
 * Be careful: any connectivity information stored in an existent cell list will be erased too.
 * The cell list will survive, but carrying no connectivity information.
 * \param[in] vertices geometry vertex structure .
 * \return false if no geometry is linked, not all vertices inserted or empty argument.
 */
bool
MimmoObject::setVertices(const bitpit::PiercedVector<bitpit::Vertex> & vertices){

    if (vertices.empty()) return false;

    getPatch()->resetVertices();

    int sizeVert = vertices.size();
    getPatch()->reserveVertices(sizeVert);

    long id;
    darray3E coords;
    bool checkTot = true;
    for (auto && val : vertices){
        id = val.getId();
        coords = val.getCoords();
        checkTot = checkTot && addVertex(coords, id);
    }

    return checkTot;
};

/*!
 *It adds one vertex to the mesh.
 * If unique-id is specified for the vertex, assign it, otherwise provide itself
 * to get a unique-id for the added vertex. The latter option is the default.
 * If the unique-id is already assigned, return with unsuccessful insertion.
 * 
 * 
 * \param[in] vertex vertex coordinates to be added 
 * \param[in] idtag  unique id associated to the vertex	
 * \return true if the vertex is successful inserted.
 */
bool
MimmoObject::addVertex(const darray3E & vertex, const long idtag){

    if(idtag != bitpit::Vertex::NULL_ID && getVertices().exists(idtag))    return false;

    bitpit::PatchKernel::VertexIterator it;
    auto patch = getPatch();
    if(idtag == bitpit::Vertex::NULL_ID){
        it = patch->addVertex(vertex);
    }else{
        it = patch->addVertex(vertex, idtag);
    }

    m_skdTreeSync = false;
    m_kdTreeSync = false;
    return true;
};

/*!
 *It adds one vertex to the mesh.
 * If unique-id is specified for the vertex, assign it, otherwise provide itself
 * to get a unique-id for the added vertex. The latter option is the default.
 * If the unique-id is already assigned, return with unsuccessful insertion.
 *
 *
 * \param[in] vertex vertex to be added
 * \param[in] idtag  unique id associated to the vertex
 * \return true if the vertex is successful inserted.
 */
bool
MimmoObject::addVertex(const bitpit::Vertex & vertex, const long idtag){

    if(idtag != bitpit::Vertex::NULL_ID && getVertices().exists(idtag))    return false;

    bitpit::PatchKernel::VertexIterator it;
    auto patch = getPatch();
    if(idtag == bitpit::Vertex::NULL_ID){
        it = patch->addVertex(vertex);
    }else{
        it = patch->addVertex(vertex, idtag);
    }

    m_skdTreeSync = false;
    m_kdTreeSync = false;
    return true;
};


/*!
 * It modifies the coordinates of a pre-existent vertex.
 * \param[in] vertex new vertex coordinates
 * \param[in] id unique-id of the vertex meant to be modified
 * \return false if no geometry is present or vertex id does not exist.
 */
bool
MimmoObject::modifyVertex(const darray3E & vertex, const long & id){

    if(!(getVertices().exists(id))) return false;
    bitpit::Vertex &vert = getPatch()->getVertex(id);
    vert.setCoords(vertex);
    m_skdTreeSync = false;
    m_kdTreeSync = false;
    return true;
};

/*!
 * Sets the cell structure of the geometry, clearing any previous cell list stored.
 * Does not do anything if class type is a point cloud (mesh type 3).
 *
 * \param[in] cells cell structure of geometry mesh.
 * \return false if no geometry is linked, not all cells are inserted or empty argument.
 */
bool
MimmoObject::setCells(const bitpit::PiercedVector<Cell> & cells){

    if (cells.empty() || !m_skdTreeSupported) return false;

    m_pidsType.clear();
    getPatch()->resetCells();

    int sizeCell = cells.size();
    getPatch()->reserveCells(sizeCell);

    long idc;
    int  nSize;
    short pid;
    bitpit::ElementType eltype;
    bool checkTot = true;
    livector1D connectivity;

    for (const auto & cell : cells){
        // get ID
        idc = cell.getId();
        //check on element type
        eltype = cell.getType();
        //check info PID 
        pid = (short)cell.getPID();
        nSize = cell.getConnectSize();
        connectivity.resize(nSize);
        auto const conn = cell.getConnect();
        for(int i=0; i<nSize; ++i){
            connectivity[i] = conn[i];
        }
        checkTot = checkTot && addConnectedCell(connectivity, eltype, pid, idc);
    }

    return checkTot;
};

/*!
 * It adds one cell with its vertex-connectivity (vertex in bitpit::PatchKernel unique id's), the type of cell to
 * be added and its own unique id. If no id is specified, teh method assigns it automatically.
 * Any kind of cell in bitpit::ElementType enum can be added according to mesh dimensionality 
 * (3D element in volume mesh, 2D in surface mesh, etc..). The method does nothing, if class type 
 * is a pointcloud one (type 3).
 * As a reminder for connectivity conn argument:
 *  - Connectivity of polygons must be defined as (nV,V1,V2,V3,V4,...) where nV is the number of vertices
 * defining the polygon and V1,V2,V3... are indices of vertices.
 *  - Connectivity of polyhedrons must be defined as (nF,nF1V, V1, V2, V3,..., nF2V, V2,V4,V5,V6,...,....)
 * where nF is the total number of faces of polyhedros, nF1V is the number of vertices composing the face 1,
 * followed by the indices of vertices which draws it. Face 2,3,..,n are defined in the same way.
 *  - Any other standard cell element is uniquely defined by its list of vertex indices.
 *
 * \param[in] conn  connectivity of target cell of geometry mesh.
 * \param[in] type  type of element to be added, according to bitpit::ElementInfo enum.
 * \param[in] idtag id of the cell
 * \return false if no geometry is linked, idtag already assigned or mismatched connectivity/element type 
 */
bool
MimmoObject::addConnectedCell(const livector1D & conn, bitpit::ElementType type, long idtag){

    if (conn.empty() || !m_skdTreeSupported) return false;
    if(idtag != bitpit::Cell::NULL_ID && getCells().exists(idtag)) return false;

    if(!checkCellConnCoherence(type, conn))  return false; 

    bitpit::PatchKernel::CellIterator it;
    auto patch = getPatch();

    if(idtag == bitpit::Cell::NULL_ID){
        it = patch->addCell(type, true, conn);
    }else{
        it = patch->addCell(type, true,conn, idtag);
    }

    m_pidsType.insert(0);		
    m_skdTreeSync = false;
    m_AdjBuilt = false;
    m_IntBuilt = false;
    return true;
};

/*!
 * It adds one cell with its vertex-connectivity (vertex in bitpit::PatchKernel unique id's), the type of cell to
 * be added and its own unique id. If no id is specified, teh method assigns it automatically.
 * Cell type and connectivity dimension of the cell are cross-checked with the mesh type of the class: if mismatch, the method
 * does not add the cell and return false. 
 * The method does nothing, if class type is a pointcloud one (type 3).
 * A part identifier PID mark can be associated to the cell.
 *
 * As a reminder for connectivity conn argument:
 *  - Connectivity of polygons must be defined as (nV,V1,V2,V3,V4,...) where nV is the number of vertices
 * defining the polygon and V1,V2,V3... are indices of vertices.
 *  - Connectivity of polyhedrons must be defined as (nF,nF1V, V1, V2, V3,..., nF2V, V2,V4,V5,V6,...,....)
 * where nF is the total number of faces of polyhedros, nF1V is the number of vertices composing the face 1,
 * followed by the indices of vertices which draws it. Face 2,3,..,n are defined in the same way.
 *  - Any other standard cell element is uniquely defined by its list of vertex indices.
 * 
 * \param[in] conn  connectivity of target cell of geometry mesh.
 * \param[in] type  type of element to be added, according to bitpit ElementInfo enum.
 * \param[in] PID   part identifier
 * \param[in] idtag id of the cell
 * \return false if no geometry is linked, idtag already assigned or mismatched connectivity/element type 
 */
bool
MimmoObject::addConnectedCell(const livector1D & conn, bitpit::ElementType type, short PID, long idtag){

    if (conn.empty() || !m_skdTreeSupported) return false;
    if(idtag != bitpit::Cell::NULL_ID && getCells().exists(idtag)) return false;

    if(!checkCellConnCoherence(type, conn))  return false; 

    bitpit::PatchKernel::CellIterator it;
    auto patch = getPatch();

    long checkedID;
    if(idtag == bitpit::Cell::NULL_ID){
        it = patch->addCell(type, true, conn);
        checkedID = it->getId();
    }else{
        it = patch->addCell(type, true,conn, idtag);
        checkedID = idtag;
    }

    setPIDCell(checkedID, PID);
    m_skdTreeSync = false;
    m_AdjBuilt = false;
    m_IntBuilt = false;
    return true;
};

/*!
 * Set PIDs for all geometry cells available. 
 * The PID list must be referred to the compact local/indexing of the cells in the class.
 * \param[in] pids PID list.
 */
void
MimmoObject::setPID(shivector1D pids){
    if((int)pids.size() != getNCells() || !m_skdTreeSupported)	return;

    m_pidsType.clear();
    int counter = 0;
    for(auto & cell: getCells()){
        m_pidsType.insert(pids[counter]);
        cell.setPID(pids[counter]);
        ++counter;
    }
};

/*!
 * Set PIDs for all geometry cells available. 
 * The PID must be provided as a map with cell unique-id as key and pid associated to it as argument.
 * \param[in] pidsMap PID amp.
 */
void
MimmoObject::setPID(std::unordered_map<long, short> pidsMap){
    if(getNCells() == 0 || !m_skdTreeSupported)	return;

    m_pidsType.clear();
    auto & cells = getCells();
    for(auto const & val: pidsMap){
        if(cells.exists(val.first)){
            cells[val.first].setPID(val.second);
            m_pidsType.insert(val.second);
        }
    }
};

/*!
 * Set the PID of a target cell
 * \param[in] id unique-id of the cell
 * \param[in] pid PID to assign on cell
 */
void
MimmoObject::setPIDCell(long id, short pid){
    auto & cells = getCells();
    if(cells.exists(id)){
        cells[id].setPID((int)pid);
        m_pidsType.insert(pid);
    }	
};

/*!
 * Update class member map m_pidsType with PIDs effectively contained 
 * in the internal/linked geometry patch.
 */
void
MimmoObject::resyncPID(){
    m_pidsType.clear();
    for(auto const & cell : getCells()){
        m_pidsType.insert( (short)cell.getPID() );
    }
};


/*!
 * Set your current class as a "hard" copy of another MimmoObject.
 * All preexistent data are destroyed and replaced by those provided by the argument, except for search Trees,
 * which are instantiated, but not filled/built/synchronized.
 * Hard means that the geometry data structure is allocated internally 
 * as an exact and stand-alone copy of the geometry data structure of the argument, 
 * indipendently on how the argument owns or links it.
 * \param[in] other MimmoObject class.
 */
void MimmoObject::setHARDCopy(const MimmoObject * other){

    m_type  = other->m_type;
    m_internalPatch = true;
    m_extpatch = NULL;
    
    switch(m_type){
        case 1:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoSurfUnstructured(2)));
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new bitpit::SurfaceSkdTree(dynamic_cast<SurfaceKernel*>(m_patch.get()))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 2:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoVolUnstructured(3)));
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new VolumeSkdTree(dynamic_cast<VolumeKernel*>(m_patch.get()))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 3:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoPointCloud));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 4:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoSurfUnstructured(1)));
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new bitpit::SurfaceSkdTree(dynamic_cast<SurfaceKernel*>(m_patch.get()))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        default:
            //never been reached
            break;
    }
    
    m_skdTreeSupported = other->m_skdTreeSupported;
    m_skdTreeSync = false;
    m_kdTreeSync = false;

    //copy data 
    const bitpit::PiercedVector<bitpit::Vertex> & pvert = other->getVertices();
    setVertices(pvert);

    if(m_skdTreeSupported){
        const bitpit::PiercedVector<bitpit::Cell> & pcell = other->getCells();
        setCells(pcell);
    }

    if(other->m_AdjBuilt)   buildAdjacencies();
    if(other->m_IntBuilt)   buildInterfaces();

};

/*!
 * Clone your MimmoObject in a new indipendent MimmoObject. All data of the current class will be "hard" copied in a new
 * MimmoObject class. Search Trees will be only instantiated but not filled/built/synchronized.
 * Hard means that the geometry data structure will be allocated internally into the new object 
 * as an exact and stand-alone copy of the geometry data structure of the current class, indipendently on how the current class
 * owns or links it.
 * \return cloned MimmoObject.
 */
std::unique_ptr<MimmoObject> MimmoObject::clone(){
    std::unique_ptr<MimmoObject> result(new MimmoObject(getType()));
    //copy data 
    result->setVertices(getVertices());
    if(m_skdTreeSupported){
        result->setCells(getCells());
    }
    if(m_AdjBuilt)   result->buildAdjacencies();
    if(m_IntBuilt)   result->buildInterfaces();
    return std::move(result);
};

/*!
 * It cleans geometry duplicated and, in case of connected tessellations, all orphan/isolated vertices.
 * \return false if the geometry member pointer is NULL.
 */
bool
MimmoObject::cleanGeometry(){
    auto patch = getPatch();
    patch->deleteCoincidentVertices();
    if(m_skdTreeSupported)  patch->deleteOrphanVertices();

    m_kdTreeSync = false;
    return true;
};

/*! 
 * Extract vertex list from an ensamble of geometry cells.
 *\param[in] cellList list of bitpit::PatchKernel ids identifying cells.
 *\return the list of bitpit::PatchKernel ids of involved vertices.
 */  
livector1D MimmoObject::getVertexFromCellList(livector1D cellList){
    if(isEmpty() || getType() == 3)   return livector1D(0);

    livector1D result;
    set<long int> ordV;
    auto patch = getPatch();
    //get conn from each cell of the list
    for(const auto id : cellList){
        if(getCells().exists(id)){
            bitpit::ConstProxyVector<long> ids = patch->getCell(id).getVertexIds();
            for(const auto & val: ids){
                ordV.insert(val);
            }
        }
    }

    result.reserve(ordV.size());
    result.insert(result.end(), ordV.begin(), ordV.end());
    return  result;
}

/*! 
 * Extract interfaces list from an ensamble of geometry cells.
 *\param[in] cellList list of bitpit::PatchKernel ids identifying cells.
 *\return the list of bitpit::PatchKernel ids of involved interfaces.
 */  
livector1D MimmoObject::getInterfaceFromCellList(livector1D cellList){
    if(isEmpty() || getType() == 3)   return livector1D(0);
    
    if(!areInterfacesBuilt())   buildInterfaces();
    livector1D result;
    set<long int> ordV;
    auto patch = getPatch();
    //get conn from each cell of the list
    for(const auto id : cellList){
        if(getCells().exists(id)){
            long * interf = patch->getCell(id).getInterfaces();
            int nIloc = patch->getCell(id).getInterfaceCount();
            for(int i=0; i<nIloc; ++i)  ordV.insert(interf[i]);
        }
    }
    
    result.reserve(ordV.size());
    result.insert(result.end(), ordV.begin(), ordV.end());
    return  result;
}


/*! 
 * Extract cell list from an ensamble of geometry vertices. Ids of all those cells whose vertex are 
 * defined inside the selection will be returned.
 * \param[in] vertexList list of bitpit::PatchKernel IDs identifying vertices.
 * \return list of bitpit::PatchKernel IDs of involved 1-Ring cells.
 */

livector1D MimmoObject::getCellFromVertexList(livector1D vertexList){
    if(isEmpty() || getType() == 3)   return livector1D(0);
   
    livector1D result;
    std::unordered_set<long int> ordV, ordC;
    ordV.insert(vertexList.begin(), vertexList.end());
    //get conn from each cell of the list
    for(auto const & cell : getPatch()->getCells()){
        bitpit::ConstProxyVector<long> vIds= cell.getVertexIds();
        bool check;
        for(const auto & id : vIds){
            check = (ordV.count(id) > 0);
            if(!check)  break ;
        }
        if(check) ordC.insert(cell.getId());
    }

    result.reserve(ordC.size());
    result.insert(result.end(), ordC.begin(), ordC.end());
    return  result;
}

/*!
 * Extract vertices at the mesh boundaries, if any. The method is meant for connected mesh only,
 * return empty list otherwise.
 * \return list of vertex unique-ids.
 */
livector1D 	MimmoObject::extractBoundaryVertexID(){

    std::unordered_map<long, std::set<int> > cellmap = extractBoundaryFaceCellID();
    if(cellmap.empty()) return livector1D(0);

    std::unordered_set<long> container;
 
    for (const auto & val : cellmap){
        bitpit::Cell & cell = getPatch()->getCell(val.first);
        for(const auto face : val.second){
            bitpit::ConstProxyVector<long> list = cell.getFaceVertexIds(face);
            for(const auto & index : list ){
                container.insert(index);
            }
        }// end loop on face
    }

    livector1D result;
    result.reserve(container.size());
    result.insert(result.end(), container.begin(), container.end());

    return result;
};

/*!
 * Extract cells who have one face at the mesh boundaries at least, if any.
 * The method is meant for connected mesh only, return empty list otherwise.
 * \return list of cell unique-ids.
 */
livector1D  MimmoObject::extractBoundaryCellID(){
    
    if(isEmpty() || m_type==3)   return livector1D(0);
    if(!areAdjacenciesBuilt())   getPatch()->buildAdjacencies();
    
    std::unordered_set<long> container;
    
    for (const auto & cell : getCells()){
        int size = cell.getFaceCount();
        
        for(int face=0; face<size; ++face){
            if(cell.isFaceBorder(face)){
                container.insert(cell.getId());
            }//endif
        }// end loop on face
    }
    
    livector1D result;
    result.reserve(container.size());
    result.insert(result.end(), container.begin(), container.end());
    
    return result;
};

/*!
 * Extract cells  who have one face at the mesh boundaries at least, if any.
 * Return the list of the local faces per cell, which lie exactly on the boundary.
 * The method is meant for connected mesh only, return empty list otherwise.
 * \return map of boundary cell unique-ids, with local boundary faces list.
 */
std::unordered_map<long, std::set<int> >  MimmoObject::extractBoundaryFaceCellID(){

    std::unordered_map<long, std::set<int> > result;
    if(isEmpty() || m_type ==3)   return result;
    if(!areAdjacenciesBuilt())   getPatch()->buildAdjacencies();

    for (const auto & cell : getCells()){
        int size = cell.getFaceCount();
        long idC = cell.getId();
        for(int face=0; face<size; ++face){
            if(cell.isFaceBorder(face)){
                result[idC].insert(face);
            }//endif
        }// end loop on face
    }
    return result;
};

/*!
 * Extract vertices at the mesh boundaries, provided the map of the cell faces at boundaries. 
 * The method is meant for connected mesh only, return empty list otherwise.
 * \param[in] cellmap map of border faces of the mesh written as cell-ID vs local cell face index.
 * \return list of vertex unique-ids.
 */
livector1D  MimmoObject::extractBoundaryVertexID(std::unordered_map<long, std::set<int> > &cellmap){
    
    if(cellmap.empty()) return livector1D(0);
    
    std::unordered_set<long> container;
    
    for (const auto & val : cellmap){
        bitpit::Cell & cell = getPatch()->getCell(val.first);
        for(const auto face : val.second){
            bitpit::ConstProxyVector<long> list = cell.getFaceVertexIds(face);
            for(const auto & index : list ){
                container.insert(index);
            }
        }// end loop on face
    }
    
    livector1D result;
    result.reserve(container.size());
    result.insert(result.end(), container.begin(), container.end());
    
    return result;
};

/*!
 * Extract all cells marked with a target PID flag.
 * \param[in]   flag    PID for extraction
 * \return  list of cells as unique-ids
 */
livector1D	MimmoObject::extractPIDCells(short flag){

    if(m_pidsType.count(flag) < 1)	return livector1D(0);

    livector1D result(getNCells());
    int counter = 0;
    for(auto const & cell : getCells()){
        if ( cell.getPID() == flag)	{
            result[counter] = cell.getId();
            ++counter;
        }
    }
    result.resize(counter);
    return  result;
};

/*!
 * Extract all cells marked with a series of target PIDs.
 * \param[in]   flag    list of PID for extraction
 * \return      list of cells as unique-ids
 */
livector1D	MimmoObject::extractPIDCells(shivector1D flag){
    livector1D result;
    result.reserve(getNCells());
    for(auto && id : flag){
        livector1D partial = extractPIDCells(id);
        result.insert(result.end(),partial.begin(), partial.end());
    }
    return(result);
};

/*!
 * Check if a given connectivity list is coherent with a bitpit::ElementType type.
 * 
 * \param[in] type cell type to check, as bitpit::ElementInfo enum
 * \param[in] conn connectivity list.
 * \return true if connectivity list is coherent with the specified cell type
 */
bool MimmoObject::checkCellConnCoherence(const bitpit::ElementType & type, const livector1D & list){

    switch(type){
        case bitpit::ElementType::VERTEX:
            return (list.size() == 1);
            break;
        case bitpit::ElementType::LINE:
            return (list.size() == 2);
            break;
        case bitpit::ElementType::TRIANGLE:
            return (list.size() == 3);
            break;
        case bitpit::ElementType::PIXEL:
            return (list.size() == 4);
            break;
        case bitpit::ElementType::QUAD:
            return (list.size() == 4);
            break;
        case bitpit::ElementType::POLYGON:
            if(list.size() < 5) return false;
            return(list.size() == std::size_t(list[0]+1));
            break;
        case bitpit::ElementType::TETRA:
            return (list.size() == 4);
            break;
        case bitpit::ElementType::VOXEL:
            return (list.size() == 8);
            break;
        case bitpit::ElementType::HEXAHEDRON:
            return (list.size() == 8);
            break;
        case bitpit::ElementType::WEDGE:
            return (list.size() == 6);
            break;
        case bitpit::ElementType::PYRAMID:
            return (list.size() == 5);
            break;
        case bitpit::ElementType::POLYHEDRON:
            if(list.size() < 9) return false;
            {
                //check if the record is a polyhedron
                long nFaces = list[0];
                std::size_t pos = 1;
                long countFaces = 0;
                while(pos < list.size() && countFaces<nFaces){
                    countFaces++;
                    pos += list[pos]+1;
                }
                return (pos == list.size() && countFaces == nFaces);
            }
            break;
        default:
            assert(false && "reached uncovered case");
            break;
    }
    return false;
};

/*!
 * Evaluate axis aligned bounding box of the current MimmoObject 
 * \param[out] pmin lowest bounding box point
 * \param[out] pmax highest bounding box point
 */
void MimmoObject::getBoundingBox(std::array<double,3> & pmin, std::array<double,3> & pmax){
    getPatch()->getBoundingBox(pmin,pmax);
    return;
}

/*!
 * Reset and build again cell skdTree of your geometry (if supports connectivity elements).
 *\param[in] value build the minimum leaf of the tree as a bounding box containing value elements at most.
 */
void MimmoObject::buildBvTree(int value){
    if(!m_skdTreeSupported || isEmpty())   return;

    if (!m_skdTreeSync){
        m_skdTree->clear();
        m_skdTree->build(value);
        m_skdTreeSync = true;
    }
    return;
}

/*!
 * Reset and build again cell skdTree of your geometry (if supports connectivity elements).
 *\param[in] value build the minimum leaf of the tree as a bounding box containing value elements at most.
 */
void MimmoObject::buildSkdTree(int value){
    if(!m_skdTreeSupported || isEmpty())   return;
    
    if (!m_skdTreeSync){
        m_skdTree->clear();
        m_skdTree->build(value);
        m_skdTreeSync = true;
    }
    return;
}

/*!
 * Reset and build again vertex kdTree of your geometry.
 */
void MimmoObject::buildKdTree(){
    if( getNVertex() == 0)  return;
    long label;
    
    if (!m_kdTreeSync){
        cleanKdTree();
        m_kdTree->nodes.resize(getNVertex() + m_kdTree->MAXSTK);

        for(auto & val : getVertices()){
            label = val.getId();
            m_kdTree->insert(&val, label);
        }
        m_kdTreeSync = true;
    }
    return;
}

/*!
 * Clean the KdTree of the class
 */
void	MimmoObject::cleanKdTree(){
    m_kdTree->n_nodes = 0;
    m_kdTree->nodes.clear();
}

/*!
 * \return true if cell-cell adjacency is built for your current mesh.
 */
bool MimmoObject::areAdjacenciesBuilt(){
    return  m_AdjBuilt;
};

/*!
 * \return true if Interfaces are built for your current mesh.
 */
bool MimmoObject::areInterfacesBuilt(){
    return  m_IntBuilt;
};

/*!
 * \return false if your mesh has open edges/faces. True otherwise
 */
bool MimmoObject::isClosedLoop(){

    if(!areAdjacenciesBuilt())	buildAdjacencies();
    bool check = true;

    auto itp = getCells().cbegin();
    auto itend = getCells().cend();
    
    while(itp != itend && check){
        
        int cAdj = itp->getAdjacencyCount();
        const long * adj = itp->getAdjacencies();

        for(int i=0; i<cAdj; ++i){
            check = check && (adj[i] != -1);
        }

        itp++;
    }

    return check;
};

/*!
 * Force the class to build cell-cell adjacency connectivity.
 */
void MimmoObject::buildAdjacencies(){
    if(m_type !=3){
        getPatch()->buildAdjacencies();
        m_AdjBuilt = true;
    }
};

/*!
 * Force the class to build Interfaces connectivity.
 * If MimmoObject does not support connectivity (as Point Clouds do) 
 * does nothing. 
 */
void MimmoObject::buildInterfaces(){
    if(m_type !=3){
        if(!areAdjacenciesBuilt()) buildAdjacencies();
        getPatch()->buildInterfaces();
        m_IntBuilt=  true;
    }
};


/*!
 * Desume Element given the vertex connectivity list associated. Polygons and Polyhedra require 
 * a special writing for their connectivity list. Please read doxy of MimmoObject(int type, dvecarr3E & vertex, livector2D * connectivity = NULL)
 * custom constructor for further information.
 * Return undefined type for unexistent or unsupported element.
 * \param[in] locConn list of vertex indices composing the cell.
 * \return cell type hold by the current connectivity argument.
 */
//TODO To review in order to implement new derived classes MimmoSurfUnstructured and MimmoVolUnstructured
bitpit::ElementType	MimmoObject::desumeElement(const livector1D & locConn){
    bitpit::ElementType result = bitpit::ElementType::UNDEFINED;

    std::size_t sizeConn = locConn.size();
    switch(m_type){
        case    1:
            if(sizeConn == 3)        result = bitpit::ElementType::TRIANGLE;
            if(sizeConn == 4)        result = bitpit::ElementType::QUAD;
            if(sizeConn > 4 && sizeConn == std::size_t(locConn[0]+1))    result= bitpit::ElementType::POLYGON;
            break;
        case    2:
            if(sizeConn == 4)        result = bitpit::ElementType::TETRA;
            if(sizeConn == 8)        result = bitpit::ElementType::HEXAHEDRON;
            if(sizeConn == 5)        result = bitpit::ElementType::PYRAMID;
            if(sizeConn == 6)        result = bitpit::ElementType::WEDGE;
            if(sizeConn > 8 ){
                //check if the record is a polyhedron
                long nFaces = locConn[0];
                std::size_t pos = 1;
                long countFaces = 0;
                while(pos < sizeConn && countFaces<nFaces){
                    countFaces++;
                    pos += locConn[pos]+1;
                }
                if(pos == sizeConn && countFaces == nFaces){
                    result= bitpit::ElementType::POLYHEDRON;
                }
            }
            break;
        case    3:
            result = bitpit::ElementType::VERTEX;
            break;
        case    4:
            result = bitpit::ElementType::LINE;
            break;
        default : 
            assert(false && "reached uncovered case");
            break;
    }

    return result;
};

/*!
 * Reset class to default constructor set-up.
 * \param[in] type type of mesh from 1 to 4; See default constructor.
 */
void MimmoObject::reset(int type){
    m_log = &bitpit::log::cout(MIMMO_LOG_FILE);
    m_type = std::max(0, type);
    if (m_type > 4){
        m_type = 1;
    }

    switch(m_type){
        case 1:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoSurfUnstructured(2)));
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new bitpit::SurfaceSkdTree(dynamic_cast<SurfaceKernel*>(m_patch.get()))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 2:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoVolUnstructured(3)));
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new VolumeSkdTree(dynamic_cast<VolumeKernel*>(m_patch.get()))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 3:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoPointCloud()));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        case 4:
            m_patch = std::move(std::unique_ptr<PatchKernel>(new MimmoSurfUnstructured(1)));
            m_skdTree = std::move(std::unique_ptr<PatchSkdTree>(new bitpit::SurfaceSkdTree(dynamic_cast<SurfaceKernel*>(m_patch.get()))));
            m_kdTree  = std::move(std::unique_ptr<KdTree<3,bitpit::Vertex,long> >(new KdTree<3,bitpit::Vertex, long>())); 
            break;
        default:
            //never been reached
            break;
    }

    m_internalPatch = true;
    m_extpatch = NULL;

    m_skdTreeSupported = (m_type != 3);
    m_skdTreeSync = false;
    m_kdTreeSync = false;
    m_AdjBuilt = false;
    m_IntBuilt = false;
}

/*!
 * Dump contents of your current MimmoObject to a stream. 
 * Search trees and adjacencies will not be dumped.
 * Write all in a binary format.
 * \param[in,out] stream to write on.
 */
void MimmoObject::dump(std::ostream & stream){
    bitpit::utils::binary::write(stream,m_type);
    getPatch()->dump(stream);
}

/*!
 * Restore contents of a dumped MimmoObject from a stream in your current class.
 * New restored data will be owned internally by the class.
 * Every data previously stored will be lost.
 * Read all in a binary format.
 * \param[in,out] stream to write on.
 */

void MimmoObject::restore(std::istream & stream){
    int type;
    bitpit::utils::binary::read(stream,type);
    reset(type);

    getPatch()->restore(stream);

    for (const auto &cell: getCells()){
        m_pidsType.insert(short(cell.getPID()));
    }
}

/*!
 * Evaluate general volume of each cell in the current mesh, 
 * according to its topology.
 * If no patch or empty patch or pointcloud patch is present 
 * in the current MimmoObject return empty map as result.
 * \param[out] volumes volume values referred to cell ids 
 */
void
MimmoObject::evalCellVolumes(bitpit::PiercedVector<double> & volumes){
    
    if(getPatch() == NULL)   return;
    if(isEmpty())       return ;

    switch (getType()){
        case 1:
        case 4:
            {
                bitpit::SurfaceKernel * p = static_cast<bitpit::SurfaceKernel *>(getPatch());
                for (const auto & cell: getCells()){
                    volumes.insert(cell.getId(), p->evalCellArea(cell.getId()));
                }
            }
            break;
        case 2:
            {
                bitpit::VolumeKernel * p = static_cast<bitpit::VolumeKernel *>(getPatch());
                for (const auto & cell: getCells()){
                    volumes.insert(cell.getId(), p->evalCellVolume(cell.getId()));
                }
            }
            break;
        default:
            //leave map empty
            break;
    }
}

/*!
 * Evaluate Aspect Ratio of each cell in the current mesh, 
 * according to its topology.
 * VERTEX and LINE elements does not support a proper definition of Aspect Ratio 
 * and will return always 0.0 as value. PointCloud and 3DCurve MimmoObject return an empty map.
 * If no patch or empty patch or pointcloud patch is present 
 * in the current MimmoObject return empty map as result.
 * \param[out] ARs aspect ratio values referred to cell id
 */
void
MimmoObject::evalCellAspectRatio(bitpit::PiercedVector<double> & ARs){

    if(getPatch() == NULL)   return;
    if(isEmpty())       return;
    
    switch (getType()){
        case 1:
        {
            bitpit::SurfaceKernel * p = static_cast<bitpit::SurfaceKernel *>(getPatch());
            int edge;
            for (const auto & cell: getCells()){
                ARs[cell.getId()] = p->evalAspectRatio(cell.getId(), edge);
            }
        }
        break;
        case 2:
        {   //Following the example of OpenFoam, the AR index is calculated as
            // the ratio S between total surface and hydraulic surface of an equilater
            //   cylinder (h=2*r) of the same volume, that is S_hyd = 1/6.0 * (V^(2/3)).
            if(!areInterfacesBuilt())   buildInterfaces();
            bitpit::VolUnstructured * p = static_cast<bitpit::VolUnstructured *>(getPatch());

            //calculate interface area
            std::unordered_map<long, double> interfaceAreas;
            for (const auto & interf: getInterfaces()){
                interfaceAreas[interf.getId()] = p->evalInterfaceArea(interf.getId());
            }

            double Svalue = 0.0;
            double sumArea;
            int size;
            for (const auto & cell: getCells()){

                sumArea = 0.0;

                size = cell.getInterfaceCount();
                const long * conInt = cell.getInterfaces();

                for(int i=0; i<size; ++i){
                    sumArea += interfaceAreas[conInt[i]];
                }

                double vol = p->evalCellVolume(cell.getId());
                if(vol <= std::numeric_limits<double>::min()){
                    Svalue = std::numeric_limits<double>::max();
                }else{
                    Svalue = sumArea/(6.0*std::pow(vol, 2.0/3.0));
                }

                ARs.insert(cell.getId(), Svalue);
            }

        }
        break;
        default:
            //leave map empty
            break;
    }
}

/*!
 * Evaluate volume of a target cell in the current mesh.
 * \param[in] id cell id.
 * \return cell volume
 */
double
MimmoObject::evalCellVolume(const long & id){
    switch (getType()){
        case 1:
        case 4:
            return static_cast<bitpit::SurfaceKernel *>(getPatch())->evalCellArea(id);
            break;
        case 2:
            return static_cast<bitpit::VolumeKernel *>(getPatch())->evalCellVolume(id);
            break;
        default:
            return 0.0;
            break;
    }
}

/*!
 * Evaluate Aspect Ratio of a target cell.
 * VERTEX and LINE elements does not support a proper definition of Aspect Ratio 
 * and will return always 0.0 as value.
 * \param[in] id of the cell
 * \return value of cell AR.
 */
double
MimmoObject::evalCellAspectRatio(const long & id){
    int edge;
    switch (getType()){
        case 1:
            return static_cast<bitpit::SurfaceKernel *>(getPatch())->evalAspectRatio(id, edge);
            break;
        case 2:
        {
            //Following the example of OpenFoam, the AR index is calculated as
            // the ratio S between total surface and hydraulic surface of an equilater
            //   cylinder (h=2*r) of the same volume, that is S_hyd = 1/6.0 * (V^(2/3)).
            if(!areInterfacesBuilt())   buildInterfaces();
            bitpit::VolUnstructured * p = static_cast<bitpit::VolUnstructured *>(getPatch());

            double Svalue = 0.0;
            double sumArea = 0.0;
            int size = p->getCell(id).getInterfaceCount();
            const long * conInt = p->getCell(id).getInterfaces();
            for(int i=0; i<size; ++i){
                sumArea += p->evalInterfaceArea(conInt[i]);
            }

            double vol = p->evalCellVolume(id);
            if(vol <= std::numeric_limits<double>::min()){
                return std::numeric_limits<double>::max();
            }else{
                return sumArea/(6.0*std::pow(vol, 2.0/3.0));
            }
        }
            break;
        default:
            return 0.0;
            break;
    }
}

}



