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
#ifndef __MIMMOOBJECT_HPP__
#define __MIMMOOBJECT_HPP__

#include "bitpit_patchkernel.hpp"
#include "bitpit_volunstructured.hpp"
#include "bitpit_surfunstructured.hpp"
#include "bitpit_SA.hpp"
#include "surface_skd_tree.hpp"
#include "volume_skd_tree.hpp"
#include "mimmoTypeDef.hpp"
#include "MimmoNamespace.hpp"

namespace mimmo{

/*!
 * \ingroup core
 */

/*!
 * \class MimmoSurfUnstructured
 * \brief Custom derivation of bitpit::SurfUnstructured class
 */
class MimmoSurfUnstructured: public bitpit::SurfUnstructured{
public:
    // Constructors
    MimmoSurfUnstructured();
    MimmoSurfUnstructured(int patch_dim);
    MimmoSurfUnstructured(const int &id, int patch_dim);
    MimmoSurfUnstructured(std::istream &stream);
    virtual ~MimmoSurfUnstructured();
    // Clone
    std::unique_ptr<bitpit::PatchKernel> clone() const;

protected:
    /*! Copy Constructor*/
    MimmoSurfUnstructured(const MimmoSurfUnstructured &) = default;
};

/*!
 * \class MimmoVolUnstructured
 * \brief Custom derivation of bitpit::VolUnstructured class
 */
class MimmoVolUnstructured: public bitpit::VolUnstructured{
public:
    // Constructors
    MimmoVolUnstructured(const int& dimension);
    MimmoVolUnstructured(const int &id, const int& dimension);
    virtual ~MimmoVolUnstructured();
    // Clone
    std::unique_ptr<bitpit::PatchKernel> clone() const;
    
protected:
    /*! Copy Constructor*/
    MimmoVolUnstructured(const MimmoVolUnstructured &) = default;
};

/*!
 * \class MimmoPointCloud
 * \brief Custom derivation of bitpit::SurfUnstructured class, for Point Cloud handling only
 */
class MimmoPointCloud: public bitpit::SurfUnstructured{
public:
    // Constructors
    MimmoPointCloud();
    MimmoPointCloud(const int &id);
    virtual ~MimmoPointCloud();
    // Clone
    std::unique_ptr<bitpit::PatchKernel> clone() const;
    
protected:
    /*! Copy Constructor*/
    MimmoPointCloud(const MimmoPointCloud &) = default;
};


/*!
* \class MimmoObject
* \brief MimmoObject is the basic geometry container for mimmo library
*
* MimmoObject is the basic container for 3D geometry unustructured mesh data structure, 
* based on bitpit::PatchKernel containers.
* MimmoObject can handle unustructured surface meshes, unustructured volume meshes, 3D point clouds and 3D tessellated curves.
* It supports interface methods to explore and handle the geometrical structure. It supports PID convention to mark subparts 
* of geometry as well as building the search-trees KdTree (3D point spatial ordering) and skdTree(Cell-AABB spatial ordering)
* to quickly retrieve vertices and cells in the data structure.
*/
class MimmoObject{

private:
    std::unique_ptr<bitpit::PatchKernel>    m_patch;           /**<Reference to INTERNAL bitpit patch handling geometry. */
    bitpit::PatchKernel *                   m_extpatch;        /**<Reference to EXTERNALLY linked patch handling geometry. */
    bool                                    m_internalPatch;   /**<True if the geometry is internally created. */

protected:
//members
    int                                                     m_type;            /**<Type of geometry (0 = undefined, 1 = surface mesh, 2 = volume mesh, 3-point cloud mesh, 4-3DCurve). */
    std::unordered_set<short>                               m_pidsType;        /**<pid type available for your geometry */
    std::unique_ptr<bitpit::PatchSkdTree>                   m_skdTree;         /**< ordered tree of geometry simplicies for fast searching purposes */
    std::unique_ptr<bitpit::KdTree<3,bitpit::Vertex,long> > m_kdTree;          /**< ordered tree of geometry vertices for fast searching purposes */
    bool                                                    m_skdTreeSync;      /**< track correct building of bvtree. Set false if any geometry modifications occur */
    bool                                                    m_kdTreeSync;     /**< track correct building of kdtree. Set false if any geometry modifications occur*/
    bool                                                    m_skdTreeSupported; /**< Flag for geometries not supporting bvTree building*/

    bool                                                    m_AdjBuilt;     /**< track correct building of adjacencies along with geometry modifications */
    bool                                                    m_IntBuilt;     /**< track correct building of interfaces  along with geometry modifications */
    bitpit::Logger*                                         m_log;          /**<Pointer to logger.*/

public:
    MimmoObject(int type = 1); 
    MimmoObject(int type, dvecarr3E & vertex, livector2D * connectivity = NULL);
    MimmoObject(int type, bitpit::PatchKernel* geometry);
    MimmoObject(int type, std::unique_ptr<bitpit::PatchKernel> & geometry);
    ~MimmoObject();

    MimmoObject(const MimmoObject & other);
    MimmoObject & operator=(MimmoObject other);

    bool                                            isEmpty();
    bool                                            isEmpty() const;
    BITPIT_DEPRECATED(bool                          isBvTreeSupported());
    bool                                            isSkdTreeSupported();
    int                                             getType();
    long                                            getNVertex()const;
    long                                            getNCells()const;
    dvecarr3E                                       getVertexCoords(liimap* mapDataInv = NULL);
    darray3E                                        getVertexCoords(long i);
    bitpit::PiercedVector<bitpit::Vertex> &         getVertices();
    const bitpit::PiercedVector<bitpit::Vertex> &   getVertices() const ;

    livector2D                                      getCompactConnectivity(liimap & mapDataInv);
    livector2D                                      getConnectivity();
    livector1D                                      getCellConnectivity(long id);
    bitpit::PiercedVector<bitpit::Cell> &           getCells();
    const bitpit::PiercedVector<bitpit::Cell> &     getCells() const;

    bitpit::PiercedVector<bitpit::Interface> &           getInterfaces();
    const bitpit::PiercedVector<bitpit::Interface> &     getInterfaces() const;
    
    livector1D                                      getCellsIds();
    bitpit::PatchKernel*                            getPatch();
    const bitpit::PatchKernel*                      getPatch() const;
    std::unordered_set<short> &                     getPIDTypeList();
    shivector1D                                     getCompactPID();
    std::unordered_map<long, short>                 getPID();

    BITPIT_DEPRECATED(bitpit::PatchSkdTree*         getBvTree());
    bitpit::PatchSkdTree*                           getSkdTree();
    bitpit::KdTree<3, bitpit::Vertex, long> *       getKdTree();
    BITPIT_DEPRECATED(bool                          isBvTreeBuilt());
    BITPIT_DEPRECATED(bool                          isKdTreeBuilt());
    BITPIT_DEPRECATED(bool                          isBvTreeSync());
    bool                          isSkdTreeSync();
    bool                          isKdTreeSync();


    bool        setVertices(const bitpit::PiercedVector<bitpit::Vertex> & vertices);
    bool        addVertex(const darray3E & vertex, const long idtag = bitpit::Vertex::NULL_ID);
    bool        addVertex(const bitpit::Vertex & vertex, const long idtag = bitpit::Vertex::NULL_ID);
    bool        modifyVertex(const darray3E & vertex, const long & id);
    bool        setCells(const bitpit::PiercedVector<bitpit::Cell> & cells);
    bool        addConnectedCell(const livector1D & locConn, bitpit::ElementType type, long idtag = bitpit::Cell::NULL_ID);
    bool        addConnectedCell(const livector1D & locConn, bitpit::ElementType type, short PID, long idtag = bitpit::Cell::NULL_ID);

    void        setPID(shivector1D ); 
    void        setPID(std::unordered_map<long, short>  ); 
    void        setPIDCell(long, short);
    void        resyncPID();

    BITPIT_DEPRECATED(void        setHARDCopy(const MimmoObject * other));
    std::unique_ptr<MimmoObject>        clone();

    bool        cleanGeometry();

    livector1D  getVertexFromCellList(livector1D cellList);
    livector1D  getCellFromVertexList(livector1D vertList, bool strict = true);
    livector1D  getInterfaceFromCellList(livector1D cellList);

    livector1D                          extractBoundaryVertexID();
    livector1D                          extractBoundaryCellID();
    std::unordered_map<long, std::set<int> > extractBoundaryFaceCellID();
    livector1D                          extractBoundaryVertexID(std::unordered_map<long, std::set<int> > & map);
    
    livector1D  extractPIDCells(short);
    livector1D  extractPIDCells(shivector1D);

    livector1D  getMapData();
    liimap      getMapDataInv();
    livector1D  getMapCell();
    liimap      getMapCellInv();

    void        getBoundingBox(std::array<double,3> & pmin, std::array<double,3> & pmax);
    BITPIT_DEPRECATED(void        buildBvTree(int value = 1));
    void        buildSkdTree(int value = 1);
    void        buildKdTree();
    void        buildAdjacencies();
    void        buildInterfaces();

    bool        areAdjacenciesBuilt();
    bool        areInterfacesBuilt();
    bool        isClosedLoop();

    bitpit::ElementType desumeElement(const livector1D &);

    void        dump(std::ostream & stream);
    void        restore(std::istream & stream);

    void   evalCellVolumes(bitpit::PiercedVector<double> &);
    void   evalCellAspectRatio(bitpit::PiercedVector<double> &);
    
    double evalCellVolume(const long & id);
    double evalCellAspectRatio(const long & id);
    
    void  getCellsNarrowBandToExtSurface(MimmoObject & surface, const double & maxdist, livector1D & idList);
    void  getCellsNarrowBandToExtSurface(MimmoObject & surface, const double & maxdist, bitpit::PiercedVector<double> & distList);
    
    void  getVerticesNarrowBandToExtSurface(MimmoObject & surface, const double & maxdist, livector1D & idList);
    void  getVerticesNarrowBandToExtSurface(MimmoObject & surface, const double & maxdist, bitpit::PiercedVector<double> & distList);
    
    std::unordered_map<long,long>   getInverseConnectivity();
    std::set<long>                  findVertexVertexOneRing(const long &, const long & );
    
protected:
    void    swap(MimmoObject & ) noexcept;
    void    reset(int type);
    
    std::unordered_set<int> elementsMap(bitpit::PatchKernel & obj);
    
private:
    bool    checkCellConnCoherence(const bitpit::ElementType & type, const livector1D & conn_);
    void    cleanKdTree();
    
};

};

#endif /* __MIMMOOBJECT_HPP__ */
