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
#include "BvTree.hpp"
#include "mimmoTypeDef.hpp"

namespace mimmo{

/*!
 * \ingroup core
 */

/*!
* \class MimmoObject
* \brief MimmoObject is the basic geometry container for mimmo library
*
* MimmoObject is the basic container for 3D geometry unustructured mesh data structure, 
* based on bitpit::PatchKernel containers.
* MimmoObject can handle unustructured surface meshes, unustructured volume meshes, 3D point clouds and 3D tessellated curves.
* It supports interface methods to explore and handle the geometrical structure. It supports PID convention to mark subparts 
* of geometry as well as building the search-trees KdTree and BvTree to quickly retrieve vertices and cells in the data structure.
*/
class MimmoObject{

protected:
//members
    int                                     m_type;            /**<Type of geometry (0 = undefined, 1 = surface mesh, 2 = volume mesh, 3-point cloud mesh, 4-3DCurve). */
    bitpit::PatchKernel*                    m_patch;           /**<Reference to bitpit patch handling geometry. */
    bool                                    m_internalPatch;   /**<True if the geometry is internally created. */
    livector1D                              m_mapData;         /**<Map of vertex ids actually set, for aligning external vertex data to bitpit::PatchKernel ordering */
    livector1D                              m_mapCell;         /**<Map of cell ids actually set, for aligning external cell data to bitpit::PatchKernel ordering*/ 
    liimap                                  m_mapDataInv;      /**<Inverse of Map of vertex ids actually set, for aligning external vertex data to bitpit::Patch ordering */
    liimap                                  m_mapCellInv;     /**<Inverse of Map of cell ids actually set, for aligning external vertex data to bitpit::Patch ordering */

    std::unordered_set<short>               m_pidsType;        /**<pid type available for your geometry */

    BvTree                                  m_bvTree;          /**< ordered tree of geometry simplicies for fast searching purposes */
    bool                                    m_bvTreeBuilt;     /**< track correct building of bvtree along with geometry modifications */
    bitpit::KdTree<3,bitpit::Vertex,long>   m_kdTree;          /**< ordered tree of geometry vertices for fast searching purposes */
    bool                                    m_kdTreeBuilt;     /**< track correct building of kdtree along eith geometry modifications */
    bool                                    m_bvTreeSupported; /**< Flag for geometries not supporting bvTree building*/

    //TODO temporary flags for modifying vertices -> ambigous effect on bvTree, kdTree
    bool                                    m_bvTreeSync;    /**< set false if Bv tree is not sync'd with geometry modifications */
    bool                                    m_kdTreeSync;    /**< set false if Bv tree is not sync'd with geometry modifications */

    bool                                    m_AdjBuilt;     /**< track correct building of adjacencies along with geometry modifications */


public:
    MimmoObject(int type = 1); 
    MimmoObject(int type, dvecarr3E & vertex, ivector2D * connectivity = NULL);
    MimmoObject(int type, bitpit::PatchKernel* geometry);
    ~MimmoObject();

    MimmoObject(const MimmoObject & other);
    MimmoObject & operator=(const MimmoObject & other);

    void                                            clear();
    bool                                            isEmpty();
    bool                                            isBvTreeSupported();
    int                                             getType();
    long                                            getNVertex()const;
    long                                            getNCells()const;
    dvecarr3E                                       getVertexCoords();
    darray3E                                        getVertexCoords(long i);
    bitpit::PiercedVector<bitpit::Vertex> &         getVertices();
    const bitpit::PiercedVector<bitpit::Vertex> &   getVertices() const ;

    ivector2D                                       getCompactConnectivity();
    livector2D                                      getConnectivity();
    livector1D                                      getCellConnectivity(long id);
    bitpit::PiercedVector<bitpit::Cell> &           getCells();
    const bitpit::PiercedVector<bitpit::Cell> &     getCells() const;

    bitpit::PatchKernel*                            getPatch();

    livector1D&                                     getMapData();
    long                                            getMapData(int i);
    liimap&                                         getMapDataInv();
    const liimap&                                   getMapDataInv()const;
    int                                             getMapDataInv(long id);
    livector1D&                                     getMapCell();
    long                                            getMapCell(int i);
    liimap&                                         getMapCellInv();
    int                                             getMapCellInv(long id);

    std::unordered_set<short> &                     getPIDTypeList();
    shivector1D                                     getCompactPID();
    std::unordered_map<long, short>                 getPID();

    bool                                            isBvTreeBuilt();
    BvTree*                                         getBvTree();
    bool                                            isKdTreeBuilt();
    bitpit::KdTree<3, bitpit::Vertex, long> *       getKdTree();
    bool                                            isBvTreeSync();
    bool                                            isKdTreeSync();

    const MimmoObject *                             getCopy();

    bool        setVertices(const bitpit::PiercedVector<bitpit::Vertex> & vertices);
    bool        addVertex(const darray3E & vertex, const long idtag = bitpit::Vertex::NULL_ID);
    bool        modifyVertex(const darray3E & vertex, long id);
    bool        setCells(const bitpit::PiercedVector<bitpit::Cell> & cells);
    bool        addConnectedCell(const livector1D & locConn, bitpit::ElementInfo::Type type, long idtag = bitpit::Cell::NULL_ID);
    bool        addConnectedCell(const livector1D & locConn, bitpit::ElementInfo::Type type, short PID, long idtag = bitpit::Cell::NULL_ID);

    bool        setPatch(int type, bitpit::PatchKernel* geometry);
    bool        setMapData();
    bool        setMapCell();

    void        setPID(shivector1D ); 
    void        setPID(std::unordered_map<long, short>  ); 
    void        setPIDCell(long, short);

    void        setSOFTCopy(const MimmoObject * other);
    void        setHARDCopy(const MimmoObject * other);

    bool        cleanGeometry();

    livector1D  getVertexFromCellList(livector1D cellList);
    livector1D  getCellFromVertexList(livector1D vertList);
    livector1D  convertLocalToVertexID(ivector1D);
    ivector1D   convertVertexIDToLocal(livector1D);
    livector1D  convertLocalToCellID(ivector1D);
    ivector1D   convertCellIDToLocal(livector1D);

    livector1D  extractBoundaryVertexID();
    livector1D  extractPIDCells(short);
    livector1D  extractPIDCells(shivector1D);

    void        getBoundingBox(std::array<double,3> & pmin, std::array<double,3> & pmax);
    void        buildBvTree(int value = 1);
    void        buildKdTree();
    void        updateBvTree();
    void        updateKdTree();
    void        buildAdjacencies();

    bool        areAdjacenciesBuilt();
    bool        isClosedLoop();
    
    bitpit::VTKElementType	desumeElement();

private:
    int     checkCellType(bitpit::ElementInfo::Type type);
    void    cleanKdTree();
};

};

#endif /* __MIMMOOBJECT_HPP__ */
