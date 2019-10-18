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

#include "mimmoTypeDef.hpp"
#include <bitpit_volunstructured.hpp>
#include <bitpit_surfunstructured.hpp>
#include <bitpit_SA.hpp>
#include <surface_skd_tree.hpp>
#include <volume_skd_tree.hpp>
#if MIMMO_ENABLE_MPI==1
#	include <mpi.h>
#endif

namespace mimmo{

/*!
 * \class MimmoSurfUnstructured
   \ingroup core
 * \brief Custom derivation of bitpit::SurfUnstructured class
 */
class MimmoSurfUnstructured: public bitpit::SurfUnstructured{
public:
    // Constructors
    MimmoSurfUnstructured();
    MimmoSurfUnstructured(int patch_dim);
    MimmoSurfUnstructured(int id, int patch_dim);
    MimmoSurfUnstructured(std::istream &stream);
    virtual ~MimmoSurfUnstructured();
    // Clone
    std::unique_ptr<bitpit::PatchKernel> clone() const override;

protected:
    /*! Copy Constructor*/
    MimmoSurfUnstructured(const MimmoSurfUnstructured & other) = default;
};

/*!
 * \class MimmoVolUnstructured
   \ingroup core
 * \brief Custom derivation of bitpit::VolUnstructured class
 */
class MimmoVolUnstructured: public bitpit::VolUnstructured{
public:
    // Constructors
    MimmoVolUnstructured(int dimension);
    MimmoVolUnstructured(int id, int dimension);
    virtual ~MimmoVolUnstructured();
    // Clone
    std::unique_ptr<bitpit::PatchKernel> clone() const override;

protected:
    /*! Copy Constructor*/
    MimmoVolUnstructured(const MimmoVolUnstructured &) = default;
};

/*!
 * \class MimmoPointCloud
 * \ingroup core
   \brief Custom derivation of bitpit::SurfUnstructured class, for Point Cloud handling only
 */
class MimmoPointCloud: public bitpit::SurfUnstructured{
public:
    // Constructors
    MimmoPointCloud();
    MimmoPointCloud(int id);
    virtual ~MimmoPointCloud();
    // Clone
    std::unique_ptr<bitpit::PatchKernel> clone() const override;


protected:
    /*! Copy Constructor*/
    MimmoPointCloud(const MimmoPointCloud &) = default;
};


/*!
* \class MimmoObject
  \ingroup core
* \brief MimmoObject is the basic geometry container for mimmo library
*
* MimmoObject is the basic container for 3D unustructured mesh data structure,
  based on bitpit::PatchKernel containers. MimmoObject can handle unustructured surface meshes,
  unustructured volume meshes, 3D point clouds and 3D tessellated curves.
  It supports interface methods to explore and handle the geometrical structure.
  It supports PID convention to mark subparts of geometry as well as building the search-trees
  KdTree (3D point spatial ordering) and skdTree(Cell-AABB spatial ordering) to quickly retrieve
  vertices and cells in the data structure.
*/
class MimmoObject{

private:
    std::unique_ptr<bitpit::PatchKernel>    m_patch;           /**<Reference to INTERNAL bitpit patch handling geometry. */
    bitpit::PatchKernel *                   m_extpatch;        /**<Reference to EXTERNALLY linked patch handling geometry. */
    bool                                    m_internalPatch;   /**<True if the geometry is internally created. */

protected:
//members
    int                                                     m_type;            /**< Type of geometry (0 = undefined, 1 = surface mesh, 2 = volume mesh, 3-point cloud mesh, 4-3DCurve). */
    std::unordered_set<long>                                m_pidsType;        /**< pid type available for your geometry */
    std::unordered_map<long, std::string>                   m_pidsTypeWNames;  /**< pid type available for your geometry, with name attached */
    std::unique_ptr<bitpit::PatchSkdTree>                   m_skdTree;         /**< ordered tree of geometry simplicies for fast searching purposes */
    std::unique_ptr<bitpit::KdTree<3,bitpit::Vertex,long> > m_kdTree;          /**< ordered tree of geometry vertices for fast searching purposes */
    bool                                                    m_skdTreeSync;     /**< track correct building of bvtree. Set false if any geometry modifications occur */
    bool                                                    m_kdTreeSync;      /**< track correct building of kdtree. Set false if any geometry modifications occur*/
    bool                                                    m_skdTreeSupported;/**< Flag for geometries not supporting skdTree building*/

    bool                                                    m_AdjBuilt;     /**< track correct building of adjacencies along with geometry modifications */
    bool                                                    m_IntBuilt;     /**< track correct building of interfaces  along with geometry modifications */
    bitpit::Logger*                                         m_log;          /**< Pointer to logger.*/

    int							m_nprocs;									/**< Total number of processors.*/
    int							m_rank;										/**< Current rank number.*/

#if MIMMO_ENABLE_MPI
	MPI_Comm 					m_communicator; 							/**< MPI communicator.*/
	long						m_ninteriorvertices = 0;					/**< Number of interior vertices.*/
	long						m_nglobalvertices = 0;						/**< Global number of vertices.*/
	std::vector<long>			m_rankinteriorvertices;						/**< Number of interior vertices for each rank.*/
	long						m_globaloffset;								/**< Points offset for local partition.*/
	std::unordered_map<int, std::vector<long>> m_pointGhostExchangeTargets;	/**< List of Ids of the local ghost points that are local points for each other processor.*/
	std::unordered_map<int, std::vector<long>> m_pointGhostExchangeSources;	/**< List of Ids of the local points that are ghost points for each other processor.*/
	std::unordered_map<int, std::vector<long>> m_pointGhostExchangeShared;	/**< List of Ids of the local points that are shared with each other processor.*/
	std::unordered_map<long, bool> m_isPointInterior;						/**< True or False if the current rank is considered the real owner (the lower rank for shared points) or not of the id-th (key) point. */
	std::unordered_map<long, long> m_pointConsecutiveId;					/**< Map id->consecutive id for vertices. */
	bool						m_pointGhostExchangeInfoSync;				/**< Track correct building of point ghost exchange info along with geometry modifications */

#endif

 	bitpit::PatchNumberingInfo	m_patchInfo;			/**< Patch Numbering Info structure for cells.*/
    bool                        m_infoSync;				/**< Track correct building of patch info along with geometry modifications */

    std::unordered_map<long, std::unordered_set<long> >	m_pointConnectivity;		/**< Point-Point connectivity. 1-Ring neighbours of each vertex.*/
    bool                        						m_pointConnectivitySync;	/**< Track correct building of points connectivity along with geometry modifications */

public:
    MimmoObject(int type = 1);
    MimmoObject(int type, dvecarr3E & vertex, livector2D * connectivity = NULL);
    MimmoObject(int type, bitpit::PatchKernel* geometry);
    MimmoObject(int type, std::unique_ptr<bitpit::PatchKernel> & geometry);
    ~MimmoObject();

    bitpit::Logger& getLog();

    bool                                            isEmpty();
    bool                                            isSkdTreeSupported();
    int                                             getType();
    long                                            getNVertices();
    long                                            getNCells()const;
    long                                            getNInternals()const;
    long                                            getNInternalVertices();
#if MIMMO_ENABLE_MPI
    long                                            getNGlobalVertices();
    long                                            getNGlobalCells();
    long                                            getPointGlobalCountOffset();
#endif
    dvecarr3E                                       getVerticesCoords(liimap* mapDataInv = NULL);
    darray3E                                        getVertexCoords(long i) const;
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
    std::unordered_set<long> &                      getPIDTypeList();
    std::unordered_map<long, std::string> &         getPIDTypeListWNames();
    livector1D                                      getCompactPID();
    std::unordered_map<long, long>                  getPID();

    bitpit::PatchSkdTree*                           getSkdTree();
    bitpit::KdTree<3, bitpit::Vertex, long> *       getKdTree();
    bitpit::PatchNumberingInfo*                     getPatchInfo();
    bool                          isSkdTreeSync();
    bool                          isKdTreeSync();
    bool                          isInfoSync();

    int getRank() const;
	int getProcessorCount() const;

    bool isPointInterior(long id);
#if MIMMO_ENABLE_MPI
	const MPI_Comm & getCommunicator() const;
    const std::unordered_map<int, std::vector<long>> & getPointGhostExchangeSources() const;
    const std::unordered_map<int, std::vector<long>> & getPointGhostExchangeTargets() const;
    bool arePointGhostExchangeInfoSync();
    void updatePointGhostExchangeInfo();
    void resetPointGhostExchangeInfo();
    bool cleanParallelSkdTreeSync();
    bool cleanParallelKdTreeSync();
    bool cleanParallelAdjacenciesSync();
    bool cleanParallelInterfacesSync();
    bool cleanParallelPointConnectivitySync();
    bool cleanParallelInfoSync();
    bool cleanParallelPointGhostExchangeInfoSync();
    void setPartitioned();
    void deleteOrphanGhostCells();
#endif

    bool        addVertex(const darray3E & vertex, long idtag = bitpit::Vertex::NULL_ID);
    bool        addVertex(const bitpit::Vertex & vertex, long idtag = bitpit::Vertex::NULL_ID);
    bool        modifyVertex(const darray3E & vertex, const long & id);

    bool        addConnectedCell(const livector1D & locConn, bitpit::ElementType type, int rank = -1);
    bool        addConnectedCell(const livector1D & locConn, bitpit::ElementType type, long idtag, int rank = -1);
    bool        addConnectedCell(const livector1D & locConn, bitpit::ElementType type, long PID, long idtag, int rank = -1);
    bool        addCell(bitpit::Cell & cell, int rank = -1);
    bool        addCell(bitpit::Cell & cell, long idtag, int rank = -1);

    void        setPID(livector1D );
    void        setPID(std::unordered_map<long, long>  );
    void        setPIDCell(long, long);
    bool        setPIDName(long, const std::string &);
    void        resyncPID();

    std::unique_ptr<MimmoObject>	clone() const ;
    void                            swap(MimmoObject & ) noexcept;

    bool        cleanGeometry();

    livector1D  getVertexFromCellList(const livector1D &cellList);
    livector1D  getCellFromVertexList(const livector1D &vertList, bool strict = true);
    livector1D  getInterfaceFromCellList(const livector1D &cellList, bool all = true);
    livector1D  getInterfaceFromVertexList(const livector1D &vertList, bool strict, bool border);

    livector1D                               extractBoundaryCellID(bool ghost = false);
    std::unordered_map<long, std::set<int> > extractBoundaryFaceCellID(bool ghost = false);
    livector1D                               extractBoundaryInterfaceID(bool ghost= false);
    livector1D                               extractBoundaryVertexID(bool ghost = false);
    livector1D                               extractBoundaryVertexID(std::unordered_map<long, std::set<int> > & map);

    livector1D  extractPIDCells(long, bool squeeze = true);
    livector1D  extractPIDCells(livector1D, bool squeeze = true);

    liimap      getMapData(bool withghosts=false);
    liimap      getMapDataInv(bool withghosts=true);
    liimap	    getMapCell(bool withghosts=true);
    liimap      getMapCellInv(bool withghosts=true);

    void        getBoundingBox(std::array<double,3> & pmin, std::array<double,3> & pmax);
    void        buildSkdTree(int value = 1);
    void        buildKdTree();
    void		buildPatchInfo();
	void        buildAdjacencies();
    void        buildInterfaces();
	void        resetAdjacencies();
    void        resetInterfaces();
    void		resetPatch();
    void    	cleanKdTree();
    void    	cleanSkdTree();

    bool        areAdjacenciesBuilt();
    bool        areInterfacesBuilt();
    bool        isClosedLoop();

    std::vector<std::vector<long> > decomposeLoop();

    bitpit::ElementType desumeElement(const livector1D &);

    void        dump(std::ostream & stream);
    void        restore(std::istream & stream);

    void   evalCellVolumes(bitpit::PiercedVector<double> &);
    void   evalCellAspectRatio(bitpit::PiercedVector<double> &);

    double evalCellVolume(const long & id);
    double evalCellAspectRatio(const long & id);

    std::array<double,3> evalCellCentroid(const long & id);
    std::array<double,3> evalInterfaceCentroid(const long & id);
    double               evalInterfaceArea(const long & id);
    std::array<double,3> evalInterfaceNormal(const long & id);

    livector1D                      getCellsNarrowBandToExtSurface(MimmoObject & surface,
                                                                   const double & maxdist,
                                                                   livector1D * seedList = nullptr);
    bitpit::PiercedVector<double>   getCellsNarrowBandToExtSurfaceWDist(MimmoObject & surface,
                                                                        const double & maxdist,
                                                                        livector1D * seedList = nullptr);


    std::unordered_map<long,long>   getInverseConnectivity();
    std::set<long>                  findVertexVertexOneRing(const long &, const long & );

    void						buildPointConnectivity();
    void						cleanPointConnectivity();
    std::unordered_set<long> &	getPointConnectivity(const long & id);
    bool						isPointConnectivitySync();

    void						triangulate();

protected:
    void    reset(int type);

    std::unordered_set<int> elementsMap(bitpit::PatchKernel & obj);

#if MIMMO_ENABLE_MPI
    void	initializeParallel();
#endif

private:
    //make copy constructor and  assignment private and not accessible (use properly clone method for copy)
    MimmoObject(const MimmoObject & other);
    MimmoObject & operator=(MimmoObject other);

    bool    checkCellConnCoherence(const bitpit::ElementType & type, const livector1D & conn_);

	/*!
        \struct VertexPositionLess
        Functional for comparing the position of two vertices.
		The comparison is made with respect to the vertex coordinates.
	*/
	struct VertexPositionLess
	{
        /*! Base constructor
            \param[in] object MimmoObject class required
        */
        VertexPositionLess(const MimmoObject &object)
			: m_object(object)
		{
		}
        /*! Base destructor */
		virtual ~VertexPositionLess() = default;

        /*!
            Compare the position coordinates of 2 MimmoObject Vertex and return true
            if the first non-coincident coordinate (x,y,z in this order) of the first vertex is
            lesser than the second vertex one.
            If all coordinates are coincident and no comparison is performed throw a runtime_error.
            \param[in] id_1 global id of first vertex
            \param[in] id_2 global id of second vertex
            \return boolean
        */
		bool operator()(const long &id_1, const long &id_2) const
		{
			std::array<double, 3> vertex_1 = m_object.getVertexCoords(id_1);
			std::array<double, 3> vertex_2 = m_object.getVertexCoords(id_2);

			for (int k = 0; k < 3; ++k) {
				if (std::abs(vertex_1[k] - vertex_2[k]) <= 1.0E-12) {
					continue;
				}
				return vertex_1[k] < vertex_2[k];
			}

			// If we are here the two vertices coincide. It's not
			// possible to define an order for the two vertices.
			std::ostringstream stream;
			stream << "It was not possible to define an order for vertices " << id_1 << " and " << id_2 << ". ";
			stream << "The two vertices have the same coordinates.";
			throw std::runtime_error (stream.str());
		}

		const MimmoObject &m_object;
	};

};

};

#endif /* __MIMMOOBJECT_HPP__ */
