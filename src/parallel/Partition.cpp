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
#include <Partition.hpp>
#include <metis.h>
#include <bitpit_operators.hpp>
#include <SkdTreeUtils.hpp>
#include <chrono>

typedef std::chrono::high_resolution_clock Clock;

namespace mimmo{

/*!Default constructor of Partition
 */
Partition::Partition(){
	m_name = "mimmo.Partition";
	m_mode = PartitionMethod::NONE;
	m_partition.clear();
	m_boundary.reset();
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
Partition::Partition(const bitpit::Config::Section & rootXML){

	m_name = "mimmo.Partition";
	m_mode = PartitionMethod::NONE;
	m_boundary.reset();
	m_partition.clear();

	std::string fallback_name = "ClassNONE";
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::string::trim(input);
	if(input == "mimmo.Partition"){
		absorbSectionXML(rootXML);
	}else{
		warningXML(m_log, m_name);
	};
}

/*!Default destructor of Partition
 */
Partition::~Partition(){
};

/*!Copy constructor of Partition.
 */
Partition::Partition(const Partition & other):BaseManipulation(other){
	m_mode = other.m_mode;
	m_partition = other.m_partition;
	m_boundary = other.m_boundary;
};

/*! It builds the input/output ports of the object
 */
void
Partition::buildPorts(){

	bool built = true;

	built = (built && createPortIn<std::unordered_map<long, int>, Partition>(this, &mimmo::Partition::setPartition, M_UMAPI));
	built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, Partition>(this, &mimmo::Partition::setGeometry, M_GEOM, true));
	built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, Partition>(this, &mimmo::Partition::setBoundaryGeometry, M_GEOM2));

	built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, Partition>(this, &mimmo::Partition::getGeometry, M_GEOM));
	built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, Partition>(this, &mimmo::Partition::getBoundaryGeometry, M_GEOM2));
	m_arePortsBuilt = built;
};

/*!
 * Get boundary geometry
 */
MimmoSharedPointer<MimmoObject>
Partition::getBoundaryGeometry(){
    return m_boundary;
}

/*!
 * Set geometry. Overload of BaseManipulation method.
 */
void
Partition::setGeometry(MimmoSharedPointer<MimmoObject> geo){
    if(geo == nullptr)    return;
    if(m_geometry == geo) return;
    m_geometry = geo;
}

/*!
 * Set boundary geometry
 */
void
Partition::setBoundaryGeometry(MimmoSharedPointer<MimmoObject> geo){
    if(geo == nullptr)    return;
    if(m_boundary == geo) return;
	m_boundary = geo;
}

/*!
 * It sets partition structure of partition block.
 * It sets the partition mode to PartitionMethod::CUSTOM.
 * \param[in] partition partition vector
 *
 */
void
Partition::setPartition(std::unordered_map<long, int> partition){
	m_partition = partition;
	m_mode = PartitionMethod::CUSTOM;
};

/*!
 * It sets partition method of partition block
 * \param[in] mode partition method
 *
 */
void
Partition::setPartitionMethod(PartitionMethod mode){
	if (mode != PartitionMethod::PARTGEOM && mode != PartitionMethod::SERIALIZE &&
	        mode != PartitionMethod::CUSTOM && mode != PartitionMethod::NONE)
	    throw std::runtime_error(m_name + " : partition method not allowed");

	m_mode = mode;
};

/*!
 * It sets partition method of partition block
 * \param[in] mode partition method
 *
 */
void
Partition::setPartitionMethod(int mode){
	if (mode != -1 && mode != 0 && mode != 1 && mode != 2)
		throw std::runtime_error(m_name + " : partition method not allowed");

	m_mode = PartitionMethod(mode);
};

/*!
 * Execution command. Clip geometry and save result in m_patch member.
 * Note: in case of seralization the older partitioned geometry is reset to empty geometry. Please consider to update pointers (use shared?) to the new serialized geometry.
 */
void
Partition::execute(){

    if(m_mode == PartitionMethod::NONE){
        // Do nothing
        return;
    }

    if( getGeometry()== nullptr){
        (*m_log)<<m_name + " : null pointer to linked geometry found."<<std::endl;
        return;
    };

    if(m_mode == PartitionMethod::PARTGEOM && getGeometry()->isDistributed()){
        (*m_log)<<m_name + " : already distributed geometry found during geometric partition."<<std::endl;
        throw std::runtime_error(m_name + " : already distributed geometry found during geometric partition.");
    };

    if (m_mode == PartitionMethod::SERIALIZE && !(getGeometry()->isDistributed())){
        (*m_log)<<m_name + " : already serialized geometry found during geometric partition."<<std::endl;
        return;
    };

    if(m_mode == PartitionMethod::CUSTOM){
        if (m_partition.size() != m_geometry->getNCells())
            throw std::runtime_error(m_name + " : partition size different from number of cells");
    }

    if (m_mode == PartitionMethod::PARTGEOM || m_mode == PartitionMethod::SERIALIZE || m_mode == PartitionMethod::CUSTOM)
    {

        // Force build adjacencies if not built
        // TODO Reset adjacencies if not already computed?
        getGeometry()->updateAdjacencies();

        // Compute partition if not custom
        if(m_mode != PartitionMethod::CUSTOM){
            computePartition();
        }

        // Compute boundary partition if boundary geometry exists
        if (getBoundaryGeometry() != nullptr){
            if (getGeometry()->getType() == 2 && getBoundaryGeometry()->getType() == 1){
                // Compute boundary partition
                computeBoundaryPartition();
            }
        }

        // Clean structures to be destroyed/reset
        getGeometry()->cleanPointConnectivity();
        getGeometry()->cleanPatchInfo();
        getGeometry()->cleanSkdTree();
        getGeometry()->cleanKdTree();
        getGeometry()->cleanBoundingBox();
#if MIMMO_ENABLE_MPI
        getGeometry()->resetPointGhostExchangeInfo();
#endif

        // Partition/Serialize the geometry
        getGeometry()->getPatch()->partition(m_partition, false, true);
        if (m_mode == PartitionMethod::SERIALIZE){
            // Sort cells and vertices with Id
            //getGeometry()->getPatch()->sortCells();
            getGeometry()->getPatch()->sortVertices();
        }

        // Adjacencies Sync and Interfaces Sync not changed by partition,
        // the bitpit partitioning maintains adjacencies and
        // interfaces if already synchronized

        // Update geometry
        getGeometry()->update();

        // Resync PID
        //TODO Add resyncPID to update method of mimmoobject with syncstatus of pid
        getGeometry()->resyncPID();

        // Partition boundary geometry
        if (getBoundaryGeometry() != nullptr){
            if (getGeometry()->getType() == 2 && getBoundaryGeometry()->getType() == 1){

                // Force build adjacencies
                getBoundaryGeometry()->updateAdjacencies();

                // Clean structures to be destroyed
                getBoundaryGeometry()->cleanPointConnectivity();
                getBoundaryGeometry()->cleanPatchInfo();
                getBoundaryGeometry()->cleanSkdTree();
                getBoundaryGeometry()->cleanKdTree();
                getBoundaryGeometry()->cleanBoundingBox();
        #if MIMMO_ENABLE_MPI
                getBoundaryGeometry()->resetPointGhostExchangeInfo();
        #endif

                // Boundary partition
                getBoundaryGeometry()->getPatch()->partition(m_boundarypartition, false, true);
                if (m_mode == PartitionMethod::SERIALIZE){
                    // Sort cells and vertices with Id
                    //getBoundaryGeometry()->getPatch()->sortCells();
                    getBoundaryGeometry()->getPatch()->sortVertices();
                }

                // Adjacencies Sync and Interfaces Sync not changed by partition,
                // the bitpit partitioning maintains adjacencies and
                // interfaces if already synchronized

                // Update ID of boundary vertices
                //NO UPDATE BECAUSE BITPIT ASSIGN THE SAME IDS DURING PARTITION (FROM SERIAL TO PARTITION IS TRUE AND VICEVERSA IF NOT MODIFIED)
                // updateBoundaryVerticesID();

                // Update boundary geometry
                getBoundaryGeometry()->update();

                // Resync PID
                getBoundaryGeometry()->resyncPID();

            }
        }
    } // end if partition mode
};

/*!
 * It computes the partition structure by using the chosen method.
 * The partition structure if empty is filled after the call.
 */
void
Partition::computePartition(){
	switch(m_mode) {
	case PartitionMethod::SERIALIZE:
		serialPartition();
		break;
	case PartitionMethod::PARTGEOM:
		parmetisPartGeom();
		break;
	default:
		break;
	}
};

/*!
 * It computes the partition structure by using the space filling method of parmetis.
 */
void
Partition::parmetisPartGeom(){

	if (!(getGeometry()->isDistributed())){

		lilimap mapcell = getGeometry()->getMapCell();
		lilimap mapcellinv = getGeometry()->getMapCellInv();

		if (m_rank == 0){

			// Use cell centers as vertices of graph

			//
			//  The number of vertices.
			//
			idx_t nvtxs = idx_t(getGeometry()->getNCells());

			//
			// Number of balancing constraints, which must be at least 1.
			//
			idx_t ncon = 1;


			// Build Adjacencies
			std::vector<idx_t> xadj(nvtxs+1);

			idx_t duenedges = 0;
			for (int i=0; i<getGeometry()->getNInternals(); i++){
				duenedges += getGeometry()->getPatch()->getCell(mapcell[i]).getAdjacencyCount();
			}

			std::vector<idx_t> adjncy(duenedges);

			idx_t i;
			idx_t j=0;
			std::vector<long> neighs;
			for (i=0; i<getGeometry()->getNInternals(); i++){
				xadj[i] = j;
				neighs.clear();
				getGeometry()->getPatch()->findCellFaceNeighs(mapcell[i], &neighs);
				for (long id : neighs){
					if (id > 0){
						adjncy[j] = idx_t(mapcellinv[id]);
						j++;
					}
				}
			}
			xadj[i] = j;

			idx_t nParts = m_nprocs;

			//
			//  On return, the edge cut volume of the partitioning solution.
			//
			idx_t objval;

			//
			//  On return, the partition vector for the graph.
			//
			std::vector<idx_t> part(nvtxs);

			int ret = METIS_PartGraphKway ( &nvtxs, &ncon, xadj.data(), adjncy.data(), nullptr, nullptr,
					nullptr, &nParts, nullptr, nullptr, nullptr, &objval, part.data() );

			m_partition.clear();
			m_partition.reserve(nvtxs);
			for (long i=0; i<nvtxs; i++){
				m_partition[mapcell[i]] = part[i];
			}
		}
	} // end if geometry is not distributed
}

/*!
 * It computes the partition structure to communicate the entire mesh to processor with rank 0.
 */
void
Partition::serialPartition(){

	if ((getGeometry()->isDistributed())){

		// Build partition to send to zero
		long ncells = getGeometry()->getPatch()->getInternalCount();
		std::unordered_map<long, int>().swap(m_partition);
		m_partition.reserve(ncells);
		for (long id : getGeometry()->getPatch()->getCells().getIds()){
			m_partition[id] = 0;
		}

	} // end if geometry is distributed
}

/*!
 * It computes the boundary path partition coherently with the volume partition.
 * The vertices are supposed with corresponding IDs between volume and surface meshes.
 */
void
Partition::computeBoundaryPartition()
{

	if (m_mode == PartitionMethod::SERIALIZE){

		if (getBoundaryGeometry()->isDistributed()){

			// Build partition to send to zero
			long ncells = getBoundaryGeometry()->getPatch()->getInternalCount();
			m_boundarypartition.clear();
			m_boundarypartition.reserve(ncells);
			for (long id : getBoundaryGeometry()->getPatch()->getCells().getIds()){
				m_boundarypartition[id] = 0;
			}

		} // end if geometry is distributed
	}
	else{

		if (!(getBoundaryGeometry()->isDistributed())){

		    std::unordered_map<long, int>().swap(m_boundarypartition);

		    // Build Skd Tree to find correspondence between surface and volume elements
		    getBoundaryGeometry()->buildSkdTree();

            if (m_rank == 0){

			    m_boundarypartition.reserve(getBoundaryGeometry()->getNCells());
                bitpit::PatchSkdTree *btree = getBoundaryGeometry()->getSkdTree();
				double tol = 1.0E-12;
                std::vector<double> distances(getBoundaryGeometry()->getNCells(), std::numeric_limits<double>::max());

                bitpit::PiercedVector<bitpit::Cell> & cells = getGeometry()->getCells();
                for (bitpit::Cell &cell : cells){
                	if (cell.isInterior()){
                		for (int iface=0; iface<cell.getFaceCount(); iface++){
                			if (cell.isFaceBorder(iface)){
                				std::array<double,3> intercenter({0.,0.,0.});
                				bitpit::ConstProxyVector<long> vertices = cell.getFaceVertexIds(iface);
                				int nV = cell.getFaceVertexCount(iface);
                				for (int iv=0; iv<nV; iv++){
                					intercenter += getGeometry()->getPatch()->getVertexCoords(vertices[iv]);
                				}
                				intercenter /= double(nV);
                				long id = bitpit::Cell::NULL_ID;
                				int maxiter(10), iter(0);
                				double r = tol;
                				double distance;
                				while(id == bitpit::Cell::NULL_ID && iter<maxiter){
                				    // Call local distance computing, the two meshes are on the same processor rank = 0
                					distance = skdTreeUtils::distance(&intercenter, btree, id, r);
                					++iter;
                					r = tol * std::pow(10, iter);
                				}
                				if(id != bitpit::Cell::NULL_ID){
                					long brawindex = getBoundaryGeometry()->getCells().getRawIndex(id);
                					if (distance < distances[brawindex]){
                						distances[brawindex] = distance;
                						m_boundarypartition[id] = m_partition.at(cell.getId());
                					}
                				}
                			}
                		}
                	}
                } // end loop on cells

				// Clean SkdTree
				getBoundaryGeometry()->cleanSkdTree();

			}

		} // end if geometry is not distributed

	} // end if mode is not serialize
}

/*!
 * Update ID of boundary vertices to be coherent with volume vertices
 */
void
Partition::updateBoundaryVerticesID()
{

    // TODO Use with the new version of bitpit with ghost vertices
    // TODO rework with update methods of mimmo

	if (m_rank != 0){

		// Compute map between old and new boundary vertices id
		getGeometry()->buildKdTree();
		bitpit::KdTree<3,bitpit::Vertex,long>  *tree = getGeometry()->getKdTree();
		std::unordered_map<long, long> old2new;
		bitpit::PiercedVector<bitpit::Vertex, long> & Vertices = getBoundaryGeometry()->getVertices();
		long idmax = 0;
		for (long id : Vertices.getIds())
			idmax = std::max(idmax, id);

		for (bitpit::Vertex & vertex : Vertices){
			long id = vertex.getId();
			std::vector<long> ids, excl;
			double tol = 1.0e-10;
			long idnew = -1;
			while (idnew < 0){
				tree->hNeighbors(&vertex,tol,&ids,nullptr);
				double distance = 1.0e+05;
				if (!ids.empty()){
					for (long id_ : ids){
						double dist = norm2(getGeometry()->getVertexCoords(id_)-vertex.getCoords());
						if (dist < distance){
							idnew = id_;
							distance = dist;
						}
					}
				}
				tol *= 2.;
			}
			long idnewpmax = idmax + idnew + 1;
			Vertices.updateId(id, idnewpmax);
			vertex.setId(idnewpmax);
//			bitpit::Vertex newVertex(idnewpmax, vertex.getCoords());
//			getBoundaryGeometry()->getPatch()->deleteVertex(id);
//			getBoundaryGeometry()->addVertex(newVertex,idnewpmax);

			old2new[id] = idnew;
		}

		Vertices.sort();
		for (bitpit::Vertex & vertex : Vertices){
			long id = vertex.getId();
			Vertices.updateId(id, id-idmax-1);
			vertex.setId(id-idmax-1);
//			bitpit::Vertex newVertex(id-idmax-1, vertex.getCoords());
//			getBoundaryGeometry()->getPatch()->deleteVertex(id);
//			getBoundaryGeometry()->addVertex(newVertex,id-idmax-1);
		}

		// Update connectivity of boundary cells
		int ic=0;
		for (bitpit::Cell & cell : getBoundaryGeometry()->getCells()){
			std::unique_ptr<long[]> connectStorage = std::unique_ptr<long[]>(new long[cell.getVertexCount()]);
			std::size_t i=0;
			for (long id : cell.getVertexIds()){
				connectStorage[i] = old2new[id];
				i++;
			}
			cell.setConnect(std::move(connectStorage));
			ic++;
		}
	}
}

/*!
 * It plots optional result of the class in execution,
 * that is the clipped geometry as standard vtk unstructured grid.
 */
void
Partition::plotOptionalResults(){

	write(getGeometry());
	std::string originalname = m_name;
	m_name = originalname + ".Boundary";
	write(getBoundaryGeometry());
	m_name = originalname;

}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
Partition::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

	/*start absorbing*/
	BITPIT_UNUSED(name);

	BaseManipulation::absorbSectionXML(slotXML,name);

	if(slotXML.hasOption("PartitionMethod")){
		std::string input = slotXML.get("PartitionMethod");
		input = bitpit::utils::string::trim(input);
		int value = 1;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setPartitionMethod(value);
	}
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
Partition::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);

	BaseManipulation::flushSectionXML(slotXML, name);

	int value = int(m_mode);
	slotXML.set("PartitionMethod", std::to_string(value));

};

}
