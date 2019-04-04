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

namespace mimmo{

/*!Default constructor of Partition
 */
Partition::Partition(){
	m_name = "mimmo.Partition";
	m_mode = PartitionMethod::PARTGEOM;
	m_partition.clear();
	m_boundary = nullptr;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
Partition::Partition(const bitpit::Config::Section & rootXML){

	m_name = "mimmo.Partition";
	m_mode = PartitionMethod::PARTGEOM;
	m_boundary = nullptr;
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
Partition::~Partition(){};

/*!Copy constructor of Partition.
 */
Partition::Partition(const Partition & other):BaseManipulation(other){
	m_mode = other.m_mode;
	m_partition = other.m_partition;
	m_boundary = other.m_boundary;
	;
};

/*! It builds the input/output ports of the object
 */
void
Partition::buildPorts(){

	bool built = true;

	built = (built && createPortIn<ivector1D, Partition>(this, &mimmo::Partition::setPartition, M_VECTORSI));
	built = (built && createPortIn<MimmoObject*, Partition>(this, &mimmo::Partition::setGeometry, M_GEOM, true));
	built = (built && createPortIn<MimmoObject*, Partition>(this, &mimmo::Partition::setBoundaryGeometry, M_GEOM2));

	built = (built && createPortOut<MimmoObject*, Partition>(this, &mimmo::Partition::getGeometry, M_GEOM));
	built = (built && createPortOut<MimmoObject*, Partition>(this, &mimmo::Partition::getBoundaryGeometry, M_GEOM2));
	m_arePortsBuilt = built;
};

/*!
 * Get geometry
 */
MimmoObject*
Partition::getGeometry(){
	return m_geometry;
}

/*!
 * Get boundary geometry
 */
MimmoObject*
Partition::getBoundaryGeometry(){
	return m_boundary;
}

/*!
 * Set boundary geometry
 */
void
Partition::setBoundaryGeometry(MimmoObject* geo){
	m_boundary = geo;
}

/*!
 * It sets partition structure of partition block
 * \param[in] partition partition vector
 *
 */
void
Partition::setPartition(ivector1D partition){
	m_partition = partition;
	if (m_partition.size() != m_geometry->getNCells())
		throw std::runtime_error(m_name + " : partition size different from number of cells");
};

/*!
 * It sets partition method of partition block
 * \param[in] mode partition method
 *
 */
void
Partition::setPartitionMethod(PartitionMethod mode){
	if (mode != PartitionMethod::PARTGEOM)
		throw std::runtime_error("Partition: partition method not allowed");

	m_mode = mode;
};

/*!
 * It sets partition method of partition block
 * \param[in] mode partition method
 *
 */
void
Partition::setPartitionMethod(int mode){
	if (mode != 1)
		throw std::runtime_error(m_name + " : partition method not allowed");

	m_mode = PartitionMethod(mode);
};

/*!
 * Execution command. Clip geometry and save result in m_patch member.
 */
void
Partition::execute(){

	if(getGeometry() == nullptr){
		(*m_log)<<m_name + " : null pointer to linked geometry found."<<std::endl;
		return;
	};
	if(getGeometry()->isEmpty()){
		(*m_log)<<m_name + " : empty linked geometry found."<<std::endl;
		return;
	};

	//Compute partition
	computePartition();

	if (getBoundaryGeometry() != nullptr){
		if (getGeometry()->getType() == 2 && getBoundaryGeometry()->getType() == 1){
			//Compute boundary partition
			computeBoundaryPartition();
		}
	}

	//partition
	std::vector<bitpit::adaption::Info> Vinfo = getGeometry()->getPatch()->partition(m_partition, true);

	//Force (Re)Build adjcencies and interfaces after partitioning
	getGeometry()->buildAdjacencies();

	//TODO Temporary reset to cancel after bitpit fixing
	getGeometry()->getPatch()->resetInterfaces();
	getGeometry()->buildInterfaces();

	//Force rebuild patch info
	getGeometry()->buildPatchInfo();

	if (getBoundaryGeometry() != nullptr){
		if (getGeometry()->getType() == 2 && getBoundaryGeometry()->getType() == 1){
			//boundary partition
			std::vector<bitpit::adaption::Info> Sinfo = getBoundaryGeometry()->getPatch()->partition(m_boundarypartition, true);

			//Update ID of boundary vertices
			updateBoundaryVerticesID();

			//Force rebuild adjacencies
			getBoundaryGeometry()->buildAdjacencies();

			//Force rebuild patch info
			getBoundaryGeometry()->buildPatchInfo();
		}
	}
};

/*!
 * It computes the partition structure by using the chosen method.
 * The partition structure if empty is filled after the call.
 */
void
Partition::computePartition(){
	switch(m_mode) {
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

	//Set communicator in case is not set
	if (!getGeometry()->getPatch()->isCommunicatorSet()){
		getGeometry()->getPatch()->setCommunicator(m_communicator);
	}

	if ((m_nprocs>1) && !(getGeometry()->getPatch()->isPartitioned())){
		if (m_rank != 0 ){
			//Reset vertices, cells and interfaces
			getGeometry()->getPatch()->reset();
		}

		getGeometry()->cleanGeometry();

		liimap mapcell = getGeometry()->getMapCellInv();

		if (m_rank == 0){

			//Renumber cell in consecutive order
			getGeometry()->getPatch()->consecutiveRenumberCells();

			//Build adjacencies to have ghost after partiioning
			for (bitpit::Cell & cell : getGeometry()->getCells())
				cell.resetAdjacencies();
			getGeometry()->buildAdjacencies();

			//Build interfaces to compute graph partitioning
			getGeometry()->getPatch()->resetInterfaces();
			getGeometry()->buildInterfaces();

			//USE CELLS CENTERS

			//
			//  The number of vertices.
			//
			idx_t nvtxs = idx_t(getGeometry()->getNCells());

			//
			// Number of balancing constraints, which must be at least 1.
			//
			idx_t ncon = 1;

			//Build Adjacencies
			std::vector<idx_t> xadj(nvtxs+1);
			long nint=0;
			for (auto inter : getGeometry()->getInterfaces()){
				if (!inter.isBorder())
					nint++;
			}
			idx_t nEdges = idx_t(nint);
			idx_t duenedges = 2*nEdges;

			std::vector<idx_t> adjncy(duenedges);

			idx_t i=0;
			idx_t j=0;
			for (bitpit::Cell & cell : getGeometry()->getCells()){
				xadj[i] = j;
				std::vector<long> neighs = getGeometry()->getPatch()->findCellFaceNeighs(cell.getId());
				for (long id : neighs){
					if (id > 0){
						adjncy[j] = idx_t(mapcell[id]);
						j++;
					}
				}
				i++;
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

			int ret = METIS_PartGraphKway ( &nvtxs, &ncon, xadj.data(), adjncy.data(), NULL, NULL,
					NULL, &nParts, NULL, NULL, NULL, &objval, part.data() );

			m_partition.resize(nvtxs);
			for (long i=0; i<nvtxs; i++){
				m_partition[i] = part[i];
			}
		}
	}

	MPI_Barrier(m_communicator);

}

/*!
 * It computes the boundary path partition coherently with the volume partition.
 * The vertices are supposed with corresponding IDs between volume and surface meshes.
 */
void
Partition::computeBoundaryPartition()
{
	if (!getBoundaryGeometry()->getPatch()->isCommunicatorSet()){
		getBoundaryGeometry()->getPatch()->setCommunicator(m_communicator);
	}

	if ((m_nprocs>1) && !(getBoundaryGeometry()->getPatch()->isPartitioned())){
		//Build adjacencies to have ghost after partiioning
		getBoundaryGeometry()->buildAdjacencies();
		m_boundarypartition.clear();
		if (m_rank == 0){
			m_boundarypartition.resize(getBoundaryGeometry()->getNCells());
			getBoundaryGeometry()->buildSkdTree();
			bitpit::PatchSkdTree *btree = getBoundaryGeometry()->getSkdTree();
			double tol = 1.0e-08;
			for (bitpit::Interface inter : getGeometry()->getInterfaces()){
				if (inter.isBorder()){
					std::array<double,3> intercenter = getGeometry()->getPatch()->evalInterfaceCentroid(inter.getId());
					long id;
					double r = tol;
					double distance = skdTreeUtils::distance(&intercenter, btree, id, r);
					if (distance <= tol){
						m_boundarypartition[getBoundaryGeometry()->getCells().getRawIndex(id)] = m_partition[getGeometry()->getCells().getRawIndex(inter.getOwner())];
					}
				}
			}
		}
		else{
			//Reset vertices, cells and interfaces
			getBoundaryGeometry()->getPatch()->reset();
		}
	}
}

/*!
 * Update ID of boundary vertices to be coherent with volume vertices
 */
void
Partition::updateBoundaryVerticesID()
{

	if (m_rank != 0){

		//Compute map between old and new boundary vertices id
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
			old2new[id] = idnew;
		}

		Vertices.sort();
		for (bitpit::Vertex & vertex : Vertices){
			long id = vertex.getId();
			Vertices.updateId(id, id-idmax-1);
			vertex.setId(id-idmax-1);
		}

		//Update connectivity of boundary cells
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
    std::string name = m_outputPlot +"/"+ m_name;
    getGeometry()->getPatch()->write(name);
	if (getBoundaryGeometry() != nullptr)
		getBoundaryGeometry()->getPatch()->write(name + "Boundary");
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
