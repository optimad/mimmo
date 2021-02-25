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
\*---------------------------------------------------------------------------*/

#include "RefineGeometry.hpp"
#if MIMMO_ENABLE_MPI
#include "mimmo_parallel.hpp"
#endif

namespace mimmo{

/*!
 * Default constructor of RefineGeometry.
 */
RefineGeometry::RefineGeometry(){
	m_name  = "mimmo.RefineGeometry";
	m_type	= RefineType(1);
	m_steps = 0;
	m_refinements = 1;
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
RefineGeometry::RefineGeometry(const bitpit::Config::Section & rootXML){

	std::string fallback_name = "ClassNONE";
	std::string input_name = rootXML.get("ClassName", fallback_name);
	input_name = bitpit::utils::string::trim(input_name);

	m_name = "mimmo.RefineGeometry";
	m_type = RefineType(1);
	m_steps = 0;
	m_refinements = 1;

	if(input_name == "mimmo.RefineGeometry"){
		absorbSectionXML(rootXML);
	}else{
		warningXML(m_log, m_name);
	};
}

/*!
 * Default destructor of RefineGeometry.
 */
RefineGeometry::~RefineGeometry(){};

/*!
 * Copy constructor of RefineGeometry.
 */
RefineGeometry::RefineGeometry(const RefineGeometry & other):BaseManipulation(other){
	m_type = other.m_type;
	m_steps = other.m_steps;
	m_refinements = other.m_refinements;
};

/*!
 * Assignement operator of RefineGeometry. Soft copy of MimmoObject
 */
RefineGeometry & RefineGeometry::operator=(RefineGeometry other){
	swap(other);
	return *this;
};


/*!
 * Swap function of RefineGeometry.
 * \param[in] x object to be swapped.
 */
void RefineGeometry::swap(RefineGeometry & x ) noexcept
{
	std::swap(m_type,x.m_type);
	std::swap(m_steps,x.m_steps);
	std::swap(m_refinements,x.m_refinements);
	BaseManipulation::swap(x);
};


/*!
 * Building ports of the class
 */
void
RefineGeometry::buildPorts(){
	bool built = true;
	built = (built && createPortIn<mimmo::MimmoSharedPointer<MimmoObject>, RefineGeometry>(this, &BaseManipulation::setGeometry, M_GEOM, true));
	built = (built && createPortOut<mimmo::MimmoSharedPointer<MimmoObject>, RefineGeometry>(this, &BaseManipulation::getGeometry, M_GEOM));
	m_arePortsBuilt = built;
}

/*!
 * Return kind of refinement set for the object.
 * \return refine type
 */
RefineType
RefineGeometry::getRefineType(){
	return m_type;
}

/*!
 * It sets refinement method of refine block
 * \param[in] type refine type
 *
 */
void
RefineGeometry::setRefineType(RefineType type){
	if (type != RefineType::TERNARY && type != RefineType::REDGREEN)
		throw std::runtime_error(m_name + " : refinement method not allowed");
	m_type = type;
};

/*!
 * It sets refinement method of refine block
 * \param[in] type refine type
 */
void
RefineGeometry::setRefineType(int type){
	if (type != 0 && type != 1)
		throw std::runtime_error(m_name + " : refinement method not allowed");

	m_type = RefineType(type);
};

/*!
 * It sets the number of refinement steps
 * \param[in] steps refinement steps
 *
 */
void
RefineGeometry::setRefineSteps(int steps){
	m_refinements = std::max(0, steps);
};

/*!
 * It sets the number of positive/negative laplacian smoothing steps
 * \param[in] steps smoothing steps
 *
 */
void
RefineGeometry::setSmoothingSteps(int steps){
	m_steps = std::max(0, steps);
};

/*!
 * Clear all stuffs in your class
 */
void
RefineGeometry::clear(){
	BaseManipulation::clear();
};


/*!Execution command.
 * It refines the target surface geometry.
 */
void
RefineGeometry::execute(){

	if(getGeometry() == nullptr){
		(*m_log)<<m_name + " : no geometry to refine found"<<std::endl;
		return;
	}

	if (getGeometry()->getType() != 1){
		(*m_log)<<m_name + " : data structure type different from Surface "<<std::endl;
		throw std::runtime_error (m_name + " : data structure type different from Surface");
	}

    if(m_type == RefineType::REDGREEN){
        bool check = checkTriangulation();
        if(!check){
            (*m_log)<<m_name + " : detected not regular triangulated surface, no redgreen refinement can be performed. Returning. "<<std::endl;
            return;
        }
    }

	// Initialize active cells with all geometry cells
    for (int i=0; i<m_refinements; i++){
        m_activecells = getGeometry()->getCellsIds();

    	if (m_type == RefineType::TERNARY){
    			ternaryRefine();
    	}
    	else if (m_type == RefineType::REDGREEN){
    			redgreenRefine();
    	}
    }
	if (m_steps>0){
		smoothing();
    }
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void RefineGeometry::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);

	std::string input;

	if(slotXML.hasOption("RefineType")){
		std::string input = slotXML.get("RefineType");
		input = bitpit::utils::string::trim(input);
		int value = 1;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setRefineType(value);
	}

	if(slotXML.hasOption("RefineSteps")){
		std::string input = slotXML.get("RefineSteps");
		input = bitpit::utils::string::trim(input);
		int value = 0;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setRefineSteps(value);
	}

	if(slotXML.hasOption("SmoothingSteps")){
		std::string input = slotXML.get("SmoothingSteps");
		input = bitpit::utils::string::trim(input);
		int value = 0;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setSmoothingSteps(value);
	}

	BaseManipulation::absorbSectionXML(slotXML, name);

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void RefineGeometry::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

	BaseManipulation::flushSectionXML(slotXML, name);
	slotXML.set("RefineType", std::to_string(int(m_type)));
	slotXML.set("RefineSteps", std::to_string(m_refinements));
	slotXML.set("SmoothingSteps", std::to_string(m_steps));

};

/*!
 * Plot stitched geometry in *.vtu file as optional result;
 */
void
RefineGeometry::plotOptionalResults(){
	write(getGeometry());
}

/*!
 * Global refinement by ternary method. A vertex is added on the barycenter of each cell.
 * One triangle for each edge is added for every cell.
 */
void
RefineGeometry::ternaryRefine(std::unordered_map<long,long> * mapping, mimmo::MimmoSharedPointer<MimmoObject> coarsepatch, mimmo::MimmoSharedPointer<MimmoObject> refinepatch)
{
	//TERNARY REFINEMENT

	// Ternary refinement force the geometry to be a triangulation by placing a new vertex
	// on the mean point of each cell (the cells have to be convex).

    // It works in parallel, supposing the cell refined in the same manner
    // if local or ghost on the owner process or on the near process.
    // The ids are not maintained unique between all the processes
#if MIMMO_ENABLE_MPI
	if (m_nprocs > 1){
		(*m_log)<< "WARNING " <<m_name <<" : uniqueness of cell/vertex ids among processes is not maintained during parallel geometry refinement."<<std::endl;
	}
#endif

	mimmo::MimmoSharedPointer<MimmoObject> geometry = getGeometry();
	std::unordered_set<long> newCells;
	std::unordered_set<long> toDelete;

	for(const long cellId : m_activecells){

		// Eval centroid
		std::array<double,3> centroid = getGeometry()->evalCellCentroid(cellId);

		// Build perimeter structure with the vertices of the cell and then a triangulation refinement
		// will be performed by placing the centroid as new vertex
		std::vector<bitpit::Vertex> perimeter;
		for (const long & id : getGeometry()->getPatch()->getCell(cellId).getVertexIds()){
			perimeter.emplace_back(getGeometry()->getPatch()->getVertex(id));
		}

		// Refine cell
		std::vector<long> generatedCells = ternaryRefineCell(cellId, perimeter, centroid);


        if (!generatedCells.empty()){

            //Add cell to todelete and refined structure
            toDelete.insert(cellId);

            // Add vertices and cell to coarse patch
            if (coarsepatch != nullptr)
            {
                bitpit::Cell & cell = getGeometry()->getPatch()->getCell(cellId);
                coarsepatch->addCell(cell, cellId);
                for (long vertexId : cell.getVertexIds()){
                    bitpit::Vertex vertex = getGeometry()->getPatch()->getVertex(vertexId);
                    coarsepatch->addVertex(vertex, vertexId);
                }
            }

            // Add entry to refine-coarse mapping and newCells structure
            for (long newCellId : generatedCells){
                if (mapping != nullptr)
                    mapping->insert({newCellId, cellId});
                newCells.insert(newCellId);
            }

        } // end if generated cells is not empty

	} // end loop on cells



	//Delete cells and clean geometries
	{
		getGeometry()->getPatch()->deleteCells(livector1D(toDelete.begin(), toDelete.end()));
        getGeometry()->setUnsyncAll();

#if MIMMO_ENABLE_MPI
        // Update adjacencies
        getGeometry()->updateAdjacencies();

        // Delete orphan ghosts
        getGeometry()->deleteOrphanGhostCells();
        getGeometry()->getPatch()->deleteOrphanVertices();
#endif

		// Update coarse patch
        getGeometry()->update();
		if (coarsepatch != nullptr){
            coarsepatch->update();
		}
	}

	// Add vertices and cells to refine patch
	if (refinepatch != nullptr){
		for (long newcellId : newCells){
			bitpit::Cell & cell = getGeometry()->getPatch()->getCell(newcellId);
			refinepatch->addCell(cell, newcellId);
			for (long vertexId : cell.getVertexIds()){
				const bitpit::Vertex & vertex = getGeometry()->getPatch()->getVertex(vertexId);
				refinepatch->addVertex(vertex, vertexId);
			}
		}

#if MIMMO_ENABLE_MPI
        // update adjacencies refine patch
        refinepatch->update();

        // Delete orphan ghosts
        refinepatch->deleteOrphanGhostCells();
        refinepatch->getPatch()->deleteOrphanVertices();
#endif

        // update refine patch
        refinepatch->update();
	}

}


/*!
 * It refines the target cell with perimeter defined by vertices and center point.
 * \param[in] cellId cell id of the cell to be refined
 * \param[in] vertices perimeter vertices of the cell
 * \param[in] center center point of the cell
 * \return Cell ids of the new refined cells insert in the geometry
 *
 * Note: the old refined cell is not deleted from the geometry in this function.
 *
 */
std::vector<long>
RefineGeometry::ternaryRefineCell(const long & cellId, const std::vector<bitpit::Vertex> & vertices, const std::array<double,3> & center)
{

    long newID, newVertID;

	std::vector<long> newCellIDs;

	bitpit::ElementType eletype;
	bitpit::ElementType eletri = bitpit::ElementType::TRIANGLE;
	livector1D connTriangle(3);

    // Take a cell copy, after the emplace back of the refined cells the reference can be changed
    const bitpit::Cell & cell = getGeometry()->getPatch()->getCell(cellId);

	eletype = cell.getType();
	long pid = cell.getPID();
    int rank = -1;

	//Number of new triangles is the number of the perimetral vertices
	std::size_t nnewTri = vertices.size();
	newCellIDs.reserve(nnewTri);

	// Add barycenter and new vertices (check for old vertices)
	newVertID = getGeometry()->addVertex(center);

	//insert new triangles from polygon subdivision
	for(std::size_t i=0; i<nnewTri; ++i){
		connTriangle[0] = newVertID;
		connTriangle[1] = vertices[ std::size_t( i % nnewTri) ].getId();
		connTriangle[2] = vertices[ std::size_t( (i+1) % nnewTri ) ].getId();
#if MIMMO_ENABLE_MPI
		// Recover cell rank
		rank = getGeometry()->getPatch()->getCellRank(cellId);
#endif
        newCellIDs.push_back(getGeometry()->addConnectedCell(connTriangle, eletri, pid, bitpit::Cell::NULL_ID, rank));
	}

	return newCellIDs;

}

/*!
 * Refinement by redgreen method. The starting elements must be all triangles.
 * All the triangle edges are splitted for every cell in case of red refinement (bulk cells),
 * only one edge is splitted in case of green refinement (boundary cells).
 * \param[out] mapping pointer to mapping container (new child cell id -> parent cell id)
 * \param[out] coarsepatch pointer to geometry containing the original parent cells
 * \param[out] refinepatch pointer to geometry containing the new children refined cells
 *
 */
void
RefineGeometry::redgreenRefine(std::unordered_map<long,long> * mapping, mimmo::MimmoSharedPointer<MimmoObject> coarsepatch, mimmo::MimmoSharedPointer<MimmoObject> refinepatch)
{

	//REDGREEN REFINEMENT
    // It works in parallel, supposing the cell refined in the same manner
    // if local or ghost on the owner process or on the near process.
    // The ids are not maintained unique between all the processes
#if MIMMO_ENABLE_MPI
    if (m_nprocs > 1){
        (*m_log)<< "WARNING " <<m_name <<" : uniqueness of cell/vertex ids among processes is not maintained during parallel geometry refinement."<<std::endl;
    }
#endif

	mimmo::MimmoSharedPointer<MimmoObject> geometry = getGeometry();
	std::unordered_set<long> newCells;
	std::unordered_set<long> toDelete;

	// Initialize tag Red, Green structure (0=no refinement, 1=green, 2=red) to 0 for all the cells
    std::unordered_map<long, int> refinementTag;
	for(const long cellId : geometry->getCellsIds()){
		refinementTag[cellId] = 0;
	}

	// Create container for edges to be splitted by using only red elements (green are included automatically)
	// In order to be unique the interface id of the patch is used
	std::set<long> edges; //edge defined as interface index

	// Create a map edge id -> new vertex id
	std::unordered_map<long,long> edgeVertexId;

	// Create map green cells -> index edge to be splitted
	std::unordered_map<long,int> greenSplitFaceIndex;

    {
	    // Build adjacencies and interfaces if not built
        geometry->updateAdjacencies();
        geometry->updateInterfaces();

        // Set active cells as reds and initialize new reds stack
        std::deque<long> newreds;
        for(const long cellId : m_activecells){
            refinementTag[cellId] = 2;
            newreds.push_back(cellId);
        }

#if MIMMO_ENABLE_MPI
		// In case of distributed mesh initialize data
		// of newreds boolean to update new reds that are
		// ghosts for other ranks
		MimmoPiercedVector<bool> isnewred;
		isnewred.initialize(geometry, MPVLocation::CELL, false);

        // Instantiate refinement tags container used only for communications
        MimmoPiercedVector<int> refinementTagCommunicated;

        // Instantiate green split face container used only for communications
        MimmoPiercedVector<int> greenSplitFaceIndexCommunicated;

        // Fill initial sources and targets values
        for (auto & source_tuple : geometry->getPatch()->getGhostCellExchangeSources()){
            int rank = source_tuple.first;
            for (long id : source_tuple.second){
                if (!refinementTagCommunicated.exists(id)){
                    refinementTagCommunicated.insert(id,-1);
                    greenSplitFaceIndexCommunicated.insert(id,-1);
                }
            }
        }
        for (auto & target_tuple : geometry->getPatch()->getGhostCellExchangeTargets()){
            int rank = target_tuple.first;
            for (long id : target_tuple.second){
                // Initialize targets (ghosts) to -1
                if (!refinementTagCommunicated.exists(id)){
                    refinementTagCommunicated.insert(id, -1);
                    greenSplitFaceIndexCommunicated.insert(id, -1);
                }
            }
        }

#endif

		// If active cells are all the cells, propagate red-green refinement
		// even if all the cells are reds, in order to build edges structure
		bool check = newreds.empty();

#if MIMMO_ENABLE_MPI
		// While global newreds stack is not empty. See below for details on local new reds.
		// Initialized to false to enter always in while loop at least one time
        bool global_check = false;

        while (!global_check){
#endif

            // While newreds stack is not empty:
            // - serch neighbors of newreds elements
            // - check in tag map if each neighbor currently it's no,green or red element
            // - if 0-> promote to green (1), if green->promote to red (2) and insert in newreds, if red->skip neighbour
            while (!check){

                long redId = newreds.front();
                newreds.pop_front();

                bitpit::Cell & redCell = geometry->getPatch()->getCell(redId);

                //Only for triangles!!!
                if (redCell.getType() == bitpit::ElementType::TRIANGLE){

                    for (int iface = 0; iface < 3; iface++){

                        if (!redCell.isFaceBorder(iface)){

                            // Find face neighbours (only one in conform case)
                            int neighscount = redCell.getAdjacencyCount(iface);
                            const long * neighs = redCell.getAdjacencies(iface);

                            for (int ineigh = 0; ineigh < neighscount; ineigh++){

                                long neighId = neighs[ineigh];

                                // Increase tag (at the end the red elements have a tag >=2)
                                refinementTag[neighId]++;

                                // Insert edge between current red and neighbor
                                // The edges structure is a set, so each edge will be not duplicated
                                // Recover interface
                                long interfaceId = geometry->getPatch()->getCell(redId).getInterface(iface);
                                bitpit::Interface & interface = geometry->getPatch()->getInterface(interfaceId);
                                edges.insert(interfaceId);

                                if (refinementTag[neighId] > 2){
                                    continue;
                                } // if neigh already red


                                if (refinementTag[neighId] == 2){

                                    // Push neigh to new reds
                                    newreds.push_back(neighId);

                                    // Destroy green entry with splitting edge index
                                    greenSplitFaceIndex.erase(neighId);

                                } // If is neigh is new red
                                else if (refinementTag[neighId] == 1){

                                    // Get if neighbor is owner of neigh
                                    bool isOwner = (interface.getOwner() == neighId);

                                    // Recover splitting face index of the neighbor
                                    int splitface;
                                    if (isOwner){
                                        splitface = interface.getOwnerFace();
                                    }
                                    else{
                                        splitface = interface.getNeighFace();
                                    }
                                    greenSplitFaceIndex[neighId] = splitface;

                                } // Else if neigh is green
                            } // end loop on neighs
                        } // end if not border face
                        else{

                            // If border face the edge is to be refined
                            // Recover border interface
                            long interfaceId = redCell.getInterface(iface);
                            edges.insert(interfaceId);

                        }// end if border face

                    }// End loop on face

                } // if triangle

                check = newreds.empty();

            } // end while stack

#if MIMMO_ENABLE_MPI
            // In case of not partitioned patch use local check
            global_check = check;

            if (geometry->isParallel()){

                // Update ghost refinement tags
                // The if needed insert new edges and put new reds in stack
                // Fill with sources and targets values
                for (auto & source_tuple : geometry->getPatch()->getGhostCellExchangeSources()){
                    int rank = source_tuple.first;
                    for (long id : source_tuple.second){
                        refinementTagCommunicated.at(id) = refinementTag.at(id);
                        if (greenSplitFaceIndex.count(id)){
                            greenSplitFaceIndexCommunicated.at(id) = greenSplitFaceIndex.at(id);
                        }
                        else{
                            greenSplitFaceIndexCommunicated.at(id) = -1;
                        }
                    }
                }
                for (auto & target_tuple : geometry->getPatch()->getGhostCellExchangeTargets()){
                    int rank = target_tuple.first;
                    for (long id : target_tuple.second){
                        // Initialize targets (ghosts) to -1
                        refinementTagCommunicated.at(id) = -1;
                        greenSplitFaceIndexCommunicated.at(id) = -1;
                    }
                }

                // Communicate tag and face index
                refinementTagCommunicated.communicateData();
                greenSplitFaceIndexCommunicated.communicateData();

                // Update refinement tags if different from communicated
                for (auto ghostIt = geometry->getPatch()->ghostCellConstBegin(); ghostIt != geometry->getPatch()->ghostCellConstEnd(); ghostIt++){
                    long ghostId = ghostIt->getId();

                    // Maximum communicated tag to 2
                    refinementTagCommunicated[ghostId] = std::min(2, refinementTagCommunicated[ghostId]);

                    // If communicated tag is different from local one update it if greater
                    if (refinementTagCommunicated[ghostId] > refinementTag[ghostId]){
                        refinementTag[ghostId] = refinementTagCommunicated[ghostId];

                        if (refinementTag[ghostId] == 1){
                            // If green refinement insert splitFaceIndex
                            greenSplitFaceIndex[ghostId] = greenSplitFaceIndexCommunicated[ghostId];
                        }
                        else if (refinementTag[ghostId] >= 2){
                            // Push ghost to new reds
                            newreds.push_back(ghostId);
                            // Destroy green entry with splitting edge index
                            greenSplitFaceIndex.erase(ghostId);
                        }
                    }
                } // end ghost iteration

                // Update global stack check
                global_check = newreds.empty();
                MPI_Allreduce(MPI_IN_PLACE, &global_check, 1, MPI_C_BOOL, MPI_LAND, m_communicator);
            }

        } // end while global stack
#endif

	} // Scope to destroy temporary containers

	// Insert new vertices
	for (long interfaceId : edges){
        // Recover vertices
		bitpit::ConstProxyVector<long> vertexIds = geometry->getPatch()->getInterface(interfaceId).getVertexIds();
		std::array<double,3> newCoordinates({0.,0.,0.});
		for (long vertexId : vertexIds){
			newCoordinates += geometry->getVertexCoords(vertexId);
		}
		newCoordinates /= double(vertexIds.size());
		long newVertexId = geometry->addVertex(newCoordinates);
		edgeVertexId[interfaceId] = newVertexId;
	}

	// Refine red and green cells
	for (auto & tupletag : refinementTag){

		long cellId = tupletag.first;
		int tag = std::min(2, tupletag.second);

		if (tag <= 0)
			continue;

		// Generated cells structure
		std::vector<long> generatedCells;

		if (tag == 2){

		    // Take reference to cell before refinement
            bitpit::Cell & cell = geometry->getPatch()->getCell(cellId);

	        // Red refinement
			// Recover new vertices
			std::vector<long> newCellVertexIds(cell.getFaceCount());
			for (int iface=0; iface<cell.getFaceCount(); iface++){
				long interfaceId = cell.getInterface(iface);
				newCellVertexIds[iface] = edgeVertexId.at(interfaceId);
			}
			// Cell refinement
			generatedCells = redRefineCell(cellId, newCellVertexIds);
		}
		else if (tag == 1){

            // Take reference to cell before refinement
            bitpit::Cell & cell = geometry->getPatch()->getCell(cellId);

            // Green refinement
			// Recover new vertex and splitted face index
			int splitFaceIndex = greenSplitFaceIndex.at(cellId);
			long interfaceId = cell.getInterface(splitFaceIndex);
			long newCellVertexId = edgeVertexId.at(interfaceId);
			// Cell refinement
			generatedCells = greenRefineCell(cellId, newCellVertexId, splitFaceIndex);
		}

		if (!generatedCells.empty()){

		    //Add cell to todelete structure
		    toDelete.insert(cellId);

		    // Add vertices and cell to coarse patch
		    if (coarsepatch != nullptr)
		    {
	            // Take reference to cell after refinement
	            // The original cell is still inside the cells structure but the new cells
	            // are already insert, so the reference has to be locally taken.
	            bitpit::Cell & cell = geometry->getPatch()->getCell(cellId);

	            coarsepatch->addCell(cell, cellId);
		        for (long vertexId : cell.getVertexIds()){
		            bitpit::Vertex & vertex = geometry->getPatch()->getVertex(vertexId);
		            coarsepatch->addVertex(vertex, vertexId);
		        }
		    }

		    // Add entry to refine-coarse mapping and newCells structure
		    for (long newCellId : generatedCells){
		        if (mapping != nullptr){
		            mapping->insert({newCellId, cellId});
		        }
		        newCells.insert(newCellId);
		    }

		} // end if generated cells is not empty

	} // end loop on refinement tag

	//Delete cells and clean geometries
	{
	    // Destroy interfaces
	    getGeometry()->destroyInterfaces();

	    // Delete original cells
        getGeometry()->getPatch()->deleteCells(livector1D(toDelete.begin(), toDelete.end()));

        // Set all sync status of mimmo to unsync/none
        getGeometry()->setUnsyncAll();

	    // Add vertices and cells to refine patch
	    if (refinepatch != nullptr){
	        for (long newcellId : newCells){
	            bitpit::Cell & cell = getGeometry()->getPatch()->getCell(newcellId);
	            refinepatch->addCell(cell, newcellId);
	            for (long vertexId : cell.getVertexIds()){
	                bitpit::Vertex & vertex = getGeometry()->getPatch()->getVertex(vertexId);
	                refinepatch->addVertex(vertex, vertexId);
	            }
	        }
	    }
#if MIMMO_ENABLE_MPI
	    // Delete orphan ghosts
        getGeometry()->deleteOrphanGhostCells();
        getGeometry()->getPatch()->deleteOrphanVertices();
#endif
	    // Update geometry
	    getGeometry()->update();

	    // Update coarse patch
	    if (coarsepatch != nullptr){
	        coarsepatch->update();
	    }

	    // Update refine patch
	    if (refinepatch != nullptr){
#if MIMMO_ENABLE_MPI
	        // Delete orphan ghosts
	        refinepatch->deleteOrphanGhostCells();
            refinepatch->getPatch()->deleteOrphanVertices();
#endif
	        // Update refine patch
	        refinepatch->update();
	    }

	} // end scope

}

/*!
 * It refines the target cell with red method.
 * \param[in] cellId Id of the target cell
 * \param[in] newVertexIds Ids of the new vertices (already in the mesh) to be used to refine the cell
 * \return Ids of the new created cells
 */
std::vector<long>
RefineGeometry::redRefineCell(const long & cellId, const std::vector<long> & newVertexIds)
{

	std::vector<long> newCellIDs;

	bitpit::ElementType eletype;
	bitpit::ElementType eletri = bitpit::ElementType::TRIANGLE;
	livector1D connTriangle(3);

    // Take a cell copy, after the emplace back of the refined cells the reference can be changed
    bitpit::Cell cell = getGeometry()->getPatch()->getCell(cellId);

	eletype = cell.getType();
	//Only for triangles
	if (eletype != eletri){
		(*m_log)<<m_name + " : red refinement allowd only for triangles. Skip element."<<std::endl;
		return newCellIDs;
	}

	long pid = cell.getPID();
	int rank = -1;
#if MIMMO_ENABLE_MPI
    // Recover cell rank
    rank = getGeometry()->getPatch()->getCellRank(cellId);
#endif

	//Number of new triangles is 4
	std::size_t sizeEle = 3;
	newCellIDs.reserve(sizeEle+1);

	//insert new triangles from red subdivision
	// Insert internal one
	// newVertexIds are supposed ordered as indices of faces (0,1,2)
	connTriangle[0] = newVertexIds[0];
	connTriangle[1] = newVertexIds[1];
	connTriangle[2] = newVertexIds[2];
	newCellIDs.push_back(getGeometry()->addConnectedCell(connTriangle, eletri, pid, bitpit::Cell::NULL_ID, rank));

	// Insert three vertex related triangles
	// Note. start from refined triangle placed on vertex index 1 of coarse triangle
	for(std::size_t i=0; i<sizeEle; ++i){
		connTriangle[0] = cell.getVertexId( int( (i+1) % sizeEle) );
		connTriangle[1] = newVertexIds[ std::size_t( (i+1) % sizeEle) ];
		connTriangle[2] = newVertexIds[ std::size_t( (i) % sizeEle ) ];
		newCellIDs.emplace_back(getGeometry()->addConnectedCell(connTriangle, eletri, pid, bitpit::Cell::NULL_ID, rank));
	}

	return newCellIDs;
}

/*!
 * It refines the target cell with green method.
 * \param[in] cellId Id of the target cell
 * \param[in] newVertexId Id of the new vertex (already in the mesh) to be used to refine the cell
 * \param[in] splitEdgeIndex Index of the edge to split
 * \return Ids of the new created cells
 */
std::vector<long>
RefineGeometry::greenRefineCell(const long & cellId, const long newVertexId, int splitEdgeIndex)
{

	std::vector<long> newCellIDs;

	bitpit::ElementType eletype;
	bitpit::ElementType eletri = bitpit::ElementType::TRIANGLE;
	livector1D connTriangle(3);

	// Take a cell copy, after the emplace back of the refined cells the reference can be changed
	bitpit::Cell cell = getGeometry()->getPatch()->getCell(cellId);

	eletype = cell.getType();
	//Only for triangles
	if (eletype != eletri){
		(*m_log)<<m_name + " : red refinement allowd only for triangles. Skip element."<<std::endl;
		return newCellIDs;
	}

	long pid = cell.getPID();
	int rank = -1;
#if MIMMO_ENABLE_MPI
    // Recover cell rank
    rank = getGeometry()->getPatch()->getCellRank(cellId);
#endif

	//Number of new triangles is 2
	std::size_t sizeEle = 3;
	newCellIDs.reserve(sizeEle-1);

	// Insert three vertex related triangles
	// Note. start from refined triangle placed on vertex index 1 of coarse triangle
	for(std::size_t i=0; i<sizeEle-1; ++i){
		connTriangle[0] = newVertexId;
		connTriangle[1] = cell.getVertexId( int( (splitEdgeIndex+1+i) % sizeEle) );
		connTriangle[2] = cell.getVertexId( int( (splitEdgeIndex+2+i) % sizeEle) );
		newCellIDs.emplace_back(getGeometry()->addConnectedCell(connTriangle, eletri, pid, bitpit::Cell::NULL_ID, rank));
	}

	return newCellIDs;

}


/*!
 * Perform the laplacian smoothing for the surface patch.
 * Each step is divided in two sub-steps, a laplacian smoothing
 * and a laplacian anti-smoothing with negative coefficient.
 * \param[in] constrainedVertices Pointer to set of constrined vertices indices
 */
void
RefineGeometry::smoothing(std::set<long> * constrainedVertices)
{

    mimmo::MimmoSharedPointer<MimmoObject> geometry = getGeometry();

    // Force build the needed structures
    if(geometry->getAdjacenciesSyncStatus() != mimmo::SyncStatus::SYNC)
        geometry->updateAdjacencies();

    if(geometry->getInterfacesSyncStatus() != mimmo::SyncStatus::SYNC)
        geometry->updateInterfaces();

    if(geometry->getInfoSyncStatus() != mimmo::SyncStatus::SYNC)
        geometry->buildPatchInfo();

    if(geometry->getPointConnectivitySyncStatus() != mimmo::SyncStatus::SYNC)
        geometry->buildPointConnectivity();

#if MIMMO_ENABLE_MPI

    if(geometry->getPointGhostExchangeInfoSyncStatus() != mimmo::SyncStatus::SYNC)
        geometry->updatePointGhostExchangeInfo();

    // Instantiate new coordinates container used only for communications
    MimmoPiercedVector<std::array<double,3>> newCoordinatesCommunicated(geometry, MPVLocation::POINT);

    // Fill initial sources and targets values
    for (auto source_tuple : geometry->getPointGhostExchangeSources()){
        int rank = source_tuple.first;
        for (long id : source_tuple.second){
            if (!newCoordinatesCommunicated.exists(id)){
                newCoordinatesCommunicated.insert(id, std::array<double,3>({{0.,0.,0.}}));
            }
        }
    }
    for (auto target_tuple : geometry->getPointGhostExchangeTargets()){
        int rank = target_tuple.first;
        for (long id : target_tuple.second){
            if (!newCoordinatesCommunicated.exists(id)){
                newCoordinatesCommunicated.insert(id,  std::array<double,3>({{0.,0.,0.}}));
            }
        }
    }

#endif


	double lambda = 0.6;
	double kappa = -0.603*lambda;

	for (int istep=0; istep < m_steps; istep++){

		{
			// First sub-step (positive) of laplacian smoothing
			std::unordered_map<long, std::array<double,3>> newcoordinates;
			std::unordered_set<long> pointconnectivity;
			std::array<double,3> newcoords, oldcoords, neighcoords;
			double weight, sumweights;
			newcoordinates.reserve(geometry->getNVertices());
			//compute new coordinates
//			for (long id : geometry->getVertices().getIds()){
            for (const bitpit::Vertex & vert : geometry->getVertices()){
                long id = vert.getId();
			    // If ghost vertex do nothing
			    if (geometry->isPointInterior(id)){

			        oldcoords = geometry->getVertexCoords(id);
			        newcoords = std::array<double,3>{{0.,0.,0.}};

			        // If constrained vertex do nothing
			        if (constrainedVertices == nullptr || !constrainedVertices->count(id)){
			            pointconnectivity = geometry->getPointConnectivity(id);
			            sumweights = 0.;
			            for (long idneigh : pointconnectivity){
			                neighcoords = geometry->getVertexCoords(idneigh);
			                std::array<double,3> distancevector = (neighcoords-oldcoords);
			                weight = 1.;
			                sumweights += weight;
			                newcoords += lambda*weight*distancevector;
			            }
			            if (sumweights > 0.)
			                newcoords /= sumweights;
			        }

			        newcoords += oldcoords;
			        newcoordinates[id] = newcoords;
			    }

			}

#if MIMMO_ENABLE_MPI
			//Communicate newcoordinates
            if (geometry->isParallel()){

                // Update ghost new coordinates

                // Fill with sources and targets values
                for (auto source_tuple : geometry->getPointGhostExchangeSources()){
                    int rank = source_tuple.first;
                    for (long id : source_tuple.second){
                        newCoordinatesCommunicated.at(id) = newcoordinates.at(id);
                    }
                }

                // Communicate new coordinates
                newCoordinatesCommunicated.communicateData();

                // Update coordinates of ghost vertices with communicated ones
                for (auto target_tuple : geometry->getPointGhostExchangeTargets()){
                    int rank = target_tuple.first;
                    for (long id : target_tuple.second){
                        newcoordinates[id] = newCoordinatesCommunicated.at(id);
                    } // end id
                } // end tuple

            } // end if geometry is parallel
#endif

			//Set new coordinates
//			for (long id : geometry->getVertices().getIds()){
            for (const bitpit::Vertex & vert : geometry->getVertices()){
                long id = vert.getId();
				geometry->modifyVertex(newcoordinates[id], id);
			}

			//Update geometry
			geometry->update();

		}

		{
			//Second sub-step (negative) of laplacian anti-smoothing
			std::unordered_map<long, std::array<double,3>> newcoordinates;
			std::unordered_set<long> pointconnectivity;
			std::array<double,3> newcoords, oldcoords, neighcoords;
			double weight, sumweights;
			newcoordinates.reserve(geometry->getNVertices());
			//compute new coordinates
//			for (long id : geometry->getVertices().getIds()){
            for (const bitpit::Vertex & vert : geometry->getVertices()){
                long id = vert.getId();
				oldcoords = geometry->getVertexCoords(id);
				newcoords = std::array<double,3>{{0.,0.,0.}};

				// If constrained vertex do nothing
				if (constrainedVertices != nullptr && !constrainedVertices->count(id)){
					pointconnectivity = geometry->getPointConnectivity(id);
					sumweights = 0.;
					for (long idneigh : pointconnectivity){
						neighcoords = geometry->getVertexCoords(idneigh);
						std::array<double,3> distancevector = (neighcoords-oldcoords);
						weight = 1.;
						sumweights += weight;
						newcoords += kappa*weight*(neighcoords-oldcoords);
					}
					if (sumweights > 0.)
						newcoords /= sumweights;
				}

				newcoords += oldcoords;
				newcoordinates[id] = newcoords;
			}

#if MIMMO_ENABLE_MPI
            //Communicate newcoordinates
            if (geometry->isParallel()){

                // Update ghost new coordinates
                // The if needed insert new edges and put new reds in stack

                // Fill with sources and targets values
                for (auto source_tuple : geometry->getPointGhostExchangeSources()){
                    int rank = source_tuple.first;
                    for (long id : source_tuple.second){
                        newCoordinatesCommunicated.at(id) = newcoordinates.at(id);
                    }
                }

                // Communicate new coordinates
                newCoordinatesCommunicated.communicateData();

                // Update coordinates of ghost vertices with communicated ones
                for (auto target_tuple : geometry->getPointGhostExchangeTargets()){
                    for (long id : target_tuple.second){
                        newcoordinates.at(id) = newCoordinatesCommunicated.at(id);
                    } // end id iteration
                } // end tuple iteration

            } // end if geometry is parallel
#endif

			//Set new coordinates
//			for (long id : geometry->getVertices().getIds()){
            for (const bitpit::Vertex & vert : geometry->getVertices()){
                long id = vert.getId();
				geometry->modifyVertex(newcoordinates[id], id);
			}

			//Update geometry
			geometry->update();

		}  // end second sub step


	} // end step iteration

}

/*!
    Check if the target geometry of the class is a regular triangulation
    \return true if a regular triangulation.
*/
bool
RefineGeometry::checkTriangulation(){
    // Check if the m_geometry is a regular triangulation
    bool check = true;
    for (const bitpit::Cell & cell : getGeometry()->getCells()){
        check = check && (cell.getType() == bitpit::ElementType::TRIANGLE);
    }
#if MIMMO_ENABLE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &check, 1, MPI_C_BOOL, MPI_LAND, m_communicator);
#endif

    return check;
}

}
