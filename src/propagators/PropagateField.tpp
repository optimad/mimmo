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

#include "customOperators.hpp"
#include "SkdTreeUtils.hpp"

namespace mimmo{

/*!
 * Constructor
 */
template<std::size_t NCOMP>
PropagateField<NCOMP>::PropagateField(){
	m_name = "mimmo.PropagateField";
	setDefaults();
};

/*!
 * Set most significant parameters to constructor defaults
 */
template<std::size_t NCOMP>
void PropagateField<NCOMP>::setDefaults(){
	this->m_isbp.clear();
	this->m_field.clear();
	this->m_bc_dir.clear();
    this->m_surface_bc_dir.clear();
	this->m_tol       = 1.0e-14;
	this->m_bsurface  = nullptr;
	this->m_dsurface  = nullptr;
	this->m_geometry  = nullptr;
	this->m_decayFactor = 1.0;
	this->m_dumping.clear();
	this->m_radius = 0.0;
	this->m_plateau = 0.0;
	this->m_dumpingActive = false;
	this->m_dumpingType = 0;
    this->m_thres = 1.E-8;
}

/*!
 * Destructor;
 */
template<std::size_t NCOMP>
PropagateField<NCOMP>::~PropagateField(){
	clear();
};

/*!
 * Copy constructor
 */
template<std::size_t NCOMP>
PropagateField<NCOMP>::PropagateField(const PropagateField<NCOMP> & other):BaseManipulation(other){
	setDefaults();
	this->m_isbp         = other.m_isbp;
	this->m_bc_dir       = other.m_bc_dir;
    this->m_surface_bc_dir = other.m_surface_bc_dir;
    this->m_tol          = other.m_tol;
	this->m_bsurface     = other.m_bsurface;
	this->m_dsurface     = other.m_dsurface;
	this->m_dumping      = other.m_dumping;
	this->m_decayFactor  = other.m_decayFactor;
	this->m_radius       = other.m_radius;
	this->m_plateau      = other.m_plateau;
	this->m_dumpingActive= other.m_dumpingActive;
	this->m_dumpingType  = other.m_dumpingType;
    this->m_thres        = other.m_thres;
};

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
template<std::size_t NCOMP>
void PropagateField<NCOMP>::swap(PropagateField<NCOMP> & x) noexcept {
	this->m_isbp.swap(x.m_isbp);
	this->m_field.swap(x.m_field);
	this->m_bc_dir.swap(x.m_bc_dir);
    this->m_surface_bc_dir.swap(x.m_surface_bc_dir);
	std::swap(this->m_tol, x.m_tol);
	std::swap(this->m_bsurface, x.m_bsurface);
	std::swap(this->m_dsurface, x.m_dsurface);
	std::swap(this->m_decayFactor, x.m_decayFactor);
	this->m_dumping.swap(x.m_dumping);
	std::swap(this->m_radius, x.m_radius);
	std::swap(this->m_plateau, x.m_plateau);
	std::swap(this->m_dumpingActive, x.m_dumpingActive);
	std::swap(this->m_dumpingType, x.m_dumpingType);
    std::swap(this->m_thres, x.m_thres);
	this->BaseManipulation::swap(x);
}


/*!
 * It builds the input/output ports of the object
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::buildPorts(){
	bool built = true;
	built = (built && createPortIn<MimmoObject*, PropagateField<NCOMP> >(this, &PropagateField<NCOMP>::setGeometry, M_GEOM, true));
	built = (built && createPortIn<MimmoObject*, PropagateField<NCOMP> >(this, &PropagateField<NCOMP>::setDirichletBoundarySurface, M_GEOM2));
	built = (built && createPortIn<MimmoObject*, PropagateField<NCOMP> >(this, &PropagateField<NCOMP>::setDumpingBoundarySurface, M_GEOM3));
	m_arePortsBuilt = built;
};



/*!
 * Set pointer to your target bulk volume geometry. Reimplemented from mimmo::BaseManipulation::setGeometry().
 * Geometry must be a of volume type (MimmoObject type = 2);
 * \param[in] geometry_ pointer to target geometry
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setGeometry(MimmoObject * geometry_){

	if (geometry_ == NULL) return;
	if (geometry_->isEmpty())   return;
	if (geometry_->getType()!= 2 ) return;

	m_geometry = geometry_;

    if(!m_geometry->areAdjacenciesBuilt()){
        m_geometry->buildAdjacencies();
    }
    //check and build the interfaces if needed
    if(!m_geometry->areInterfacesBuilt()){
        m_geometry->buildInterfaces();
    }
}

/*!
 * Sets the portion of boundary mesh relative to geometry target
 * that must be constrained with Dirichlet conditions.
 * \param[in] bsurface Boundary patch.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setDirichletBoundarySurface(MimmoObject* bsurface){
	if (bsurface == nullptr)       return;
	if (bsurface->isEmpty())    return;
	if (bsurface->getType()!= 1 ) return;

	m_bsurface = bsurface;
}

/*!
 * It sets the sub-portion of boundary mesh to be used for dumping calculation.
 * Sub-boundary and bulk meshes must be Vertex-Id coherent.
 * Must be a valid surface type mesh (MimmoObject type =1).
 * \param[in] bdumping Boundary sub-portion.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setDumpingBoundarySurface(MimmoObject* bdumping){
	if (bdumping == nullptr)       return;
	if (bdumping->isEmpty())    return;
	if (bdumping->getType()!= 1 ) return;

	m_dsurface = bdumping;
}

/*!
 * Activate dumping control of artificial diffusivity(see class doc).
 * \param[in] flag boolean true activate, false deactivate.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setDumping(bool flag){
	m_dumpingActive = flag;
}


/*!
 * Set the inner dumping radius p (see class doc).
 * \param[in] plateau inner distance.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setDumpingInnerDistance(double plateau){
	m_plateau = plateau;
}

/*!
 * Set the outer dumping distance r (see class doc).
 * \param[in] radius outer distance.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setDumpingOuterDistance(double radius){
	m_radius = radius;
}

/*!
 * If dumping is active, set the type of dumping, if distance control based (0) or volume cell control based (1).
 * \param[in] type of control for dumping [0,1].
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setDumpingType(int type){
	m_dumpingType = std::max(0, std::min(1, type));
}

/*!
 * Set the dumping factor.
 * \param[in] dump Exponential of dumping function.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setDecayFactor(double decay){
	m_decayFactor = decay;
}


/*!
 * It sets the tolerance on residuals for solver convergence on both smoothing and laplacian solver.
 * \param[in] tol Convergence tolerance.
 */
template <std::size_t NCOMP>
void PropagateField<NCOMP>::setTolerance(double tol){
	m_tol = tol;
}

/*!
 * It sets a lower threshold such that all the cells who have a norm of solution field
 * above it will be marked for update purposes.
 * \param[in] thres update lower threshold
 */
template <std::size_t NCOMP>
void PropagateField<NCOMP>::setUpdateThreshold(double thres){
	m_thres = std::max(std::numeric_limits<double>::min(), thres);
}


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
template <std::size_t NCOMP>
void PropagateField<NCOMP>::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){


	BITPIT_UNUSED(name);

	//start absorbing
	BaseManipulation::absorbSectionXML(slotXML, name);

	if(slotXML.hasOption("Tolerance")){
		std::string input = slotXML.get("Tolerance");
		input = bitpit::utils::string::trim(input);
		double value = 1.0E-12;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
			value = std::fmax(0.0, value);
		}
		setTolerance(value);
	}

    if(slotXML.hasOption("UpdateThres")){
		std::string input = slotXML.get("UpdateThres");
		input = bitpit::utils::string::trim(input);
		double value = 1.0E-8;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setUpdateThreshold(value);
	}

	if(slotXML.hasOption("Dumping")){
		std::string input = slotXML.get("Dumping");
		input = bitpit::utils::string::trim(input);
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setDumping(value);
	}

	if(slotXML.hasOption("DecayFactor")){
		std::string input = slotXML.get("DecayFactor");
		input = bitpit::utils::string::trim(input);
		double value = 1.0;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
			value = std::fmax(0.0, value);
		}
		setDecayFactor(value);
	}

	if(m_dumpingActive){

		if(slotXML.hasOption("DumpingInnerDistance")){
			std::string input = slotXML.get("DumpingInnerDistance");
			input = bitpit::utils::string::trim(input);
			double value = 0.0;
			if(!input.empty()){
				std::stringstream ss(input);
				ss >> value;
				value = std::fmax(0.0, value);
			}
			setDumpingInnerDistance(value);
		}

		if(slotXML.hasOption("DumpingOuterDistance")){
			std::string input = slotXML.get("DumpingOuterDistance");
			input = bitpit::utils::string::trim(input);
			double value = 1.0e+18;
			if(!input.empty()){
				std::stringstream ss(input);
				ss >> value;
				value = std::fmax(0.0, value);
			}
			setDumpingOuterDistance(value);
		}

		if(slotXML.hasOption("DumpingType")){
			std::string input = slotXML.get("DumpingType");
			input = bitpit::utils::string::trim(input);
			int value = 0;
			if(!input.empty()){
				std::stringstream ss(input);
				ss >> value;
			}
			setDumpingType(value);
		}
	}
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
template <std::size_t NCOMP>
void PropagateField<NCOMP>::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);

	BaseManipulation::flushSectionXML(slotXML, name);

	slotXML.set("Tolerance",std::to_string(m_tol));
    if(m_thres > 0)    slotXML.set("UpdateThres",std::to_string(m_thres));
    slotXML.set("Dumping", std::to_string(int(m_dumpingActive)));
	if(m_dumpingActive){
		slotXML.set("DumpingInnerDistance",std::to_string(m_plateau));
		slotXML.set("DumpingOuterDistance",std::to_string(m_radius));
		slotXML.set("DumpingType",std::to_string(m_dumpingType));
	}
	slotXML.set("DecayFactor",std::to_string(m_decayFactor));

};

/*!
 * Restore data as in class default construction.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::clear(){
	BaseManipulation::clear();
	setDefaults();
};

/*!
 * Check coherence of the input data of the class, in particular:
 * - check if boundary patches shares nodes with the bulk geometry
 * - check if data defined on boundary patches are coherent with them
 * if check is good, append the type of BC info on m_isbp.
 * \return true if coherence is satisfied, false otherwise.
 */
template <std::size_t NCOMP>
bool
PropagateField<NCOMP>::checkBoundariesCoherence(){

    //Clean the old m_isbp and initialize it again.
    initializeBoundaryInfo();

    //1st step verify coherence of the Dirichlet point field on boundary surface
    // with the dirichlet boundary surface provided
    if (m_surface_bc_dir.getGeometry() != m_bsurface || !m_surface_bc_dir.completeMissingData({0.0})){
        return false;
    }

    //1st step check if points-ID of m_bsurface get a boundary-interface-patch on the bulk geometry;
    livector1D interfaces = getGeometry()->getInterfaceFromVertexList(m_bsurface->getVertices().getIds(), true, true);
    if(interfaces.empty()){
        return false;
    }

    //update the m_isbp marking dirichlet interfaces
    for(long id: interfaces){
        m_isbp.at(id) = 1;
    }
    return true;
}

/*!
 * Distribute the input bc values on the border interfaces of the bulk mesh.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::distributeBCOnBoundaryInterfaces(){
    //transfer point field info of boundary dirichlet on the volume mesh interfaces.
    MimmoPiercedVector<std::array<double,NCOMP>> temp(m_geometry, MPVLocation::POINT);
    temp.reserve(m_surface_bc_dir.size());
    for(auto it=m_surface_bc_dir.begin(); it!=m_surface_bc_dir.end(); ++it){
        temp.insert(it.getId(), *it );
    }

    //interpolate now point data to interface data
    m_bc_dir.clear();
    m_bc_dir = temp.pointDataToBoundaryInterfaceData();
}

/*!
 * It computes the dumping function (artificial diffusivity)
 * used to modulate laplacian solution over the target mesh.
 * The dumping function is now a scalar field data on CELLS (ghost included)
 * and it will be stored in internal member m_dumping.
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::computeDumpingFunction(){

	bitpit::PatchKernel * patch_ = getGeometry()->getPatch();
	double dist;
	long ID;

	const double maxd(m_radius);
	m_dumping.clear();
    m_dumping.reserve(getGeometry()->getNCells());
    m_dumping.setGeometry(getGeometry());
    m_dumping.setDataLocation(MPVLocation::CELL);

	for (auto it = patch_->cellBegin(); it!=patch_->cellEnd(); ++it){
		m_dumping.insert(it.getId(), 1.0);
	}

	if (!m_dumpingActive || m_decayFactor <= 1.E-12) return;

	//MODULATING DUMPING WITH DISTANCE
	MimmoObject * dumptarget= m_dsurface;
	if(m_dsurface == nullptr)  dumptarget = m_bsurface;

	bitpit::PiercedVector<double> distFactor = getGeometry()->getCellsNarrowBandToExtSurfaceWDist(*dumptarget, maxd);

	double distanceMax = std::pow((maxd/m_plateau), m_decayFactor);
	for(auto it = distFactor.begin(); it !=distFactor.end(); ++it){
		if(*it < m_plateau){
			(*it) = 1.0;
		}else{
			(*it) = (std::pow(maxd/(*it), m_decayFactor) -1.0) / (distanceMax -1.0);
		}
	}

    if(m_dumpingType != 1) { // no volume stuff, fill dumping with distance part only
        for(auto it=distFactor.begin(); it!=distFactor.end(); ++it){
            m_dumping.at(it.getId()) = (distanceMax - 1.0)*(*it) + 1.0;;
        }
        return; //you can return.
    }

    // Evaluating the volume part.
	bitpit::PiercedVectorStorage<double> volFactor;
	volFactor.setStaticKernel(&distFactor); //, bitpit::PiercedSyncMaster::SyncMode::SYNC_MODE_DISABLED);
	volFactor.fill(1.0);

	//evaluating cell volumes
    double locvol;
	double volmax = 0.0, volmin=1.0E18;
	for(auto it=distFactor.begin(); it!=distFactor.end(); ++it){
		locvol = getGeometry()->evalCellVolume(it.getId());
		if(locvol <= std::numeric_limits<float>::min()){
			(*m_log)<<"Warning in "<<m_name<<". Detected cells with almost zero or negative volume"<<std::endl;
            locvol = std::numeric_limits<float>::min(); //to assess myself around a 1.E-38 as minimum.
		}
        volFactor.rawAt(it.getRawIndex()) = locvol;
		volmin = std::min(volmin,locvol);
		volmax = std::max(volmax,locvol);
	}

	//evaluate the volume normalized function and store it in dumping.
	for(auto it = distFactor.begin(); it !=distFactor.end(); ++it){
        m_dumping.at(it.getId()) = std::pow(1.0 + (volmax -volmin)/volFactor.rawAt(it.getRawIndex()), *it);
	}

    //all done if you are here, you computed also the volume part and put together all the stuffs.
}

/*!
 * It updates an existent dumping function (artificial diffusivity) stored in m_dumping.
 * used before to modulate laplacian solution over the target mesh.
 * It will extract the pool of cells with diffusivity > 1.0. Thus, it will update
 * the narrow band inside the radius of influence, and recalculate distance and volumes eventually.
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::updateDumpingFunction(){
    if(!m_dumpingActive) return;

	bitpit::PatchKernel * patch_ = getGeometry()->getPatch();
	double dist;
	long ID;

    const double maxd(m_radius);
    //get the list of elements in m_dumping with diffusivity > 1.0;
    // at the same time reset the dumping function values to 1.0;
    livector1D seedlist;
    seedlist.reserve(m_dumping.size());
    for(auto it= m_dumping.begin(); it!=m_dumping.end(); ++it){
        if(*it > 1.0){
            seedlist.push_back(it.getId());
            *it = 1.0;
        }
    }

	//MODULATING DUMPING WITH DISTANCE
	MimmoObject * dumptarget= m_dsurface;
	if(m_dsurface == nullptr)  dumptarget = m_bsurface;

	bitpit::PiercedVector<double> distFactor = getGeometry()->getCellsNarrowBandToExtSurfaceWDist(*dumptarget, maxd, &seedlist);
    seedlist.clear();

	double distanceMax = std::pow((maxd/m_plateau), m_decayFactor);
	for(auto it = distFactor.begin(); it !=distFactor.end(); ++it){
		if(*it < m_plateau){
			(*it) = 1.0;
		}else{
			(*it) = (std::pow(maxd/(*it), m_decayFactor) -1.0) / (distanceMax -1.0);
		}
	}

    if(m_dumpingType != 1) { // no volume stuff, fill dumping with distance part only
        for(auto it=distFactor.begin(); it!=distFactor.end(); ++it){
            m_dumping.at(it.getId()) = (distanceMax - 1.0)*(*it) + 1.0;;
        }
        return; //you can return.
    }

    // Evaluating the volume part.
    bitpit::PiercedVectorStorage<double> volFactor;
    volFactor.setStaticKernel(&distFactor); //, bitpit::PiercedSyncMaster::SyncMode::SYNC_MODE_DISABLED);
    volFactor.fill(1.0);

    //evaluating cell volumes
    double locvol;
    double volmax = 0.0, volmin=1.0E18;
    for(auto it=distFactor.begin(); it!=distFactor.end(); ++it){
    	locvol = getGeometry()->evalCellVolume(it.getId());
    	if(locvol <= std::numeric_limits<float>::min()){
    		(*m_log)<<"Warning in "<<m_name<<". Detected cells with almost zero or negative volume"<<std::endl;
            locvol = std::numeric_limits<float>::min(); //to assess myself around a 1.E-38 as minimum.
    	}
        volFactor.rawAt(it.getRawIndex()) = locvol;
    	volmin = std::min(volmin,locvol);
    	volmax = std::max(volmax,locvol);
    }

    //evaluate the volume normalized function and store it in dumping.
    for(auto it = distFactor.begin(); it !=distFactor.end(); ++it){
        m_dumping.at(it.getId()) = std::pow(1.0 + (volmax -volmin)/volFactor.rawAt(it.getRawIndex()), *it);
    }
    //all done if you are here, you computed also the volume part and put together all the stuffs.
}

/*!
 * Prepare your system solver, feeding the laplacian stencils you previosly calculated
 * with FVolStencil::computeFVLaplacianStencil method.
 * Provide the map that get consecutive Index from Global Pierced vector Index system for CELLS
 * The stencil will be renumerated with the consecutiveIdIndexing provided.
 *
 * param[in] laplacianStencils pointer to MPV structure of laplacian stencils.
 * param[in] map of consecutive cell ID from Global PV indexing (typically get from MimmoObject::getMapcellInv)
 */

template<std::size_t NCOMP>
void
PropagateField<NCOMP>::initializeLaplaceSolver(FVolStencil::MPVDivergence * laplacianStencils, const liimap & maplocals){

    bitpit::KSPOptions &solverOptions = m_solver->getKSPOptions();

	solverOptions.rtol      = m_tol;
	solverOptions.subrtol   = m_tol;
    // total number of local DOFS, determines size of matrix
    long nDOFs = laplacianStencils->size();

    // total number of non-zero elements in the stencils.
    long nNZ(0);
    for(auto it=laplacianStencils->begin(); it!=laplacianStencils->end(); ++it){
        nNZ += it->size();
    }

#if MIMMO_ENABLE_MPI==1
    //instantiate the SparseMatrix
    bitpit::SparseMatrix matrix(m_communicator, getGeometry()->getPatch()->isPartitioned(), nDOFs, nDOFs, nNZ);
#else
    //instantiate the SparseMatrix
    bitpit::SparseMatrix matrix(nDOFs, nDOFs, nNZ);
#endif

    std::vector<long> mapsort(laplacianStencils->size());
    long id, ind;
    for(auto it=laplacianStencils->begin(); it!=laplacianStencils->end(); ++it){
    	id = it.getId();
    	ind = maplocals.at(id);
#if MIMMO_ENABLE_MPI
    	ind -= getGeometry()->getPatchInfo()->getCellGlobalCountOffset();
#endif
    	mapsort[ind] = id;
    }

    //Add ordered rows
    for(auto id : mapsort){
    	bitpit::StencilScalar item = laplacianStencils->at(id);
        //renumber values on the fly.
        item.renumber(maplocals);
    	matrix.addRow(item.size(), item.patternData(), item.weightData());
    }

    //assembly the matrix;
    matrix.assembly();

    //clean up the previous stuff in the solver.
    m_solver->clear();
    // now you can initialize the m_solver with this matrix.
    m_solver->initialize(matrix);

}

/*!
 * Update your system solver, feeding the cell based laplacian stencils you want to update in the matrix.
 * This method works with any valid subset of stencils in the mesh, but require the solver matrix to be initialized
 * and to have the new stencils with the same id pattern as they had at the time of the matrix initialization.
 * Provide the map that get consecutive Index from Global Pierced vector Index system for CELLS
 * The stencil will be renumerated with the consecutiveIdIndexing provided.
 *
 * param[in] laplacianStencils pointer to MPV structure of laplacian stencils subset to feed as update.
 * param[in] map of consecutive cell ID from Global PV indexing (typically get from MimmoObject::getMapcellInv)
 */

template<std::size_t NCOMP>
void
PropagateField<NCOMP>::updateLaplaceSolver(FVolStencil::MPVDivergence * laplacianStencils, const liimap & maplocals){

    // total number of local DOFS, determines size of matrix

// #if MIMMO_ENABLE_MPI==1
//     long nDOFs = m_solver->getColGlobalCount();
// #else
    long nDOFs = m_solver->getColCount();
// #endif
    long nupdate = laplacianStencils->size();

    // total number of non-zero elements in the stencils.
    long nNZ(0);
    for(auto it=laplacianStencils->begin(); it!=laplacianStencils->end(); ++it){
        nNZ += it->size();
    }

#if MIMMO_ENABLE_MPI==1
    //instantiate the SparseMatrix
    bitpit::SparseMatrix upelements(m_communicator, getGeometry()->getPatch()->isPartitioned(), nupdate, nDOFs, nNZ);
#else
    //instantiate the SparseMatrix
    bitpit::SparseMatrix upelements(nupdate, nDOFs, nNZ);
#endif


    // store the local ind of the rows involved, while filling the matrix of update values.
    std::vector<long> rows_involved;
    rows_involved.reserve(laplacianStencils->size());

    long id, ind;
    for(auto it=laplacianStencils->begin(); it != laplacianStencils->end(); ++it){
        id = it.getId();
    	ind = maplocals.at(id);
#if MIMMO_ENABLE_MPI
	    ind -= getGeometry()->getPatchInfo()->getCellGlobalCountOffset();
#endif
    	rows_involved.push_back(ind);
        bitpit::StencilScalar item(*it);
        item.renumber(maplocals);
        upelements.addRow(item.size(), item.patternData(), item.weightData());
    }
    //assembly the update matrix;
    upelements.assembly();

    //call the solver update;
    m_solver->update(rows_involved, upelements);
}

/*!
 * This method evaluate the bc corrections for a singular run of the system solver,
 * update the system matrix in m_solver and evaluate the rhs part due to bc.
 * After you call this method, you are typically ready to solve the laplacian system.
 * The type of bc @ interface are directly desumed from class map member m_isbp.
 * The method requires the Laplacian m_solver to be initialized. No ghost are taken into account.
 *
 * Moreover it needs as input :
 * \param[in] comp target component of bc conditions.
 * \param[in] unused boolean
 * \param[in] borderLaplacianStencil list of laplacian Stencil on border cells, where the bc is temporarely imposed as homogeneous Neumann
 * \param[in] borderCCGradientStencil list of Center cell gradient stencils defined on border cells.
 * \param[in] maplocals map from global id numbering to local system solver numbering.
 * \param[in,out] rhs, vector of right-hand-side's to append constant data from bc corrections.
   */

template<std::size_t NCOMP>
void
PropagateField<NCOMP>::assignBCAndEvaluateRHS(std::size_t comp, bool unused,
                                              FVolStencil::MPVDivergence * borderLaplacianStencil,
                                              FVolStencil::MPVGradient * borderCCGradientStencil,
                                              const liimap & maplocals,
                                              dvector1D & rhs)
{
    BITPIT_UNUSED(unused);
    //resize rhs to the number of internal cells
    MimmoObject * geo = getGeometry();
    rhs.resize(geo->getPatch()->getInternalCount(), 0.0);

    if (!m_solver->isInitialized()) {
        (*m_log)<<"Warning in "<<m_name<<". Unable to assign BC to the system. The solver is not yet initialized."<<std::endl;
        return;
    }

    if (!borderLaplacianStencil || !borderCCGradientStencil) {
        (*m_log)<<"Warning in "<<m_name<<". Unable to reach border cells stencils data. Nothing to do."<<std::endl;
        return;
    }


    //extract all interfaces from border cells.
    livector1D interfaceList = geo->getInterfaceFromCellList(borderLaplacianStencil->getIds());

    //copy laplacian stencils in a work mpv .
    FVolStencil::MPVDivergenceUPtr lapwork (new FVolStencil::MPVDivergence(*borderLaplacianStencil));

    //correct the original border laplacian stencils applying the Dirichlet/Neumann corrections.
    //renumber it and update the laplacian matrix and fill the rhs.
    bitpit::StencilVector correction;
    bitpit::PiercedVector<bitpit::Interface> & mesh_interfaces = geo->getInterfaces();
    bitpit::PiercedVector<bitpit::Cell> & mesh_cells = geo->getCells();

    long idOwner, idInterface;
    std::array<double,3> interfaceNormal, interfaceCentroid, ownerCentroid;
    double volume, iarea;
    for(long idInterface : interfaceList){

        correction.clear(false);

        if(!mesh_interfaces.at(idInterface).isBorder()) continue; //skip non-border interfaces;

        idOwner = mesh_interfaces.at(idInterface).getOwner();
        interfaceNormal = geo->evalInterfaceNormal(idInterface);
        volume = geo->evalCellVolume(idOwner);
        iarea = geo->evalInterfaceArea(idInterface);

        //apply the correction relative to bc @ interface.
        switch(m_isbp.at(idInterface)){
            case 1: //DIRICHLET
                    ownerCentroid = geo->evalCellCentroid(idOwner);
                    interfaceCentroid = geo->evalInterfaceCentroid(idInterface);
                    //get the correction.
                    correction = FVolStencil::correctionDirichletBCFaceGradient(m_bc_dir.at(idInterface)[comp],
                                                idOwner,ownerCentroid, interfaceCentroid, interfaceNormal,
                                                dotProduct(interfaceCentroid-ownerCentroid, interfaceNormal),
                                                borderCCGradientStencil->at(idOwner) );
                    break;
            default: //NEUMANN for non zero flux, change the first entry.
                    correction = FVolStencil::correctionNeumannBCFaceGradient(0.0, interfaceNormal);
                    break;
        }

        //calculate the laplacian correction and push it in work laplacian
        lapwork->at(idOwner) += (m_dumping.at(idOwner) * iarea / volume) * dotProduct(correction, interfaceNormal);

    }

    // now its time to update the solver matrix and to extract the rhs contributes.
    updateLaplaceSolver(lapwork.get(), maplocals);

    // now get the rhs
    for(auto it = lapwork->begin(); it != lapwork->end();++it){
        auto index = maplocals.at(it.getId());
        #if MIMMO_ENABLE_MPI
                //correct index if in parallel
            	index -= getGeometry()->getPatchInfo()->getCellGlobalCountOffset();
        #endif
        rhs[index] -= it->getConstant();
    }
}

/*!
 * It solves the laplacian problem, one field at a time.
 * Basically it is a wrapper to m_solver method solve();
 * Before calling this method be sure to have initialized the Laplacian linear system
 *
 * \param[in] rhs on internal cells;
 * \param[in,out] result on internal cells. In input is an initial solution, in output is the reusl of computation.
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::solveLaplace(const dvector1D &rhs, dvector1D &result){

    //just to be sure, resize the vector result
	result.resize(getGeometry()->getPatch()->getInternalCount(), 0.);

    // Check if the internal solver is initialized
    if (!m_solver->isInitialized()) {
        (*m_log)<<"Warning in "<<m_name<<". Unable to solve the system. The solver is not yet initialized."<<std::endl;
        return;
    }

    // Solve the system
	m_solver->solve(rhs, &result);

    //think I've done my job.
}

/*!
 * Utility to put laplacian solution into a Cell based MPV and after point interpolation
 * directly in m_field (cleared and refreshed).
 * Ghost communication is already taken into account in case of mpi run.
 * \param[in] results data of laplacian solutions collected in raw vectors.
 * \param[in] mapglobals map to retrieve cell global id from locals raw laplacian indexing
 * \param[out] markedcells (OPTIONAL) list of cells whose field norm is greater then m_tolerance.
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::reconstructResults(const dvector2D & results, const liimap & mapglobals, livector1D * markedcells)
{
    if(results.size() != NCOMP){
        (*m_log) << "WARNING in "<<m_name<<" . A field with dimension different from" <<NCOMP<<" is feeded to reconstructResults. m_field is not touched"<<std::endl;
        return;
    }
    // push result in a mpv linked to target mesh and on cell location.
    MimmoObject * geo = getGeometry();
    std::unique_ptr<MimmoPiercedVector<std::array<double,NCOMP> > > mpvres(new MimmoPiercedVector<std::array<double,NCOMP> >(geo, MPVLocation::CELL));
    mpvres->reserve(geo->getNCells());
    long id,counter;
    std::array<double, NCOMP> temp;
    for(auto & pair : mapglobals){
        id = pair.second;
        if (geo->getPatch()->getCell(id).isInterior()){
            counter = pair.first;
#if MIMMO_ENABLE_MPI
            counter -= geo->getPatchInfo()->getCellGlobalCountOffset();
#endif
            for(int i=0; i<NCOMP; ++i){
                temp[i] = results[i][counter];
            }
            mpvres->insert(id, temp);
        }
        else{
            for(int i=0; i<NCOMP; ++i){
                temp[i] = 0.0;
            }
            mpvres->insert(id, temp);
        }
    }

#if MIMMO_ENABLE_MPI
    // Creating ghost communications for exchanging solved values
    if (geo->getPatch()->isPartitioned()) {
        m_ghostStreamer = std::unique_ptr<MimmoDataBufferStreamer<NCOMP>>(new MimmoDataBufferStreamer<NCOMP>(mpvres.get()));
        m_ghostTag = createGhostCommunicator(true);
        m_ghostCommunicator->addData(m_ghostStreamer.get());
    }

    if (geo->getPatch()->isPartitioned()) {
        // Send data
        m_ghostCommunicator->startAllExchanges();
        // Receive data
        m_ghostCommunicator->completeAllExchanges();
    }
#endif

    if(markedcells){
        markedcells->clear();
        markedcells->reserve(mpvres->size());
        int counter = 0;
        for(auto it=mpvres->begin(); it != mpvres->end(); ++it){
            if(norm2(*it) > m_thres){
                markedcells->push_back(it.getId());
                counter++;
            }
        }
        markedcells->resize(counter);
    }

    // interpolate result to POINT location.
    m_field.clear();
    m_field = mpvres->cellDataToPointData();


#if MIMMO_ENABLE_MPI
    // Creating point ghost communications for exchanging interpolated values
    if (geo->getPatch()->isPartitioned()) {
    	//Force update exchange info
    	getGeometry()->updatePointGhostExchangeInfo();
        m_pointGhostStreamer = std::unique_ptr<MimmoPointDataBufferStreamer<NCOMP>>(new MimmoPointDataBufferStreamer<NCOMP>(mpvres.get()));
        m_pointGhostTag = createPointGhostCommunicator(true);
        m_pointGhostCommunicator->addData(m_pointGhostStreamer.get());
    }

    if (geo->getPatch()->isPartitioned()) {
        // Send data
        m_pointGhostCommunicator->startAllExchanges();
        // Receive data
        m_pointGhostCommunicator->completeAllExchanges();
    }
#endif

}



/*!
 * Initialize m_isbp member getting all the real border interfaces (no ghost) and
 * assigning to them a condition of type 0 (Neummann)
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::initializeBoundaryInfo(){
    //prepare the m_isbp;
    m_isbp.clear();
    // get all the real boundary interfaces (no ghost)
    livector1D list = m_geometry->extractBoundaryInterfaceID(false);
    m_isbp.reserve(list.size());
    //prepare m_isbp pushing neumann condition( 0 type ) to all interfaces of the list.
    for(const auto & val: list){
        m_isbp.insert(val, 0);
    }
}


#if MIMMO_ENABLE_MPI
/*!
    Creates a new ghost communicator.

    \param continuous defines if the communicator will be set in continuous mode
    \return The tag associated to the newly created communicator.
*/
template<std::size_t NCOMP>
int PropagateField<NCOMP>::createGhostCommunicator(bool continuous)
{
    // Create communicator
    GhostCommunicator *ghostCommunicator = new GhostCommunicator(getGeometry()->getPatch());
    ghostCommunicator->resetExchangeLists();
    ghostCommunicator->setRecvsContinuous(continuous);

    // Communicator tag
    int tag = ghostCommunicator->getTag();

    // Add ghost communicator
    m_ghostCommunicator = std::unique_ptr<GhostCommunicator>(ghostCommunicator);

    // Return communicator tag
    return tag;
}

/*!
    Creates a new point ghost communicator.

    \param continuous defines if the communicator will be set in continuous mode
    \return The tag associated to the newly created communicator.
*/
template<std::size_t NCOMP>
int PropagateField<NCOMP>::createPointGhostCommunicator(bool continuous)
{
    // Create communicator
    PointGhostCommunicator *ghostCommunicator = new PointGhostCommunicator(getGeometry());
    ghostCommunicator->resetExchangeLists();
    ghostCommunicator->setRecvsContinuous(continuous);

    // Communicator tag
    int tag = ghostCommunicator->getTag();

    // Add ghost communicator
    m_pointGhostCommunicator = std::unique_ptr<PointGhostCommunicator>(ghostCommunicator);

    // Return communicator tag
    return tag;
}

#endif


} //end of mimmo namespace
