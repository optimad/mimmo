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
	this->m_tol       = 1.0e-12;
	this->m_bsurface  = nullptr;
	this->m_dsurface  = nullptr;
	this->m_geometry  = nullptr;
	this->m_decayFactor = 1.0;
	this->m_dumping.clear();
	this->m_radius = 0.0;
	this->m_plateau = 0.0;
	this->m_dumpingActive = false;
	this->m_dumpingType = 0;
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
	this->m_dumpingType = other.m_dumpingType;
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
		double value = 1.0e-12;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
			value = std::fmax(0.0, value);
		}
		setTolerance(value);
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

	bitpit::PiercedVector<double> distFactor;
	getGeometry()->getCellsNarrowBandToExtSurface(*dumptarget, maxd, distFactor);

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
		if(locvol <= 0.0){
			throw std::runtime_error("Detected cells with zero or negative volume");
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
 * The method evaluates the Gradient stencils at each interface of the mesh. If
 * Interfaces are not built for the current mesh, it builds them.
 *
 * These interface gradient stencils are used to assembly the laplacian stencils.
 * Dirichlet and Neumann Boundary conditions are directly provided on their corrispective
 * border gradient stencil with a neutral value of 1.0 (Dirichlet field value or normal Neumann flux value).
 * Since this part does not affect the proper stencil pattern and is stored into it
 * as Constant, you can put unitary boundary conditions and fix the correct Dirichlet
 * or Neumann flux value later, by multiplying it to the right stencil Constant.
 * \param[in] alternativeData unused
 */

template<std::size_t NCOMP>
FVolStencil::MPVGradientUPtr
PropagateField<NCOMP>::computeGradientStencilsWithNeutralBC(FVolStencil::MPVGradientUPtr alternativeData){

    BITPIT_UNUSED(alternativeData);
    //start to evaluate the Gradient Stencil @ Interfaces.
    FVolStencil::MPVGradientUPtr gradients = FVolStencil::computeFVFaceGradientStencil(*(getGeometry()));

    //The gradient stencil @ boundary interfaces is the CCell Gradient stencil of the owner up to now.
    // You need to modify it to take into account of boundary conditions.
    // We do not push for now the proper bc condition, but impose a Neutral value 1.0.
    // Once the boundary interface Gradient stencil is built we can impose the right value
    // multiplying it to the constant part of the stencil (the only one influenced by it).

    //So loop on m_isbp : for each boundary interface get its type;

    long id, ownerID;
    std::array<double,3> interfaceNormal, interfaceCentroid, ownerCentroid;
    MimmoObject * geo = getGeometry();
    bitpit::PiercedVector<bitpit::Interface> & interfaces = geo->getInterfaces();

    for(auto it = m_isbp.begin(); it!=m_isbp.end(); ++it){
        id = it.getId();
        interfaceNormal = geo->evalInterfaceNormal(id);
        switch(*it){
            case 1: //dirichlet
                ownerID = interfaces.at(id).getOwner();
                ownerCentroid = geo->evalCellCentroid(ownerID);
                interfaceCentroid = geo->evalInterfaceCentroid(id);

                gradients->at(id) = FVolStencil::computeDirichletBCFaceGradient(1.0, ownerID,
                                   ownerCentroid, interfaceCentroid, interfaceNormal,
                                   dotProduct(interfaceCentroid-ownerCentroid, interfaceNormal),
                                   gradients->at(id) );
                break;
            default: //neumann
                gradients->at(id) = FVolStencil::computeNeumannBCFaceGradient(1.0, interfaceNormal, gradients->at(id));
                break;
        }
    }

    return gradients;
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

    //fill in the stencil rows, after renumbering the stencil.
    for(auto it=laplacianStencils->begin(); it!=laplacianStencils->end(); ++it){
        it->renumber(maplocals);
    }

    //Build an aux map instead of sorting laplacianStencils
//    laplacianStencils->sort();
    std::vector<long> mapsort(laplacianStencils->size());
    for(auto it=laplacianStencils->begin(); it!=laplacianStencils->end(); ++it){
    	long id = it.getId();
    	long ind = maplocals.at(id);
    	mapsort[ind - getGeometry()->getPatchInfo()->getCellGlobalCountOffset()] = id;
    }

//    //Add ordered rows
//    for(auto it=laplacianStencils->begin(); it!=laplacianStencils->end(); ++it){
//    	matrix.addRow(it->size(), it->patternData(), it->weightData());
//    }
    //Add ordered rows
    for(auto id : mapsort){
    	auto item = laplacianStencils->at(id);
    	matrix.addRow(item.size(), item.patternData(), item.weightData());
    }

    //assembly the matrix;
    matrix.assembly();

    //clean up the previous stuff in the solver.
    m_solver->clear();
    // now you can initialize the m_solver with this matrix.
    m_solver->initialize(matrix);

//    //dump matrix
//    m_solver->dump(".", "system");

}

/*!
 * Append to RHS vector the part of due to boundary conditions of the field with component comp.
 * The method does the following:
 * - Start from the faceGradient with Neutral conditions 1.0;
 * - Save their constant part;
 * - Get the real bc's (on the target component comp) and multiply them to faceGradient stencil constant part
 * - Recalculate the laplacian stencils on all the cell influenced by the border with
 *   method FVolStencil::computeFVBorderLaplacianStencil
 * - get the constant part only of these laplacian stencils and sum them to the consecutive ordered vector rhs position.
 * - restore faceGradient old neutral constant.
 *
 * In this way the rhs is ready to be used by the solver
 *
 * \param[in] comp target component of bc conditions.
 * \param[in] faceGradientStencils pointer to list of faceGradientStencils with bc Neutral
 * \param[in] maplocals map from GlobalIndex PV to local consecutive index
 * \param[in,out] rhs, vector of rhs to append data
   */

template<std::size_t NCOMP>
void
PropagateField<NCOMP>::appendToRHSFromBorderFluxes(std::size_t comp,
                                                   FVolStencil::MPVGradient * faceGradientStencils,
                                                   const liimap & maplocals,
                                                   dvector1D & rhs)
{

    if(comp >= NCOMP){
        std::string message = "Error in PropagateField::appendToRHSFromBorderFlux comp can be at most " + std::to_string(NCOMP-1);
        throw std::runtime_error(message.c_str());
    }

    //resize rhs to the number of internal cells
    rhs.resize(getGeometry()->getPatch()->getInternalCount(), 0.0);

    //store the constant part of stencils relative to border interfaces.
    // and update it with the right boundary condition inside m_bc_dir.
    std::vector<std::array<double,3> > oldConstantStore;

    oldConstantStore.reserve(m_isbp.size());

    for(auto it = m_isbp.begin(); it != m_isbp.end(); ++it){
        std::array<double,3> &temp = faceGradientStencils->at(it.getId()).getConstant();
        oldConstantStore.push_back(temp);
        //update this unitary constant value multiplying the real bc.
        switch(*it){
            case 1: //dirichlet
                temp *= m_bc_dir.at(it.getId())[comp];
                break;
            default : //neumann homogeneous
                temp *= 0.0;
                break;
        }
    }

    // now recalculate the boundary stencils of laplacian:
    FVolStencil::MPVDivergenceUPtr borderlap =
                FVolStencil::computeFVBorderLaplacianStencil(*faceGradientStencils, &m_dumping);

    //Subtract their constant part from rhs.
    for(auto it= borderlap->begin(); it!=borderlap->end(); ++it){
        auto index = maplocals.at(it.getId()) - getGeometry()->getPatchInfo()->getCellGlobalCountOffset();
        rhs[index] += -1.0 * it->getConstant();
    }

    //Return to the old neutral value on border stencils.
    std::vector<std::array<double,3>>::iterator itold = oldConstantStore.begin();
    for(auto it = m_isbp.begin(); it != m_isbp.end(); ++it){
        faceGradientStencils->at(it.getId()).getConstant() = *itold;
        ++itold;
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

    // Check if the internal solver is initialized
    if (!m_solver->isInitialized()) {
        throw std::runtime_error("Unable to solve the system. The solver is not yet initialized.");
    }

    //just to be sure, resize the vector result
	result.resize(getGeometry()->getPatch()->getInternalCount(), 0.);
    // Solve the system
	m_solver->solve(rhs, &result);

    //think I've done my job.
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

#endif


} //end of mimmo namespace
