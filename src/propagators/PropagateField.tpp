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
    this->m_np        = 0;
    this->m_gamma        = 1.0;
    this->m_laplace   = true;
    this->m_sstep     = 10;
    this->m_convergence = false;
    this->m_tol = 1.0e-12;
    this->m_bsurface  = NULL;
    this->m_dsurface  = NULL;
    this->m_geometry  = NULL;
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
    this->m_np           = other.m_np;
    this->m_gamma        = other.m_gamma;
    this->m_laplace      = other.m_laplace;
    this->m_sstep        = other.m_sstep;
    this->m_convergence  = other.m_convergence;
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
    std::swap(this->m_np, x.m_np);
    std::swap(this->m_gamma, x.m_gamma);
    std::swap(this->m_laplace, x.m_laplace);
    std::swap(this->m_sstep, x.m_sstep);
    std::swap(this->m_convergence, x.m_convergence);
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
    m_np = m_geometry->getNVertex();
}

/*! 
 * Sets the portion of boundary mesh relative to geometry target
 * that must be constrained with Dirichlet conditions.
 * \param[in] bsurface Boundary patch.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setDirichletBoundarySurface(MimmoObject* bsurface){
    if (bsurface == NULL)       return;
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
    if (bdumping == NULL)       return;
    if (bdumping->isEmpty())    return;
    if (bdumping->getType()!= 1 ) return;

    m_dsurface = bdumping;
}

/*! 
 * It sets the weight constant used in stencil computing. The constant (gamma) is used
 * to scale an exponential distance function taken as weight function in the stencil computing
 * for each point.
 * \param[in] gamma Weight constant.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setWeightConstant(double gamma){
    m_gamma = gamma;
}

/*!
 * It sets the number of steps of smoothing to propagate the boundary conditions over
 * the points cloud.
 * \param[in] ns Number of smoothing steps.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setSmoothingSteps(int ns){
    m_sstep = ns;
}

/*! 
 * It sets the solver used during the propagation of the surface constraint field.
 * \param[in] solveLaplacian true for Laplacian linear system solver, false for iterative smoothing technique.
 */
template <std::size_t NCOMP>
void PropagateField<NCOMP>::setSolver(bool solveLaplacian){
    m_laplace = solveLaplacian;
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
 * It sets if the smoothing solver must reach the convergence with prescribed tolerance set by setTolerance method.
 * \param[in] convergence Convergence flag.
 */
template <std::size_t NCOMP>
void PropagateField<NCOMP>::setConvergence(bool convergence){
    m_convergence = convergence;
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

    if(slotXML.hasOption("Solver")){
        std::string input = slotXML.get("Solver");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setSolver(value);
    }

    if(slotXML.hasOption("WeightConstant")){
        std::string input = slotXML.get("WeightConstant");
        input = bitpit::utils::string::trim(input);
        double value = 1.0;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
            value = std::fmax(0.0, value);
        }
        setWeightConstant(value);
    }

    if(slotXML.hasOption("SmoothingSteps")){
        std::string input = slotXML.get("SmoothingSteps");
        int value =1;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss>>value;
        }
        setSmoothingSteps(value);
    };

    if(slotXML.hasOption("Convergence")){
        std::string input = slotXML.get("Convergence");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setConvergence(value);
    }

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

    slotXML.set("WeightConstant",std::to_string(m_gamma));
    slotXML.set("Solver", std::to_string(int(m_laplace)));
    slotXML.set("SmoothingSteps",std::to_string(m_sstep));
    slotXML.set("Convergence",std::to_string(int(m_convergence)));
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
 * It computes the dumping function used in weights computation.
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::computeDumpingFunction(){

    bitpit::PatchKernel * patch_ = getGeometry()->getPatch();
    double dist;
    long ID;

    /* Maxdist should be the maximum distance between
     * boundaries with zero values and
     * boundaries with values different from zero.
     */
    //TODO compute it
    const double maxd(m_radius);
    m_dumping.clear();

    for (auto const & vertex : patch_->getVertices()){
        ID = vertex.getId();
        m_dumping.insert(ID, 1.0);
    }

    if (m_dumpingActive && m_decayFactor > 1.0e-12){

        //MODULATING DUMPING WITH DISTANCE
        MimmoObject * dumptarget= m_dsurface;
        if(m_dsurface == NULL)  dumptarget = m_bsurface;

        bitpit::PiercedVector<double> distFactor;
        getGeometry()->getVerticesNarrowBandToExtSurface(*dumptarget, maxd, distFactor);
        
        double distanceMax = std::pow((maxd/m_plateau), m_decayFactor);
        for(auto it = distFactor.begin(); it !=distFactor.end(); ++it){
            if(*it < m_plateau){
                (*it) = 1.0;
            }else{
                (*it) = (std::pow(maxd/(*it), m_decayFactor) -1.0) / (distanceMax -1.0);
            }
        }

        //evaluating volume on each vertex
        bitpit::PiercedVectorStorage<double> volFactor;
        volFactor.setStaticKernel(&distFactor); //, bitpit::PiercedSyncMaster::SyncMode::SYNC_MODE_DISABLED);
        volFactor.fill(1.0);

        if(m_dumpingType == 1){
            //evaluating cell volumes
            double volmax = 0.0, volmin=1.0E18;
            livector1D cellList = getGeometry()->getCellFromVertexList(distFactor.getIds(), false);
            bitpit::PiercedVector<double> volumes;
            volumes.reserve(cellList.size());
            //evaluate volumes on each cell and save it
            
            for(const auto & idC: cellList){
                auto it = volumes.insert(idC, getGeometry()->evalCellVolume(idC)); 
                if(*it <= 0.0){
                    throw std::runtime_error("Detected cells with negative volume");
                }
                volmin = std::min(volmin,*it);
                volmax = std::max(volmax,*it);
            }
            
            //pass on vertices assigning the min value of volume for each vertices.
            for(auto & id : cellList){
                bitpit::ConstProxyVector<long> vids = patch_->getCell(id).getVertexIds();
                for(auto & idV : vids){
                    if(distFactor.exists(idV))   volFactor[idV] = std::min(volFactor[idV], volumes[id]);
                }
            }
            
            //evaluate the volume normalized function
            for(auto it = distFactor.begin(); it !=distFactor.end(); ++it){
                std::size_t pos = it.getRawIndex();
                volFactor.rawAt(pos) = std::pow(1.0 + (volmax -volmin)/volFactor.rawAt(pos), distFactor.rawAt(pos));
            }
        }

        
        //get an average of distance and volume functions.
        for(auto it=distFactor.begin(); it!=distFactor.end(); ++it){
            long id = it.getId();
            double val;
            if(m_dumpingType == 1){
                val = volFactor[id];
            }else{
                val = (distanceMax - 1.0)*distFactor[id] + 1.0;
            }
            m_dumping[id] = val;
        }
    }
}


/*!
 * It applies a smoothing filter for a defined number of step.
 * \param[in] nstep desired number of smoothing steps
 * \param[in] stencils stencil-ids of laplace operator on target mesh nodes
 * \param[in] weights  associated to stencils
 * \param[in] dataInv map of local node indexing of stencils vs global mesh node indexing.
 *                    If map is empty, the method itself provides its evaluation from class target mesh.
 * \param[out] field resulting field after smoothing
 */
template<std::size_t NCOMP>
void 
PropagateField<NCOMP>::solveSmoothing(int nstep,
                               ivector2D &stencils,
                               dvector2D &weights,
                               dvector1D &rhs,
                               liimap & dataInv,
                               MimmoPiercedVector<std::array<double, NCOMP> > & field)
{
    long ID;
    int ind;
    field.clear();

    //initialize field
    field.clear();
    for (auto vertex : getGeometry()->getVertices()){
        ID = vertex.getId();
        field.insert(ID, std::array<double, NCOMP>({}));
    }
    
    //found maxval in rhs for normalization purposes
    double maxval = 0.0;
    for (const auto & val : rhs){
        maxval = std::max(maxval, val);
    }
    
    dvector1D guess = rhs;
    
    (*m_log)<< m_name <<" starts field propagation."<<std::endl;
    
    for (int istep = 0; istep < nstep; istep++){

        if (!m_convergence) (*m_log)<<m_name << " smoothing step : " << istep+1 << " / " << nstep <<std::endl;

        double maxdiff = 0.0;
        for (int locId=0; locId<NCOMP*m_np; ++locId){

            double delta = 0.0;
            std::size_t sizeStencil = stencils[locId].size();
            for (std::size_t i=0; i<sizeStencil; ++i){
                delta += -1.0*guess[stencils[locId][i]]*weights[locId][i];
            }
            delta += rhs[locId];
            guess[locId] += delta;
            
            if (m_convergence){
                maxdiff = std::max(maxdiff, delta);
            }
        }//end for vertex

        if (m_convergence) (*m_log)<< m_name<<" residual : " << maxdiff <<std::endl;

        //convergence
        if (m_convergence){
            if (maxdiff <= m_tol)
                istep = nstep;
            else{
                nstep = istep+2;
            }
        }
    }// end step

    (*m_log)<< m_name<<" ends field propagation."<<std::endl;

    for (auto vertex : getGeometry()->getVertices()){
        for (int icomp=0; icomp<NCOMP; ++icomp ){ 
            ID = vertex.getId();
            ind = dataInv[ID];
            field[ID][icomp] = guess[ind + icomp*m_np];
        }
    }
    
    field.setDataLocation(MPVLocation::POINT);
    field.setGeometry(getGeometry());
}

/*!
 * It solves the laplacian problem. Stencils, weights and rhs must be already corrected to 
 * account of boundary condition of the problem. See calculateStencilsLaplace and calculateRHSLaplace method.
 * 
 * \param[in] stencils stencil-ids of laplace operator on target mesh nodes
 * \param[in] weights  associated to stencils
 * \param[in] rhs right-hand-side of laplacian linear system 
 * \param[in] dataInv map of local node indexing of stencils vs global mesh node indexing.
 *                    If map is empty, the method itself provides its evaluation from class target mesh.
 * \param[out] field target solution
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::solveLaplace(ivector2D &stencils,
                             dvector2D &weights,
                             dvector1D &rhs,
                             liimap &dataInv,
                             MimmoPiercedVector<std::array<double, NCOMP> > & field)
{
    //initialize field
    field.clear();
    for (auto vertex : getGeometry()->getVertices()){
        long int ID = vertex.getId();
        field.insert(ID, std::array<double, NCOMP>({}));
    }

    // Create the system for solving the pressure
    m_solver = std::unique_ptr<mimmo::SystemSolver>(new mimmo::SystemSolver(false));

    // Initialize the system
    KSPOptions &solverOptions = m_solver->getKSPOptions();
    solverOptions.nullspace = false;
    solverOptions.rtol      = m_tol;
    solverOptions.subrtol   = m_tol;

#if ENABLE_MPI==1
        m_solver->initialize(stencils, weights, rhs, ghosts);
#else
        m_solver->initialize(stencils, weights, rhs);
#endif

    // Solve the system
    m_solver->solve();

    // Get the solution
    const double *solution = m_solver->getSolutionRawReadPtr();

    long ID;
    int ind;
    for (auto vertex : getGeometry()->getVertices()){
        for (int icomp=0; icomp<NCOMP; ++icomp ){ 
            ID = vertex.getId();
            ind = dataInv[ID];
            field[ID][icomp] = solution[ind + icomp*m_np];
        }
    }

    // Clear the solver
    m_solver->clear();
    field.setDataLocation(MPVLocation::POINT);
    field.setGeometry(getGeometry());
}

/*!
 * Given the target geometry mesh, evaluate the stencils and the weights of laplacian operator.
 * Boundary conditions corrections on the laplacian operator are not applied.
 * Please note the diagonal element is always on the end of each stencils/weights subvector
 * 
 * \param[in] dataInv map of local node indexing of stencils vs global mesh node indexing.
 *                    If map is empty, the method itself provides its evaluation from class target mesh.
 * \param[out] stencils stencil-ids of laplace operator on target mesh nodes
 * \param[out] weights  associated to stencils
 * 
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::computeStencils(liimap &dataInv,
                                  ivector2D &stencils,
                                  dvector2D &weights)
{
    MimmoPiercedVector<livector1D> conn; 
    MimmoPiercedVector<dvector1D>  wgt;

    computeRingConnectivity(conn);
    computeRingWeights(conn,wgt);
    
    stencils.clear();
    weights.clear();
    stencils.resize(NCOMP*m_np);
    weights.resize(NCOMP*m_np);

    long ID;
    int ind;
    //Create stencils 
    for (auto vertex : getGeometry()->getVertices()){
        ID = vertex.getId();
        ind = dataInv[ID];
        //bulk evaluation
        for(int comp=0; comp<NCOMP; ++comp){
            for (const long & IDN : conn[ID]){
                stencils[ind+comp*m_np].push_back(dataInv[IDN]+comp*m_np);
            }
            stencils[ind+comp*m_np].push_back(ind+comp*m_np);
            weights[ind+comp*m_np] = -1.0*wgt[ID];
            weights[ind+comp*m_np].push_back(1.0);
        }
    }

}

/*!
 * It computes the vertex ring connectivity structure for each mesh vertex.
 * \param[out] conn ring connectivity
 */
template<std::size_t NCOMP>
void PropagateField<NCOMP>::computeRingConnectivity(MimmoPiercedVector<livector1D> & conn){
    //DOESN'T WORK FOR POINTS CLOUD OR CURVES
    
    bitpit::PatchKernel * patch_ = getGeometry()->getPatch();

    conn.reserve(m_np);

    //map for edges already visited
    std::map<std::pair<long, long>, bool> visitedge;
    
    //compute connectivity between points
    for (const auto & cell : patch_->getCells()){
        int edgecount;
        //         if (type == 1 || patch_->getDimension() == 2){
        //             edgecount = cell.getFaceCount();
        //         }
        //         else{
        edgecount = cell.getEdgeCount();
        //         }
        for (int ie=0; ie<edgecount; ++ie){
            bitpit::ConstProxyVector<long> econn;
            //             if (type == 1 || patch_->getDimension() == 2){
            //                 econn = cell.getFaceConnect(ie);
            //             }
            //             else{
            econn = cell.getEdgeConnect(ie);
            //             }
            if (!visitedge[std::pair<long, long>(econn[0],econn[1])]){
                
                visitedge[std::pair<long, long>(econn[0],econn[1])] = true;
                visitedge[std::pair<long, long>(econn[1],econn[0])] = true;
                
                if (!conn.exists(econn[0])){
                    conn.insert(econn[0], livector1D(1, econn[1]));
                }else{
                    conn[econn[0]].push_back(econn[1]);
                }
                
                if (!conn.exists(econn[1])){
                    conn.insert(econn[1], livector1D(1, econn[0]));
                }else{
                    conn[econn[1]].push_back(econn[0]);
                }
            }
        }
    }
}

/*! 
 * It computes the weights associated to the vertex ring connectivity conn.
 * \param[in] conn vertex ring connectivity
 * \param[out] wgt weights associated.
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::computeRingWeights(MimmoPiercedVector<livector1D> & conn, 
                                   MimmoPiercedVector<dvector1D> & wgt)
{
    bitpit::PatchKernel * patch_ = getGeometry()->getPatch();
    
    int nsize;
    livector1D ids;
    double dist;
    double sumdist;
    
    wgt.reserve(conn.size());
    darray3E point;
    dvector1D lweights;
    for (long ID : conn.getIds()){
        ids = conn[ID];
        point = patch_->getVertex(ID).getCoords();
        nsize = ids.size();
        lweights.resize(nsize);
        sumdist = 0.0;
        for (int j=0; j<nsize; j++){
            dist = norm2(point-patch_->getVertex(ids[j]).getCoords());
            lweights[j] = m_dumping[ids[j]] / (std::pow(dist, m_gamma));
            sumdist += lweights[j];
        }
        lweights /= sumdist;
        wgt.insert(ID, lweights);
    }
    if (!m_execPlot) m_dumping.clear();
}

}
