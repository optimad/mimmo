
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

#include "PropagateField.hpp"
#include "customOperators.hpp"
#include "SkdTreeUtils.hpp"

namespace mimmo{

/*!
 * Constructor
 */
PropagateField::PropagateField(){
    m_name = "mimmo.PropagateField";
    setDefaults();
};

/*!
 * Set most significant parameters to constructor defaults
 */
void PropagateField::setDefaults(){
    m_isbp.clear();
    m_np        = 0;
    m_nbp        = 0;
    m_conn.clear();
    m_gamma        = 1.0;
    m_weights.clear();
    m_laplace   = true;
    m_sstep     = 10;
    m_convergence = false;
    m_tol = 1.0e-05;
    m_bsurface  = NULL;
    m_dsurface  = NULL;
    m_geometry  = NULL;
    m_decayFactor = 1.0;
    m_dumping.clear();
    m_radius = 0.0;
    m_plateau = 0.0;
    m_dumpingActive = false;
}

/*!
 * Destructor;
 */
PropagateField::~PropagateField(){
    clear();
};

/*!
 * Copy constructor
 */
PropagateField::PropagateField(const PropagateField & other):BaseManipulation(other){
    setDefaults();
    m_isbp         = other.m_isbp;
    m_np           = other.m_np;
    m_nbp          = other.m_nbp;
    m_conn         = other.m_conn;
    m_gamma        = other.m_gamma;
    m_weights      = other.m_weights;
    m_laplace      = other.m_laplace;
    m_sstep        = other.m_sstep;
    m_convergence  = other.m_convergence;
    m_tol          = other.m_tol;
    m_bsurface     = other.m_bsurface;
    m_dsurface     = other.m_dsurface;
    m_dumping      = other.m_dumping;
    m_decayFactor  = other.m_decayFactor;
    m_radius       = other.m_radius;
    m_plateau      = other.m_plateau;
    m_dumpingActive= other.m_dumpingActive;

};

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void PropagateField::swap(PropagateField & x) noexcept {
    m_isbp.swap(x.m_isbp);
    std::swap(m_np, x.m_np);
    std::swap(m_nbp, x.m_nbp);
    m_conn.swap(x.m_conn);
    std::swap(m_gamma, x.m_gamma);
    m_weights.swap(x.m_weights);
    std::swap(m_laplace, x.m_laplace);
    std::swap(m_sstep, x.m_sstep);
    std::swap(m_convergence, x.m_convergence);
    std::swap(m_tol, x.m_tol);
    std::swap(m_bsurface, x.m_bsurface);
    std::swap(m_dsurface, x.m_dsurface);
    std::swap(m_decayFactor, x.m_decayFactor);
    m_dumping.swap(x.m_dumping);
    std::swap(m_radius, x.m_radius);
    std::swap(m_plateau, x.m_plateau);
    std::swap(m_dumpingActive, x.m_dumpingActive);
    BaseManipulation::swap(x);
}

/*! 
 * It builds the input/output ports of the object
 */
void
PropagateField::buildPorts(){
    bool built = true;
    built = (built && createPortIn<MimmoObject*, PropagateField>(this, &PropagateField::setGeometry, M_GEOM, true));
    built = (built && createPortIn<MimmoObject*, PropagateField>(this, &PropagateField::setDirichletBoundarySurface, M_GEOM2));
    built = (built && createPortIn<MimmoObject*, PropagateField>(this, &PropagateField::setDumpingBoundarySurface, M_GEOM3));
    m_arePortsBuilt = built;
};

/*!
 * It gets the number of points in the cloud.
 * \return Number of points.
 */
int
PropagateField::getNPoints(){
    return(m_np);
}

/*! 
 * It gets the number of boundary points in the cloud.
 * \return Number of boundary points.
 */
int
PropagateField::getNBoundaryPoints(){
    return(m_nbp);
}

/*!
 * Set pointer to your target bulk volume geometry. Reimplemented from mimmo::BaseManipulation::setGeometry().
 * Geometry must be a of volume type (MimmoObject type = 2);
 * \param[in] geometry_ pointer to target geometry
 */
void
PropagateField::setGeometry(MimmoObject * geometry_){

    if (geometry_ == NULL) return;
    if (geometry_->isEmpty())   return;
    if ( geometry_->getType()!= 2 ) return;

    m_geometry = geometry_;
}

/*! 
 * Sets the portion of boundary mesh relative to geometry target
 * that must be constrained with Dirichlet conditions.
 * \param[in] bsurface Boundary patch.
 */
void
PropagateField::setDirichletBoundarySurface(MimmoObject* bsurface){
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
void
PropagateField::setDumpingBoundarySurface(MimmoObject* bdumping){
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
void
PropagateField::setWeightConstant(double gamma){
    m_gamma = gamma;
}

/*!
 * It sets the number of steps of smoothing to propagate the boundary conditions over
 * the points cloud.
 * \param[in] ns Number of smoothing steps.
 */
void
PropagateField::setSmoothingSteps(int ns){
    m_sstep = ns;
}

/*! 
 * It sets the solver used during the propagation of the surface constraint field.
 * \param[in] solveLaplacian true for Laplacian linear system solver, false for iterative smoothing technique.
 */
void PropagateField::setSolver(bool solveLaplacian){
    m_laplace = solveLaplacian;
}

/*!
 * Activate dumping control of artificial diffusivity(see class doc).
 * \param[in] flag boolean true activate, false deactivate.
 */
void
PropagateField::setDumping(bool flag){
    m_dumpingActive = flag;
}


/*!
 * Set the inner dumping radius p (see class doc).
 * \param[in] plateau inner distance.
 */
void
PropagateField::setDumpingInnerDistance(double plateau){
    m_plateau = plateau;
}

/*!
 * Set the outer dumping distance r (see class doc).
 * \param[in] radius outer distance.
 */
void
PropagateField::setDumpingOuterDistance(double radius){
    m_radius = radius;
}

/*!
 * Set the dumping factor.
 * \param[in] dump Exponential of dumping function.
 */
void
PropagateField::setDecayFactor(double decay){
    m_decayFactor = decay;
}

/*! 
 * It sets if the solver must reach the convergence with prescribed tolerance set by setTolerance method.
 * \param[in] convergence Convergence flag.
 */
void PropagateField::setConvergence(bool convergence){
    m_convergence = convergence;
}

/*!
 * It sets the tolerance on residuals for solver convergence.
 * \param[in] tol Convergence tolerance.
 */
void PropagateField::setTolerance(double tol){
    m_tol = tol;
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){


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
        double value = 1.0e-05;
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
    }
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

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
    }
};



/*!
 * Restore data as in class default construction.
 */
void
PropagateField::clear(){
    BaseManipulation::clear();
    setDefaults();
};

/*!
 * Execution command. After the execution the result constraint field is stored in the class.
 */
void
PropagateField::execute(){

    if(getGeometry() == NULL){
        (*m_log)<<"Error in "<<m_name<<" .No target volume mesh linked"<<std::endl;
        throw std::runtime_error("Error in PropagateField execute. No target volume mesh linked");
    }

    if(m_bsurface == NULL ){
        (*m_log)<<"Error in "<<m_name<<" .No Dirichlet Boundary patch linked"<<std::endl;
        throw std::runtime_error("Error in PropagateField execute. No Dirichlet Boundary patch linked");
    }

    if(!checkBoundariesCoherence()){
        (*m_log)<<"Error in "<<m_name<<" .Boundary patches linked are uncoherent with target bulk geometry"
                "or bc-fields not coherent with boundary patches"<<std::endl;
        throw std::runtime_error("Error in PropagateField execute. Boundary patches linked are uncoherent" 
                "with target bulk geometry or bc-fields not coherent with boundary patches");
    }

    computeConnectivity();
    computeDumpingFunction();
    computeWeights();

    if (m_laplace){
        solveLaplace();
    }else{
        solveSmoothing(m_sstep);
    }
}

/*! 
 * It computes the connectivity structure between input points.
 */
void
PropagateField::computeConnectivity(){

    //DOESN'T WORK FOR POINTS CLOUD OR CURVES

    bitpit::PatchKernel * patch_ = getGeometry()->getPatch();
    //     int type = getGeometry()->getType();

    m_conn.resize(m_np);
    m_weights.resize(m_np);

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

                if (!m_conn.exists(econn[0])){
                    m_conn.insert(econn[0], livector1D(1, econn[1]));
                }else{
                    m_conn[econn[0]].push_back(econn[1]);
                }

                if (!m_conn.exists(econn[1])){
                    m_conn.insert(econn[1], livector1D(1, econn[0]));
                }else{
                    m_conn[econn[1]].push_back(econn[0]);
                }
            }
        }
    }
}

/*! 
 * It computes the weight values used in stencil computation.
 */
void
PropagateField::computeWeights(){

    bitpit::PatchKernel * patch_ = getGeometry()->getPatch();

    int nsize;
    livector1D ids;
    double dist;
    double sumdist;

    m_weights.clear();
    darray3E point;
    dvector1D lweights;
    for (long ID : m_conn.getIds()){
        ids = m_conn[ID];
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
        m_weights.insert(ID, lweights);
    }
    if (!m_execPlot) m_dumping.clear();
}

/*!
 * It computes the dumping function used in weights computation.
 */
void
PropagateField::computeDumpingFunction(){

    bitpit::PatchKernel * patch_ = getGeometry()->getPatch();
    double dist;
    long ID;

    /* Maxdist should be the maximum distance between
     * boundaries with zero values and
     * boundaries with values different from zero.
     */
    //TODO compute it
    double maxd = m_radius;

    m_dumping.clear();

    for (auto const & vertex : patch_->getVertices()){
        ID = vertex.getId();
        m_dumping.insert(ID, 1.0);
    }

    if (m_dumpingActive && m_decayFactor > 1.0e-12){

        //MODULATING DUMPING WITH DISTANCE
        MimmoObject * dumptarget= m_dsurface;
        if(m_dsurface == NULL)  dumptarget = m_bsurface;

        if(!dumptarget->isSkdTreeSync()){
            dumptarget->buildSkdTree();
        }
        auto tree = dumptarget->getSkdTree();

        double valmax1 = 1./std::pow((maxd/m_plateau), m_decayFactor);
        long idsupp;
        darray3E point;
        double val;
        double dist;
        std::map<long, double> eta;
        bitpit::PiercedVector<double> arguments;
        std::vector<long> vertices;
        vertices.reserve(patch_->getVertexCount());
        for(auto const & vertex: patch_->getVertices()){
            ID = vertex.getId();
            point = vertex.getCoords();
            dist = skdTreeUtils::distance(&point, tree, idsupp, maxd);
            if (dist < maxd){
                if(dist <= m_plateau){
                    eta[ID] = 1.0;
                }else{
                    eta[ID] = (std::pow(maxd/dist, m_decayFactor) -1.0)*valmax1;
                }
                vertices.push_back(ID);
                arguments.insert(ID, 0.0);
            }
        }
        vertices.resize(arguments.size());
        std::vector<long> cellids = getGeometry()->getCellFromVertexList(vertices);

        //prepare information on volumes and aspectRatio.
         bitpit::PiercedVector<double> volumes, AR;
         for (const long & idc : cellids){
             volumes.insert(idc, getGeometry()->evalCellVolume(idc));
             AR.insert(idc, getGeometry()->evalCellAspectRatio(idc));
         }

         double volmax = 0.0;
         double volmin = 1.0E18;
         double armax = 0.0;
         double armin = 1.0E18;
         for(const auto & val : volumes){
             volmax = std::max(volmax, val);
             volmin = std::min(volmin, val);
         }
         for(const auto & val : AR){
             armax = std::max(armax, val);
             armin = std::min(armin, val);
         }

         double volmaxmin = (volmax - volmin);
         double armin_max = (1./armin - 1./armax);

         for(auto it=volumes.begin(); it!=volumes.end(); ++it){
             *it = (1.0 + volmaxmin/(*it)*armin_max*AR[it.getId()]); //*AR[it.getId()];
         }

        for(auto & id : cellids){
            bitpit::Cell cell = patch_->getCell(id);
            bitpit::ConstProxyVector<long> vids = cell.getVertexIds();
            for(auto & idV : vids){
                if(eta.count(idV)){
                    arguments[idV] = std::max(arguments[idV], volumes[id]);
                }
            }
        }
        for(auto it=eta.begin(); it!=eta.end(); ++it){
            m_dumping[it->first] = std::pow(arguments[it->first], it->second);
        }
    }
}

//--------------------------------------
//--------------------------------------
// SCALARFIELD (USUALLY FILTER DISPLACEMENTS)
//--------------------------------------
//--------------------------------------

/*!
 * Constructor
 */
PropagateScalarField::PropagateScalarField():PropagateField(){
    m_name = "mimmo.PropagateScalarField";
    setDefaults();
};

/*!
 * Set most significant parameters to constructor defaults
 */
void PropagateScalarField::setDefaults(){
    PropagateField::setDefaults();
    m_bc_dir.clear();
    m_field.clear();

}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
PropagateScalarField::PropagateScalarField(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.PropagateScalarField";
    setDefaults();

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.PropagateScalarField"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Destructor;
 */
PropagateScalarField::~PropagateScalarField(){
    clear();
};

/*!
 * Copy constructor
 */
PropagateScalarField::PropagateScalarField(const PropagateScalarField & other):PropagateField(other){
    m_bc_dir    = other.m_bc_dir;
    m_field     = other.m_field;
};

/*!
 * Assignment operator of the class
 */
PropagateScalarField & PropagateScalarField::operator=(PropagateScalarField other){
    swap(other);
    return *this;
};

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void PropagateScalarField::swap(PropagateScalarField & x) noexcept {
    m_bc_dir.swap(x.m_bc_dir);
    m_field.swap(x.m_field);
    PropagateField::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
PropagateScalarField::buildPorts(){
    bool built = true;
    built = (built && createPortIn<dmpvector1D, PropagateScalarField>(this, &PropagateScalarField::setDirichletConditions, M_FILTER));
    built = (built && createPortOut<dmpvector1D, PropagateScalarField>(this, &PropagateScalarField::getPropagatedField, M_FILTER));
    PropagateField::buildPorts();
    m_arePortsBuilt = built;
};

/*!
 * It gets the resulting propagated field on the whole bulk mesh.
 * \return Deformation field.
 */
dmpvector1D
PropagateScalarField::getPropagatedField(){
    return(m_field);
}

/*!
 * It sets the Dirichlet conditions for scalar field on the previously linked
 * Dirichlet Boundary patch.
 * \param[in] bc dirichlet conditions
 */
void
PropagateScalarField::setDirichletConditions(dmpvector1D bc){
    if (bc.isEmpty()) return;
    m_bc_dir = bc;
}


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateScalarField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
    BITPIT_UNUSED(name);
    //start absorbing
    PropagateField::absorbSectionXML(slotXML, name);
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateScalarField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    PropagateField::flushSectionXML(slotXML, name);
};

/*!
 * Clear all data actually stored in the class
 */
void
PropagateScalarField::clear(){
    PropagateField::clear();
    setDefaults();
};


/*!
 * Check coherence of the input data of the class, in particular:
 * - check if boundary patches are referred to target bulk mesh 
 * - check if boundary conditions vectors are coherent with relative boundary patches.
 * In case of successfull check, the internal member m_ibsp is filled out.
 * \return true if coherence is satisfied, false otherwise.
 */
bool PropagateScalarField::checkBoundariesCoherence(){

    //initialize m_isbp
    bitpit::PiercedVector<bitpit::Vertex> & pVtarget = m_geometry->getVertices();

    for(const auto & vert: pVtarget){
        m_isbp.insert(vert.getId(), std::make_pair(false, 0));
    }

    //1st step: verify boundary IDs of Dirichlet boundary patch and target are coherent
    // and fill m_isbp with flag true and mark 1 for Dirichlet condition.
    long id;
    for(const auto & vert: m_bsurface->getVertices()){
        id= vert.getId();
        if(!pVtarget.exists(id)) {
            m_isbp.clear();
            return false;
        }
        m_isbp[id].first = true;
        m_isbp[id].second = 1;
    }

    //2nd step verify coherence of the Dirichlet field with boundary surface
    if(m_bc_dir.getGeometry() != m_bsurface || !m_bc_dir.completeMissingData(0.0)){
        m_isbp.clear();
        return false;
    }

    return true;
}


/*!
 * It applies a smoothing filter for a defined number of step.
 * \param[in] nstep desired number of smoothing steps
 */
void
PropagateScalarField::solveSmoothing(int nstep){


    int nsize;
    livector1D ids;
    long ID;

    {
        m_field.clear();
        double maxval = 0.0;
        for (auto vertex : getGeometry()->getVertices()){
            ID = vertex.getId();
            if (m_isbp[ID].first){
                m_field.insert(ID, m_bc_dir[ID]);
                maxval = std::max(maxval, std::abs(m_bc_dir[ID]));
            }
            else{
                m_field.insert(ID, 0.0);
            }
        }

        (*m_log)<< m_name <<" starts field propagation."<<std::endl;
        for (int istep = 0; istep < nstep; istep++){

            if (!m_convergence) (*m_log)<<"Smoothing step : " << istep+1 << " / " << nstep <<std::endl;

            double maxdiff;
            maxdiff = 0.0;
            for (auto vertex : getGeometry()->getVertices()){

                ID = vertex.getId();
                double value = m_field[ID];

                if (m_isbp[ID].first)    continue;

                m_field[ID] = 0.0;
                nsize = m_conn[ID].size();
                ids = m_conn[ID];
                for (int j=0; j<nsize; j++){
                    m_field[ID] += m_field[ids[j]]*m_weights[ID][j];
                }

                if (m_convergence){
                    maxdiff = std::max(maxdiff, std::abs(value - m_field[ID])/maxval);
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

        m_field.setDataLocation(MPVLocation::POINT);
        m_field.setGeometry(getGeometry());
    }

    m_conn.clear();
    m_weights.clear();
    //m_bc_dir.clear();
}

/*!
 * It solves the laplacian problem.
 */
void
PropagateScalarField::solveLaplace(){


    m_field.clear();
    //initialization
    for (auto vertex : getGeometry()->getVertices()){
        long int ID = vertex.getId();
        if (m_isbp[ID].first){
            m_field.insert(ID, m_bc_dir[ID]);
        }else{
            m_field.insert(ID, 0.0);
        }
    }

    liimap  dataInv = m_geometry->getMapDataInv();

    // Create the system for solving the pressure
    bool debug = false;
    m_solver = std::unique_ptr<mimmo::SystemSolver>(new mimmo::SystemSolver(debug));

    // Initialize the system
    KSPOptions &solverOptions = m_solver->getKSPOptions();
    solverOptions.nullspace = false;
    solverOptions.rtol      = 1e-12;
    solverOptions.subrtol   = 1e-12;

    {
        localivector2D stencils(m_conn.size());
        localdvector2D weights(m_conn.size());
        localdvector1D rhs(m_conn.size());

        //Create stencils for petsc and prepare RHS
        for (auto vertex : getGeometry()->getVertices()){
            long int ID = vertex.getId();
            int ind = dataInv[ID];
            if (m_isbp[ID].first ){

                stencils[ind] = ivector1D(1, ind);
                weights[ind] = dvector1D(1, 1.0);

                rhs[ind] = m_bc_dir[ID];
            }else{
                for (long IDN : m_conn[ID]){
                    stencils[ind].push_back(dataInv[IDN]);
                }
                stencils[ind].push_back(ind);
                weights[ind] = -1.0*m_weights[ID];
                weights[ind].push_back(1.0);
                rhs[ind] = 0.0;
            }
        }

        m_weights.clear();
        m_conn.clear();
        //WARNING releasing of boundary conditions! Done for memory saving
        //m_bc_dir.clear();


#if ENABLE_MPI==1
        m_solver->initialize(stencils, weights, rhs, ghosts);
#else
        m_solver->initialize(stencils, weights, rhs);
#endif

        // Solve the system
        m_solver->solve();

        // Get the solution
        const double *solution = m_solver->getSolutionRawReadPtr();

        for (auto vertex : getGeometry()->getVertices()){
            long int ID = vertex.getId();
            int ind = dataInv[ID];
            m_field[ID] = solution[ind];
        }

        // Clear the solver
        m_solver->clear();

    }

    m_field.setDataLocation(MPVLocation::POINT);
    m_field.setGeometry(getGeometry());

}


/*!
 * Plot optional results on vtu unstructured grid file
 */
void
PropagateScalarField::plotOptionalResults(){


    if(getGeometry() == NULL || getGeometry()->isEmpty())    return;

    bitpit::VTKUnstructuredGrid& vtk = getGeometry()->getPatch()->getVTK();
    dvector1D data;
    for (auto val : m_field){
        data.push_back(val);
    }
    vtk.addData("field", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, data);

    dvector1D datad(m_dumping.size());
    int count = 0;
    for (auto val : m_dumping){
        datad[count] = val;
        count++;
    }
    vtk.addData("dumping", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, datad);

    vtk.setCounter(getId());
    getGeometry()->getPatch()->write(m_name +"_field");
    vtk.removeData("field");
    vtk.removeData("dumping");
    vtk.unsetCounter();

};





//--------------------------------------
//--------------------------------------
// VECTORFIELD (USUALLY GEOMETRY DISPLACEMENTS)
//--------------------------------------
//--------------------------------------

/*!
 * Constructor
 */
PropagateVectorField::PropagateVectorField():PropagateField(){
    m_name = "mimmo.PropagateVectorField";
    setDefaults();
};

/*!
 * Set most significant parameters to constructor defaults
 */
void PropagateVectorField::setDefaults(){
    PropagateField::setDefaults();
    m_bc_dir.clear();
    m_field.clear();
    m_slipsurface = NULL;
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
PropagateVectorField::PropagateVectorField(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.PropagateVectorField";
    setDefaults();

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.PropagateVectorField"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Destructor;
 */
PropagateVectorField::~PropagateVectorField(){
    clear();
};

/*!
 * Copy constructor
 */
PropagateVectorField::PropagateVectorField(const PropagateVectorField & other):PropagateField(other){
    m_bc_dir    = other.m_bc_dir;
    m_slipsurface = other.m_slipsurface;
    m_field     = other.m_field;
};

/*!
 * Assignment operator of the class
 */
PropagateVectorField & PropagateVectorField::operator=(PropagateVectorField other){
    swap(other);
    return *this;
};

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void PropagateVectorField::swap(PropagateVectorField & x) noexcept {
    m_bc_dir.swap(x.m_bc_dir);
    m_field.swap(x.m_field);
    std::swap(m_slipsurface, x.m_slipsurface);
    PropagateField::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
PropagateVectorField::buildPorts(){
    bool built = true;
    built = (built && createPortIn<dmpvecarr3E, PropagateVectorField>(this, &PropagateVectorField::setDirichletConditions, M_GDISPLS));
    built = (built && createPortIn<MimmoObject *, PropagateVectorField>(this, &PropagateVectorField::setSlipBoundarySurface, M_GEOM4));
    built = (built && createPortOut<dmpvecarr3E, PropagateVectorField>(this, &PropagateVectorField::getPropagatedField, M_GDISPLS));
    PropagateField::buildPorts();
    m_arePortsBuilt = built;
};

/*!
 * It gets the resulting deformation field on points cloud.
 * \return Deformation field.
 */
dmpvecarr3E
PropagateVectorField::getPropagatedField(){
    return(m_field);
}

/*! 
 * Sets the portion of boundary mesh relative to geometry target
 * that must be constrained component-wise with zero normal field throughout boundary surface.
 * This patch is optional. If nothing is linked, the relative boundary is 
 * solved free of any conditions.
 * \param[in] surface Boundary patch.
 */
void
PropagateVectorField::setSlipBoundarySurface(MimmoObject* surface){
    if (surface == NULL)       return;
    if (surface->isEmpty())    return;
    if (surface->getType()!= 1 ) return;

    m_slipsurface = surface;

}

/*!
 * It sets the Dirichlet conditions for each component of the vector field on the previously linked
 * Dirichlet Boundary patch.
 * \param[in] bc dirichlet conditions
 */
void
PropagateVectorField::setDirichletConditions(dmpvecarr3E bc){
    if (bc.isEmpty()) return;
    m_bc_dir = bc;
}


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateVectorField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
    BITPIT_UNUSED(name);
    //start absorbing
    PropagateField::absorbSectionXML(slotXML, name);
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateVectorField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    PropagateField::flushSectionXML(slotXML, name);
};


/*!
 * Clear all data actually stored in the class
 */
void
PropagateVectorField::clear(){
    PropagateField::clear();
    setDefaults();
};

/*!
 * Check coherence of the input data of the class, in particular:
 * - check if boundary patches are referred to target bulk mesh 
 * - check if boundary conditions vectors are coherent with relative boundary patches.
 * In case of successfull check, the internal member m_ibsp is filled out.
 * \return true if coherence is satisfied, false otherwise.
 */
bool PropagateVectorField::checkBoundariesCoherence(){

    //initialize m_isbp
    bitpit::PiercedVector<bitpit::Vertex> & pVtarget = m_geometry->getVertices();

    for(const auto & vert: pVtarget){
        m_isbp.insert(vert.getId(), std::make_pair(false, 0));
    }

    //1st step: verify boundary IDs of Dirichlet boundary patch and target are coherent
    // and fill m_isbp with flag true and mark 1 for Dirichlet condition.
    long id;
    for(const auto & vert: m_bsurface->getVertices()){
        id= vert.getId();
        if(!pVtarget.exists(id)) {
            m_isbp.clear();
            return false;
        }
        m_isbp[id].first = true;
        m_isbp[id].second = 1;
    }

    //2nd step verify coherence of the Dirichlet field with boundary surface
    if(m_bc_dir.getGeometry() != m_bsurface || !m_bc_dir.completeMissingData({{0.0, 0.0,0.0}})){
        m_isbp.clear();
        return false;
    }

    // verify if optional Neumann condition are set. If not, exit with true condition
    if(m_slipsurface == NULL)    return true; //ok

    //neumann surface is not null, double step check as m_bsurface
    //3rd step: verify boundary IDs of Neumann boundary patch and target are coherent
    // and fill m_isbp with flag true and mark 2 for Neumann condition.
    for(const auto & vert: m_slipsurface->getVertices()){
        id= vert.getId();
        if(!pVtarget.exists(id)) {
            m_isbp.clear();
            return false;
        }
        if(!m_isbp.exists(id)){
            m_isbp[id].first = true;
            m_isbp[id].second = 2;
        }
    }


    m_vNormals.clear();
    m_vNormals.reserve(m_slipsurface->getNVertex());

    bitpit::ConstProxyVector<long> verts;
    std::size_t size;
    long idN;
    //save the vertex Normals using the boundary surface m_neu_surface;
    for(const auto & cell: m_slipsurface->getCells()){
        verts= cell.getVertexIds();
        size = verts.size();
        for(std::size_t i=0; i<size; ++i){
            idN = verts[i];
            if(m_isbp[idN].second ==2 && !m_vNormals.exists(idN)){
                m_vNormals.insert(idN, static_cast<bitpit::SurfaceKernel*>(m_slipsurface->getPatch())->evalVertexNormal(cell.getId(), i));
            }
        }
    }

    //if it is survived, then it's all ok.
    return true;
}


/*!
 * It applies a smoothing filter for a defined number of step.
 * \param[in] nstep desired number of smoothing steps
 */
void
PropagateVectorField::solveSmoothing(int nstep){

    int nsize;
    livector1D ids;
    livector1D noids;
    long ID;

    {

        m_field.clear();
        double maxval = 0.0;
        for (const auto & vertex : getGeometry()->getVertices()){
            ID = vertex.getId();
            if (m_isbp[ID].first && m_isbp[ID].second ==1 ){
                m_field.insert(ID, m_bc_dir[ID]);
                maxval = std::max(maxval, norm2(m_bc_dir[ID]));
            }
            else{
                m_field.insert(ID, darray3E({0.0, 0.0, 0.0}));
            }
        }

        (*m_log)<< m_name <<" starts field propagation."<<std::endl;
        for (int istep = 0; istep < nstep; istep++){

            if (!m_convergence) (*m_log)<<m_name << " smoothing step : " << istep+1 << " / " << nstep <<std::endl;

            double maxdiff;
            maxdiff = 0.0;
            for (const auto & vertex : getGeometry()->getVertices()){

                ID = vertex.getId();
                darray3E value = m_field[ID];

                if (m_isbp[ID].first && m_isbp[ID].second == 1) continue;

                m_field[ID] = {{0.0,0.0,0.0}};
                nsize = m_conn[ID].size();
                ids = m_conn[ID];
                for (int j=0; j<nsize; j++){
                    m_field[ID] += m_field[ids[j]]*m_weights[ID][j];
                }

                if(m_isbp[ID].first && m_isbp[ID].second == 2){
                    int candidate = 0;
                    if(std::abs(m_vNormals[ID][1]) > std::abs(m_vNormals[ID][candidate])) candidate = 1;
                    if(std::abs(m_vNormals[ID][2]) > std::abs(m_vNormals[ID][candidate])) candidate = 2;

                    m_field[ID][candidate] = -1.0*(m_vNormals[ID][(candidate+1)%3] * m_field[ID][(candidate+1)%3] +
                            m_vNormals[ID][(candidate+2)%3] * m_field[ID][(candidate+2)%3]) / m_vNormals[ID][candidate];

                }
                if (m_convergence){
                    maxdiff = std::max(maxdiff, norm2(value - m_field[ID])/maxval);
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

        m_field.setDataLocation(MPVLocation::POINT);
        m_field.setGeometry(getGeometry());

    }

    m_conn.clear();
    m_weights.clear();
    m_vNormals.clear();
}

/*!
 * It solves the laplacian problem.
 */
void
PropagateVectorField::solveLaplace(){

    m_field.clear();
    for (auto vertex : getGeometry()->getVertices()){
        long int ID = vertex.getId();
        if (m_isbp[ID].first && m_isbp[ID].second ==1){
            m_field.insert(ID, m_bc_dir[ID]);
        }else{
            m_field.insert(ID, darray3E({0.0, 0.0, 0.0}));
        }
    }

    liimap  dataInv = m_geometry->getMapDataInv();

    // Create the system for solving the pressure
    bool debug = false;
    m_solver = std::unique_ptr<mimmo::SystemSolver>(new mimmo::SystemSolver(debug));

    // Initialize the system
    KSPOptions &solverOptions = m_solver->getKSPOptions();
    solverOptions.nullspace = false;
    solverOptions.rtol      = 1e-12;
    solverOptions.subrtol   = 1e-12;

    {
        int connSize = m_conn.size();
        localivector2D stencils(3*connSize);
        localdvector2D weights(3*connSize);
        localdvector1D rhs(3*connSize,0.0);

        //Create stencils for petsc
        for (auto vertex : getGeometry()->getVertices()){
            long int ID = vertex.getId();
            int ind = dataInv[ID];
            if (m_isbp[ID].first && m_isbp[ID].second ==1){
                for(int comp=0; comp<3; ++comp){
                    stencils[ind+comp*connSize] = ivector1D(1, ind+comp*connSize);
                    weights[ind+comp*connSize] = dvector1D(1, 1.0);
                    rhs[ind+comp*connSize] = m_bc_dir[ID][comp]; 
                }
            }else{
                for(int comp=0; comp<3; ++comp){
                    for (const long & IDN : m_conn[ID]){
                        stencils[ind+comp*connSize].push_back(dataInv[IDN]+comp*connSize);
                    }
                    stencils[ind+comp*connSize].push_back(ind+comp*connSize);
                    weights[ind+comp*connSize] = -1.0*m_weights[ID];
                    weights[ind+comp*connSize].push_back(1.0);
                }
            }

            if(m_isbp[ID].first && m_isbp[ID].second == 2){
                int comp = 0;
                if(std::abs(m_vNormals[ID][1]) > std::abs(m_vNormals[ID][comp])) comp = 1;
                if(std::abs(m_vNormals[ID][2]) > std::abs(m_vNormals[ID][comp])) comp = 2;

                stencils[ind+comp*connSize].resize(3); 
                stencils[ind+comp*connSize] = {{ind, ind+connSize, ind+2*connSize}};
                weights[ind+comp*connSize].resize(3);
                for(int i=0; i<3; ++i)
                    weights[ind +comp*connSize][i] = m_vNormals[ID][i]/m_vNormals[ID][comp];
            }
        }
        m_weights.clear();
        m_conn.clear();
        m_vNormals.clear();

#if ENABLE_MPI==1
        m_solver->initialize(stencils, weights, rhs, ghosts);
#else
        m_solver->initialize(stencils, weights, rhs);
#endif

        // Solve the system
        m_solver->solve();

        // Get the solution
        const double *solution = m_solver->getSolutionRawReadPtr();

        for (auto vertex : getGeometry()->getVertices()){
            for (int icomp=0; icomp<3; ++icomp ){ 
                long int ID = vertex.getId();
                int ind = dataInv[ID];
                m_field[ID][icomp] = solution[ind + icomp*connSize];
            }
        }

        // Clear the solver
        m_solver->clear();
    }

    m_field.setDataLocation(MPVLocation::POINT);
    m_field.setGeometry(getGeometry());

}

/*!
 * Plot optional results on vtu unstructured grid file
 */
void
PropagateVectorField::plotOptionalResults(){


    if(getGeometry() == NULL || getGeometry()->isEmpty())    return;

    bitpit::VTKUnstructuredGrid& vtk = getGeometry()->getPatch()->getVTK();
    dvecarr3E data(m_field.size());
    int count = 0;
    for (auto val : m_field){
        data[count] = val;
        count++;
    }
    vtk.addData("field", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, data);

    dvector1D datad(m_dumping.size());
    count = 0;
    for (auto val : m_dumping){
        datad[count] = val;
        count++;
    }
    vtk.addData("dumping", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, datad);

    vtk.setCounter(getId());
    getGeometry()->getPatch()->write(m_name +"_field");
    vtk.removeData("field");
    vtk.removeData("dumping");
    vtk.unsetCounter();

};

/*!
 * Directly apply deformation field to target geometry.
 */
void
PropagateVectorField::apply(){
    if (getGeometry() == NULL) return;
    if (getGeometry()->isEmpty() || m_field.isEmpty()) return;
    darray3E vertexcoords;
    long int ID;
    for (const auto & vertex : m_geometry->getVertices()){
        vertexcoords = vertex.getCoords();
        ID = vertex.getId();
        vertexcoords += m_field[ID];
        getGeometry()->modifyVertex(vertexcoords, ID);
    }
}





}
