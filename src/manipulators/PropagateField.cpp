/*----------------------------------------------------------------------------*\
 *
 *  mimic
 *
 *  Optimad Engineering S.r.l. ("COMPANY") CONFIDENTIAL
 *  Copyright (c) 2015-2017 Optimad Engineering S.r.l., All Rights Reserved.
 *
 *  --------------------------------------------------------------------------
 *
 *  NOTICE:  All information contained herein is, and remains the property
 *  of COMPANY. The intellectual and technical concepts contained herein are
 *  proprietary to COMPANY and may be covered by Italian and Foreign Patents,
 *  patents in process, and are protected by trade secret or copyright law.
 *  Dissemination of this information or reproduction of this material is
 *  strictly forbidden unless prior written permission is obtained from
 *  COMPANY. Access to the source code contained herein is hereby forbidden
 *  to anyone except current COMPANY employees, managers or contractors who
 *  have executed Confidentiality and Non-disclosure agreements explicitly
 *  covering such access.
 *
 *  The copyright notice above does not evidence any actual or intended
 *  publication or disclosure of this source code, which includes information
 *  that is confidential and/or proprietary, and is a trade secret, of
 *  COMPANY. ANY REPRODUCTION, MODIFICATION, DISTRIBUTION, PUBLIC PERFORMANCE,
 *  OR PUBLIC DISPLAY OF OR THROUGH USE  OF THIS  SOURCE CODE  WITHOUT THE
 *  EXPRESS WRITTEN CONSENT OF COMPANY IS STRICTLY PROHIBITED, AND IN
 *  VIOLATION OF APPLICABLE LAWS AND INTERNATIONAL TREATIES. THE RECEIPT OR
 *  POSSESSION OF THIS SOURCE CODE AND/OR RELATED INFORMATION DOES NOT CONVEY
 *  OR IMPLY ANY RIGHTS TO REPRODUCE, DISCLOSE OR DISTRIBUTE ITS CONTENTS, OR
 *  TO MANUFACTURE, USE, OR SELL ANYTHING THAT IT  MAY DESCRIBE, IN WHOLE OR
 *  IN PART.
 *
\*----------------------------------------------------------------------------*/

#include "PropagateField.hpp"
#include "Operators.hpp"

/*!
 * Custom Operator - for bitpit::Vertex
 */
std::array<double, 3> operator-(const bitpit::Vertex &v1, const bitpit::Vertex &v2){
    std::array<double, 3> coords2 = v1.getCoords() - v2.getCoords();
    return ( coords2 );
};


namespace mimmo{

/*!
 * Constructor
 */
PropagateField::PropagateField(){
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
    m_bsurface  = NULL;
    m_laplace   = false;
    m_sstep     = 10;
    m_bPointsToCompute = true;
    m_dumping.clear();
    m_dumpingFactor = 0.0;
    m_convergence = false;
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
    m_isbp        = other.m_isbp;
    m_np        = other.m_np;
    m_nbp        = other.m_nbp;
    m_conn        = other.m_conn;
    m_gamma        = other.m_gamma;
    m_weights    = other.m_weights;
    m_bsurface  = other.m_bsurface;
    m_laplace   = other.m_laplace;
    m_sstep     = other.m_sstep;
    m_dumping     = other.m_dumping;
    m_dumpingFactor  = other.m_dumpingFactor;
    m_convergence     = other.m_convergence;
};

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void PropagateField::swap(PropagateField & x) noexcept {
    std::swap(m_isbp, x.m_isbp);
    std::swap(m_np, x.m_np);
    std::swap(m_nbp, x.m_nbp);
    std::swap(m_conn, x.m_conn);
    std::swap(m_gamma, x.m_gamma);
    std::swap(m_weights, x.m_weights);
    std::swap(m_bsurface, x.m_bsurface);
    std::swap(m_laplace, x.m_laplace);
    std::swap(m_sstep, x.m_sstep);
    std::swap(m_dumping, x.m_dumping);
    std::swap(m_dumpingFactor, x.m_dumpingFactor);
    BaseManipulation::swap(x);
}

/*! 
 * It builds the input/output ports of the object
 */
void
PropagateField::buildPorts(){
    bool built = true;
    built = (built && createPortIn<MimmoObject*, PropagateField>(this, &PropagateField::setGeometry, M_GEOM, true));
    built = (built && createPortIn<MimmoObject*, PropagateField>(this, &PropagateField::setBoundarySurface, M_GEOM2));
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
 * Set pointer to your target geometry. Reimplemented from mimmo::BaseManipulation::setGeometry() 
 * \param[in] geometry_ pointer to target geometry
 */
void
PropagateField::setGeometry(MimmoObject * geometry_){

    if (geometry_ == NULL) return;
    int type = geometry_->getType();
    if ( type != 1 && type != 2 ) return;

    m_geometry = geometry_;
    m_np = getGeometry()->getNVertex();
    m_isbp.clear();
    for (auto vertex : getGeometry()->getVertices()){
        m_isbp.insert(vertex.getId(), false);
    }
    if (m_bsurface != NULL && m_bPointsToCompute == true){
        setBoundarySurface(m_bsurface);
    }
}

/*! 
 * It sets the boundary points of the cloud starting from an input boundary surface.
 * \param[in] bsurface Boundary surface object (mimmo object with coinciding ID for coinciding vertex with volume patch of points cloud).
 */
void
PropagateField::setBoundarySurface(MimmoObject* bsurface){
    if (bsurface == NULL) return;
    int type = bsurface->getType();
    if ( type != 1 && type != 2 && type != 3 && type != 4) return;

    m_bsurface = bsurface;
    m_bPointsToCompute = true;
    if (m_geometry == NULL) return;
    m_nbp = bsurface->getNVertex();
    long ID;
    for (auto vertex : bsurface->getVertices()){
        ID = vertex.getId();
        m_isbp[ID] = true;
    }
    m_bPointsToCompute = false;
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
 * \param[in] solveLaplacian true for Laplacian problem solver (not available); false for iterative smoothing technique.
 */
void PropagateField::setSolver(bool solveLaplacian){
    m_laplace = solveLaplacian;
    //TODO REMOVE IT WHEN LAPLACIAN SOLVER IMPLEMENTED !!!
    //Forced to smoothing
    m_laplace = false;
}

/*!
 * Set the dumping factor.
 * \param[in] dump Exponential of dumping function.
 */
void
PropagateField::setDumpingFactor(double dump){
    m_dumpingFactor = dump;
}

/*!
 * Set the dumping radius.
 * \param[in] radius Support radius of dumping function.
 */
void
PropagateField::setDumpingRadius(double radius){
    m_radius = radius;
}

/*! 
 * It sets if the solver has to reach the convergence.
 * \param[in] convergence Convergence flag.
 */
void PropagateField::setConvergence(bool convergence){
    m_convergence = convergence;
}

/*!
 * Clear all data actually stored in the class
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

    //Old algorithm to review....
    if (m_laplace)
    {
        //TODO insert here your laplacian solver for connected meshes
    }
    else{
        smoothing();
    }

    //    m_conn.clear();
    m_weights.clear();

}

/*!
 * Smoothing execution command
 */
void
PropagateField::smoothing(){
    if (m_geometry == NULL) return;
    computeConnectivitySmoothing();
    computeDumpingFunction();
    computeWeightsSmoothing();
    solveSmoothing(m_sstep);
}

/*! 
 * It computes the connectivity structure between input points.
 */
void
PropagateField::computeConnectivitySmoothing(){

    //DOESN'T WORK FOR POINTS CLOUD OR CURVES

    if (m_geometry == NULL) return;

    bitpit::PatchKernel * patch_ = getGeometry()->getPatch();

    int type = getGeometry()->getType();

    double h;
    int nsize;
    livector1D ids;
    livector1D noids;
    long ID;
    double dist;

    m_conn.resize(m_np);
    m_weights.resize(m_np);

    //map for edges already visited
    std::map<std::pair<long, long>, bool> visitedge;

    //compute connectivity between points
    for (auto cell : patch_->getCells()){
        int edgecount;
        if (type == 1 || patch_->getDimension() == 2){
            edgecount = cell.getFaceCount();
        }
        else{
            edgecount = cell.getEdgeCount();
        }
        for (int ie=0; ie<edgecount; ie++){
            livector1D econn;
            if (type == 1 || patch_->getDimension() == 2){
                econn = cell.getFaceConnect(ie);
            }
            else{
                econn = cell.getEdgeConnect(ie);
            }
            if (!visitedge[std::pair<long, long>(econn[0],econn[1])]){
                visitedge[std::pair<long, long>(econn[0],econn[1])] = true;
                visitedge[std::pair<long, long>(econn[1],econn[0])] = true;
                if (!m_conn.exists(econn[0])){
                    m_conn.insert(econn[0], livector1D(1, econn[1]));
                }
                else{
                    m_conn[econn[0]].push_back(econn[1]);
                }
                if (!m_conn.exists(econn[1])){
                    m_conn.insert(econn[1], livector1D(1, econn[0]));
                }
                else{
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
PropagateField::computeWeightsSmoothing(){

    if (m_geometry == NULL) return;

    bitpit::PatchKernel * patch_ = getGeometry()->getPatch();
    if (m_conn.size() != patch_->getVertexCount()) return;

    int nsize;
    livector1D ids;
    double dist;
    double sumdist;
    long ID;

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
            //            m_weights[m_mapID[ID]][j] = 1/((std::abs(m_gamma))*dist);
            //            sumdist += 1/((std::abs(m_gamma))*dist);
            //            //m_weights[m_mapID[ID]][j] = 1/dist;
            lweights[j] = 1/(std::pow(dist, m_gamma))*m_dumping[ids[j]];
            sumdist += lweights[j];
        }
        lweights /= sumdist;
        m_weights.insert(ID, lweights);
    }
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    /*
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

    if(slotXML.hasOption("minSearchRadius")){
        std::string input = slotXML.get("minSearchRadius");
        input = bitpit::utils::string::trim(input);
        double value = 0.0;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
            value = std::fmax(0.0, value);
        }
        setRmin(value);
    }

    if(slotXML.hasOption("maxSearchRadius")){
        std::string input = slotXML.get("maxSearchRadius");
        input = bitpit::utils::string::trim(input);
        double value = 1.0;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
            value = std::fmax(0.0, value);
        }
        setRmax(value);
    }

    if(slotXML.hasOption("SearchSteps")){
        std::string input = slotXML.get("SearchSteps");
        int value =10000;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss>>value;
        }
        setNR(value);
    };

    if(slotXML.hasOption("MinimumNeighs")){
        std::string input = slotXML.get("MinimumNeighs");
        int value =6;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss>>value;
        }
        setMinimumNeighbors(value);
    };

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

    if(slotXML.hasOption("ComputeKdTree")){
        std::string input = slotXML.get("ComputeKdTree");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setComputeKDTree(value);
    }

    if(slotXML.hasOption("ComputeConnectivity")){
        std::string input = slotXML.get("ComputeConnectivity");
        input = bitpit::utils::string::trim(input);
        bool value = true;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setComputeConnectivity(value);
    }

    if(slotXML.hasOption("RBFSupportRadius")){
        std::string input = slotXML.get("RBFSupportRadius");
        input = bitpit::utils::string::trim(input);
        double value = 1.0e+18;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
            value = std::fmax(0.0, value);
        }
        setSupportRadius(value);
    }
     */
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    /*
    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("Solver", std::to_string(int(m_laplace)));

    slotXML.set("minSearchRadius",std::to_string(m_rmin));
    slotXML.set("maxSearchRadius",std::to_string(m_rmax));
    slotXML.set("SearchSteps",std::to_string(m_nr));
    slotXML.set("MinimumNeighs",std::to_string(m_minNeigh));
    slotXML.set("WeightConstant",std::to_string(m_gamma));
    slotXML.set("SmoothingSteps",std::to_string(m_sstep));
    slotXML.set("ComputeKdTree",std::to_string(int(m_kdtree)));
    slotXML.set("ComputeConnectivity",std::to_string(int(m_neigh)));
    slotXML.set("RBFSupportRadius",std::to_string(m_suppR));
     */
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
    m_bc.clear();
    m_field.clear();
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
    m_bc        = other.m_bc;
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
    std::swap(m_bc, x.m_bc);
    std::swap(m_field, x.m_field);
    PropagateField::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
PropagateVectorField::buildPorts(){
    bool built = true;
    built = (built && createPortIn<dmpvecarr3E, PropagateVectorField>(this, &PropagateVectorField::setBoundaryConditions, M_GDISPLS));
    built = (built && createPortOut<dmpvecarr3E, PropagateVectorField>(this, &PropagateVectorField::getField, M_GDISPLS));
    PropagateField::buildPorts();
    m_arePortsBuilt = built;
};

/*!
 * It gets the resulting deformation field on points cloud.
 * \return Deformation field.
 */
dmpvecarr3E
PropagateVectorField::getField(){
    return(m_field);
}

/*!
 * It sets the boundary condition on boundary points of the cloud.
 * \param[in] bc Value of Dirchlet boundary conditions on boundary points.
 */
void
PropagateVectorField::setBoundaryConditions(dmpvecarr3E bc){
    m_bc = bc;
}

/*!
 * Clear all data actually stored in the class
 */
void
PropagateVectorField::clear(){
    PropagateField::clear();
    setDefaults();
};

/*!
 * It computes the dumping function used in weights computation.
 */
void
PropagateVectorField::computeDumpingFunction(){

    if (m_geometry == NULL ) return;

    bitpit::PatchKernel * patch_ = getGeometry()->getPatch();
    if (m_conn.size() != patch_->getVertexCount()) return;

    int nsize;
    livector1D ids;
    double dist;
    long ID;

    /*Maxdist should be the maximum distance between boundaries with zero values and
     * boundaries with values different from zero.
     */
    //TODO compute it
    double maxd = m_radius;

    m_dumping.clear();
    darray3E point;
    double val;
    if (m_dumpingFactor <= 1.0e-12){
        for (auto const & vertex : patch_->getVertices()){
            ID = vertex.getId();
            m_dumping.insert(ID, 1.0);
        }
    }
    else{


        bitpit::SurfUnstructured boundDef(999, 1, 2);
        boundDef.setExpert(true);
        for (auto const & vertex : m_bsurface->getVertices()){
            ID = vertex.getId();
            if (norm2(m_bc[ID]) >= 1.0e-12){
                boundDef.addVertex(vertex, ID);
            }
        }
        for (auto const & cell : m_bsurface->getCells()){
            ID = cell.getId();
            auto conn = cell.getConnect();
            std::vector<long> vconn(cell.getVertexCount());
            for (int i=0; i<cell.getVertexCount(); i++)
                vconn[i] = conn[i];
            for (int i=0; i<cell.getVertexCount(); i++){
                if (norm2(m_bc[vconn[i]]) >= 1.0e-12){
                    boundDef.addCell(cell.getType(), cell.isInterior(), vconn);
                    break;
                }
            }
        }
        boundDef.buildAdjacencies();
        bitpit::SurfaceSkdTree tree(&boundDef);
        tree.build();

        for (auto const & vertex : patch_->getVertices()){
            ID = vertex.getId();
            point = vertex.getCoords();
            dist = max(1.0e-08, tree.evalPointDistance(point));
            //            val = 1.0 + m_dumpingFactor*bitpit::rbf::wendlandc2( dist/maxd );
            val = std::max(1.0, std::pow((maxd/dist), m_dumpingFactor));
            m_dumping.insert(ID, val);
        }
    }
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

    bitpit::PatchKernel * patch_ = getGeometry()->getPatch();

    {

        m_field.clear();
        for (auto vertex : getGeometry()->getVertices()){
            ID = vertex.getId();
            if (m_isbp[ID]){
                m_field.insert(ID, m_bc[ID]);
            }
            else{
                m_field.insert(ID, darray3E({0.0, 0.0, 0.0}));
            }
        }

        (*m_log)<< m_name <<" starts field propagation."<<std::endl;
        for (int istep = 0; istep < nstep; istep++){

            (*m_log)<<"Smoothing step : " << istep+1 << " / " << nstep <<std::endl;

            double maxdiff;
            maxdiff = 0.0;
            for (auto vertex : getGeometry()->getVertices()){

                ID = vertex.getId();
                darray3E value = m_field[ID];

                if (!m_isbp[ID]){

                    m_field[ID] = darray3E({0.0, 0.0, 0.0});
                    nsize = m_conn[ID].size();
                    ids = m_conn[ID];
                    for (int j=0; j<nsize; j++){
                        m_field[ID] += m_field[ids[j]]*m_weights[ID][j];
                    }
                }
                if (m_convergence)
                    maxdiff = std::max(maxdiff, norm2(value - m_field[ID]));
            }//end for vertex

            //convergence
            if (m_convergence){
                if (maxdiff <= 1.0e-08)
                    istep = nstep;
                else
                    istep = 0;
            }

        }// end step
        (*m_log)<< m_name<<" ends field propagation."<<std::endl;
    }
}

/*!
 * Plot optional results on vtu unstructured grid file
 */
void
PropagateVectorField::plotOptionalResults(){


    if(getGeometry() == NULL || getGeometry()->isEmpty())    return;

    bitpit::VTKUnstructuredGrid& vtk = getGeometry()->getPatch()->getVTK();
    dvecarr3E data;
    for (auto val : m_field){
        data.push_back(val);
    }
    vtk.addData("field", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, data);

    dvector1D datad(m_dumping.size());
    int count = 0;
    for (auto val : m_dumping){
        datad[count] = val;
        count++;
    }
    vtk.addData("dumping", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, datad);

    vtk.setCounter(getId());
    vtk.setName(m_name +"_field");
    getGeometry()->getPatch()->write();
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


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateVectorField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

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

    if(slotXML.hasOption("WeightFactor")){
        std::string input = slotXML.get("WeightFactor");
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

    if(slotXML.hasOption("DumpingFactor")){
        std::string input = slotXML.get("DumpingFactor");
        input = bitpit::utils::string::trim(input);
        double value = 1.0;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
            value = std::fmax(0.0, value);
        }
        setDumpingFactor(value);
    }

    if(slotXML.hasOption("DumpingRadius")){
        std::string input = slotXML.get("DumpingRadius");
        input = bitpit::utils::string::trim(input);
        double value = 1.0e+18;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
            value = std::fmax(0.0, value);
        }
        setDumpingRadius(value);
    }
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateVectorField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("Solver", std::to_string(int(m_laplace)));
    slotXML.set("WeightFactor",std::to_string(m_gamma));
    slotXML.set("SmoothingSteps",std::to_string(m_sstep));
    slotXML.set("Convergence",std::to_string(int(m_convergence)));
    slotXML.set("DumpingFactor",std::to_string(m_dumpingFactor));
    slotXML.set("DumpingRadius",std::to_string(m_radius));
};



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
    m_bc.clear();
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
    m_bc        = other.m_bc;
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
    std::swap(m_bc, x.m_bc);
    std::swap(m_field, x.m_field);
    PropagateField::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
PropagateScalarField::buildPorts(){
    bool built = true;
    built = (built && createPortIn<dmpvector1D, PropagateScalarField>(this, &PropagateScalarField::setBoundaryConditions, M_FILTER));
    built = (built && createPortOut<dmpvector1D, PropagateScalarField>(this, &PropagateScalarField::getField, M_FILTER));
    PropagateField::buildPorts();
    m_arePortsBuilt = built;
};

/*!
 * It gets the resulting deformation field on points cloud.
 * \return Deformation field.
 */
dmpvector1D
PropagateScalarField::getField(){
    return(m_field);
}

/*!
 * It sets the boundary condition on boundary points of the cloud.
 * \param[in] bc Value of Dirchlet boundary conditions on boundary points.
 */
void
PropagateScalarField::setBoundaryConditions(dmpvector1D bc){
    m_bc = bc;
}

/*!
 * Clear all data actually stored in the class
 */
void
PropagateScalarField::clear(){
    PropagateField::clear();
    setDefaults();
};

/*!
 * It computes the dumping function used in weights computation.
 */
void
PropagateScalarField::computeDumpingFunction(){

    if (m_geometry == NULL ) return;

    bitpit::PatchKernel * patch_ = getGeometry()->getPatch();
    if (m_conn.size() != patch_->getVertexCount()) return;

    int nsize;
    livector1D ids;
    double dist;
    long ID;

    /*Maxdist should be the maximum distance between boundaries with zero values and
     * boundaries with values different from zero.
     */
    //TODO compute it
    double maxd = m_radius;

    m_dumping.clear();
    darray3E point;
    double val;
    if (m_dumpingFactor <= 1.0e-12){
        for (auto const & vertex : patch_->getVertices()){
            ID = vertex.getId();
            m_dumping.insert(ID, 1.0);
        }
    }
    else{


        bitpit::SurfUnstructured boundDef(999, 1, 2);
        boundDef.setExpert(true);
        for (auto const & vertex : m_bsurface->getVertices()){
            ID = vertex.getId();
            if (std::abs(m_bc[ID]) >= 1.0e-12){
                boundDef.addVertex(vertex, ID);
            }
        }
        for (auto const & cell : m_bsurface->getCells()){
            ID = cell.getId();
            auto conn = cell.getConnect();
            std::vector<long> vconn(cell.getVertexCount());
            for (int i=0; i<cell.getVertexCount(); i++)
                vconn[i] = conn[i];
            for (int i=0; i<cell.getVertexCount(); i++){
                if (std::abs(m_bc[vconn[i]]) >= 1.0e-12){
                    boundDef.addCell(cell.getType(), cell.isInterior(), vconn);
                    break;
                }
            }
        }
        boundDef.buildAdjacencies();
        bitpit::SurfaceSkdTree tree(&boundDef);
        tree.build();

        for (auto const & vertex : patch_->getVertices()){
            ID = vertex.getId();
            point = vertex.getCoords();
            dist = max(1.0e-08, tree.evalPointDistance(point));
            val = std::max(1.0, std::pow((maxd/dist), m_dumpingFactor));
            m_dumping.insert(ID, val);
        }
    }
}

/*!
 * It applies a smoothing filter for a defined number of step.
 * \param[in] nstep desired number of smoothing steps
 */
void
PropagateScalarField::solveSmoothing(int nstep){


    int nsize;
    livector1D ids;
    livector1D noids;
    long ID;

    bitpit::PatchKernel * patch_ = getGeometry()->getPatch();

    {

        m_field.clear();
        for (auto vertex : getGeometry()->getVertices()){
            ID = vertex.getId();
            if (m_isbp[ID]){
                m_field.insert(ID, m_bc[ID]);
            }
            else{
                m_field.insert(ID, 0.0);
            }
        }

        (*m_log)<< m_name <<" starts field propagation."<< std::endl;
        nstep = std::max(1, nstep);
        for (int istep = 0; istep < nstep; istep++){

            (*m_log)<<"Smoothing step : " << istep+1 << " / " << nstep <<std::endl;

            double maxdiff;
            maxdiff = 0.0;
            for (auto vertex : getGeometry()->getVertices()){

                ID = vertex.getId();
                double value = m_field[ID];

                if (!m_isbp[ID]){

                    m_field[ID] = 0.0;
                    nsize = m_conn[ID].size();
                    ids = m_conn[ID];
                    for (int j=0; j<nsize; j++){
                        m_field[ID] += m_field[ids[j]]*m_weights[ID][j];
                    }
                }
                if (m_convergence)
                    maxdiff = std::max(maxdiff, std::abs(value - m_field[ID]));
            }//end for vertex

            //convergence
            if (m_convergence){
                if (maxdiff <= 1.0e-08)
                    istep = nstep;
                else
                    istep = 0;

                std::cout << maxdiff << std::endl;
            }

        }// end step
        (*m_log)<< m_name<<" ends field propagation."<<std::endl;
    }
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
    vtk.setName(m_name +"_field");
    getGeometry()->getPatch()->write();
    vtk.removeData("field");
    vtk.removeData("dumping");
    vtk.unsetCounter();

};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateScalarField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){


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

    if(slotXML.hasOption("WeightFactor")){
        std::string input = slotXML.get("WeightFactor");
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

    if(slotXML.hasOption("DumpingFactor")){
        std::string input = slotXML.get("DumpingFactor");
        input = bitpit::utils::string::trim(input);
        double value = 1.0;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
            value = std::fmax(0.0, value);
        }
        setDumpingFactor(value);
    }

    if(slotXML.hasOption("DumpingRadius")){
        std::string input = slotXML.get("DumpingRadius");
        input = bitpit::utils::string::trim(input);
        double value = 1.0e+18;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
            value = std::fmax(0.0, value);
        }
        setDumpingRadius(value);
    }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateScalarField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("Solver", std::to_string(int(m_laplace)));
    slotXML.set("WeightFactor",std::to_string(m_gamma));
    slotXML.set("SmoothingSteps",std::to_string(m_sstep));
    slotXML.set("Convergence",std::to_string(int(m_convergence)));
    slotXML.set("DumpingFactor",std::to_string(m_dumpingFactor));
    slotXML.set("DumpingRadius",std::to_string(m_radius));
};


}
