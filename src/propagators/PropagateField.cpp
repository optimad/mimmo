
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
#include <CG.hpp>
namespace mimmo {
//--------------------------------------
//--------------------------------------
// SCALARFIELD (USUALLY FILTER DISPLACEMENTS)
//--------------------------------------
//--------------------------------------

/*!
 * Constructor
 */
PropagateScalarField::PropagateScalarField():PropagateField<1>(){
    m_name = "mimmo.PropagateScalarField";
    setDefaults();
};

/*!
 * Set most significant parameters to constructor defaults
 */
void PropagateScalarField::setDefaults(){
    PropagateField<1>::setDefaults();
    m_thres = -1.0;
    m_nstep = 1;
}
/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
PropagateScalarField::PropagateScalarField(const bitpit::Config::Section & rootXML):PropagateField<1>(){

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
PropagateScalarField::~PropagateScalarField(){};

/*!
 * Copy constructor
 */
PropagateScalarField::PropagateScalarField(const PropagateScalarField & other):PropagateField<1>(other){
    m_nstep = other.m_nstep;
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
    std::swap(m_nstep, x.m_nstep);
    PropagateField<1>::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
PropagateScalarField::buildPorts(){

    PropagateField<1>::buildPorts();

    bool built = m_arePortsBuilt;
    built = (built && createPortIn<dmpvector1D*, PropagateScalarField>(this, &PropagateScalarField::addDirichletConditions, M_FILTER, true));
    built = (built && createPortOut<dmpvector1D*, PropagateScalarField>(this, &PropagateScalarField::getPropagatedField, M_FILTER));

    m_arePortsBuilt = built;
};

/*!
 * It gets the resulting propagated field on the whole bulk mesh.
 * \return Deformation field.
 */
dmpvector1D*
PropagateScalarField::getPropagatedField(){
    m_tempfield.clear();
    m_tempfield.reserve(m_field.size());
    m_tempfield.setDataLocation(m_field.getDataLocation());
    m_tempfield.setGeometry(m_field.getGeometry());
    for(auto it = m_field.begin(); it != m_field.end(); ++it){
        m_tempfield.insert(it.getId(), (*it)[0]);
    }

    return &m_tempfield;
}

/*!
  Add Dirichlet conditions for scalar field on the previously linked
  Dirichlet Boundary patches list.
  \param[in] bc dirichlet conditions
 */
void
PropagateScalarField::addDirichletConditions(dmpvector1D * bc){
    //avoid linking null field or field with null geometry inside.
    if (!bc) return;
    if(!bc->getGeometry()) return;
    if(bc->getGeometry()->getType() != 1) return;

    //store it in temporary structure for dirichlet
    m_tempDirichletBcs.push_back(MimmoPiercedVector<std::array<double,1>> (bc->getGeometry(), bc->getDataLocation()));
    m_tempDirichletBcs.back().reserve(bc->size());
    for(auto it = bc->begin(); it != bc->end(); ++it){
        m_tempDirichletBcs.back().insert(it.getId(), std::array<double,1>({*it}));
    }
    //insert the pointer of the temp structure in the official list of bcs.
    m_dirichletBcs.insert(&m_tempDirichletBcs.back());
    //save also bc geometry in the list of Dirichlet surfaces.
    m_dirichletSurfaces.insert(m_tempDirichletBcs.back().getGeometry());
}

/*!
 * Force solver to get deformation in a finite number of substep
 * \param[in] sstep number of substep. Default is 1.
 */
void
PropagateScalarField::setSolverMultiStep(unsigned int sstep){
    unsigned int loc(1);
    m_nstep = std::max(loc,sstep);
}


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateScalarField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
    BITPIT_UNUSED(name);
    //start absorbing
    PropagateField<1>::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("MultiStep")){
        std::string input = slotXML.get("MultiStep");
        input = bitpit::utils::string::trim(input);
        double value = 1.0;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        unsigned int value2 = 1;
        if(value >= 1.0) value2 = (unsigned int)value;
        setSolverMultiStep(value2);
    }
}

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateScalarField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    PropagateField<1>::flushSectionXML(slotXML, name);
    slotXML.set("MultiStep",std::to_string(m_nstep));
};

/*!
 * Clear all data actually stored in the class
 */
void
PropagateScalarField::clear(){
    PropagateField<1>::clear();
    m_tempfield.clear();
    m_tempDirichletBcs.clear();
    setDefaults();
};

/*!
 * Plot optional results on vtu unstructured grid file
 */
void
PropagateScalarField::plotOptionalResults(){

    if(getGeometry() == nullptr)    return;
    //Recover the scalar field
    dmpvector1D * field = getPropagatedField();
    //Recover the narrowBandDistances.
    MimmoPiercedVector<double> nbc(getGeometry(), MPVLocation::POINT);
    nbc = m_banddistances; //recover data of the raw PiercedVector;
    nbc.completeMissingData(1.E+18);
    //Set names
    field->setName("field");
    m_damping.setName("damping");
    nbc.setName("narrowband");
    //Write
    write(getGeometry(), *field, m_damping, nbc);

};

/*!
 * Execution command. After the execution the result constraint field is stored in the class.
 */
void
PropagateScalarField::execute(){

    MimmoSharedPointer<MimmoObject> geo = getGeometry();
    if(!geo){
        (*m_log)<<"Error in "<<m_name<<" .No target volume mesh linked"<<std::endl;
        throw std::runtime_error("Error in "+m_name+" .No target volume mesh linked");
    }

    if(m_dirichletSurfaces.empty()){
        (*m_log)<<"Warning in "<<m_name<<" .No Dirichlet Boundary patch linked"<<std::endl;
    }

#if MIMMO_ENABLE_MPI
    //be sure ghost info are available
    if(!geo->isInfoSync()) geo->buildPatchInfo();
    if(!geo->arePointGhostExchangeInfoSync()) geo->updatePointGhostExchangeInfo();
#endif

    if(!checkBoundariesCoherence()){
        (*m_log)<<"Warning in "<<m_name<<" .Boundary patches linked are uncoherent with target bulk geometry"<<std::endl;
    }

    (*m_log) << bitpit::log::priority(bitpit::log::NORMAL);
    (*m_log) << bitpit::log::context("mimmo");

    //allocate the solver;
    m_solver = std::unique_ptr<bitpit::SystemSolver>(new bitpit::SystemSolver(m_print));
    //get this inverse map -> you will need it to compact the stencils.
    liimap dataInv;

    //check if damping or narrow band control are active,
    //initialize their reference surfaces and compute them
    if(m_dampingActive){
        if(m_dampingSurfaces.empty())   m_dampingSurfaces = m_dirichletSurfaces;
        initializeUniqueSurface(m_dampingSurfaces, m_dampingUniSurface);
    }

    if(m_bandActive){
        if(m_bandSurfaces.empty())   m_bandSurfaces = m_dirichletSurfaces;
        initializeUniqueSurface(m_bandSurfaces, m_bandUniSurface);
    }

    //initialize damping and Narrow band Control. If UniSurfaces are null the methods set:
    // - unitary m_damping field.
    // - empty m_banddistances member.
    //
    initializeDampingFunction();
    updateNarrowBand();

    // Graph Laplace method on points

    //store the id of the border nodes only;
    livector1D borderPointsID = geo->extractBoundaryVertexID(false);

    //get this inverse map -> you will need it to compact the stencils.
    dataInv = geo->getMapDataInv(true);

    //pass dirichlet bc point information to bulk m_bc_dir member.
    distributeBCOnBoundaryPoints();

    // compute the laplacian stencils
    GraphLaplStencil::MPVStencilUPtr laplaceStencils = GraphLaplStencil::computeLaplacianStencils(geo, m_tol, &m_damping);

    //modify stencils if Narrow band is active i.e. m_banddistances is not empty.
    //This is directly managed in the method.
    modifyStencilsForNarrowBand(laplaceStencils);

    // initialize the laplacian Matrix in solver and squeeze out the laplace stencils and save border cells only.
    initializeLaplaceSolver(laplaceStencils.get(), dataInv);
    laplaceStencils->squeezeOutExcept(borderPointsID);
    borderPointsID.clear();

    //create step bc in case of multistep.
    MimmoPiercedVector<std::array<double, 1> > stepBCdir(geo, MPVLocation::POINT);
    stepBCdir.reserve(m_bc_dir.size());
    for(auto it = m_bc_dir.begin(); it != m_bc_dir.end(); ++it){
        stepBCdir.insert(it.getId(), *it/double(m_nstep));
    }

    //solve
    std::vector<std::vector<double>> result(1);
    // multistep subiteration. Grid does not change, boundaries are forced each step with a constant increment, so:
    for(int istep=0; istep < m_nstep; ++istep){

        for(auto it = m_bc_dir.begin(); it != m_bc_dir.end(); ++it){
            *it = double(istep + 1) * stepBCdir.at(it.getId());
        }
        //update the solver matrix applying bc and evaluate the rhs;
        dvector1D rhs(geo->getNInternalVertices(), 0.0);
        //USELESS FOR EACH STEP?...
        assignBCAndEvaluateRHS(0, false, laplaceStencils.get(), dataInv, rhs);
        solveLaplace(rhs, result[0]);
        (*m_log)<<m_name<<" solved step "<<istep+1<<" out of total steps "<<m_nstep<<std::endl;
    }

    dataInv.clear();
    //reconstruct getting the direct node map -> you will need it uncompact the system solution in global id.
    liimap mapdata = geo->getMapData(true);
    reconstructResults(result, mapdata);
    // now data are direcly pushed in m_field.

    //clear UniSurface for Damping and NarrowBand if any
    m_bandUniSurface = nullptr;
    m_dampingUniSurface = nullptr;
    //clear temp bc;
    m_bc_dir.clear();
    //clear the solver;
    m_solver->clear();
    (*m_log) << bitpit::log::priority(bitpit::log::DEBUG);
}


//--------------------------------------
//--------------------------------------
// VECTORFIELD (USUALLY GEOMETRY DISPLACEMENTS)
//--------------------------------------
//--------------------------------------

/*!
 * Constructor
 */
PropagateVectorField::PropagateVectorField():PropagateField<3>(){
    m_name = "mimmo.PropagateVectorField";
    setDefaults();
};

/*!
 * Set most significant parameters to constructor defaults
 */
void PropagateVectorField::setDefaults(){
    PropagateField<3>::setDefaults();
    m_nstep = 1;
    m_forcePlanarSlip = false;
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
PropagateVectorField::PropagateVectorField(const bitpit::Config::Section & rootXML):PropagateField<3>(){

    m_name = "mimmo.PropagateVectorField";
    setDefaults();

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.PropagateVectorField"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    }
}

/*!
 * Destructor;
 */
PropagateVectorField::~PropagateVectorField(){};

/*!
 * Copy constructor
 */
PropagateVectorField::PropagateVectorField(const PropagateVectorField & other):PropagateField<3>(other){
    m_nstep = other.m_nstep;
    m_forcePlanarSlip = other.m_forcePlanarSlip;
    m_slipSurfaces = other.m_slipSurfaces;
    m_slipReferenceSurfaces = other.m_slipReferenceSurfaces;
    m_periodicSurfaces = other.m_periodicSurfaces;

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
    std::swap(m_nstep, x.m_nstep);
    std::swap(m_forcePlanarSlip, x.m_forcePlanarSlip);
    std::swap(m_slipSurfaces, x.m_slipSurfaces);
    std::swap(m_slipReferenceSurfaces,x.m_slipReferenceSurfaces);
    std::swap(m_periodicSurfaces, x.m_periodicSurfaces);
    PropagateField<3>::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
PropagateVectorField::buildPorts(){

    PropagateField<3>::buildPorts();
	bool built = m_arePortsBuilt;
	built = (built && createPortIn<dmpvecarr3E *, PropagateVectorField>(this, &PropagateVectorField::addDirichletConditions, M_GDISPLS, true));
	built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, PropagateVectorField>(this, &PropagateVectorField::addSlipBoundarySurface, M_GEOM4));
	built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, PropagateVectorField>(this, &PropagateVectorField::addSlipReferenceSurface, M_GEOM6));
	built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, PropagateVectorField>(this, &PropagateVectorField::addPeriodicBoundarySurface, M_GEOM5));
    built = (built && createPortOut<dmpvecarr3E *, PropagateVectorField>(this, &PropagateVectorField::getPropagatedField, M_GDISPLS));
	m_arePortsBuilt = built;
};

/*!
 * It gets the resulting deformation field on points.
 * \return Deformation field.
 */
dmpvecarr3E*
PropagateVectorField::getPropagatedField(){
	return &m_field;
}

/*!
 * \return true if the class is forcing the slip surface to be treated as plane.
    See forcePlanarSlip method docs.
 */
bool
PropagateVectorField::isForcingPlanarSlip(){
    return m_forcePlanarSlip;
}

/*!
 * Add the portion of boundary mesh to identify zone of the bulk volume target
 * where the field is reprojected onto a reference slip geometry.
 * This patch is optional. If nothing is linked, the relative boundary is
 * solved free of any conditions.
 * WARNING: in case of internal holes inside the slip surface, quality of the slip computation
 * might be affected. Slip surface holes are NOT filled automatically right now.
 * \param[in] surface Boundary patch.
 */
void
PropagateVectorField::addSlipBoundarySurface(MimmoSharedPointer<MimmoObject> surface){
    if (!surface)       return;
    if (surface->getType()!= 1 ) return;
    m_slipSurfaces.insert(surface);
}

/*!
 * Add a portion of surface mesh against which the field on slip boundaries is
   reprojected. The ensemble of linked patches define the final slip reference surfaces.
 * This patch is optional. If nothing is linked, slip boundaries linked for identification
   will be used as slip reference geometries.
 * \param[in] surface Slip reference surface.
 */
void
PropagateVectorField::addSlipReferenceSurface(MimmoSharedPointer<MimmoObject> surface){
    if (!surface)       return;
    if (surface->getType()!= 1 ) return;
    m_slipReferenceSurfaces.insert(surface);
}

/*!
 * Add a portion of boundary mesh relative to geometry target
 * that must be constrained as periodic boundary surface.
   As stated in the general doc of the class, this is a sort of special
   slip condition, which assign to the boundary patch border nodes a zero
   Dirichlet condition.
 * This patch is optional. If nothing is linked, the relative boundary is
 * solved free of any conditions.
 * \param[in] surface Boundary patch.
 */
void
PropagateVectorField::addPeriodicBoundarySurface(MimmoSharedPointer<MimmoObject> surface){
    if (!surface)       return;
    if (surface->getType()!= 1 ) return;
    m_periodicSurfaces.insert(surface);
    //save also surface in the list of slip surfaces.
    m_slipSurfaces.insert(surface);

}


/*!
 * Instead of using the external reference slip surfaces to reproject the final
   field on slip boundaries, this method forces the class to use an implicit plane,
   computed as an average plane of bulk volume boundaries identified as slip.
 * If boundary surfaces are quasi-planar, this trick allow the User to
   overcome the effect of internal holes, without affecting too much quality of the
 * slip conditions computed.
 * Results are not guaranteed in case of highly non-planar slip surfaces.
 * \param[in] planar true, use the average plane of the slip surface to calculate slip conditions, false use the slip refernce surface externally linked.
 */
void
PropagateVectorField::forcePlanarSlip(bool planar){
    m_forcePlanarSlip = planar;
}


/*!
 * Add a Dirichlet condition field for each patch linked as Dirichlet (see addDirichletBoundarySurface).
 * \param[in] bc dirichlet conditions
 */
void
PropagateVectorField::addDirichletConditions(dmpvecarr3E * bc){
    if(!bc) return;
    if(!bc->getGeometry()) return;
    if(bc->getGeometry()->getType() != 1) return;
    //insert the field pointer in the official list of bcs.
    m_dirichletBcs.insert(bc);
    //save also bc geometry in the list of Dirichlet surfaces.
    m_dirichletSurfaces.insert(bc->getGeometry());
}

/*!
 * Force solver to get deformation in a finite number of substep
 * \param[in] sstep number of substep. Default is 1.
 */
void
PropagateVectorField::setSolverMultiStep(unsigned int sstep){
	unsigned int loc(1);
	m_nstep = std::max(loc,sstep);
}

/*!
 * Clear all data actually stored in the class
 */
void
PropagateVectorField::clear(){
    PropagateField<3>::clear();
    m_slipSurfaces.clear();
    m_slipReferenceSurfaces.clear();
    m_slipUniSurface = nullptr;

    m_periodicSurfaces.clear();
    m_periodicBoundaryPoints.clear();
    setDefaults();
};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateVectorField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
    BITPIT_UNUSED(name);
    //start absorbing
    PropagateField<3>::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("MultiStep")){
        std::string input = slotXML.get("MultiStep");
        input = bitpit::utils::string::trim(input);
        double value = 1.0;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        unsigned int value2 = 1;
        if(value >= 1.0) value2 = value;
        setSolverMultiStep(value2);
    }

    if(slotXML.hasOption("ForcePlanarSlip")){
        std::string input = slotXML.get("ForcePlanarSlip");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        forcePlanarSlip(value);
    }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateVectorField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    PropagateField<3>::flushSectionXML(slotXML, name);
    slotXML.set("MultiStep", std::to_string(int(m_nstep)));
    slotXML.set("ForcePlanarSlip", std::to_string(int(m_forcePlanarSlip)));
};


/*!
 * Check coherence of the input data of the class, in particular:
   - check if Dirichlet surfaces belongs to the bulk mesh.
   - check if NarrowBand surfaces belongs to the bulk mesh.
   - check if Damping surfaces belongs to the bulk mesh.
   - check if Slip surfaces belongs to bulk mesh (this include periodic surfaces if any)
   - set the internal members m_slip_bc_dir (init to zero) and m_periodicBoundaryPoints
 * \return true if coherence is satisfied, false otherwise.
 */
bool PropagateVectorField::checkBoundariesCoherence(){

    // this check the part related to the dirichlet/damping/narrowband information
    if(!PropagateField<3>::checkBoundariesCoherence()){
        return false;
    }

    //check slip surfaces.(this include alse the periodic surfaces if any. See method addPeriodicBoundarySurface)
    std::unordered_set<long> nodeList;
    for(MimmoSharedPointer<MimmoObject> obj : m_slipSurfaces){
        std::vector<long> temp = obj->getVerticesIds(true); //only on internals
        nodeList.insert(temp.begin(), temp.end());
    }
    bitpit::PiercedVector<bitpit::Vertex> meshVertices = m_geometry->getVertices();
    for(long id : nodeList){
        if(!meshVertices.exists(id)){
            nodeList.clear();
            return false;
        }
    }

    //initialize to zero the m_slip_bc_dir structure
    m_slip_bc_dir.clear();
    std::array<double,3> zero({{0.0,0.0,0.0}});
    m_slip_bc_dir.setGeometry(m_geometry);
    m_slip_bc_dir.setDataLocation(MPVLocation::POINT);
    m_slip_bc_dir.reserve(nodeList.size());
    for(long id: nodeList){
        m_slip_bc_dir.insert(id, zero);
    }

    // run over periodic surfaces and retain border patch vertices common with the volume mesh
    m_periodicBoundaryPoints.clear();
    for (MimmoSharedPointer<MimmoObject> obj : m_periodicSurfaces){
        std::vector<long> tempboundary = obj->extractBoundaryVertexID(false); //no ghost, only internals
        //I have all internals and operations before guarantees me that
        //all internal points of periodic surfaces are present in the bulk mesh
        m_periodicBoundaryPoints.insert(tempboundary.begin(), tempboundary.end());
    }

    // all done.
	return true;
}

/*!
 * Subdivide m_bc_dir for multi step purposes.
 */
void
PropagateVectorField::subdivideBC(){
    for (auto & val : m_bc_dir){
        val /= double(m_nstep);
    }
}

/*!
 * restore m_bc_dir, after subdivision
 */
void
PropagateVectorField::restoreBC(){
    for (auto & val : m_bc_dir){
        val *= double(m_nstep);
    }
}

/*!
 * Apply deformation field to target bulk geometry.
   At the same time it deforms:
   - m_dampingUniSurface if any
   - m_bandUniSurface if any
   m_slipUniSurface is not meant to be deformed, no matter what.
   Input List surfaces are for recog/identification purpose only, and must not be deformed.
 */
void
PropagateVectorField::apply(){
    MimmoSharedPointer<MimmoObject> target = getGeometry();
    bitpit::PiercedVector<bitpit::Vertex> & verts = target->getVertices();

    //deform bulk.
    for (auto it= verts.begin(); it != verts.end(); ++it){
        target->modifyVertex( it->getCoords() + m_field.at(it.getId()), it.getId() );
    }

    dmpvecarr3E serialized_bf;

    //Check the damping and deform m_dampingUniSurface
    if(m_dampingActive && m_dampingUniSurface.get() != nullptr){

        if(serialized_bf.size() == 0) serialized_bf = getBoundaryPropagatedField();

        for(auto it= m_dampingUniSurface->getPatch()->vertexBegin(); it!= m_dampingUniSurface->getPatch()->vertexEnd(); ++it){
            m_dampingUniSurface->modifyVertex(it->getCoords() + serialized_bf.at(it.getId()), it.getId());
        }
    }

    //Check the narrowband and deform m_bandUniSurface
    if(m_bandActive && m_bandUniSurface.get() != nullptr){
        if(serialized_bf.size() == 0) serialized_bf = getBoundaryPropagatedField();

        for(auto it= m_bandUniSurface->getPatch()->vertexBegin(); it!= m_bandUniSurface->getPatch()->vertexEnd(); ++it){
            m_bandUniSurface->modifyVertex(it->getCoords() + serialized_bf.at(it.getId()), it.getId());
        }
    }
}
/*!
 * restore geometry to target vertices and re-evaluate m_field as whole
 * \param[in] vertices list  to be restored
 */
void
PropagateVectorField::restoreGeometry(bitpit::PiercedVector<bitpit::Vertex> & vertices){
    MimmoSharedPointer<MimmoObject> target = getGeometry();
    bitpit::PiercedVector<bitpit::Vertex> &currentmesh = target->getVertices();

    //restore bulk.
    long ID;
    for (auto it= vertices.begin(); it!=vertices.end(); ++it){
        ID = it.getId();
        const std::array<double,3> &coords = it->getCoords();
        m_field.at(ID) = currentmesh.at(ID).getCoords() - coords;
        target->modifyVertex(coords, ID);
    }

    //damping and narrowband if active have their UniSurface morphed.
    //But since the UniSurfaces are not inputs, but internal temp resources, I see no utility
    //in restoring them
}

/*!
 * Plot optional results on vtu unstructured grid file
 */
void
PropagateVectorField::plotOptionalResults(){

    if(getGeometry() == nullptr)    return;
    //Recover the scalar field
    dmpvecarr3E * field = getPropagatedField();
    //Recover the narrowBandDistances.
    MimmoPiercedVector<double> nbc(getGeometry(), MPVLocation::POINT);
    nbc = m_banddistances; //recover data of the raw PiercedVector;
    nbc.completeMissingData(1.E+18);
    //Set names
    field->setName("field");
    m_damping.setName("damping");
    nbc.setName("narrowband");
    //Write
    write(getGeometry(), *field, m_damping, nbc);

};


/*!
 * Given an ensemble of points marked as moving, enrich this list adding the first neighbours of
 * these vertices.
 *
 * \param[in,out] vertexList pool of moving points in input, 1Ring augmented in output.
 */
void PropagateVectorField::propagateMaskMovingPoints(livector1D & vertexList) {

    std::unordered_set<long> core(vertexList.begin(), vertexList.end());

    if (!getGeometry()->isPointConnectivitySync()){
        getGeometry()->buildPointConnectivity();
    }
    std::unordered_set<long> tempV1;
    for(long id: vertexList){
        tempV1 = getGeometry()->getPointConnectivity(id);
        core.insert(tempV1.begin(), tempV1.end());
    }
    vertexList.clear();
    vertexList.insert(vertexList.end(), core.begin(), core.end());
}

/*!
 * OVERRIDE Base class:
 * This method evaluate the bc corrections for a singular run of the system solver,
 * update the system matrix in m_solver and evaluate the rhs part due to bc.
 * After you call this method, you are typically ready to solve the laplacian system.
 * The type of bc @ nodes are directly desumed from class nodes member m_bc_dir and m_slip_bc_dir.
 * The method requires the Laplacian m_solver to be initialized. No ghost are taken into account.
 *
 * Slip conditions are calculated in two call of the method, using a predictor-corrector scheme.
 * First call creates a set of slip Neumann homogeneous bc condition. Equations are solved and
   a guess solution is predicted on candidate slip surfaces.
   At this point, computeSlipBCCorrector method is invoked. The guess solution on slip surfaces
   is corrected through a reprojection on reference slip surface and stored in the member m_slip_bc_dir.
 * On second call, the corrected m_slip_bc_dir is applied on nodes of slip surfaces, as a Dirichlet condition.
   The slipCorrect boolean rules the switch between slip predictor/corrector mode.
 *
 * \param[in] comp target component of bc conditions (0,1,2).
 * \param[in] slipCorrect true to apply the bc slip corrector step(Dirichlet), false for bc slip predictor(Neumann).
 * \param[in] borderLaplacianStencil list of laplacian Stencil on border nodes, where the bc is temporarely imposed as homogeneous Neumann
 * \param[in] maplocals map from global id numbering to local system solver numbering.
 * \param[in,out] rhs vector of right-hand-side's to append constant data from bc corrections.
 */
void
PropagateVectorField::assignBCAndEvaluateRHS(std::size_t comp, bool slipCorrect,
                                            GraphLaplStencil::MPVStencil * borderLaplacianStencil,
                                            const liimap & maplocals, dvector1D & rhs)
{

    //resize rhs to the number of internal cells
    MimmoSharedPointer<MimmoObject> geo = getGeometry();
    rhs.resize(geo->getNInternalVertices(), 0.0);

    if (!m_solver->isAssembled()) {
        (*m_log)<<"Warning in "<<m_name<<". Unable to assign BC to the system. The solver is not yet initialized."<<std::endl;
        return;
    }

    if (!borderLaplacianStencil) {
        (*m_log)<<"Warning in "<<m_name<<". Unable to reach border cells stencils data. Nothing to do."<<std::endl;
        return;
    }

    //copy laplacian stencils in a work mpv .
    GraphLaplStencil::MPVStencilUPtr lapwork(new GraphLaplStencil::MPVStencil(*borderLaplacianStencil));

    //correct the original border laplacian stencils applying the Dirichlet
    //conditions and slip conditions.
    //Neumann are implicitely imposed by graph-laplacian scheme.
    //renumber it and update the laplacian matrix and fill the rhs.
    bitpit::StencilScalar correction;

    //loop on all slip boundary nodes first.
    //Correct if it is the correction step. If not the neumann condition are automatically imposed by Graph-Laplace scheme.
    if(slipCorrect && !m_slipSurfaces.empty()){
        for(auto it = m_slip_bc_dir.begin(); it!=m_slip_bc_dir.end(); ++it){
            long id = it.getId();
            if (!geo->isPointInterior(id)) continue;

            //apply the correction relative to bc @ dirichlet node.
            correction.clear(true);
            correction.appendItem(id, 1.);
            correction.sumConstant(-1.0 * ((*it)[comp]) );
            //Fix to zero the old stencil (the update of system solver doesn't substitute but modify or append new pattern item and weights)
            lapwork->at(id) *= 0.;
            lapwork->at(id).setConstant(0.);
            //sum the correction
            lapwork->at(id) += correction;
        }
    }

    //add zero dirichlet for all periodic points if any
    //Loop on all periodic boundary points and force to be fixed (dirichlet 0)
    for (long id : m_periodicBoundaryPoints){

        if (!geo->isPointInterior(id)) continue;

        correction.clear(true);
        correction.appendItem(id, 1.);
        correction.sumConstant(0.);
        //Fix to zero the old stencil (the update of system solver doesn't substitute but modify or append new pattern item and weights)
        lapwork->at(id) *= 0.;
        lapwork->at(id).setConstant(0.);
        //sum the correction
        lapwork->at(id) += correction;
    }

    //loop on all dirichlet boundary nodes -> They have priority on all other conditions.
    for(auto it = m_bc_dir.begin(); it!=m_bc_dir.end(); ++it){
        long id = it.getId();
        if (!geo->isPointInterior(id)) continue;

        //apply the correction relative to bc @ dirichlet node.
        correction.clear(true);
        correction.appendItem(id, 1.);
        correction.sumConstant(-1.0* ( (*it)[comp]) );
        //Fix to zero the old stencil (the update of system solver doesn't substitute but modify or append new pattern item and weights)
        lapwork->at(id) *= 0.;
        lapwork->at(id).setConstant(0.);
        //sum the correction
        lapwork->at(id) += correction;
    }

    // now its time to update the solver matrix and to extract the rhs contributes.
    updateLaplaceSolver(lapwork.get(), maplocals);

    // now get the rhs
    for(auto it = lapwork->begin(); it != lapwork->end();++it){
        auto index = maplocals.at(it.getId());
#if MIMMO_ENABLE_MPI
        //correct index if in parallel
        index -= getGeometry()->getPointGlobalCountOffset();
#endif
        rhs[index] -= it->getConstant();
    }
}

/*!
 * Starting from a guess solution on mesh points, find a set of Dirichlet
   conditions on slip boundaries so that the guess deformation is reprojected
   onto the slip reference surface.
 * The reprojected value will be stored as bc in m_slip_bc_dir and it is used
   in assignBCAndEvaluateRHS forcing the boolean variable slipCorrect to true.
   See current class method assignBCAndEvaluateRHS doxy
 *
 * \param[in] guessSolutionOnPoint guess laplacian solution of the vector field on mesh POINTS
 */
void
PropagateVectorField::computeSlipBCCorrector(const MimmoPiercedVector<std::array<double,3> > & guessSolutionOnPoint){

    //first step: extract solutions on twin border nodes of m_slip_bc_dir;
    // I'm sure from checkBoundariesCoherence that m_slip_bc_dir share the same id
    // of guessSolutionOnPoint
    for(auto it=m_slip_bc_dir.begin(); it!=m_slip_bc_dir.end(); ++it){
        *it = guessSolutionOnPoint.at(it.getId());
    }
    //calculate projection of deformed mesh node under guessSolution deformation;
    // the projection deformation will be added to the guess Solution
    std::unordered_map<long, std::array<double,3> > projectionVector;

    if(m_forcePlanarSlip){

        // VERSION USING THE AVERAGE NORMAL!
        // average slip plane features (normal and point) must be already available
        // since initializeSlipSurfaceAsPlane method invoke.
        //loop on surface points.
        for(auto it=m_slip_bc_dir.begin(); it!=m_slip_bc_dir.end(); ++it){
            long idV = it.getId();
            std::array<double,3> point = m_geometry->getVertexCoords(idV) + *it;
            projectionVector[idV] = bitpit::CGElem::projectPointPlane(point, m_AVGslipCenter,m_AVGslipNormal) - point;
        }
    }else{

        //VERSION USING REFERENCE EXTERNAL SURFACE, that must be stored in m_slipUniSurface
        //at this point this surface must be allocated and not null. --> invoke initializeUniqueSurface
        //applied to list m_slipReferenceSurfaces
        bitpit::PatchSkdTree *tree = m_slipUniSurface->getSkdTree(); //(method directly build skdtree if not built)
        //loop on surface points.
        for(auto it=m_slip_bc_dir.begin(); it!=m_slip_bc_dir.end(); ++it){
            long idV = it.getId();
            std::array<double,3> point = m_geometry->getVertexCoords(idV) + *it;
            double r = std::max(1.0E-05, norm2(*it));
            projectionVector[idV] = skdTreeUtils::projectPoint(&point, tree, r) - point;
        }
    }

    // add the projectionVector correction on m_slip_bc_dir
    for(auto it=m_slip_bc_dir.begin(); it!=m_slip_bc_dir.end(); ++it){
        (*it) += projectionVector.at(it.getId());
    }
    //correction done.
}


/*!
 * Calculate Average Plane normal and point, from provided list of m_slipReferenceSurfaces,
   and store them in internal members m_AVGslipNormal and m_AVGslipCenter.
   m_slipUniSurface is ignored.
   This method must be called once and for all in execution, and only in case
   m_forcePlanarSlip is true and a non empty list of m_slipReferenceSurfaces is provided.
 */
void PropagateVectorField::initializeSlipSurfaceAsPlane(){

    m_AVGslipCenter.fill(0.0);
    m_AVGslipNormal.fill(0.0);
    long countV = 0;
    long countC = 0;

    //loop on m_slipReferenceSurfaces list
    for(MimmoSharedPointer<MimmoObject> obj : m_slipReferenceSurfaces){

#if MIMMO_ENABLE_MPI
        // be sure ghost information on point/cell are synchronized
        if(!obj->isInfoSync())                    obj->buildPatchInfo();
        if(!obj->arePointGhostExchangeInfoSync()) obj->updatePointGhostExchangeInfo();
#endif

        bitpit::SurfaceKernel * surfkernss = dynamic_cast<bitpit::SurfaceKernel*>(obj->getPatch());
        if(surfkernss == nullptr){
            throw std::runtime_error("PropagateVectorField::initializeSlipSurface -> SurfaceKernel dynamic cast failed!");
        }
        //start evaluating barycenter of interior points.
        for(auto itV = surfkernss->vertexBegin(); itV != surfkernss->vertexEnd(); ++itV){
            if(obj->isPointInterior(itV.getId())){
                m_AVGslipCenter += itV->getCoords();
                ++countV;
            }
        }

        //evaluating average facet normal of interior cells.
        for(auto itC = surfkernss->internalBegin(); itC != surfkernss->internalEnd(); ++itC){
            m_AVGslipNormal += surfkernss->evalFacetNormal(itC.getId());
            ++countC;
        }
    }

#if MIMMO_ENABLE_MPI
    //communicate with other processors summing barycenters, normals
    //number of vertices and number of cells and then average it.
    MPI_Allreduce(MPI_IN_PLACE, &countV, 1, MPI_LONG, MPI_SUM, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, &countC, 1, MPI_LONG, MPI_SUM, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, m_AVGslipCenter.data(), 3, MPI_DOUBLE, MPI_SUM, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, m_AVGslipNormal.data(), 3, MPI_DOUBLE, MPI_SUM, m_communicator);
#endif

    // perform the average
    if(countV > 0)  m_AVGslipCenter /= double(countV);
    if(countC > 0)  m_AVGslipNormal /= double(countC);

}

/*!
    Return values of m_field associated to boundary nodes of the bulk volume mesh.
    In case of MPI Version, returned structure is serialized, i.e. the boundary values
    collected are not only those owned by the local rank, but
    also those coming from all other ranks. In the end all ranks will had a copy
    of exactly the same structure.
 */
dmpvecarr3E PropagateVectorField::getBoundaryPropagatedField(){

    MimmoPiercedVector<std::array<double,3> > mpvres(m_geometry, MPVLocation::POINT);
    mpvres.reserve(m_field.size());

    //fill mpvres with the value of m_field on boundary nodes
    //FOR MPI version retain only values on interior boundary nodes.
    std::vector<long> bIds =  m_geometry->extractBoundaryVertexID(false);
    for(long id : bIds){
        mpvres.insert(id, m_field.at(id));
    }

#if MIMMO_ENABLE_MPI
    //MPI stuffs
    //Send my own and receive mpvres from others.

    //prepare my own send buffer.
    mimmo::OBinaryStream myrankDataBuffer;
    myrankDataBuffer << (std::size_t)mpvres.size();
    auto itE = mpvres.cend();
    for (auto it=mpvres.cbegin(); it!=itE; it++){
        myrankDataBuffer << it.getId();
        myrankDataBuffer << *it;
    }

    long myrankDataBufferSize = myrankDataBuffer.getSize();

    for (int sendRank=0; sendRank<m_nprocs; sendRank++){

        if (m_rank != sendRank){
            // receive data from other ranks.
            long defBufferSize;
            MPI_Recv(&defBufferSize, 1, MPI_LONG, sendRank, 900, m_communicator, MPI_STATUS_IGNORE);
            mimmo::IBinaryStream defBuffer(defBufferSize);
            MPI_Recv(defBuffer.data(), defBuffer.getSize(), MPI_CHAR, sendRank, 910, m_communicator, MPI_STATUS_IGNORE);

            MimmoPiercedVector<std::array<double,3>> temp;
            std::size_t nP;
            defBuffer >> nP;
            std::array<double,3> val;
            long int id;
            for (std::size_t i = 0; i < nP; ++i) {
                defBuffer >> id;
                defBuffer >> val;
                temp.insert(id, val);
            }

            // insert this part in mpvres.
            for (auto it = temp.begin(); it!=temp.end(); ++it) {
                if(!mpvres.exists(it.getId())){
                    mpvres.insert(it.getId(), *it);
                }
            }
        }else{
            //send to all other except me the def data.
            for (int recvRank=0; recvRank<m_nprocs; recvRank++){
                if (m_rank != recvRank){
                    MPI_Send(&myrankDataBufferSize, 1, MPI_LONG, recvRank, 900, m_communicator);
                    MPI_Send(myrankDataBuffer.data(), myrankDataBuffer.getSize(), MPI_CHAR, recvRank, 910, m_communicator);
                }
            }
        }
    }// end external sendrank loop

    MPI_Barrier(m_communicator);

#endif

    // shrink structure
    mpvres.shrinkToFit();
    //and return
    return mpvres;
}

/*!
 * Execution command. After the execution the result constraint field is stored in the class.
 */
void
PropagateVectorField::execute(){

    MimmoSharedPointer<MimmoObject> geo = getGeometry();
    if(!geo){
        (*m_log)<<"Error in "<<m_name<<" .No target volume mesh linked"<<std::endl;
        throw std::runtime_error("Error in "+m_name+" .No target volume mesh linked");
    }

    if(m_dirichletSurfaces.empty()){
        (*m_log)<<"Warning in "<<m_name<<" .No Dirichlet Boundary patch linked"<<std::endl;
    }

#if MIMMO_ENABLE_MPI
    //be sure ghost info are available
    if(!geo->isInfoSync()) geo->buildPatchInfo();
    if(!geo->arePointGhostExchangeInfoSync()) geo->updatePointGhostExchangeInfo();
#endif

    if(!checkBoundariesCoherence()){
        (*m_log)<<"Warning in "<<m_name<<" .Boundary patches linked are uncoherent with target bulk geometry"<<std::endl;
    }
    //after this call m_slip_bc_dir is initialized and periodic points stored, in case.


    (*m_log) << bitpit::log::priority(bitpit::log::NORMAL);
    (*m_log) << bitpit::log::context("mimmo");
    //allocate the solver;
    m_solver = std::unique_ptr<bitpit::SystemSolver>(new bitpit::SystemSolver(m_print));

    //declare inverse and direct map
    liimap dataInv, data;

    //check if damping or narrow band control are active,
    //initialize their reference surfaces and compute them
    if(m_dampingActive){
        if(m_dampingSurfaces.empty())   m_dampingSurfaces = m_dirichletSurfaces;
        initializeUniqueSurface(m_dampingSurfaces, m_dampingUniSurface);
    }

    if(m_bandActive){
        if(m_bandSurfaces.empty())   m_bandSurfaces = m_dirichletSurfaces;
        initializeUniqueSurface(m_bandSurfaces, m_bandUniSurface);
    }

    //check the slip part
    if(!m_slipSurfaces.empty()){
        if(m_slipReferenceSurfaces.empty())   m_slipReferenceSurfaces = m_slipSurfaces;
        if(m_forcePlanarSlip){
            initializeSlipSurfaceAsPlane(); //-> this is useful for the implicit plane reproj
        }else{
            initializeUniqueSurface(m_slipReferenceSurfaces, m_slipUniSurface);  //-> this is needed for the external surface reproj
        }
    }

    //initialize damping and Narrow band Control. If UniSurfaces are null the methods set:
    // - unitary m_damping field.
    // - empty m_banddistances member.
    //
    initializeDampingFunction();
    updateNarrowBand();

    // Graph Laplace method on points

    //store the id of the border nodes only;
    livector1D borderPointsID = geo->extractBoundaryVertexID(false);

    //get this inverse map -> you will need it to compact the stencils.
    dataInv = geo->getMapDataInv(true);
    //get this direct map -> you will need it to deflate compact solution of the system.
    data = geo->getMapData(true);

    // compute the laplacian stencils
    GraphLaplStencil::MPVStencilUPtr laplaceStencils = GraphLaplStencil::computeLaplacianStencils(geo, m_tol, &m_damping);

    //modify stencils if Narrow band is active i.e. m_banddistances is not empty.
    //This is directly managed in the method.
    modifyStencilsForNarrowBand(laplaceStencils);

    // initialize the laplacian Matrix in solver and squeeze out the laplace stencils and save border cells only.
    initializeLaplaceSolver(laplaceStencils.get(), dataInv);
    laplaceStencils->squeezeOutExcept(borderPointsID);
    borderPointsID.clear();

    //declare results here and keep it during the loop to re-use the older steps.
    std::vector<std::vector<double>> results(3);

    //since bc is constant, even in case of multistep, once and for all
    //pass dirichlet bc point information to bulk m_bc_dir internal member.
    distributeBCOnBoundaryPoints();
    //PREPARE THE MULTISTEP;
    bitpit::PiercedVector<bitpit::Vertex> undeformedTargetVertices;
    std::unique_ptr<livector1D> movingElementList = nullptr;
    if(m_nstep > 1){
        subdivideBC(); //this subdivide m_bc_dir just prepared in distributeBCOnBoundaryPoints
        undeformedTargetVertices = geo->getVertices();
        movingElementList = std::unique_ptr<livector1D>(new livector1D());
    }

    //loop on multistep
    for(int istep=0; istep < m_nstep; ++istep){

        //3-COMPONENT SYSTEM SOLVING ---> ///////////////////////////////////////////////////////////////////////
        // solve the field component by component
        //first loop -> if slip is enforced in some walls, this is the PREDICTOR
        //stage of guess solution with 0-Neumann on slip walls
        for(int comp = 0; comp<3; ++comp){
            //prepare the right hand side ;
            dvector1D rhs(geo->getNInternalVertices(), 0.0);
            assignBCAndEvaluateRHS(comp, false, laplaceStencils.get(), dataInv, rhs);
            //solve
            results[comp].resize(rhs.size(), 0.0);
            solveLaplace(rhs, results[comp]);
        }

        //if I have slip walls active, it needs a corrector stage for slip boundaries;
        if(!m_slipSurfaces.empty()){
            //reconstruct result on mesh points (ghost included)-> stored in m_field.
            reconstructResults(results, data);
            //compute the correction/reprojection @ slip walls
            computeSlipBCCorrector(m_field);
            //now you have a set of BC Dirichlet condition m_slip_bc_dir internal.
            // so loop again on the components, reusing the previous result as starting guess, and setting
            // the boolean of slipCorrect to true (corrector stage of slip, read Dirichlet from m_slip_bc_dir)
            for(int comp = 0; comp<3; ++comp){
                //prepare the right hand side ;
                dvector1D rhs(geo->getNInternalVertices(), 0.0);
                assignBCAndEvaluateRHS(comp, true, laplaceStencils.get(), dataInv, rhs);
                //solve
                solveLaplace(rhs, results[comp]);
            }
        }

        //RECONSTRUCT M_FIELD --> //////////////////////////////////////////////
        // get the list of moving nodes also
        reconstructResults(results, data, movingElementList.get());

        // if you are in multistep stage apply the deformation to the mesh.
        // here is morphed also the m_dampingUniSurface if active dumping and
        // m_bandUniSurface if active narrow band. Slip surfaces are untouched.
        if(m_nstep > 1){
            apply();
        }

        //UPDATING LAPLACIAN, DAMPING, NARROW BAND etc...
        // if in multistep stage continue to update the other laplacian stuff up to "second-to-last" step.
        if(istep < m_nstep-1){

            //if damping active, update the damping function using m_dumpingUniSurface deformed.
            updateDampingFunction();
            //if narrowband control active, update the narrowband, with m_bandUnisurface deformed
            updateNarrowBand();

            //enlarge the list of marked moving nodes adding the 1st vertex ring.
            propagateMaskMovingPoints(*(movingElementList.get()));

            // update the laplacian stencils
            GraphLaplStencil::MPVStencilUPtr updateLaplaceStencils = GraphLaplStencil::computeLaplacianStencils(geo, movingElementList.get(), m_tol, &m_damping);
            movingElementList->clear();

            //apply modification to the interested stencils if narrow band control is active
            modifyStencilsForNarrowBand(updateLaplaceStencils);

            //store the update boundary stencils in laplaceStencils structure
            laplaceStencils->getDataFrom(*(updateLaplaceStencils.get()), true); //only common elements are updated.

            // update the laplacian Matrix in solver free the updateLaplaceStencils
            updateLaplaceSolver(updateLaplaceStencils.get(), dataInv);
            //clear updateLaplaceStencils
            updateLaplaceStencils = nullptr;
        }

        (*m_log)<<m_name<<" solved step "<<istep+1<<" out of total steps "<<m_nstep<<std::endl;

    } //end of multistep loop;


    if(m_nstep > 1){
        //this take the geometry to the original state and update the deformation field as the current
        //deformed grid minus the undeformed state (this directly on POINTS).
        restoreGeometry(undeformedTargetVertices);
    }

    //free some memory from internal stuff

    //clear UniSurface for Damping and NarrowBand if any
    m_bandUniSurface = nullptr;
    m_dampingUniSurface = nullptr;

    //clear slipUnisurface and m_periodicBoundaryPoints
    m_slipUniSurface = nullptr;
    m_periodicBoundaryPoints.clear();

    //clear temp bc;
    m_bc_dir.clear();
    m_slip_bc_dir.clear();

    //clear the solver;
    m_solver->clear();

    (*m_log) << bitpit::log::priority(bitpit::log::DEBUG);
}

} //end of mimmo namespace
