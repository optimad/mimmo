
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
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
PropagateScalarField::PropagateScalarField(const bitpit::Config::Section & rootXML):PropagateField<1>(){

    m_name = "mimmo.PropagateScalarField";

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
PropagateScalarField::PropagateScalarField(const PropagateScalarField & other):PropagateField<1>(other){};

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
    PropagateField<1>::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
PropagateScalarField::buildPorts(){
    bool built = true;
    built = (built && createPortIn<dmpvector1D, PropagateScalarField>(this, &PropagateScalarField::setDirichletConditions, M_FILTER));
    built = (built && createPortOut<dmpvector1D, PropagateScalarField>(this, &PropagateScalarField::getPropagatedField, M_FILTER));
    PropagateField<1>::buildPorts();
    m_arePortsBuilt = built;
};

/*!
 * It gets the resulting propagated field on the whole bulk mesh.
 * \return Deformation field.
 */
dmpvector1D
PropagateScalarField::getPropagatedField(){
    dmpvector1D field;
    field.reserve(m_field.size());
    field.setDataLocation(m_field.getDataLocation());
    field.setGeometry(m_field.getGeometry());
    for(auto it = m_field.begin(); it != m_field.end(); ++it){
        field.insert(it.getId(), (*it)[0]);
    }

    return(field);
}

/*!
 * It sets the Dirichlet conditions for scalar field on the previously linked
 * Dirichlet Boundary patch.
 * \param[in] bc dirichlet conditions
 */
void
PropagateScalarField::setDirichletConditions(dmpvector1D bc){
    if (bc.isEmpty()) return;
    m_surface_bc_dir.reserve(bc.size());
    m_surface_bc_dir.setDataLocation(bc.getDataLocation());
    m_surface_bc_dir.setGeometry(bc.getGeometry());
    for(auto it = bc.begin(); it != bc.end(); ++it){
        m_surface_bc_dir.insert(it.getId(), std::array<double,1>({*it}));
    }
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
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateScalarField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    PropagateField<1>::flushSectionXML(slotXML, name);
};

/*!
 * Clear all data actually stored in the class
 */
void
PropagateScalarField::clear(){
    PropagateField<1>::clear();
};


/*!
 * Check coherence of the input data of the class, in particular:
 * - check if boundary patches are referred to target bulk mesh
 * - check if boundary conditions vectors are coherent with relative boundary patches.
 * In case of successfull check, the internal member m_ibsp is filled out.
 * \return true if coherence is satisfied, false otherwise.
 */
bool PropagateScalarField::checkBoundariesCoherence(){

    //Clean the old m_isbp and initialize it again.
    initializeBoundaryInfo();
    //1st step verify coherence of the Dirichlet point field on boundary surface
    // with the dirichlet boundary surface provided
    if(m_surface_bc_dir.getGeometry() != m_bsurface || !m_surface_bc_dir.completeMissingData({0.0})){
        return false;
    }

    //transfer point field info of boundary dirichlet on the volume mesh interfaces.
    MimmoPiercedVector<std::array<double,1>> temp;
    temp.setGeometry(m_geometry);
    temp.setDataLocation(MPVLocation::POINT);
    temp.reserve(m_surface_bc_dir.size());
    for(auto it=m_surface_bc_dir.begin(); it!=m_surface_bc_dir.end(); ++it){
        temp.insert(it.getId(), *it );
    }
    //interpolate now point data to interface data
    m_bc_dir.clear();
    m_bc_dir = temp.pointDataToBoundaryInterfaceData();
    if(m_bc_dir.isEmpty()){
        return false;
    }

    //update the m_isbp marking dirichlet interfaces
    for(auto it=m_bc_dir.begin(); it != m_bc_dir.end(); ++it){
        m_isbp.at(it.getId()) = 1;
    }

    return true;
}

/*!
 * Plot optional results on vtu unstructured grid file
 */
void
PropagateScalarField::plotOptionalResults(){


    if(getGeometry() == NULL || getGeometry()->isEmpty())    return;

    bitpit::VTKUnstructuredGrid& vtk = getGeometry()->getPatch()->getVTK();

    std::vector<std::array<double,1> > dataraw = m_field.getDataAsVector();
    dvector1D data;
    data.reserve(dataraw.size());
    for (auto val : dataraw){
        data.push_back(val[0]);
    }
    vtk.addData("field", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, data);

    dvector1D datad = m_dumping.getDataAsVector() ;
    vtk.addData("dumping", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, datad);

    vtk.setCounter(getId());
    getGeometry()->getPatch()->write(m_name +"_field");
    vtk.removeData("field");
    vtk.removeData("dumping");
    vtk.unsetCounter();

};

/*!
 * Execution command. After the execution the result constraint field is stored in the class.
 */
void
PropagateScalarField::execute(){

    if(getGeometry() == NULL){
        (*m_log)<<"Error in "<<m_name<<" .No target volume mesh linked"<<std::endl;
        throw std::runtime_error("Error in PropagateScalarField execute. No target volume mesh linked");
    }

    if(m_bsurface == NULL ){
        (*m_log)<<"Error in "<<m_name<<" .No Dirichlet Boundary patch linked"<<std::endl;
        throw std::runtime_error("Error in PropagateScalarField execute. No Dirichlet Boundary patch linked");
    }

    if(!checkBoundariesCoherence()){
        (*m_log)<<"Error in "<<m_name<<" .Boundary patches linked are uncoherent with target bulk geometry"
        "or bc-fields not coherent with boundary patches"<<std::endl;
        throw std::runtime_error("Error in PropagateScalarField execute. Boundary patches linked are uncoherent"
        "with target bulk geometry or bc-fields not coherent with boundary patches");
    }
    //allocate the solver;
    m_solver = std::unique_ptr<bitpit::SystemSolver>(new bitpit::SystemSolver());

    //get this inverse map -> you will need it to compact the stencils.
    liimap dataInv = getGeometry()->getMapCellInv();
    //compute the dumping.
    computeDumpingFunction();
    //compute the gradient stencils @ interface with Neutral conditions.
    FVolStencil::MPVGradientUPtr faceGradients  = computeGradientStencilsWithNeutralBC();
    // compute the laplacian stencils ;
    FVolStencil::MPVDivergenceUPtr laplaceStencils = FVolStencil::computeFVLaplacianStencil(*(faceGradients.get()), &m_dumping);
    // initialize the laplacian Matrix in solver and release the laplace stencils.
    initializeLaplaceSolver(laplaceStencils.get(), dataInv);
    laplaceStencils = nullptr;

    //prepare the right hand side and release the faceGradients;
    dvector1D rhs(getGeometry()->getPatch()->getInternalCount(), 0.0);
    appendToRHSFromBorderFluxes(0, faceGradients.get(), dataInv, rhs);
    faceGradients = nullptr;
    dataInv.clear();

    //solve
    dvector1D result(rhs.size(), 0.0);
    solveLaplace(rhs, result);

    // push result in a mpv linked to target mesh and on cell location.
    MimmoPiercedVector<std::array<double,1> > tempres;
    tempres.setGeometry(m_geometry);
    tempres.setDataLocation(MPVLocation::CELL);
    tempres.reserve(result.size());
    livector1D cellmap = getGeometry()->getMapCell();
    int counter = 0;
    for(const double & val: result){
        tempres.insert(cellmap[counter], std::array<double,1>({val}));
        ++counter;
    }
    cellmap.clear();
    // interpolate result to POINT location.
    m_field = tempres.cellDataToPointData();

    //clear the solver;
    m_solver->clear();
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
    m_nstep = 1;
    m_slipsurface = NULL;
};

/*!
 * Set most significant parameters to constructor defaults
 */
void PropagateVectorField::setDefaults(){
    PropagateField<3>::setDefaults();
    m_nstep = 1;
    m_slipsurface = NULL;
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
PropagateVectorField::PropagateVectorField(const bitpit::Config::Section & rootXML):PropagateField<3>(){

    m_name = "mimmo.PropagateVectorField";
    m_nstep = 1;
    m_slipsurface = NULL;

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
PropagateVectorField::PropagateVectorField(const PropagateVectorField & other):PropagateField<3>(other){
    m_slipsurface = other.m_slipsurface;
    m_nstep = other.m_nstep;
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
    std::swap(m_slipsurface, x.m_slipsurface);
    std::swap(m_nstep, x.m_nstep);
    PropagateField<3>::swap(x);
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
    PropagateField<3>::buildPorts();
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

    //TODO turn on when it's ready
    //m_slipsurface = surface;
    m_slipsurface = nullptr;

}

/*!
 * It sets the Dirichlet conditions for each component of the vector field on the previously linked
 * Dirichlet Boundary patch.
 * \param[in] bc dirichlet conditions
 */
void
PropagateVectorField::setDirichletConditions(dmpvecarr3E bc){
    if (bc.isEmpty()) return;
    m_surface_bc_dir = bc;
}

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
        if(value >= 0.0) value2 = value;
        setSolverMultiStep(value2);
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
};


/*!
 * Clear all data actually stored in the class
 */
void
PropagateVectorField::clear(){
    PropagateField<3>::clear();
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

    //Clean the old m_isbp and initialize it again (set all to Neumann).
    initializeBoundaryInfo();


    //2nd step verify coherence of the Dirichlet point field on boundary surface
    // with the dirichlet boundary surface provided
    if(m_surface_bc_dir.getGeometry() != m_bsurface || !m_surface_bc_dir.completeMissingData({{0.0, 0.0,0.0}})){
        return false;
    }

    //transfer point field info of boundary dirichlet on the volume mesh interfaces.
    MimmoPiercedVector<std::array<double,3>> temp;
    temp.setGeometry(m_geometry);
    temp.setDataLocation(MPVLocation::POINT);
    temp.reserve(m_surface_bc_dir.size());
    for(auto it=m_surface_bc_dir.begin(); it!=m_surface_bc_dir.end(); ++it){
        temp.insert(it.getId(), *it );
    }
    //interpolate now point data to interface data
    m_bc_dir.clear();
    m_bc_dir = temp.pointDataToBoundaryInterfaceData();

    //check the part of slip surface condition.
    if(m_slipsurface){
        std::vector<long> slipInterfaceList = getGeometry()->getInterfaceFromVertexList(m_slipsurface->getVertices().getIds(), true, true);

        //update the m_isbp marking slip interfaces (2)
        for(long id : slipInterfaceList){
            m_isbp.at(id) = 2;
        }
    }
    //update the m_isbp marking dirichlet interfaces
    for(auto it=m_bc_dir.begin(); it != m_bc_dir.end(); ++it){
        m_isbp.at(it.getId()) = 1;
    }

    return true;
}


/*!
 * Plot optional results on vtu unstructured grid file
 */
void
PropagateVectorField::plotOptionalResults(){


    if(getGeometry() == NULL || getGeometry()->isEmpty())    return;

    bitpit::VTKUnstructuredGrid& vtk = getGeometry()->getPatch()->getVTK();
    dvecarr3E data = m_field.getDataAsVector();
    vtk.addData("field", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, data);

    dvector1D datad = m_dumping.getDataAsVector();
    vtk.addData("dumping", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, datad);

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
    if (!getGeometry()) return;
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
 * Force solver to get deformation in a finite number of substep
 * \param[in] sstep number of substep. Default is 1.
 */
void
PropagateVectorField::setSolverMultiStep(unsigned int sstep){
    unsigned int loc(1);
    m_nstep = std::max(loc,sstep);
}

/*!
 * subdivide dirichlet boundary conditions  for multi step purposes
 */
void
PropagateVectorField::subdivideBC(){
    for (auto & val : m_bc_dir){
        val /= double(m_nstep);
    }
}

/*!
 * restore dirichlet boundary conditions for multi step purposes
 */
void
PropagateVectorField::restoreBC(){
    for (auto & val : m_bc_dir){
        val *= double(m_nstep);
    }
}

/*!
 * restore geometry to target vertices
 * \param[in] list of vertices to be restored
 */

void
PropagateVectorField::restoreGeometry(bitpit::PiercedVector<bitpit::Vertex> & vertices){

    for (bitpit::Vertex & vertex : vertices){
        long ID = vertex.getId();
        m_field[ID] = getGeometry()->getVertexCoords(ID) - vertex.getCoords();
        getGeometry()->modifyVertex(vertex.getCoords(), ID);
    }

}

/*!
 * Execution command. After the execution the result constraint field is stored in the class.
 */
void
PropagateVectorField::execute(){

    if(getGeometry() == NULL){
        (*m_log)<<"Error in "<<m_name<<" .No target volume mesh linked"<<std::endl;
        throw std::runtime_error("Error in PropagateVectorField execute. No target volume mesh linked");
    }

    if(m_bsurface == NULL ){
        (*m_log)<<"Error in "<<m_name<<" .No Dirichlet Boundary patch linked"<<std::endl;
        throw std::runtime_error("Error in PropagateVectorField execute. No Dirichlet Boundary patch linked");
    }

    if(!checkBoundariesCoherence()){
        (*m_log)<<"Error in "<<m_name<<" .Boundary patches linked are uncoherent with target bulk geometry"
        "or bc-fields not coherent with boundary patches"<<std::endl;
        throw std::runtime_error("Error in PropagateVectorField execute. Boundary patches linked are uncoherent"
        "with target bulk geometry or bc-fields not coherent with boundary patches");
    }
    //allocate the solver;
    m_solver = std::unique_ptr<bitpit::SystemSolver>(new bitpit::SystemSolver());

    //get this inverse map -> you will need it to compact the stencils.
    liimap dataInv = getGeometry()->getMapCellInv();

    //compute the dumping.
    computeDumpingFunction();

    //compute the gradient stencils @ interface with Neutral conditions.
    FVolStencil::MPVGradientUPtr faceGradients  = computeGradientStencilsWithNeutralBC();

    // compute the laplacian stencils ;
    FVolStencil::MPVDivergenceUPtr laplaceStencils = FVolStencil::computeFVLaplacianStencil(*(faceGradients.get()), &m_dumping);

    // initialize the laplacian Matrix in solver and release the laplace stencils.
    initializeLaplaceSolver(laplaceStencils.get(), dataInv);
    laplaceStencils = nullptr;

    // solve the field component by component
    std::array<std::vector<double>, 3> results;
    for(int comp = 0; comp<3; ++comp){
        //prepare the right hand side and release the faceGradients;
        dvector1D rhs(getGeometry()->getPatch()->getInternalCount(), 0.0);
        appendToRHSFromBorderFluxes(comp, faceGradients.get(), dataInv, rhs);
        //solve
        results[comp].resize(rhs.size(), 0.0);
        solveLaplace(rhs, results[comp]);
    }
    faceGradients = nullptr;
    dataInv.clear();

    // push result in a mpv linked to target mesh and on cell location.
    MimmoPiercedVector<std::array<double,3> > tempres;
    std::size_t locsize = results[0].size();
    tempres.setGeometry(m_geometry);
    tempres.setDataLocation(MPVLocation::CELL);
    tempres.reserve(locsize);
    livector1D cellmap = getGeometry()->getMapCell();
    int counter = 0;
    for(int counter = 0; counter<locsize; ++counter){
        tempres.insert(cellmap[counter], std::array<double,3>({results[0][counter],results[1][counter],results[2][counter]}));
    }
    cellmap.clear();
    // interpolate result to POINT location.
    m_field = tempres.cellDataToPointData();

    //clear the solver;
    m_solver->clear();


//     livector2D stencils;
//     dvector2D weights;
//     dvector1D rhs;
//     liimap dataInv = getGeometry()->getMapDataInv();
//
//     if (m_laplace){
//         bitpit::PiercedVector<bitpit::Vertex> vertices0;
//         if (m_nstep > 1){
//             subdivideBC();
//             vertices0      = getGeometry()->getVertices();
//         }
//
//
//         for(int istep=0; istep<m_nstep; istep++){
//             //TODO need to provide an updater for stencils weights and dumping function to speed up the multi-step stage.
//             computeDumpingFunction();
//             computeStencils(dataInv, stencils, weights);
//             correctStencils(dataInv, stencils, weights);
//             computeRHS(m_bc_dir, dataInv, rhs);
//
//             solveLaplace(stencils, weights, rhs, dataInv, m_field);
//
//             if (m_nstep > 1){
//                 apply();
//                 if(m_dumpingActive){
//                     //update vertices of candidate dumping surface
//                     MimmoObject * dumptarget = m_dsurface;
//                     if(dumptarget == NULL) dumptarget = m_bsurface;
//                     for(const auto & vert: dumptarget->getVertices()){
//                         dumptarget->modifyVertex(getGeometry()->getVertexCoords(vert.getId()), vert.getId());
//                     }
//                 }
//                 (*m_log)<<"                        "<<m_name<<" performing substep :"<<std::to_string(istep+1)<<std::endl;
//
//             }
//         }//end loop step
//
//         if (m_nstep > 1){
//             restoreGeometry(vertices0);
//             if(m_dumpingActive){
//                 //update vertices of candidate dumping surface
//                 MimmoObject * dumptarget = m_dsurface;
//                 if(dumptarget == NULL) dumptarget = m_bsurface;
//                 for(const auto & vert: dumptarget->getVertices()){
//                     dumptarget->modifyVertex(getGeometry()->getVertexCoords(vert.getId()), vert.getId());
//                 }
//             }
//         }
//
//         restoreBC();
//
//     }else{
//         computeDumpingFunction();
//         computeStencils(dataInv, stencils, weights);
//         correctStencils(dataInv, stencils, weights);
//         computeRHS(m_bc_dir, dataInv, rhs);
//
//         solveSmoothing(m_sstep, stencils, weights, rhs, dataInv, m_field);
//     }
}

} //end of mimmo namespace
