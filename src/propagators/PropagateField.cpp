
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
	m_thres = -1.0;
	m_nstep = 1;
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
	m_nstep = 1;
	m_thres = -1.0;

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
//	if (bc.isEmpty()) return;
	m_surface_bc_dir.reserve(bc.size());
	m_surface_bc_dir.setDataLocation(bc.getDataLocation());
	m_surface_bc_dir.setGeometry(bc.getGeometry());
	for(auto it = bc.begin(); it != bc.end(); ++it){
		m_surface_bc_dir.insert(it.getId(), std::array<double,1>({*it}));
	}
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
		if(value >= 1.0) value2 = value;
		setSolverMultiStep(value2);
	}
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void PropagateScalarField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);
	slotXML.set("MultiStep",std::to_string(m_nstep));
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
 * Plot optional results on vtu unstructured grid file
 */
void
PropagateScalarField::plotOptionalResults(){

	if(getGeometry() == NULL)    return;

	bitpit::VTKUnstructuredGrid& vtk = getGeometry()->getPatch()->getVTK();
#if MIMMO_ENABLE_MPI
	getGeometry()->getPatch()->setVTKWriteTarget(bitpit::PatchKernel::WriteTarget::WRITE_TARGET_CELLS_INTERNAL);
#endif

	std::vector<std::array<double,1> > dataraw = m_field.getInternalDataAsVector();
	dvector1D data;
	data.reserve(dataraw.size());
	for (auto val : dataraw){
		data.push_back(val[0]);
	}
	vtk.addData("field", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, data);

	dvector1D datad = m_dumping.getInternalDataAsVector();
	vtk.addData("dumping", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, datad);

	vtk.setCounter(getId());
    getGeometry()->getPatch()->getVTK().setDirectory(m_outputPlot+"/");
    getGeometry()->getPatch()->getVTK().setName(m_name+"_field");
	getGeometry()->getPatch()->write();
	vtk.removeData("field");
	vtk.removeData("dumping");
	vtk.unsetCounter();

};

/*!
 * Execution command. After the execution the result constraint field is stored in the class.
 */
void
PropagateScalarField::execute(){

	MimmoObject * geo = getGeometry();
	if(!geo){
		(*m_log)<<"Warning in "<<m_name<<" .No target volume mesh linked"<<std::endl;
	}

	if(!m_bsurface ){
		(*m_log)<<"Warning in "<<m_name<<" .No Dirichlet Boundary patch linked"<<std::endl;
	}

	if(!checkBoundariesCoherence()){
		(*m_log)<<"Warning in "<<m_name<<" .Boundary patches linked are uncoherent with target bulk geometry"
				"or bc-fields not coherent with boundary patches"<<std::endl;
	}

	(*m_log) << bitpit::log::priority(bitpit::log::NORMAL);
	(*m_log) << bitpit::log::context("mimmo");

	//allocate the solver;
	m_solver = std::unique_ptr<bitpit::SystemSolver>(new bitpit::SystemSolver(m_print));

	//get this inverse map -> you will need it to compact the stencils.
	liimap dataInv;


	//compute the dumping.
    initializeDumpingSurface();
    computeDumpingFunction();

	//Switch on solver method
	if (m_method == PropagatorMethod::FINITEVOLUMES){

		// Finite Volumes method

		//store the id of the border cells only;
		livector1D borderCellsID = geo->extractBoundaryCellID(false);

		//get this inverse map -> you will need it to compact the stencils.
		dataInv = geo->getMapCellInv(true);

		//pass bc point information to bulk interfaces.
		distributeBCOnBoundaryInterfaces();

		//compute the center cell gradients.
		FVolStencil::MPVGradientUPtr ccellGradients = FVolStencil::computeFVCellGradientStencil(*geo);
		//compute the gradient stencils @ interface with homogeneous Neumann.
		FVolStencil::MPVGradientUPtr faceGradients  = FVolStencil::computeFVFaceGradientStencil(*geo, ccellGradients.get());
		//and squeeze out cell gradients and save the border cells only.
		ccellGradients->squeezeOutExcept(borderCellsID);

		// compute the laplacian stencils and free faceGradients;
		FVolStencil::MPVDivergenceUPtr laplaceStencils = FVolStencil::computeFVLaplacianStencil(*(faceGradients.get()), m_tol, &m_dumping);
		faceGradients  = nullptr;

		// initialize the laplacian Matrix in solver and squeeze out the laplace stencils and save border cells only.
		initializeLaplaceSolver(laplaceStencils.get(), dataInv);
		laplaceStencils->squeezeOutExcept(borderCellsID);
		borderCellsID.clear();

		MimmoPiercedVector<std::array<double, 1> > stepBCdir(geo, MPVLocation::INTERFACE);
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
			dvector1D rhs(geo->getPatch()->getInternalCount(), 0.0);
			assignBCAndEvaluateRHS(0, false, laplaceStencils.get(), ccellGradients.get(), dataInv, rhs);
			solveLaplace(rhs, result[0]);
			(*m_log)<<m_name<<" solved step "<<istep+1<<" out of total steps "<<m_nstep<<std::endl;
		}

		dataInv.clear();
		//reconstruct getting the direct cell map -> you will need it uncompact the system solution in global id.
		reconstructResults(result, geo->getMapCell(true));
		// now data are direcly pushed in m_field.

		if (m_forceDirichletConditions){
			//force boundary Dirichlet on POINTS on m_field;
			for(auto it=m_surface_bc_dir.begin(); it != m_surface_bc_dir.end(); ++it){
				m_field.at(it.getId()) = *it;
			}
		}

	}
	else if (m_method == PropagatorMethod::GRAPHLAPLACE){

		// Graph Laplace method on points

		//store the id of the border nodes only;
		livector1D borderPointsID = geo->extractBoundaryVertexID(false);

		//get this inverse map -> you will need it to compact the stencils.
		dataInv = geo->getMapDataInv(true);

		//pass bc point information to bulk interfaces.
		distributeBCOnBoundaryPoints();

		// compute the laplacian stencils
		GraphLaplStencil::MPVStencilUPtr laplaceStencils = GraphLaplStencil::computeLaplacianStencils(*geo, m_tol, &m_dumping);

		// initialize the laplacian Matrix in solver and squeeze out the laplace stencils and save border cells only.
		initializeLaplaceSolver(laplaceStencils.get(), dataInv);
		laplaceStencils->squeezeOutExcept(borderPointsID);
		borderPointsID.clear();

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

	}// end if solver method

//#if MIMMO_ENABLE_MPI
//	communicatePointGhostData(&m_field);
//#endif

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
	m_nstep = 1;
	m_slipsurface = nullptr;
    m_forcePlanarSlip = false;
};

/*!
 * Set most significant parameters to constructor defaults
 */
void PropagateVectorField::setDefaults(){
	PropagateField<3>::setDefaults();
	m_nstep = 1;
	m_slipsurface = nullptr;
    m_forcePlanarSlip = false;
	m_slip_bc_dir.clear();
	m_surface_slip_bc_dir.clear();
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
PropagateVectorField::PropagateVectorField(const bitpit::Config::Section & rootXML):PropagateField<3>(){

	m_name = "mimmo.PropagateVectorField";
	m_nstep = 1;
	m_slipsurface = nullptr;
    m_forcePlanarSlip = false;

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
	m_slip_bc_dir = other.m_slip_bc_dir;
	m_surface_slip_bc_dir = other.m_surface_slip_bc_dir;
    m_forcePlanarSlip = other.m_forcePlanarSlip;
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
	std::swap(m_originalslipsurface, x.m_originalslipsurface);
	std::swap(m_nstep, x.m_nstep);
	std::swap(m_slip_bc_dir, x.m_slip_bc_dir);
	std::swap(m_surface_slip_bc_dir, x.m_surface_slip_bc_dir);
    std::swap(m_forcePlanarSlip, x.m_forcePlanarSlip);
    std::swap(m_AVGslipNormal, x.m_AVGslipNormal);
    std::swap(m_AVGslipCenter, x.m_AVGslipCenter);

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
 * \return true if the class is forcing the slip surface to be treated as plane. See forcePlanarSlip method docs.
 */
bool
PropagateVectorField::isForcingPlanarSlip(){
	return m_forcePlanarSlip;
}

/*!
 * Sets the portion of boundary mesh relative to geometry target
 * that must be constrained component-wise with zero normal field throughout boundary surface.
 * This patch is optional. If nothing is linked, the relative boundary is
 * solved free of any conditions.
 * WARNING: in case of internal holes inside the slip surface, quality of the slip computation
 * might be affected. Slip surface holes are NOT filled automatically right now.
 * \param[in] surface Boundary patch.
 */
void
PropagateVectorField::setSlipBoundarySurface(MimmoObject* surface){
	if (!surface)       return;
	if (surface->getType()!= 1 ) return;
	m_slipsurface = surface;
}

/*!
 * If an almost planar surface is provided as slip surface,
 * this method forces the class to treat the surface as planar, computing an average plane.
 * This trick can manage quasi-planar surface with internal holes, without affecting quality of the
 * slip conditions computations.
 * Results are not guaranteed in case of highly non-planar slip surface.
 * \param[in] planar true, use the average plane of the slip surface to calculate slip conditions, false use the slip surface as it is.
 */
void
PropagateVectorField::forcePlanarSlip(bool planar){
	m_forcePlanarSlip = planar;
}
/*!
 * It sets the Dirichlet conditions for each component of the vector field on the previously linked
 * Dirichlet Boundary patch.
 * \param[in] bc dirichlet conditions
 */
void
PropagateVectorField::setDirichletConditions(dmpvecarr3E bc){
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
 * Clear all data actually stored in the class
 */
void
PropagateVectorField::clear(){
	PropagateField<3>::clear();
	setDefaults();
};

/*!
 * Check coherence of the input data of the class, in particular:
 * - check if boundary patches shares nodes with the bulk geometry
 * - check if data defined on boundary patches are coherent with them
 * if check is good, append the type of BC info on m_isbp.
 * \return true if coherence is satisfied, false otherwise.
 */
bool PropagateVectorField::checkBoundariesCoherence(){

	if (m_method == PropagatorMethod::FINITEVOLUMES){
		//Clean the old m_isbp and initialize it again (set all to Neumann - 0).
		initializeBoundaryInfo();
	}

	// this check the part related to the dirichlet information
	if(!PropagateField<3>::checkBoundariesCoherence()){
		return false;
	}

	//check the part of slip surface condition if present.
	if(m_slipsurface){
		livector1D slipoints = m_slipsurface->getVertices().getIds();

		//create the surface_slip point based dirichlet vector, for internal use only.
		m_surface_slip_bc_dir.reserve(slipoints.size());
		for(long id : slipoints){
			if (!m_surface_slip_bc_dir.exists(id))
				m_surface_slip_bc_dir.insert(id, {{0.0,0.0,0.0}});
			else
				m_surface_slip_bc_dir[id] = std::array<double,3>({{0.0,0.0,0.0}});
		}
		m_surface_slip_bc_dir.setGeometry(m_slipsurface);
		m_surface_slip_bc_dir.setDataLocation(MPVLocation::POINT);

		if (m_method == PropagatorMethod::FINITEVOLUMES){

			livector1D slipInterfaceList = getGeometry()->getInterfaceFromVertexList(slipoints, true, true);

			if(slipInterfaceList.empty()){
				return false;
			}

			//update the m_isbp marking slip interfaces (2)
			for(long id : slipInterfaceList){
				if (m_isbp.at(id) ==1) continue; //keep the hands off from dirichlet.
				m_isbp.at(id) = 2;
			}
		}
	}

	return true;
}

/*!
 * Instantiate m_originalslipsurface, that is the surface where you need to evaluate
 * slip conditions. In SERIAL version m_originalslipsurface is the m_slipsurface passed as boundary.
 * In MPI version m_originalslipsurface is the whole slip surface rebuilt from every portion "m_slipsurface"
 * owned by each rank. Once recollected the whole surface is sent as it is to all ranks.
 * This surface needs to remain untouched and undeformed during Multistep iterations.
 */
void PropagateVectorField::initializeSlipSurface(){

#if MIMMO_ENABLE_MPI

	if (getTotalProcs() > 1){

		//MPI version
		m_originalslipsurface = std::move(std::unique_ptr<MimmoObject>(new MimmoObject(1)));

		//fill serialized geometry
		for (bitpit::Vertex &vertex : m_slipsurface->getVertices()){
			long vertexId = vertex.getId();
			if (m_slipsurface->isPointInterior(vertexId)){
				m_originalslipsurface->addVertex(vertex, vertexId);
			}
		}
		for (bitpit::Cell &cell : m_slipsurface->getCells()){
			if (cell.isInterior()){
				m_originalslipsurface->addCell(cell, cell.getId());
			}
		}


		//Receive vertices and cells
		for (int sendRank=0; sendRank<m_nprocs; sendRank++){

			if (m_rank != sendRank){

				// Vertex data
				long vertexBufferSize;
				MPI_Recv(&vertexBufferSize, 1, MPI_LONG, sendRank, 100, m_communicator, MPI_STATUS_IGNORE);

				bitpit::IBinaryStream vertexBuffer(vertexBufferSize);
				MPI_Recv(vertexBuffer.data(), vertexBuffer.getSize(), MPI_CHAR, sendRank, 110, m_communicator, MPI_STATUS_IGNORE);

				// Cell data
				long cellBufferSize;
				MPI_Recv(&cellBufferSize, 1, MPI_LONG, sendRank, 200, m_communicator, MPI_STATUS_IGNORE);

				bitpit::IBinaryStream cellBuffer(cellBufferSize);
				MPI_Recv(cellBuffer.data(), cellBuffer.getSize(), MPI_CHAR, sendRank, 210, m_communicator, MPI_STATUS_IGNORE);

				// There are no duplicate in the received vertices, but some of them may
				// be already a local vertex of a interface cell.
				//TODO GENERALIZE IT
				//NOTE! THE COHINCIDENT VERTICES ARE SUPPOSED TO HAVE THE SAME ID!!!!
				long nRecvVertices;
				vertexBuffer >> nRecvVertices;
				m_originalslipsurface->getPatch()->reserveVertices(m_originalslipsurface->getPatch()->getVertexCount() + nRecvVertices);

				// Do not add the vertices with Id already in serialized geometry
				for (long i = 0; i < nRecvVertices; ++i) {
					bitpit::Vertex vertex;
					vertexBuffer >> vertex;
					long vertexId = vertex.getId();

					if (!m_originalslipsurface->getVertices().exists(vertexId)){
						m_originalslipsurface->addVertex(vertex, vertexId);
					}
				}

				//Receive and add all Cells
				long nReceivedCells;
				cellBuffer >> nReceivedCells;
				m_originalslipsurface->getPatch()->reserveCells(m_originalslipsurface->getPatch()->getCellCount() + nReceivedCells);

				for (long i = 0; i < nReceivedCells; ++i) {
					// Cell data
					bitpit::Cell cell;
					cellBuffer >> cell;

					long cellId = cell.getId();

					// Add cell
					m_originalslipsurface->addCell(cell, cellId);

				}
			}
			else{

				//Send local vertices and local cells to ranks

				//
				// Send vertex data
				//
				bitpit::OBinaryStream vertexBuffer;
				long vertexBufferSize = 0;
				long nVerticesToCommunicate = 0;

				// Fill buffer with vertex data
				vertexBufferSize += sizeof(long);
				for (long vertexId : m_slipsurface->getVertices().getIds()){
					if (m_slipsurface->isPointInterior(vertexId)){
						vertexBufferSize += m_slipsurface->getVertices()[vertexId].getBinarySize();
						nVerticesToCommunicate++;
					}
				}
				vertexBuffer.setSize(vertexBufferSize);

				vertexBuffer << nVerticesToCommunicate;
				for (long vertexId : m_slipsurface->getVertices().getIds()){
					if (m_slipsurface->isPointInterior(vertexId)){
						vertexBuffer << m_slipsurface->getVertices()[vertexId];
					}
				}

				for (int recvRank=0; recvRank<m_nprocs; recvRank++){

					if (m_rank != recvRank){

						// Communication
						MPI_Send(&vertexBufferSize, 1, MPI_LONG, recvRank, 100, m_communicator);
						MPI_Send(vertexBuffer.data(), vertexBuffer.getSize(), MPI_CHAR, recvRank, 110, m_communicator);

					}
				}

				//
				// Send cell data
				//
				bitpit::OBinaryStream cellBuffer;
				long cellBufferSize = 0;
				long nCellsToCommunicate = 0;

				// Fill the buffer with cell data
				cellBufferSize += sizeof(long);
				for (const long cellId : m_slipsurface->getCellsIds()) {
					if (m_slipsurface->getCells()[cellId].isInterior()){
						cellBufferSize += sizeof(int) + sizeof(int) + m_slipsurface->getCells()[cellId].getBinarySize();
						nCellsToCommunicate++;
					}
				}
				cellBuffer.setSize(cellBufferSize);

				cellBuffer << nCellsToCommunicate;
				for (const long cellId : m_slipsurface->getCellsIds()) {
					if (m_slipsurface->getCells()[cellId].isInterior()){
						const bitpit::Cell &cell = m_slipsurface->getCells()[cellId];
						// Cell data
						cellBuffer << cell;
					}
				}

				for (int recvRank=0; recvRank<m_nprocs; recvRank++){

					if (m_rank != recvRank){
						// Communication
						MPI_Send(&cellBufferSize, 1, MPI_LONG, recvRank, 200, m_communicator);
						MPI_Send(cellBuffer.data(), cellBuffer.getSize(), MPI_CHAR, recvRank, 210, m_communicator);

					}
				}

			}

		}// end external sendrank loop
	}
	else
#endif
	{
		//serial version
		m_originalslipsurface = m_slipsurface->clone();
	}
}

/*!
 * Instead of instantiating m_originalslipsurface, calculate the average normal and
 * barycenter of m_slipsurface, and store them in internal members m_AVGslipNormal and m_AVGslipCenter.
 */
void PropagateVectorField::initializeSlipSurfaceAsPlane(){

    m_AVGslipCenter.fill(0.0);
    m_AVGslipNormal.fill(0.0);
    long countV = 0;
    long countC = 0;

#if MIMMO_ENABLE_MPI
    if(!m_slipsurface->arePointGhostExchangeInfoSync()) m_slipsurface->updatePointGhostExchangeInfo();
#endif

bitpit::SurfaceKernel * surfkernss = dynamic_cast<bitpit::SurfaceKernel*>(m_slipsurface->getPatch());
    //start evaluating barycenter of interior points.
    for(auto itV = surfkernss->vertexBegin(); itV != surfkernss->vertexEnd(); ++itV){
#if MIMMO_ENABLE_MPI
        if(m_slipsurface->isPointInterior(itV.getId()))
#endif
        {
            m_AVGslipCenter += itV->getCoords();
            ++countV;
        }
    }

    //evaluating average facet normal of interior cells.
    for(auto itC = surfkernss->internalBegin(); itC != surfkernss->internalEnd(); ++itC){
        m_AVGslipNormal += surfkernss->evalFacetNormal(itC.getId());
        ++countC;
    }


#if MIMMO_ENABLE_MPI
    //communicate with other processors summing barycenters, normals
    //number of vertices and number of cells and then average it.
    MPI_Allreduce(MPI_IN_PLACE, &countV, 1, MPI_LONG, MPI_SUM, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, &countC, 1, MPI_LONG, MPI_SUM, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, m_AVGslipCenter.data(), 3, MPI_DOUBLE, MPI_SUM, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, m_AVGslipNormal.data(), 3, MPI_DOUBLE, MPI_SUM, m_communicator);

#endif

    if(countV > 0)  m_AVGslipCenter /= double(countV);
    if(countC > 0)  m_AVGslipNormal /= double(countC);

}


/*!
 * If dumping is active deform with value of m_field the internal dumping geometry
 * m_originalDumpingSurface, which is a reconstruction of all dumping surfaces throughout the procs
 * All procs has the same m_originalDumpingSurface.
 * This method ask to all procs to recover the part of deformation field m_field referred to nodes
 * of m_originaDumpingSurface. Procs Mutually exchange their own parts to have all of them a unique deformation field for
 * their own m_originalDumpingSurface. Then the deformation is applied.
 */
void PropagateVectorField::deformDumpingSurface(){


std::unique_ptr<MimmoPiercedVector<std::array<double,3> > >
        mpvres(new MimmoPiercedVector<std::array<double,3> >(m_originalDumpingSurface.get(), MPVLocation::POINT) );

    mpvres->reserve(m_originalDumpingSurface->getNVertices());

    MimmoObject * dumptarget = m_dsurface;
    if(!dumptarget) dumptarget = m_bsurface;

    for(long id: dumptarget->getVertices().getIds()){
#if MIMMO_ENABLE_MPI
        if(dumptarget->isPointInterior(id))
#endif
        {
            mpvres->insert(id, m_field.at(id));
        }
    }

#if MIMMO_ENABLE_MPI
	//Send my own and receive deformations by other.

    //prepare my own send buffer.
    bitpit::OBinaryStream myrankDataBuffer;
    myrankDataBuffer << *(mpvres.get());
    long myrankDataBufferSize = myrankDataBuffer.getSize();

    for (int sendRank=0; sendRank<m_nprocs; sendRank++){

		if (m_rank != sendRank){
            // receive data from other ranks.
			long defBufferSize;
			MPI_Recv(&defBufferSize, 1, MPI_LONG, sendRank, 900, m_communicator, MPI_STATUS_IGNORE);
			bitpit::IBinaryStream defBuffer(defBufferSize);
			MPI_Recv(defBuffer.data(), defBuffer.getSize(), MPI_CHAR, sendRank, 910, m_communicator, MPI_STATUS_IGNORE);

            MimmoPiercedVector<std::array<double,3>> temp;
            defBuffer >> temp;

			// insert this part in mpvres.
			for (auto it = temp.begin(); it!=temp.end(); ++it) {
                if(!mpvres->exists(it.getId())){
                    mpvres->insert(it.getId(), *it);
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

    // you can deform
    bitpit::PiercedVector<bitpit::Vertex> & verts = m_originalDumpingSurface->getVertices();
	for (auto it= verts.begin(); it != verts.end(); ++it){
		m_originalDumpingSurface->modifyVertex(it->getCoords() + mpvres->at(it.getId()), it.getId());
	}


}

/*!
 * Distribute the input m_surface_slip_bc_dir on the border interfaces of the bulk mesh.
 */
void PropagateVectorField::distributeSlipBCOnBoundaryInterfaces(){

	if(!m_slipsurface)  {
		return;
	}
	//transfer slip point field info of boundary "dirichlet" on the volume mesh interfaces.
	MimmoPiercedVector<std::array<double,3>> temp(m_geometry, MPVLocation::POINT);
	temp.reserve(m_surface_slip_bc_dir.size());
	for(auto it=m_surface_slip_bc_dir.begin(); it!=m_surface_slip_bc_dir.end(); ++it){
		temp.insert(it.getId(), *it );
	}

	//interpolate now point data to interface data
	m_slip_bc_dir.clear();
	m_slip_bc_dir = temp.pointDataToBoundaryInterfaceData();

}

/*!
 * Plot optional results on vtu unstructured grid file
 */
void
PropagateVectorField::plotOptionalResults(){

	if(!getGeometry())    return;

	bitpit::VTKUnstructuredGrid& vtk = getGeometry()->getPatch()->getVTK();
#if MIMMO_ENABLE_MPI
	getGeometry()->getPatch()->setVTKWriteTarget(bitpit::PatchKernel::WriteTarget::WRITE_TARGET_CELLS_INTERNAL);
#endif

	dvecarr3E data = m_field.getInternalDataAsVector();
	vtk.addData("field", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, data);

	dvector1D datad = m_dumping.getInternalDataAsVector();
	vtk.addData("dumping", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, datad);

	vtk.setCounter(getId());
    getGeometry()->getPatch()->getVTK().setDirectory(m_outputPlot+"/");
    getGeometry()->getPatch()->getVTK().setName(m_name+"_field");
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
	MimmoObject * target = getGeometry();
	if (!target) return;

	bitpit::PiercedVector<bitpit::Vertex> & verts = target->getVertices();

	for (auto it= verts.begin(); it != verts.end(); ++it){
		target->modifyVertex(it->getCoords()+m_field.at(it.getId()), it.getId());
	}

	//Modify m boundary surface
	MimmoObject * bsurf = m_bsurface;
	bitpit::PiercedVector<bitpit::Vertex> & vertbound = bsurf->getVertices();
	for (auto it= vertbound.begin(); it != vertbound.end(); ++it){
		bsurf->modifyVertex(it->getCoords()+m_field.at(it.getId()), it.getId());
	}

	//if the dumping is active morph the internal member m_originaDumpingSurface. It is enough.
	if(m_dumpingActive && m_originalDumpingSurface.get() != nullptr){
	       deformDumpingSurface();
	}

	//if the slipsurface is active apply this deformation also to the slipsurface(re-evaluation of normals)
	if(m_slipsurface){
		bitpit::PiercedVector<bitpit::Vertex> & vertslip = m_slipsurface->getVertices();
		for (auto it= vertslip.begin(); it != vertslip.end(); ++it){
			m_slipsurface->modifyVertex(it->getCoords()+m_field.at(it.getId()), it.getId());
		}
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
	for (auto & val : m_surface_bc_dir){
		val /= double(m_nstep);
	}
}

/*!
 * restore dirichlet boundary conditions for multi step purposes
 */
void
PropagateVectorField::restoreBC(){
	for (auto & val : m_surface_bc_dir){
		val *= double(m_nstep);
	}

}

/*!
 * restore geometry to target vertices and reevaluate m_field as whole
 * \param[in] list of vertices to be restored
 */
void
PropagateVectorField::restoreGeometry(bitpit::PiercedVector<bitpit::Vertex> & vertices){
	MimmoObject * target = getGeometry();
	bitpit::PiercedVector<bitpit::Vertex> &currentmesh = target->getVertices();
	long ID;
	for (auto it= vertices.begin(); it!=vertices.end(); ++it){
		ID = it.getId();
		const std::array<double,3> &coords = it->getCoords();
		m_field.at(ID) = currentmesh.at(ID).getCoords() - coords;
		target->modifyVertex(coords, ID);
	}

	//Restore m boundary surface
	MimmoObject * bsurf = m_bsurface;
	bitpit::PiercedVector<bitpit::Vertex> & vertbound = bsurf->getVertices();
	for (auto it= vertbound.begin(); it != vertbound.end(); ++it){
		*it = vertices.at(it.getId());
	}

    //m_dumping partitioned members are not deformed. m_originalDumpingSurface, their
    // reconstruction, is deformed. By I have no interest in restoring it,
    // since is not a User feeded input.

	//restore the slip surface
	if(m_slipsurface){
		bitpit::PiercedVector<bitpit::Vertex> &vertslip = m_slipsurface->getVertices();
		for(auto it = vertslip.begin(); it != vertslip.end(); ++it){
			*it = vertices.at(it.getId());
		}
	}
}

/*!
 * Given an ensemble of cell marked as moving, enrich this list adding the first ring of
 * their vertex neighbours V1Ring, and the second ring of face neighbours F2Ring.
 *
 * \param[in,out] pool of moving cells in input, V1Ring+F2Ring augmented in output.
 */
void PropagateVectorField::propagateMaskMovingCells(livector1D & cellList) {

	std::unordered_set<long> core(cellList.begin(), cellList.end());
	auto * patch = getGeometry()->getPatch();

	livector1D tempV1, tempF2;
	for(long id: cellList){
		//get the V1Ring.
		tempV1 = patch->findCellVertexNeighs(id, true);

		//run the list element by element. Any new element, evaluate its face neighs in F2Ring and push it in core.
		for(long candidate : tempV1){
			if(core.count(candidate) > 0 ) continue;

			core.insert(candidate);
			tempF2 = patch->findCellFaceNeighs(candidate);
			core.insert(tempF2.begin(), tempF2.end());
		}
	}

	cellList.clear();
	cellList.insert(cellList.end(), core.begin(), core.end());
}

/*!
 * Given an ensemble of points marked as moving, enrich this list adding the first neighbours of
 * these vertices.
 *
 * \param[in,out] pool of moving points in input, 1Ring augmented in output.
 */
void PropagateVectorField::propagateMaskMovingPoints(livector1D & vertexList) {

	std::unordered_set<long> core(vertexList.begin(), vertexList.end());
	auto * patch = getGeometry()->getPatch();

    if (!getGeometry()->isPointConnectivitySync())
    	getGeometry()->buildPointConnectivity();

	std::unordered_set<long> tempV1;
	for(long id: vertexList){
		tempV1 = getGeometry()->getPointConnectivity(id);

		//run the list element by element. Any new element, evaluate its face neighs in F2Ring and push it in core.
		for(long candidate : tempV1){
			if(core.count(candidate) > 0 ) continue;
			core.insert(candidate);
		}
	}

	vertexList.clear();
	vertexList.insert(vertexList.end(), core.begin(), core.end());
}

/*!
 * OVERRIDE Base class: This method evaluate the bc corrections for a singular run of the system solver,
 * update the system matrix in m_solver and evaluate the rhs part due to bc.
 * After you call this method, you are typically ready to solve the laplacian system.
 * The type of bc @ interface are directly desumed from class map member m_isbp.
 * The method requires the Laplacian m_solver to be initialized. No ghost are taken into account.
 *
 * It manages also the slip condition using a predictor-corrector scheme:
 * First create a set of slip Neumann homogeneous bc condition to predict a guess solution.
 * Then when the correct bc is defined in m_slip_bc_dir from guess solution, create a set
 * of slip Dirichlet bc condition, reading it form m_slip_bc_dir.
 *
 * Moreover it needs as input :
 * \param[in] comp target component of bc conditions.
 * \param[in] slipCorrect true to apply the bc slip corrector step(Dirichlet), false for bc slip predictor(Neumann).
 * \param[in] borderLaplacianStencil list of laplacian Stencil on border cells, where the bc is temporarely imposed as homogeneous Neumann
 * \param[in] borderCCGradientStencil list of Center cell gradient stencils defined on border cells.
 * \param[in] maplocals map from global id numbering to local system solver numbering.
 * \param[in,out] rhs, vector of right-hand-side's to append constant data from bc corrections.
 */
void
PropagateVectorField::assignBCAndEvaluateRHS(std::size_t comp, bool slipCorrect,
		FVolStencil::MPVDivergence * borderLaplacianStencil,
		FVolStencil::MPVGradient * borderCCGradientStencil,
		const liimap & maplocals,
		dvector1D & rhs)
{
	MimmoObject * geo = getGeometry();

	//resize rhs to the number of internal cells
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
	FVolStencil::MPVDivergenceUPtr lapwork(new FVolStencil::MPVDivergence(*borderLaplacianStencil));

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
		case 2: //SLIP/IMPERMEABILITY

			if(!slipCorrect){ //predictor stage - push neumann condition.
				correction = FVolStencil::correctionNeumannBCFaceGradient(0.0, interfaceNormal);
			}else{ //corrector stage - push dirichlet condition from m_slip_bc_dir you have fixed.
				ownerCentroid = geo->evalCellCentroid(idOwner);
				interfaceCentroid = geo->evalInterfaceCentroid(idInterface);
				//get the correction.
				correction = FVolStencil::correctionDirichletBCFaceGradient(m_slip_bc_dir.at(idInterface)[comp],
						idOwner,ownerCentroid, interfaceCentroid, interfaceNormal,
						dotProduct(interfaceCentroid-ownerCentroid, interfaceNormal),
						borderCCGradientStencil->at(idOwner) );
			}
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
 * OVERRIDE Base class: This method evaluate the bc corrections for a singular run of the system solver,
 * update the system matrix in m_solver and evaluate the rhs part due to bc.
 * After you call this method, you are typically ready to solve the laplacian system.
 * The type of bc @ nodes are directly desumed from class nodes member m_bc_dir and m_surface_slip_bc_dir.
 * The method requires the Laplacian m_solver to be initialized. No ghost are taken into account.
 *
 * It manages also the slip condition using a predictor-corrector scheme:
 * First create a set of slip Neumann homogeneous bc condition to predict a guess solution.
 * Then when the correct bc is defined in m_slip_bc_dir from guess solution, create a set
 * of slip Dirichlet bc condition, reading it form m_slip_bc_dir.
 *
 * Moreover it needs as input :
 * \param[in] comp target component of bc conditions.
 * \param[in] slipCorrect true to apply the bc slip corrector step(Dirichlet), false for bc slip predictor(Neumann).
 * \param[in] borderLaplacianStencil list of laplacian Stencil on border nodes, where the bc is temporarely imposed as homogeneous Neumann
 * \param[in] maplocals map from global id numbering to local system solver numbering.
 * \param[in,out] rhs, vector of right-hand-side's to append constant data from bc corrections.
 */
void
PropagateVectorField::assignBCAndEvaluateRHS(std::size_t comp, bool slipCorrect,
		GraphLaplStencil::MPVStencil * borderLaplacianStencil,
		const liimap & maplocals,
		dvector1D & rhs)
{
	MimmoObject * geo = getGeometry();

	//resize rhs to the number of internal cells
	rhs.resize(geo->getNInternalVertices(), 0.0);

	if (!m_solver->isInitialized()) {
		(*m_log)<<"Warning in "<<m_name<<". Unable to assign BC to the system. The solver is not yet initialized."<<std::endl;
		return;
	}

	if (!borderLaplacianStencil) {
		(*m_log)<<"Warning in "<<m_name<<". Unable to reach border cells stencils data. Nothing to do."<<std::endl;
		return;
	}

	//copy laplacian stencils in a work mpv .
	GraphLaplStencil::MPVStencilUPtr lapwork(new GraphLaplStencil::MPVStencil(*borderLaplacianStencil));

	//correct the original border laplacian stencils applying the Dirichlet conditions and slip conditions.
	//Nuemann are implicitely imposed by graph-laplacian scheme.
	//renumber it and update the laplacian matrix and fill the rhs.
	bitpit::StencilScalar correction;
	bitpit::PiercedVector<bitpit::Vertex> & mesh_vertices = geo->getVertices();

	//loop on all dirichlet boundary nodes.
	for(long id : m_bc_dir.getIds()){
#if MIMMO_ENABLE_MPI
		if (geo->isPointInterior(id))
#endif
		{
			//apply the correction relative to bc @ dirichlet node.
			correction.clear(true);
			correction.appendItem(id, 1.);
			correction.sumConstant(-m_bc_dir[id][comp]);
			//Fix to zero the old stencil (the update of system solver doesn't substitute but modify or append new pattern item and weights)
			lapwork->at(id) *= 0.;
			lapwork->at(id) += correction;
		}
	}


	//loop on all slip boundary nodes.
	//Correct if it is the correction step. If not the neumann condition are automatically imposed by Graph-Laplace scheme.
	if(slipCorrect){
		for(long id : m_surface_slip_bc_dir.getIds()){
#if MIMMO_ENABLE_MPI
			if (geo->isPointInterior(id))
#endif
			{
				//apply the correction relative to bc @ dirichlet node.
				correction.clear(true);
				correction.appendItem(id, 1.);
				correction.sumConstant(-m_surface_slip_bc_dir[id][comp]);
				//Fix to zero the old stencil (the update of system solver doesn't substitute but modify or append new pattern item and weights)
				lapwork->at(id) *= 0.;
				lapwork->at(id).setConstant(0.);
				lapwork->at(id) += correction;
			}
		}
	}

	// now its time to update the solver matrix and to extract the rhs contributes.
	//NOTE: USE THE FINITE VOLUMES METHOD, THE STRUCTURES ARE THE SAME!
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
 * Starting from a guess solution on mesh points, find a set of dirichlet condition on slip
 * boundary so that the dot product of vector solution at slip wall and slip wall local normal is zero.
 * The temporary dirichlet bc will be stored in m_slip_bc_dir and can be used in assignBCAndEvaluateRHS
 * forcing the boolean variable slipCorrect to true.
 *
 * \param[in] guessSolutionOnPoint guess laplacian solution of the vector field on mesh POINTS
 */
void
PropagateVectorField::computeSlipBCCorrector(const MimmoPiercedVector<std::array<double,3> > & guessSolutionOnPoint)
{
	if(!m_slipsurface) return;
	bitpit::SurfaceKernel * slipsurf = dynamic_cast<bitpit::SurfaceKernel *> (m_slipsurface->getPatch());

//	//If not initialized initialize original slip surface
//	if (!m_originalslipsurface)
//		m_originalslipsurface = m_slipsurface->clone();


	//first step: extract solutions on twin border nodes of m_surface_slip_bc_dir;
	long id;
	for(auto it=m_surface_slip_bc_dir.begin(); it!=m_surface_slip_bc_dir.end(); ++it){
		id = it.getId();
		if(guessSolutionOnPoint.exists(id))
			*it = guessSolutionOnPoint.at(id);
	}

//	//precalculate vertex normals on m_slipsurface;
//	std::unordered_map<long, std::array<double,3> > mapNormals;
//	//loop on surface cells.
//	for( const bitpit::Cell & cell: m_slipsurface->getCells()){
//		bitpit::ConstProxyVector<long> vertexList = cell.getVertexIds();
//		std::size_t local = 0;
//		for(long idV: vertexList){
//			if(mapNormals.count(idV) < 1){
//				mapNormals[idV] = slipsurf->evalVertexNormal(cell.getId(), local);
//			}
//			++local;
//		}
//	}
//
//	// force the correction on m_surface_slip_bc_dir
//	std::array<int, 3> pord;
//	for(auto it=m_surface_slip_bc_dir.begin(); it!=m_surface_slip_bc_dir.end(); ++it){
//		id = it.getId();
//		std::array<double,3> & normal = mapNormals.at(id);
//		std::array<double,3> & bcval  = *it;
//
//		bcval -= dotProduct(bcval, normal) * normal; //remove its normal component.
//	}

    //precalculate projection of new point on m_slipsurface;
    std::unordered_map<long, std::array<double,3> > projectionVector;

    if(m_forcePlanarSlip){
    	  // VERSION USING THE AVERAGE NORMAL!
    	//loop on surface points.
    	for(auto it=m_surface_slip_bc_dir.begin(); it!=m_surface_slip_bc_dir.end(); ++it){
    		long idV = it.getId();
    		std::array<double,3> point = m_slipsurface->getVertexCoords(idV) + m_surface_slip_bc_dir[idV];
    		std::array<double,3> projpoint = bitpit::CGElem::projectPointPlane(point, m_AVGslipCenter,m_AVGslipNormal);
    		projectionVector[idV] = projpoint - point;
    	}
    }else{
	//precalculate projection of new point on m_slipsurface;
    	if (!m_originalslipsurface->isSkdTreeSync())
    		m_originalslipsurface->buildSkdTree();

    	bitpit::PatchSkdTree *tree = m_originalslipsurface->getSkdTree();
    	//loop on surface points.
    	for(auto it=m_surface_slip_bc_dir.begin(); it!=m_surface_slip_bc_dir.end(); ++it){
    		long idV = it.getId();
    		bitpit::Vertex &vertex = m_slipsurface->getPatch()->getVertex(idV);
    		std::array<double,3> point = vertex.getCoords() + m_surface_slip_bc_dir[idV];
    		 double r = std::max(1.0e-02, norm2(m_surface_slip_bc_dir[idV])*1.25);
    		std::array<double,3> projpoint = skdTreeUtils::projectPoint(&point, tree, r);
    		projectionVector[idV] = projpoint - point;
    	}
    }



	// force the correction on m_surface_slip_bc_dir
	for(auto it=m_surface_slip_bc_dir.begin(); it!=m_surface_slip_bc_dir.end(); ++it){
		id = it.getId();
        (*it) += projectionVector.at(id);
	}

//
//	//TODO DEBUG
//	{
//		bitpit::VTKUnstructuredGrid& vtk = m_slipsurface->getPatch()->getVTK();
//#if MIMMO_ENABLE_MPI
//		m_slipsurface->getPatch()->setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_ALL);
//#endif
//
//		dvecarr3E data(m_slipsurface->getPatch()->getInternalCount());
//		data = m_surface_slip_bc_dir.getDataAsVector();
//		vtk.addData("slip-bc", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, data);
//
//		vtk.setCounter(vtk.getCounter()+1);
//		m_slipsurface->getPatch()->write(m_name);
//		vtk.removeData("slip-bc");
//	}


	// done.
}

/*!
 * Execution command. After the execution the result constraint field is stored in the class.
 */
void
PropagateVectorField::execute(){

	MimmoObject * geo = getGeometry();
	if(!geo){
		(*m_log)<<"Warning in "<<m_name<<" .No target volume mesh linked"<<std::endl;
	}

	if(!m_bsurface){
		(*m_log)<<"Warning in "<<m_name<<" .No Dirichlet Boundary patch linked"<<std::endl;
	}

	if(!checkBoundariesCoherence()){
		(*m_log)<<"Warning in "<<m_name<<" .Boundary patches linked are uncoherent with target bulk geometry"
				"or bc-fields not coherent with boundary patches"<<std::endl;
	}

	if(m_slipsurface){
        if(m_forcePlanarSlip){
            initializeSlipSurfaceAsPlane();
        }else{
            initializeSlipSurface();
        }
	}

	(*m_log) << bitpit::log::priority(bitpit::log::NORMAL);
	(*m_log) << bitpit::log::context("mimmo");
	//INITIALIZATION --->////////////////////////////////////////////////////////////////////////////////////
	//allocate the solver;
	m_solver = std::unique_ptr<bitpit::SystemSolver>(new bitpit::SystemSolver(m_print));

	//get the inverse and the direct map -> you will need it to compact the stencil/ and recover the results respectively.
	liimap dataInv;
	liimap data;

	//compute the dumping.
    initializeDumpingSurface();
    computeDumpingFunction();

	//PREPARE THE MULTISTEP;
	bitpit::PiercedVector<bitpit::Vertex> undeformedTargetVertices;
	std::unique_ptr<livector1D> movingElementList = nullptr;
	if(m_nstep > 1){
		subdivideBC();
		undeformedTargetVertices = geo->getVertices();
		movingElementList = std::unique_ptr<livector1D>(new livector1D());
	}

	//Switch on solver method
	if (m_method == PropagatorMethod::FINITEVOLUMES){

		// Finite Volumes method

		//store the id of the border cells only;
		livector1D borderCellsID = geo->extractBoundaryCellID(false);

		//get the inverse and the direct map -> you will need it to compact the stencil/ and recover the results respectively.
		dataInv = geo->getMapCellInv(true);
		data = geo->getMapCell(true);

		//compute all the center cell gradients.
		FVolStencil::MPVGradientUPtr ccellGradients = FVolStencil::computeFVCellGradientStencil(*geo);

		//compute the gradient stencils @ interface with homogeneous Neumann.
		FVolStencil::MPVGradientUPtr faceGradients  = FVolStencil::computeFVFaceGradientStencil(*geo, ccellGradients.get());

		//and squeeze out cell gradients and save the border cells only.
		ccellGradients->squeezeOutExcept(borderCellsID);

		// compute the laplacian stencils and free faceGradients;
		FVolStencil::MPVDivergenceUPtr laplaceStencils = FVolStencil::computeFVLaplacianStencil(*(faceGradients.get()), m_tol, &m_dumping);
		faceGradients  = nullptr;

		// initialize the laplacian Matrix in solver and squeeze out the laplace stencils and save border cells only.
		initializeLaplaceSolver(laplaceStencils.get(), dataInv);
		laplaceStencils->squeezeOutExcept(borderCellsID);
		borderCellsID.clear();

		//declare results here and keep it during the loop to re-use the older steps.
		std::vector<std::vector<double>> results(3);

		//loop on multistep
		for(int istep=0; istep < m_nstep; ++istep){

			distributeBCOnBoundaryInterfaces(); //distribute the original dirichlet condition on bulk interfaces.
			//PLEASE NOTE, this redistribution is essential, since the shape of the boundary changes during
			// multistep, and the interpolation weights from point to interfaces changes accordingly.

			//3-COMPONENT SYSTEM SOLVING ---> ///////////////////////////////////////////////////////////////////////
			// solve the field component by component
			//first loop -> if slip is enforced in some walls, this is the predictor stage of guess solution with 0-Neumann on slip walls
			for(int comp = 0; comp<3; ++comp){
				//prepare the right hand side ;
				dvector1D rhs(geo->getPatch()->getInternalCount(), 0.0);
				assignBCAndEvaluateRHS(comp, false, laplaceStencils.get(), ccellGradients.get(), dataInv, rhs);
				//solve
				results[comp].resize(rhs.size(), 0.0);
				solveLaplace(rhs, results[comp]);
			}

			//if I have a slip wall active, it needs a corrector stage for slip boundaries;
			if(m_slipsurface){
				//reconstruct result on mesh points (ghost included)-> stored in m_field.
				reconstructResults(results, data);
				//correct the BC of the slip walls
				computeSlipBCCorrector(m_field);
				// distribute this slip dirichlet on the bulk interfaces
				distributeSlipBCOnBoundaryInterfaces();
				//now you have a set of dirichlet BC in m_slip_bc_dir pierced vector internal.

				// so loop again on the components, reusing the previous result as starting guess, and setting
				// the boolean of assignBC as true (corrector stage of slip, read Dirichlet from m_slip_bc_dir)
				for(int comp = 0; comp<3; ++comp){
					//prepare the right hand side ;
					dvector1D rhs(geo->getPatch()->getInternalCount(), 0.0);
					assignBCAndEvaluateRHS(comp, true, laplaceStencils.get(), ccellGradients.get(), dataInv, rhs);
					//solve
					solveLaplace(rhs, results[comp]);
				}
			}
			//RECONSTRUCT STAGE --> /////////////////////////////////////////////////////////////////////////////////
			reconstructResults(results, data, movingElementList.get());

			//force boundary on slip
			if(m_slipsurface){
				//if slip is active, force SLIP CORRECTION on POINTS (stored in m_surface_slip_bc_dir) on m_field.
				for(auto it=m_surface_slip_bc_dir.begin(); it != m_surface_slip_bc_dir.end(); ++it){
					m_field.at(it.getId()) = *it;
				}
			}

			if (m_forceDirichletConditions){
				//at last force boundary Dirichlet on POINTS on m_field; Because Dirichlet RULEZ
				for(auto it=m_surface_bc_dir.begin(); it != m_surface_bc_dir.end(); ++it){
					m_field.at(it.getId()) = *it;
				}
			}

#if MIMMO_ENABLE_MPI
			communicatePointGhostData(&m_field);
#endif

			// if you are in multistep stage apply the deformation to the mesh.
            // here is morphed also the m_originalDumpingSurface if active dumping.
			if(m_nstep > 1){
				apply();
			}

			// if in multistep stage continue to update the other laplacian stuff up to "second-to-last" step.
			if(istep < m_nstep-1){
				//update the dumping function.
//				updateDumpingFunction();
				computeDumpingFunction();

				//enlarge the moving cell list taking its first vertex neighs and its second face neighs.
				propagateMaskMovingCells(*(movingElementList.get()));
				//update the center cell gradients and clear moving cells..
				FVolStencil::MPVGradientUPtr updateCcellGradients = FVolStencil::computeFVCellGradientStencil(*geo, movingElementList.get());
				movingElementList->clear();
				//and update the border ccellgradients with this new updated values
				ccellGradients->getDataFrom(*(updateCcellGradients.get()), true); //only common elements are updated
				//compute the gradient stencils @ interface with homogeneous Neumann.
				FVolStencil::MPVGradientUPtr updateFaceGradients  = FVolStencil::updateFVFaceGradientStencil(*geo, *(updateCcellGradients.get()));
				//we can destroy the updateCcellGradients
				updateCcellGradients = nullptr;

				// compute the update laplacian stencils, free faceGradients, update the old border laplacian stencils;
				FVolStencil::MPVDivergenceUPtr updateLaplaceStencils = FVolStencil::computeFVLaplacianStencil(*(updateFaceGradients.get()), m_tol, &m_dumping);
				faceGradients  = nullptr;
				laplaceStencils->getDataFrom(*(updateLaplaceStencils.get()), true); //only common elements are updated.

				// update the laplacian Matrix in solver free the updateLaplaceStencils
				updateLaplaceSolver(updateLaplaceStencils.get(), dataInv);
				updateLaplaceStencils = nullptr;

			}
		    (*m_log)<<m_name<<" solved step "<<istep<<" out of total steps "<<m_nstep<<std::endl;
	} //end of multistep loop;

	}
	else if (m_method == PropagatorMethod::GRAPHLAPLACE){

		// Graph Laplace method on points

		//store the id of the border nodes only;
		livector1D borderPointsID = geo->extractBoundaryVertexID(false);

		//get this inverse map -> you will need it to compact the stencils.
		dataInv = geo->getMapDataInv(true);
		data = geo->getMapData(true);

		// compute the laplacian stencils
		GraphLaplStencil::MPVStencilUPtr laplaceStencils = GraphLaplStencil::computeLaplacianStencils(*geo, m_tol, &m_dumping);

		// initialize the laplacian Matrix in solver and squeeze out the laplace stencils and save border cells only.
		initializeLaplaceSolver(laplaceStencils.get(), dataInv);
		laplaceStencils->squeezeOutExcept(borderPointsID);
		borderPointsID.clear();

		//declare results here and keep it during the loop to re-use the older steps.
		std::vector<std::vector<double>> results(3);

		//loop on multistep
		for(int istep=0; istep < m_nstep; ++istep){

			//Useless?
			distributeBCOnBoundaryPoints(); //distribute the original dirichlet condition on bulk interfaces.
			//PLEASE NOTE, this redistribution is essential, since the shape of the boundary changes during
			// multistep, and the interpolation weights from point to interfaces changes accordingly.

			//3-COMPONENT SYSTEM SOLVING ---> ///////////////////////////////////////////////////////////////////////
			// solve the field component by component
			//first loop -> if slip is enforced in some walls, this is the predictor stage of guess solution with 0-Neumann on slip walls
			for(int comp = 0; comp<3; ++comp){
				//prepare the right hand side ;
				dvector1D rhs(geo->getNInternalVertices(), 0.0);
				assignBCAndEvaluateRHS(comp, false, laplaceStencils.get(), dataInv, rhs);

				//solve
				results[comp].resize(rhs.size(), 0.0);
				solveLaplace(rhs, results[comp]);

			}

			//if I have a slip wall active, it needs a corrector stage for slip boundaries;
			if(m_slipsurface){
				//reconstruct result on mesh points (ghost included)-> stored in m_field.
				reconstructResults(results, data);
				//correct the BC of the slip walls
				computeSlipBCCorrector(m_field);
				// this slip dirichlet on the bulk nodes don't need to be distribute (they are already on points)
				//now you have a set of dirichlet BC in m_slip_bc_dir pierced vector internal.

				// so loop again on the components, reusing the previous result as starting guess, and setting
				// the boolean of assignBC as true (corrector stage of slip, read Dirichlet from m_slip_bc_dir)
				for(int comp = 0; comp<3; ++comp){
					//prepare the right hand side ;
					dvector1D rhs(geo->getNInternalVertices(), 0.0);
					assignBCAndEvaluateRHS(comp, true, laplaceStencils.get(), dataInv, rhs);
					//solve
					solveLaplace(rhs, results[comp]);
				}
			}
			//RECONSTRUCT STAGE --> /////////////////////////////////////////////////////////////////////////////////
			reconstructResults(results, data, movingElementList.get());

			// if you are in multistep stage apply the deformation to the mesh.
            // here is morphed also the m_originalDumpingSurface if active dumping.
			if(m_nstep > 1){
				apply();
			}


//    		//TODO DEBUG
//    		{
//    			bitpit::VTKUnstructuredGrid& vtk = getGeometry()->getPatch()->getVTK();
//#if MIMMO_ENABLE_MPI
//    			getGeometry()->getPatch()->setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_ALL);
//#endif
//
//    			dvecarr3E data(m_bsurface->getPatch()->getInternalCount());
//    			data = m_field.getDataAsVector();
//    			vtk.addData("field", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, data);
//
//    			vtk.setCounter(istep);
//    			getGeometry()->getPatch()->write(m_name + "_step");
//    			vtk.removeData("field");
//    			vtk.unsetCounter();
//    		}

			// if in multistep stage continue to update the other laplacian stuff up to "second-to-last" step.
			if(istep < m_nstep-1){

				//update the dumping function. using m_originalDumpingSurface deformed.
//				updateDumpingFunction();
				computeDumpingFunction();

				//enlarge the moving cell list taking its first vertex neighs and its second face neighs.
				propagateMaskMovingPoints(*(movingElementList.get()));

				// update the laplacian stencils
				GraphLaplStencil::MPVStencilUPtr updateLaplaceStencils = GraphLaplStencil::computeLaplacianStencils(*geo, movingElementList.get(), m_tol, &m_dumping);
				movingElementList->clear();

				laplaceStencils->getDataFrom(*(updateLaplaceStencils.get()), true); //only common elements are updated.

				// update the laplacian Matrix in solver free the updateLaplaceStencils
				updateLaplaceSolver(updateLaplaceStencils.get(), dataInv);
				updateLaplaceStencils = nullptr;

			}

			(*m_log)<<m_name<<" solved step "<<istep<<" out of total steps "<<m_nstep<<std::endl;

		} //end of multistep loop;


	} //end if on method


	if(m_nstep > 1){
		//this take the geometry to the original state and update the deformation field as the current
		//deformed grid minus the undeformed state (this directly on POINTS).
		restoreGeometry(undeformedTargetVertices);
		restoreBC();
	}

	//clear the solver;
	m_solver->clear();
	(*m_log) << bitpit::log::priority(bitpit::log::DEBUG);
}

} //end of mimmo namespace
