/*---------------------------------------------------------------------------*\
 * 
 *  CAMILO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License Commercial (//TODO Temporary header of license)
 *  This file is part of CAMILO.
 *
 *  CAMILO is a commercial software: you do not own rights to redistribute it 
 * 	and/or modify it both in source or pre-build formats
 *  Please contact Optimad offices for any further informations				
 *
 *  You should have received a copy of the Camilo Commercial License
 *  along with CAMILO, as well as the key to unlock the software.
 *
 \ *----------------*-----------------------------------------------------------*/
#include "MeshSelection.hpp"
#include "levelSet.hpp"
#include <cstddef>
using namespace mimmino;

//------------------------------------------------------------------------
//GENERIC	SELECTION class **********************************************
//------------------------------------------------------------------------

/*!
 * Basic Constructor
 */
GenericSelection::GenericSelection(){
	m_type = SelectionType::UNDEFINED;
	m_topo = 1; //default to surface geometry
	m_dual = false; //default to exact selection
};

/*!
 * Basic Destructor
 */
GenericSelection::~GenericSelection(){
	m_subpatch.reset(nullptr);
};

/*!
 * Copy Constructor
 */
GenericSelection::GenericSelection(const GenericSelection & other){
	*this = other;
};

/*!
 * Copy operator
 */
GenericSelection & GenericSelection::operator=(const GenericSelection & other){
	*(static_cast<BaseManipulation *>(this)) = *(static_cast<const BaseManipulation * >(&other));
	m_type = other.m_type;
	m_topo = other.m_topo;
	m_dual = other.m_dual;
	m_subpatch.reset(nullptr);
	std::unique_ptr<MimmoObject> dum(new MimmoObject(*(other.m_subpatch.get())));
	m_subpatch = std::move(dum);
	return *this;
};

/*!
 * It builds the input/output ports of the object
 */
void GenericSelection::buildPorts(){

	bool built = true;

	//input
	built = (built && createPortIn<MimmoObject *, GenericSelection>(this, &GenericSelection::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	built = (built && createPortIn<bool, GenericSelection>(this, &GenericSelection::setDual, PortType::M_VALUEB, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::BOOL));

	//output
	built = (built && createPortOut<MimmoObject *, GenericSelection>(this, &GenericSelection::getPatch, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	built = (built && createPortOut<livector1D, GenericSelection>(this, &GenericSelection::constrainedBoundary, PortType::M_VECTORLI, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::LONG));
	m_arePortsBuilt = built;
};


/*!
 * Return type of method for selection. See SelectionType enum for available methods.
 */
SelectionType	GenericSelection::whichMethod(){
	return	m_type;
};

/*!
 * Return pointer by copy to subpatch extracted by the class
 */
MimmoObject*	GenericSelection::getPatch(){
	return	m_subpatch.get();
};

/*!
 * Return pointer by copy to subpatch extracted by the class. Const overloading
 */
const MimmoObject*	GenericSelection::getPatch() const{
	return	m_subpatch.get();
};

/*!
 * Set link to target geometry for your selection. Reimplementation of 
 * mimmo::BaseManipulation::setGeometry();
 */
void	GenericSelection::setGeometry( MimmoObject * target){
	if(target->isEmpty()) return;
	m_geometry = target;
	//set topology informations
	m_topo = target->getType();
	
};


/*!
 * Set your class behavior selecting a portion of a target geoemetry. Given a initial set up, gets the dual
 * result (its negative) of current selection. For instance, in a class extracting geometry inside the volume of an elemental shape,
 * gets all other parts of it not included in the shape. For a class extracting portions of geometry
 * matching an external tessellation, gets all the other parts not matching it. For a class extracting portion of 
 * geometry identified by a PID list, get all other parts not marked by such list.
 * \param[in] flag to activate/deactivate "dual" feature set true/false
 */
void 	GenericSelection::setDual(bool flag ){
	m_dual = flag;
}

/*!
 * Return actual status of "dual" feature of the class. See setDual method.
 * \return  true/false for "dual" feature activated or not
 */
bool 	GenericSelection::isDual(){
	return m_dual;
};

/*!
 * Return list of constrained boundary nodes (all those boundary nodes of
 * the subpatch extracted which are not part of the boundary of the mother geometry)
 * \return list of boundary nodes ID, according to bitpit::PatchKernel indexing
 */
livector1D	GenericSelection::constrainedBoundary(){
	
	if(getGeometry()->isEmpty() || getPatch()->isEmpty())	return livector1D(0);
	livector1D sonV, fatherV;
	
	sonV	= getPatch()->extractBoundaryVertexID();
	bitpit::PiercedVector<bitpit::Cell> & existentCells = getPatch()->getCells();
	
	//going to create the dual surface to the current extraction.
	//create your subpatch.
	std::unique_ptr<MimmoObject> temp(new MimmoObject(m_topo));
	bitpit::PatchKernel * tri = getGeometry()->getPatch();
	
	bitpit::PiercedVector<bitpit::Vertex> & mapV = temp->getPatch()->getVertices();
	
	livector1D TT;
	
	long idV,idC;
	int sizeCC;
	int counter=0;
	bitpit::ElementInfo::Type eltype;
	
	for(auto && cell : tri->getCells()){
		
		idC = cell.getId();
		if (!(existentCells.exists(idC))){
			eltype = cell.getType();
			sizeCC = cell.getVertexCount();
			TT.resize(sizeCC);
		
			for(int i=0; i<sizeCC; ++i){
				idV = cell.getVertex(i);
				TT[i] = idV;
				if(!mapV.exists(idV))	temp->addVertex(tri->getVertexCoords(idV),idV);
			}
		temp->addConnectedCell(TT,eltype,idC);
		TT.clear();
		counter++;
		}
	}
	
	if(counter == 0)	return sonV;
		
	livector1D dualpV = temp->extractBoundaryVertexID();
	if(dualpV.empty())	return sonV;
	
	std::sort(sonV.begin(), sonV.end());
	std::sort(dualpV.begin(), dualpV.end());
	
	livector1D::iterator itS, itD;
	livector1D matching(sonV.size());
	
	counter = 0;
	itS =  sonV.begin();
	itD =  dualpV.begin();
	
	 	while (itS != sonV.end()) {
	 		if (itD != dualpV.end()){
				long valS = *itS;
				long valD = *itD;
				
				if(valD == valS){
					matching[counter] = valS;
					++counter;
					++itS;
					++itD;
				}else if(valD > valS){	
					++itS;
				}else{
					++itD;
				}
			}else{
				itS = sonV.end();
			}
		}
// 	livector1D::iterator found, itB, itE = sonV.end();
// 	livector1D toErase(sonV.size());
	
// 	for( itB = sonV.begin();  itB != itE; itB++ ){
// 		found = std::find(fatherV.begin(), fatherV.end(), *itB);
// 		if( found != fatherV.end() ){	
// 			toErase[counter] = (*itB);
// 			++counter;
// 		}
// 	}
// 	toErase.resize(counter);
// 	itE = toErase.begin();
// 	itB = sonV.begin();
// 	
//     livector1D result(sonV.size());
// 	counter= 0;
// 	while (itB != sonV.end()) {
// 		long val = *itB;
// 		if (itE == toErase.end() || val != *itE) {
// 			result[counter] = val;
// 			++counter;
// 		} else {
// 			++itE;
// 		}
// 		++itB;
// 	}
	matching.resize(counter);
	return matching;
};


/*! 
 * Execute your object. A selection is extracted and trasferred in
 * an indipendent MimmoObject structure pointed by m_subpatch member
 */
void GenericSelection::execute(){
	
	if(getGeometry()->isEmpty()) return;
	
	m_subpatch.reset(nullptr);
	
    livector1D extracted = extractSelection();

    if(extracted.empty()) return;
	
	//create your subpatch.
	std::unique_ptr<MimmoObject> temp(new MimmoObject(m_topo));
	bitpit::PatchKernel * tri = getGeometry()->getPatch();
	
	bitpit::PiercedVector<bitpit::Vertex> & mapV = temp->getPatch()->getVertices();
	
	livector1D TT;

	long idV;
	int sizeCC;
	bitpit::ElementInfo::Type eltype;
	
	if (m_topo != 3){
	    for(auto && idCell : extracted){

	        bitpit::Cell & cell = tri->getCell(idCell);
	        eltype = cell.getType();
	        sizeCC = cell.getVertexCount();
	        TT.resize(sizeCC);

	        for(int i=0; i<sizeCC; ++i){
	            idV = cell.getVertex(i);
	            TT[i] = idV;

	            if(!mapV.exists(idV))	temp->addVertex(tri->getVertexCoords(idV),idV);
	        }
	        temp->addConnectedCell(TT,eltype,idCell);
	        TT.clear();
	    }
	}
	else{
        for(auto && idV : extracted){
            temp->addVertex(tri->getVertexCoords(idV),idV);
        }
	}

	m_subpatch = std::move(temp);
	tri = NULL;
};

/*!
 * Plot optional result of the class in execution, that is the selected patch
 * as standard vtk unstructured grid.
 */
void GenericSelection::plotOptionalResults(){
	if(getPatch()->isEmpty()) return;
	std::string name = m_name + "_" + std::to_string(getClassCounter()) +  "_Patch";
	getPatch()->getPatch()->write(name);
}

//------------------------------------------------------------------------
//SELECTION	BY BOX class 	**********************************************
//------------------------------------------------------------------------

/*!
 * Basic Constructor
 */
SelectionByBox::SelectionByBox(){
	m_name = "MiMMiNO.SelectionByBox";
	m_type = SelectionType::BOX;
};

/*!
 * Custom Constructor
 * \param[in] origin of the box->baricenter
 * \param[in] span	 of the box, width/height/depth
 * \param[in] target	pointer to a target geometry
 */
SelectionByBox::SelectionByBox(darray3E origin, darray3E span, MimmoObject * target){
	m_name = "MiMMiNO.SelectionByBox";
	m_type = SelectionType::BOX;
	setGeometry(target);
	setOrigin(origin);
	setSpan(span[0],span[1],span[2]);
};

/*!
 * Destructor
 */
SelectionByBox::~SelectionByBox(){};

/*!
 * Copy Constructor
 */
SelectionByBox::SelectionByBox(const SelectionByBox & other){
	*this = other;
};

/*!
 * Copy operator
 */
SelectionByBox & SelectionByBox::operator=(const SelectionByBox & other){
	*(static_cast<GenericSelection * >(this)) = *(static_cast<const GenericSelection *>(&other));
	*(static_cast<Cube * >(this)) = *(static_cast<const Cube *>(&other));
	return *this;
};

/*!
 * It builds the input/output ports of the object
 */
void SelectionByBox::buildPorts(){

	bool built = true;

	//inheritance
	GenericSelection::buildPorts();

	//input
	built = (built && createPortIn<darray3E, SelectionByBox>(this, &SelectionByBox::setOrigin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<darray3E, SelectionByBox>(this, &SelectionByBox::setSpan, PortType::M_SPAN, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dmatrix33E, SelectionByBox>(this, &SelectionByBox::setRefSystem, PortType::M_AXES, mimmo::pin::containerTAG::ARR3ARR3, mimmo::pin::dataTAG::FLOAT));

	m_arePortsBuilt = built;
};

/*!
 * Clear content of the class
 */
void SelectionByBox::clear(){
	m_subpatch.reset(nullptr);
	BaseManipulation::clear();
};



/*!
 * Extract portion of target geometry which are enough near to the external geometry provided
 * \return ids of cell of target tesselation extracted 
 */
livector1D SelectionByBox::extractSelection(){
	switch(m_topo){
		case 3:
			if(m_dual)	return  excludeCloudPoints(getGeometry());
			else		return	includeCloudPoints(getGeometry());
			break;
		default:
			if(m_dual)	return  excludeGeometry(getGeometry());
			else		return	includeGeometry(getGeometry());
			break;
	}
	
	return livector1D(0); //never been reached
};


/*!
 * Get infos from a XML bitpit::Config::section. The parameters available are
 * 
 * 1) Dual     - boolean to get straight what given by selection method or its exact dual
 * 2) Origin   - array of 3 doubles identifying origin
 * 3) Span	   - span of the box (width, height, depth)
 * 4) RefSystem - reference system of the box;
 * 
 * Geometry is mandatorily passed through ports. 
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionByBox::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	//start absorbing
	if(slotXML.hasOption("Dual")){
		std::string input = slotXML.get("Dual");
		input = bitpit::utils::trim(input);
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setDual(value);
	}
	
	if(slotXML.hasOption("Origin")){
		std::string input = slotXML.get("Origin");
		input = bitpit::utils::trim(input);
		darray3E temp = {{0.0,0.0,0.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp[0]>>temp[1]>>temp[2];
			setOrigin(temp);
		}else{
			setOrigin(temp);
		}	
	}
	
	if(slotXML.hasOption("Span")){
		std::string input = slotXML.get("Span");
		input = bitpit::utils::trim(input);
		darray3E temp = {{1.0,1.0,1.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp[0]>>temp[1]>>temp[2];
			setSpan(temp);
		}else{
			setSpan(temp);
		}	
	}	

	if(slotXML.hasSection("RefSystem")){
		
		bitpit::Config::Section & axesXML = slotXML.getSection("RefSystem");
		dmatrix33E axes;
		for(int i=0; i<3; ++i){
			axes[i].fill(0.0);
			axes[i][i] = 1.0;
		}
		
		if(axesXML.hasOption("axis0")){
			std::string input = axesXML.get("axis0");
			input = bitpit::utils::trim(input);
			if(!input.empty()){
				std::stringstream ss(input);
				ss>>axes[0][0]>>axes[0][1]>>axes[0][2];
			}
		}
		
		if(axesXML.hasOption("axis1")){
			std::string input = axesXML.get("axis1");
			input = bitpit::utils::trim(input);
			if(!input.empty()){
				std::stringstream ss(input);
				ss>>axes[1][0]>>axes[1][1]>>axes[1][2];
			}
		}

		if(axesXML.hasOption("axis2")){
			std::string input = axesXML.get("axis2");
			input = bitpit::utils::trim(input);
			if(!input.empty()){
				std::stringstream ss(input);
				ss>>axes[2][0]>>axes[2][1]>>axes[2][2];
			}
		}
		
		setRefSystem(axes);
	}	
	
		if(slotXML.hasOption("PlotInExecution")){
			std::string input = slotXML.get("PlotInExecution");
			input = bitpit::utils::trim(input);
			bool value = false;
			if(!input.empty()){
				std::stringstream ss(input);
				ss >> value;
			}
			setPlotInExecution(value);
		}
		
		if(slotXML.hasOption("OutputPlot")){
			std::string input = slotXML.get("OutputPlot");
			input = bitpit::utils::trim(input);
			std::string temp = ".";
			if(!input.empty())	setOutputPlot(input);
			else			  	setOutputPlot(temp);
		}
	
	return;	
};

/*!
 * Plot infos from a XML bitpit::Config::section. The parameters available are
 * 
 * 1) Dual     - boolean to get straight what given by selection method or its exact dual
 * 2) Origin   - array of 3 doubles identifying origin
 * 3) Span	   - span of the box (width, height, depth)
 * 4) RefSystem - reference system of the box;
 * 
 * Geometry is mandatorily passed through ports. 
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionByBox::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
	int value = m_dual;
	slotXML.set("Dual", std::to_string(value));
	
	{
		darray3E org = getOrigin();
		std::stringstream ss;
		ss<<std::scientific<<org[0]<<'\t'<<org[1]<<'\t'<<org[2];
		slotXML.set("Origin",ss.str());
	}
	
	{
		darray3E span = getSpan();
		std::stringstream ss;
		ss<<std::scientific<<span[0]<<'\t'<<span[1]<<'\t'<<span[2];
		slotXML.set("Span",ss.str());
	}
	
	{
		dmatrix33E axes = getRefSystem();
		bitpit::Config::Section & axesXML = slotXML.addSection("RefSystem");
		
		for(int i=0; i<3; ++i){
			std::string name = "axis"+std::to_string(i);
			std::stringstream ss;
			ss<<std::scientific<<axes[i][0]<<'\t'<<axes[i][1]<<'\t'<<axes[i][2];
			axesXML.set(name, ss.str());
		}
	}
	
		if(isPlotInExecution()){
			slotXML.set("PlotInExecution", std::to_string(1));
		}
		
		if(m_outputPlot != "."){
			slotXML.set("OutputPlot", m_outputPlot);
		}
	
	return;
};




//------------------------------------------------------------------------
//SELECTION	BY CYLINDER class 	**********************************************
//------------------------------------------------------------------------

/*!
 * Basic Constructor
 */
SelectionByCylinder::SelectionByCylinder(){
	m_name = "MiMMiNO.SelectionByCylinder";
	m_type = SelectionType::CYLINDER;
};

/*!
 * Custom Constructor. Pay attention span of angular coordinate must be at most 2*pi.
 * \param[in] origin of the cylinder->baricenter
 * \param[in] span	 of the cylinder, basis radius/ span of angular coord in radians/height
 * \param[in] infLimTheta	starting origin of the angular coordinate. default is 0 radians.
 * \param[in] mainAxis	orientation of the cylinder height axis 
 * \param[in] target	pointer to a target geometry 
 */
SelectionByCylinder::SelectionByCylinder(darray3E origin, darray3E span, double infLimTheta, darray3E mainAxis, MimmoObject * target){
	m_name = "MiMMiNO.SelectionByCylinder";
	m_type = SelectionType::CYLINDER;
	setGeometry(target);
	setOrigin(origin);
	setSpan(span[0],span[1],span[2]);
	setInfLimits(infLimTheta,1);
	setRefSystem(2, mainAxis);
};

/*!
 * Destructor
 */
SelectionByCylinder::~SelectionByCylinder(){};

/*!
 * Copy Constructor
 */
SelectionByCylinder::SelectionByCylinder(const SelectionByCylinder & other){
	*this = other;
};

/*!
 * Copy operator
 */
SelectionByCylinder & SelectionByCylinder::operator=(const SelectionByCylinder & other){
	*(static_cast<GenericSelection * >(this)) = *(static_cast<const GenericSelection *>(&other));
	*(static_cast<Cylinder * >(this)) = *(static_cast<const Cylinder *>(&other));
	return *this;
};

/*!
 * It builds the input/output ports of the object
 */
void SelectionByCylinder::buildPorts(){

	bool built = true;

	//inheritance
	GenericSelection::buildPorts();

	//input
	built = (built && createPortIn<darray3E, SelectionByCylinder>(this, &SelectionByCylinder::setOrigin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<darray3E, SelectionByCylinder>(this, &SelectionByCylinder::setSpan, PortType::M_SPAN, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dmatrix33E, SelectionByCylinder>(this, &SelectionByCylinder::setRefSystem, PortType::M_AXES, mimmo::pin::containerTAG::ARR3ARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<darray3E, SelectionByCylinder>(this, &SelectionByCylinder::setInfLimits, PortType::M_INFLIMITS, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));

	m_arePortsBuilt = built;
};

/*!
 * Clear content of the class
 */
void SelectionByCylinder::clear(){
	m_subpatch.reset(nullptr);
	BaseManipulation::clear();
};

/*!
 * Extract portion of target geometry which are enough near to the external geometry provided
 * \return ids of cell of target tesselation extracted 
 */
livector1D SelectionByCylinder::extractSelection(){
	switch(m_topo){
		case 3:
			if(m_dual)	return  excludeCloudPoints(getGeometry());
			else		return	includeCloudPoints(getGeometry());
			break;
		default:
			if(m_dual)	return  excludeGeometry(getGeometry());
			else		return	includeGeometry(getGeometry());
			break;
	}
	
	return livector1D(0); //never been reached
};



/*!
 * Get infos from a XML bitpit::Config::section. The parameters available are
 * 
 * 1) Dual     - boolean to get straight what given by selection method or its exact dual
 * 2) Origin   - array of 3 doubles identifying origin
 * 3) Span	   - span of the cylinder (base radius, azimuthal span(in radians), height)
 * 4) RefSystem - reference system of the cylinder (z/2 is the cylinder height axis);
 * 5) InfLimits - inferior limits for span dimensioning
 * 
 * Geometry is mandatorily passed through ports. 
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionByCylinder::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	//start absorbing
	if(slotXML.hasOption("Dual")){
		std::string input = slotXML.get("Dual");
		input = bitpit::utils::trim(input);
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setDual(value);
	}
	
	if(slotXML.hasOption("Origin")){
		std::string input = slotXML.get("Origin");
		input = bitpit::utils::trim(input);
		darray3E temp = {{0.0,0.0,0.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp[0]>>temp[1]>>temp[2];
			setOrigin(temp);
		}else{
			setOrigin(temp);
		}	
	}
	
	if(slotXML.hasOption("Span")){
		std::string input = slotXML.get("Span");
		input = bitpit::utils::trim(input);
		darray3E temp = {{1.0,2.0*M_PI,1.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp[0]>>temp[1]>>temp[2];
			setSpan(temp);
		}else{
			setSpan(temp);
		}	
	}	
	
	if(slotXML.hasSection("RefSystem")){
		
		bitpit::Config::Section & axesXML = slotXML.getSection("RefSystem");
		dmatrix33E axes;
		for(int i=0; i<3; ++i){
			axes[i].fill(0.0);
			axes[i][i] = 1.0;
		}
		
		if(axesXML.hasOption("axis0")){
			std::string input = axesXML.get("axis0");
			input = bitpit::utils::trim(input);
			if(!input.empty()){
				std::stringstream ss(input);
				ss>>axes[0][0]>>axes[0][1]>>axes[0][2];
			}
		}
		
		if(axesXML.hasOption("axis1")){
			std::string input = axesXML.get("axis1");
			input = bitpit::utils::trim(input);
			if(!input.empty()){
				std::stringstream ss(input);
				ss>>axes[1][0]>>axes[1][1]>>axes[1][2];
			}
		}
		
		if(axesXML.hasOption("axis2")){
			std::string input = axesXML.get("axis2");
			input = bitpit::utils::trim(input);
			if(!input.empty()){
				std::stringstream ss(input);
				ss>>axes[2][0]>>axes[2][1]>>axes[2][2];
			}
		}
		setRefSystem(axes);
	}		
	
	if(slotXML.hasOption("InfLimits")){
		std::string input = slotXML.get("InfLimits");
		input = bitpit::utils::trim(input);
		darray3E temp = {{0.0,0.0,0.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp[0]>>temp[1]>>temp[2];
			setInfLimits(temp);
		}else{
			setInfLimits(temp);
		}	
	}	
	
	
		if(slotXML.hasOption("PlotInExecution")){
			std::string input = slotXML.get("PlotInExecution");
			input = bitpit::utils::trim(input);
			bool value = false;
			if(!input.empty()){
				std::stringstream ss(input);
				ss >> value;
			}
			setPlotInExecution(value);
		}
		
		if(slotXML.hasOption("OutputPlot")){
			std::string input = slotXML.get("OutputPlot");
			input = bitpit::utils::trim(input);
			std::string temp = ".";
			if(!input.empty())	setOutputPlot(input);
			else			  	setOutputPlot(temp);
		}
	
	return;	
};

/*!
 * Plot infos from a XML bitpit::Config::section. The parameters available are
 * 
 * 1) Dual     - boolean to get straight what given by selection method or its exact dual
 * 2) Origin   - array of 3 doubles identifying origin
 * 3) Span	   - span of the cylinder (base radius, azimuthal span(in radians), height)
 * 4) RefSystem - reference system of the cylinder (z/2 is the cylinder height axis);
 * 5) InfLimits - inferior limits for span dimensioning
 * 
 * Geometry is mandatorily passed through ports. 
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionByCylinder::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
	int value = m_dual;
	slotXML.set("Dual", std::to_string(value));
	
	
	{
		darray3E org = getOrigin();
		std::stringstream ss;
		ss<<std::scientific<<org[0]<<'\t'<<org[1]<<'\t'<<org[2];
		slotXML.set("Origin",ss.str());
	}
	
	{
		darray3E span = getSpan();
		std::stringstream ss;
		ss<<std::scientific<<span[0]<<'\t'<<span[1]<<'\t'<<span[2];
		slotXML.set("Span",ss.str());
	}
	
	{
		dmatrix33E axes = getRefSystem();
		bitpit::Config::Section & axesXML = slotXML.addSection("RefSystem");
		
		for(int i=0; i<3; ++i){
			std::string name = "axis"+std::to_string(i);
			std::stringstream ss;
			ss<<std::scientific<<axes[i][0]<<'\t'<<axes[i][1]<<'\t'<<axes[i][2];
			axesXML.set(name, ss.str());
		}
	}
	
	{
		darray3E inflim = getInfLimits();
		std::stringstream ss;
		ss<<std::scientific<<inflim[0]<<'\t'<<inflim[1]<<'\t'<<inflim[2];
		slotXML.set("InfLimits",ss.str());
	}
	
		if(isPlotInExecution()){
			slotXML.set("PlotInExecution", std::to_string(1));
		}
		
		if(m_outputPlot != "."){
			slotXML.set("OutputPlot", m_outputPlot);
		}
	
	return;
};



//------------------------------------------------------------------------
//SELECTION	BY SPHERE class 	******************************************
//------------------------------------------------------------------------

/*!
 * Basic Constructor
 */
SelectionBySphere::SelectionBySphere(){
	m_name = "MiMMiNO.SelectionBySphere";
	m_type = SelectionType::SPHERE;
};

/*!
 * Custom Constructor. Pay attention span of angular and polar coords are at most 2*pi and pi respectively. 
 *  Inf limits of polar coordinate must be > 0 and < pi.
 * \param[in] origin of the sphere->baricenter
 * \param[in] span	 of the cylinder, main radius/ span of angular coord in radians/span of the polar coord in radians
 * \param[in] infLimTheta	starting origin of the angular coordinate. default is 0 radians.
 * \param[in] infLimPhi	starting origin of the polar coordinate. default is 0 radians.
 * \param[in] target	pointer to a target geometry 
 */
SelectionBySphere::SelectionBySphere(darray3E origin, darray3E span, double infLimTheta, double infLimPhi, MimmoObject * target){
	m_name = "MiMMiNO.SelectionBySphere";
	m_type = SelectionType::SPHERE;
	setGeometry(target);
	setOrigin(origin);
	setSpan(span[0],span[1],span[2]);
	setInfLimits(infLimTheta,1);
	setInfLimits(infLimPhi,2);
};

/*!
 * Destructor
 */
SelectionBySphere::~SelectionBySphere(){};

/*!
 * Copy Constructor
 */
SelectionBySphere::SelectionBySphere(const SelectionBySphere & other){
	*this = other;
};

/*!
 * Copy operator
 */
SelectionBySphere & SelectionBySphere::operator=(const SelectionBySphere & other){
	*(static_cast<GenericSelection * >(this)) = *(static_cast<const GenericSelection *>(&other));
	*(static_cast<Sphere * >(this)) = *(static_cast<const Sphere *>(&other));
	return *this;
};

/*!
 * It builds the input/output ports of the object
 */
void SelectionBySphere::buildPorts(){

	bool built = true;

	//inheritance
	GenericSelection::buildPorts();

	//input
	built = (built && createPortIn<darray3E, SelectionBySphere>(this, &SelectionBySphere::setOrigin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<darray3E, SelectionBySphere>(this, &SelectionBySphere::setSpan, PortType::M_SPAN, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dmatrix33E, SelectionBySphere>(this, &SelectionBySphere::setRefSystem, PortType::M_AXES, mimmo::pin::containerTAG::ARR3ARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<darray3E, SelectionBySphere>(this, &SelectionBySphere::setInfLimits, PortType::M_INFLIMITS, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));

	m_arePortsBuilt = built;
};

/*!
 * Clear content of the class
 */
void SelectionBySphere::clear(){
	m_subpatch.release();
	BaseManipulation::clear();
};

/*!
 * Extract portion of target geometry which are enough near to the external geometry provided
 * \return ids of cell of target tesselation extracted 
 */
livector1D SelectionBySphere::extractSelection(){
	switch(m_topo){
		case 3:
			if(m_dual)	return  excludeCloudPoints(getGeometry());
			else		return	includeCloudPoints(getGeometry());
			break;
		default:
			if(m_dual)	return  excludeGeometry(getGeometry());
			else		return	includeGeometry(getGeometry());
			break;
	}
	
	return livector1D(0); //never been reached
};

/*!
 * Get infos from a XML bitpit::Config::section. The parameters available are
 * 
 * 1) Dual     - boolean to get straight what given by selection method or its exact dual
 * 2) Origin   - array of 3 doubles identifying origin
 * 3) Span	   - span of the sphere (radius, azimuthal span(in radians), polar span(in radians))
 * 4) RefSystem - reference system of the sphere (z/2 is the pole/pole axis);
 * 5) InfLimits - inferior limits for span dimensioning
 * 
 * Geometry is mandatorily passed through ports. 
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionBySphere::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	//start absorbing
	if(slotXML.hasOption("Dual")){
		std::string input = slotXML.get("Dual");
		input = bitpit::utils::trim(input);
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setDual(value);
	}
	
	if(slotXML.hasOption("Origin")){
		std::string input = slotXML.get("Origin");
		input = bitpit::utils::trim(input);
		darray3E temp = {{0.0,0.0,0.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp[0]>>temp[1]>>temp[2];
			setOrigin(temp);
		}else{
			setOrigin(temp);
		}	
	}
	
	if(slotXML.hasOption("Span")){
		std::string input = slotXML.get("Span");
		input = bitpit::utils::trim(input);
		darray3E temp = {{1.0,2.0*M_PI,M_PI}};
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp[0]>>temp[1]>>temp[2];
			setSpan(temp);
		}else{
			setSpan(temp);
		}	
	}	
	
	if(slotXML.hasSection("RefSystem")){
		
		bitpit::Config::Section & axesXML = slotXML.getSection("RefSystem");
		dmatrix33E axes;
		for(int i=0; i<3; ++i){
			axes[i].fill(0.0);
			axes[i][i] = 1.0;
		}
		
		if(axesXML.hasOption("axis0")){
			std::string input = axesXML.get("axis0");
			input = bitpit::utils::trim(input);
			if(!input.empty()){
				std::stringstream ss(input);
				ss>>axes[0][0]>>axes[0][1]>>axes[0][2];
			}
		}
		
		if(axesXML.hasOption("axis1")){
			std::string input = axesXML.get("axis1");
			input = bitpit::utils::trim(input);
			if(!input.empty()){
				std::stringstream ss(input);
				ss>>axes[1][0]>>axes[1][1]>>axes[1][2];
			}
		}
		
		if(axesXML.hasOption("axis2")){
			std::string input = axesXML.get("axis2");
			input = bitpit::utils::trim(input);
			if(!input.empty()){
				std::stringstream ss(input);
				ss>>axes[2][0]>>axes[2][1]>>axes[2][2];
			}
		}
		setRefSystem(axes);
	}		
	
	if(slotXML.hasOption("InfLimits")){
		std::string input = slotXML.get("InfLimits");
		input = bitpit::utils::trim(input);
		darray3E temp = {{0.0,0.0,0.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp[0]>>temp[1]>>temp[2];
			setInfLimits(temp);
		}else{
			setInfLimits(temp);
		}	
	}	
	
	
		if(slotXML.hasOption("PlotInExecution")){
			std::string input = slotXML.get("PlotInExecution");
			input = bitpit::utils::trim(input);
			bool value = false;
			if(!input.empty()){
				std::stringstream ss(input);
				ss >> value;
			}
			setPlotInExecution(value);
		}
		
		if(slotXML.hasOption("OutputPlot")){
			std::string input = slotXML.get("OutputPlot");
			input = bitpit::utils::trim(input);
			std::string temp = ".";
			if(!input.empty())	setOutputPlot(input);
			else			  	setOutputPlot(temp);
		}
	
	return;	
};

/*!
 * Plot infos from a XML bitpit::Config::section. The parameters available are
 * 
 * 1) Dual     - boolean to get straight what given by selection method or its exact dual
 * 2) Origin   - array of 3 doubles identifying origin
 * 3) Span	   - span of the sphere (radius, azimuthal span(in radians), polar span(in radians))
 * 4) RefSystem - reference system of the sphere (z/2 is the pole/pole axis);
 * 5) InfLimits - inferior limits for span dimensioning
 * 
 * 
 * Geometry is mandatorily passed through ports. 
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionBySphere::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
	int value = m_dual;
	slotXML.set("Dual", std::to_string(value));
	
	
	{
		darray3E org = getOrigin();
		std::stringstream ss;
		ss<<std::scientific<<org[0]<<'\t'<<org[1]<<'\t'<<org[2];
		slotXML.set("Origin",ss.str());
	}
	
	{
		darray3E span = getSpan();
		std::stringstream ss;
		ss<<std::scientific<<span[0]<<'\t'<<span[1]<<'\t'<<span[2];
		slotXML.set("Span",ss.str());
	}
	
	{
		dmatrix33E axes = getRefSystem();
		bitpit::Config::Section & axesXML = slotXML.addSection("RefSystem");
		
		for(int i=0; i<3; ++i){
			std::string name = "axis"+std::to_string(i);
			std::stringstream ss;
			ss<<std::scientific<<axes[i][0]<<'\t'<<axes[i][1]<<'\t'<<axes[i][2];
			axesXML.set(name, ss.str());
		}
	}
	
	{
		darray3E inflim = getInfLimits();
		std::stringstream ss;
		ss<<std::scientific<<inflim[0]<<'\t'<<inflim[1]<<'\t'<<inflim[2];
		slotXML.set("InfLimits",ss.str());
	}
	
		if(isPlotInExecution()){
			slotXML.set("PlotInExecution", std::to_string(1));
		}
		
		if(m_outputPlot != "."){
			slotXML.set("OutputPlot", m_outputPlot);
		}
	
	return;
};


//------------------------------------------------------------------------
//SELECTION	BY MAPPING class 	******************************************
//------------------------------------------------------------------------

/*!
 * Basic Constructor. Need to know kind of topology chosen 1-3D surface, 2-VolumeMesh
 * . Other options are not available, if forced trigger default value of 1.
 */
SelectionByMapping::SelectionByMapping(int topo){
	m_name = "MiMMiNO.SelectionByMapping";
	m_type = SelectionType::MAPPING;
	m_tolerance = 1.E-08;
	
	topo = std::min(1,topo);
	if(topo > 2)	topo=1;
	m_topo = topo;
	
	
	m_allowedType.resize(3);
	m_allowedType[1].insert(FileType::STL);
	m_allowedType[1].insert(FileType::STVTU);
	m_allowedType[1].insert(FileType::SQVTU);
	m_allowedType[1].insert(FileType::NAS);
	m_allowedType[2].insert(FileType::VTVTU);
	m_allowedType[2].insert(FileType::VHVTU);
	
};

/*!
 * Custom Constructor. Non null geometry set also the topology allowed bu the class.
 * \param[in]	geolist	list of external files to read from, for comparison purposes 
 * \param[in]	target	pointer to target geometry
 * \param[in]	tolerance	proximity criterium tolerance 
 */
SelectionByMapping::SelectionByMapping(std::unordered_map<std::string, int> & geolist, MimmoObject * target, double tolerance){
	m_name = "MiMMiNO.SelectionByMapping";
	m_type = SelectionType::MAPPING;
	m_tolerance = 1.E-08;
	
	m_allowedType.resize(3);
	m_allowedType[1].insert(FileType::STL);
	m_allowedType[1].insert(FileType::STVTU);
	m_allowedType[1].insert(FileType::SQVTU);
	m_allowedType[1].insert(FileType::NAS);
	m_allowedType[2].insert(FileType::VTVTU);
	m_allowedType[2].insert(FileType::VHVTU);
	
	if(target != NULL){
		m_topo = target->getType();
		m_topo = std::min(1, m_topo);
		if(m_topo > 2)	m_topo = 1;
		setGeometry(target);
		setFiles(geolist);
		m_tolerance = tolerance;
	}
};

/*!
 * Destructor
 */
SelectionByMapping::~SelectionByMapping(){};

/*!
 * copy Constructor
 */
SelectionByMapping::SelectionByMapping(const SelectionByMapping & other){
	*this = other;
};

/*!
 * Copy Operator
 */
SelectionByMapping & SelectionByMapping::operator=(const SelectionByMapping & other){
	*(static_cast<GenericSelection * >(this)) = *(static_cast<const GenericSelection * >(&other));
	m_tolerance = other.m_tolerance;
	m_geolist = other.m_geolist;
	m_allowedType = other.m_allowedType;
	return *this;
};

/*!
 * It builds the input/output ports of the object
 */
void SelectionByMapping::buildPorts(){

	bool built = true;

	//inheritance
	GenericSelection::buildPorts();

	m_arePortsBuilt = built;
};

/*!
 * Return tolerance actually set in your class
 */
double	SelectionByMapping::getTolerance(){
	return	m_tolerance;
};

/*!
 * Set Proximity tolerance. Under this threshold compared geometries are considered coincident
 */
void 	SelectionByMapping::setTolerance(double tol){
	if(tol == 0.0){
		tol = 1.e-8;
	}
	m_tolerance = tol;
};

/*!
 * Set link to target geometry for your selection. Reimplementation of 
 * GenericSelection::setGeometry();
 */
void	SelectionByMapping::setGeometry( MimmoObject * target){
	
	if(target->getType() != m_topo){
		std::cout<<"SelectionMapping Cannot support current geometry. Topology not supported."<<std::endl;
		return;
	}
	m_geometry = target;
	
};

/*!
 * Return the actual list of external geometry files to read from and compare your target geometry for mapping purposes
 */
const std::unordered_map<std::string, int> & SelectionByMapping::getFiles() const{
	return	m_geolist;
}

/*!
 * Set a list of external geometry files to read from and compare your target geometry for mapping purposes
 * \param[in] files list external geometries to be read
 */
void	SelectionByMapping::setFiles(std::unordered_map<std::string, int>  files){
	for(auto && val : files){
		addFile(val);
	}	
};

/*!
 * Add a new file to a list of external geometry files
 *\param[in] file of external geometry to be read
 */
void 	SelectionByMapping::addFile(std::pair<std::string, int> file){
	int type = m_topo;
	if(m_allowedType[type].find(file.second) != m_allowedType[type].end()){									
		m_geolist.insert(file);
	}	
};

/*!
 * Remove an existent file to a list of external geometry files. If not in the list, do nothing 
 * \param[in] file to be removed from the list 
 */
void 	SelectionByMapping::removeFile(std::string file){
	 if(m_geolist.find(file) != m_geolist.end())	m_geolist.erase(file);
};

/*!
 * Empty your list of file for mapping 
 */
void 	SelectionByMapping::removeFiles(){
	m_geolist.clear();
};

/*!
 * Clear your class
 */
void SelectionByMapping::clear(){
	m_subpatch.reset(nullptr);
	removeFiles();
	m_topo = 0;
	BaseManipulation::clear();
};

/*!
 * Extract portion of target geometry which are enough near to the external geometry provided
 * \return ids of cell of target tesselation extracted 
 */
livector1D SelectionByMapping::extractSelection(){
	
	if(!(getGeometry()->isBvTreeBuilt()))	getGeometry()->buildBvTree();
	std::set<long> cellList;
	
	for (auto && file : m_geolist){
		livector1D list = getProximity(file);
		cellList.insert(list.begin(), list.end());
	}
	
	//check if dual selection is triggered
	livector1D result;
	if(m_dual){
		livector1D totID = getGeometry()->getMapCell();
		result.resize(totID.size() - cellList.size());
		if(result.size() == 0) return result;
		
		std::sort(totID.begin(), totID.end());
		int counter = 0;
		auto tot_it  = totID.begin();
		auto cand_it = cellList.begin();
		while(tot_it != totID.end()){
			long val = *tot_it;
			if (cand_it == cellList.end() || val != *cand_it) {
				result[counter] = val;
				++counter;
			} else {
				++cand_it;
			}
			++tot_it;
		}
	}else{
		result.insert(result.end(), cellList.begin(), cellList.end());
	}	
	return	result;
};

/*!
 * Return portion of target geometry near to an external geometry
 * \param[in] file of external geometry to be compared.
 * \param[in] nbins	number of total raw points for level set evaluation
 */
livector1D SelectionByMapping::getProximity(std::pair<std::string, int> val){
	
	svector1D info = extractInfo(val.first);
	
	MimmoGeometry * geo = new MimmoGeometry();
	geo->setRead(true);
	geo->setWrite(false);
	geo->setReadDir(info[0]);
	geo->setReadFilename(info[1]);
	geo->setReadFileType(val.second);
	geo->setBuildBvTree(true);
	geo->execute();

	if(geo->getGeometry()->getNVertex() == 0 || geo->getGeometry()->getNCells() == 0 ){ 
		std::cout<<"failed to read geometry in SelectionByMapping::getProximity"<<std::endl;
		return livector1D();
	}	
	livector1D result = mimmo::bvTreeUtils::selectByPatch(geo->getGeometry()->getBvTree(), getGeometry()->getBvTree(), m_tolerance);
	delete geo; 
	geo=NULL;
	
	return	result;
};

/*!
 * Extract root dir/filename/tag from an absolute file pattern
 */
svector1D SelectionByMapping::extractInfo(std::string file){
	
	std::string root, name, tag,temp;
	std::string key1=".", key2="/\\";
	
	std::size_t found = file.find_last_of(key2); 
	root = file.substr(0, found);
	temp = file.substr(found+1);
	
	found = temp.find_last_of(key1);
	name = temp.substr(0,found);
	tag = temp.substr(found+1);
	
	svector1D result(3);
	result[0] = root;
	result[1] = name;
	result[2] = tag;
	
	return 	result;
}


/*!
 * Get infos from a XML bitpit::Config::section. The parameters available are
 * 
 * 1) Dual       - boolean to get straight what given by selection method or its exact dual
 * 2) Topology   - set Topology for your mapping 1- 3D surface, 2- Volume mesh, 0 none;
 * 3) Tolerance  - tolerance for detect proximity volume in which perform mapping 
 * 4) Files	     - list by their filepaths of external geometries to be mapped on target one 
 * 
 * Geometry is mandatorily passed through ports. 
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionByMapping::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){

	//checking topology
	if(slotXML.hasOption("Topology")){
		std::string input = slotXML.get("Topology");
		input = bitpit::utils::trim(input);
		int temp = -1;
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp;
		}
		if(m_topo != temp)	return;
	}	
	
	//start absorbing
	if(slotXML.hasOption("Dual")){
		std::string input = slotXML.get("Dual");
		input = bitpit::utils::trim(input);
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setDual(value);
	}

	if(slotXML.hasOption("Tolerance")){
		std::string input = slotXML.get("Tolerance");
		input = bitpit::utils::trim(input);
		double temp = 1.0E-8;
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp;
		}
			setTolerance(temp);
			
	}
	
	std::unordered_map<std::string, int> mapp;
	if(slotXML.hasSection("Files")){
		
		bitpit::Config::Section & filesXML = slotXML.getSection("Files");
		
		for(auto & subfile : filesXML.getSections()){
			std::string path;
			std::string tag;
			
			if(subfile.second->hasOption("fullpath"))	{
				path = subfile.second->get("fullpath");
				path = bitpit::utils::trim(path);
			}
		    if(subfile.second->hasOption("tag")){
				tag = subfile.second->get("tag");
				tag = bitpit::utils::trim(tag);
				//check tag;
				auto maybe_tag = FileType::_from_string_nothrow(tag.c_str());
				if(!maybe_tag)	tag.clear();
				else	tag = maybe_tag->_to_string();
			}	
			
			if(!path.empty() && !tag.empty()){
				mapp[path] = (int) FileType::_from_string(tag.c_str());
			}	
		}

		setFiles(mapp);
		
	}	
	
		if(slotXML.hasOption("PlotInExecution")){
			std::string input = slotXML.get("PlotInExecution");
			input = bitpit::utils::trim(input);
			bool value = false;
			if(!input.empty()){
				std::stringstream ss(input);
				ss >> value;
			}
			setPlotInExecution(value);
		}
		
		if(slotXML.hasOption("OutputPlot")){
			std::string input = slotXML.get("OutputPlot");
			input = bitpit::utils::trim(input);
			std::string temp = ".";
			if(!input.empty())	setOutputPlot(input);
			else			  	setOutputPlot(temp);
		}
	
	return;	
};

/*!
 * Plot infos from a XML bitpit::Config::section. The parameters available are
 * 
 * 1) Dual       - boolean to get straight what given by selection method or its exact dual
 * 2) Tolerance  - tolerance for detect proximity volume in which perform mapping 
 * 3) Files	     - list by their filepaths of external geometries to be mapped on target one 
 * 
 * 
 * Geometry is mandatorily passed through ports. 
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionByMapping::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));

	slotXML.set("Topology", m_topo);
	
	int value = m_dual;
	slotXML.set("Dual", std::to_string(value));
	
	
	if(m_tolerance != 1.E-08){
		std::stringstream ss;
		ss<<std::scientific<<m_tolerance;
		slotXML.set("Tolerance", ss.str());
	}
	
	bitpit::Config::Section & filesXML = slotXML.addSection("Files");
	
	int counter = 0;
	for(auto & file : m_geolist){
		std::string name = "file"+std::to_string(counter);
		bitpit::Config::Section & local = filesXML.addSection(name);
		local.set("fullpath", file.first);
		std::string typetag = (FileType::_from_integral(file.second))._to_string(); 
		local.set("tag", typetag);
		++counter;
	}
	
		if(isPlotInExecution()){
			slotXML.set("PlotInExecution", std::to_string(1));
		}
		
		if(m_outputPlot != "."){
			slotXML.set("OutputPlot", m_outputPlot);
		}
	
	return;
};



//------------------------------------------------------------------------
//SELECTION	BY PID class 	******************************************
//------------------------------------------------------------------------

/*!
 * Basic Constructor
 */
SelectionByPID::SelectionByPID(){
	m_name = "MiMMiNO.SelectionByPID";
};

/*!
 * Custom Constructor
 * \param[in] pidlist list of pid to be included into selection
 * \param[in] target	pointer to taret geometry
 */
SelectionByPID::SelectionByPID(shivector1D & pidlist, MimmoObject * target){
	m_name = "MiMMiNO.SelectionByPID";
	setGeometry(target);
	setPID(pidlist);
};

/*!
 * Destructor
 */
SelectionByPID::~SelectionByPID(){};

/*!
 * Copy constructor
 */
SelectionByPID::SelectionByPID(const SelectionByPID & other){
	*this = other;
};

/*!
 * Copy Operator
 */
SelectionByPID & SelectionByPID::operator=(const SelectionByPID & other){
	*(static_cast<GenericSelection * >(this)) = *(static_cast<const GenericSelection * >(&other));
	m_activePID = other.m_activePID;
	return *this;
};

/*!
 * It builds the input/output ports of the object
 */
void SelectionByPID::buildPorts(){

	bool built = true;

	//inheritance
	GenericSelection::buildPorts();

	//input
	built = (built && createPortIn<short, SelectionByPID>(this, &SelectionByPID::setPID, PortType::M_VALUESI, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::SHORT));
	built = (built && createPortIn<std::vector<short>, SelectionByPID>(this, &SelectionByPID::setPID, PortType::M_VECTORSI, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::SHORT));

	m_arePortsBuilt = built;
};


/*!
 * Return all available PID for current geometry linked. Void list means not linked target geometry
 */
shivector1D 	SelectionByPID::getPID(){
	shivector1D result(m_activePID.size());
	int counter;
	for(auto && val : m_activePID){
		result[counter] = val.first;
		++counter;
	}
	return 	result;
};
/*!
 * Return active /not active for selection PIDs for current geometry linked. 
 * \param[in] active boolean to get active-true, not active false PID for selection
 */
shivector1D	SelectionByPID::getActivePID(bool active){
	shivector1D result;
	for(auto && val : m_activePID){
		if(val.second == active)	result.push_back(val.first);
	}
	return(result);
};


/*!
 * Set pointer to your target geometry. Reimplemented from mimmo::BaseManipulation::setGeometry()
 */

void 	SelectionByPID::setGeometry(MimmoObject * target ){
	if(target == NULL) return;
	m_geometry = target;
	
	std::unordered_set<short> & pids = 	target->getPIDTypeList();
	std::unordered_set<short>::iterator it, itE = pids.end();
	
	for(it = pids.begin(); it!=itE; ++it){
		m_activePID.insert(std::make_pair(*it,false));
	}
};

/*!
 * Activate flagged PID i. If i<0, activates all PIDs available.
 * Modifications are available after applying them with syncPIDList method; 
 * \param[in] i PID to be activated.
 */
void 	SelectionByPID::setPID(short i){
	if(m_setPID.count(-1) >0 || (!m_setPID.empty() && i==-1))	m_setPID.clear();
	m_setPID.insert(i);
	
};

/*!
 * Activate a list of flagged PID i. SelectionByPID::setPID(short i).
 * Modifications are available after applying them with syncPIDList method; 
 * \param[in] list	list of PID to be activated
 */
void 	SelectionByPID::setPID(shivector1D list){
	for(auto && index : list){
		setPID(index);
	}
};
/*!
 * Deactivate flagged PID i. If i<0, deactivates all PIDs available.
 * if i > 0 but does not exist in PID available list, do nothing
 * Modifications are available after applying them with syncPIDList method; 
 * \param[in] i PID to be deactivated.
 */

void 	SelectionByPID::removePID(short i){
	if(i>0){
			if(m_setPID.count(i) >0) m_setPID.erase(i);
	}else{
		m_setPID.clear();
	}
};

/*!
 * Deactivate a list of flagged PID i. SelectionByPID::removePID(short i = -1).
 * Modifications are available after applying them with syncPIDList method; 
 * \param[in] list	list of PID to be deactivated
 */
void 	SelectionByPID::removePID(shivector1D list){
	for(auto && index : list){
		removePID(index);
	}
};

/*!
 * Clear your class
 */
void SelectionByPID::clear(){
	m_subpatch.reset(nullptr);
	m_activePID.clear();
	BaseManipulation::clear();
};

/*!
 * Extract portion of target geometry which are enough near to the external geometry provided
 * \return ids of cell of target tesselation extracted 
 */
livector1D SelectionByPID::extractSelection(){
	livector1D result;
	std::set<long> extraction;
	
	syncPIDList();
	shivector1D pids = getActivePID();
	
	for(auto && pid : pids){
		livector1D local = getGeometry()->extractPIDCells(pid);
		extraction.insert(local.begin(),local.end());
	}
	
	//check if dual selection is triggered
	
	if(m_dual){
		livector1D totID = getGeometry()->getMapCell();
		result.resize(totID.size() - extraction.size());
		if(result.size() == 0) return result;
		
		std::sort(totID.begin(), totID.end());
		int counter = 0;
		auto tot_it  = totID.begin();
		auto cand_it = extraction.begin();
		while(tot_it != totID.end()){
			long val = *tot_it;
			if (cand_it == extraction.end() || val != *cand_it) {
				result[counter] = val;
				++counter;
			} else {
				++cand_it;
			}
			++tot_it;
		}
	}else{
		result.insert(result.end(), extraction.begin(), extraction.end());
	}	
	return	result;
};

/*!
 * Checks & synchronizes user given pid list w/ current pid list available by linked geometry.
 * Activate positive matching PIDs and erase negative matches from user list. 
 * if geometry is not currently linked does nothing.
 */
void SelectionByPID::syncPIDList(){
	if(getGeometry() == NULL)	return;
	
	if(m_setPID.count(-1) == 0){
		
		std::unordered_set<short>::iterator itU;
		shivector1D negative;
		for(itU = m_setPID.begin(); itU != m_setPID.end(); ++itU){
			short value = *itU;
			if(m_activePID.count(value) > 0)	m_activePID[value] = true;
			else	negative.push_back(value);
		}
	
		for(auto &val: negative)	m_setPID.erase(val);

		
	}else{
		m_setPID.clear();
		for(auto &val: m_activePID ){	
			m_setPID.insert(val.first);
			val.second = true;
			
		}
	}	
}


/*!
* Get infos from a XML bitpit::Config::section. The parameters available are
* 
* 1) Dual       - boolean to get straight what given by selection method or its exact dual
* 2) nPID		- number of PID to be selected
* 3) PID   		- set PID to select. -1 select all available PID
* 
* Geometry is mandatorily passed through ports. 
* 
* \param[in] slotXML 	bitpit::Config::Section of XML file
* \param[in] name   name associated to the slot
*/
void SelectionByPID::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	//start absorbing
	if(slotXML.hasOption("Dual")){
		std::string input = slotXML.get("Dual");
		input = bitpit::utils::trim(input);
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setDual(value);
	}
	
	int nPID = 0;
	if(slotXML.hasOption("nPID")){
		std::string input = slotXML.get("nPID");
		input = bitpit::utils::trim(input);
		nPID = 0;
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>nPID;
			std::max(0, nPID);
		}
	}
	shivector1D pidlist(nPID, -1);
	
	if(slotXML.hasOption("PID")){
		std::string input = slotXML.get("PID");
		input = bitpit::utils::trim(input);
		if(!input.empty()){
			std::stringstream ss(input);
			for(int i=0; i<nPID; ++i){
				ss>>pidlist[i];
			}	
		}
	}
	
	setPID(pidlist);
	
	
		if(slotXML.hasOption("PlotInExecution")){
			std::string input = slotXML.get("PlotInExecution");
			input = bitpit::utils::trim(input);
			bool value = false;
			if(!input.empty()){
				std::stringstream ss(input);
				ss >> value;
			}
			setPlotInExecution(value);
		}
		
		if(slotXML.hasOption("OutputPlot")){
			std::string input = slotXML.get("OutputPlot");
			input = bitpit::utils::trim(input);
			std::string temp = ".";
			if(!input.empty())	setOutputPlot(input);
			else			  	setOutputPlot(temp);
		}
	
	return;	
};

/*!
 * Plot infos from a XML bitpit::Config::section. The parameters available are
 * 
 * 1) Dual       - boolean to get straight what given by selection method or its exact dual
 * 2) nPID		- number of PID to be selected
 * 3) PID   		- set PID to select. -1 select all available PID
 * 
 * 
 * Geometry is mandatorily passed through ports. 
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionByPID::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
	int value = m_dual;
	slotXML.set("Dual", std::to_string(value));
	
	
	shivector1D selected = getActivePID(true);
	int size = selected.size();
	
	
	if(size != 0){
		slotXML.set("nPID", std::to_string(size));
		
		std::stringstream ss;
		for(int i=0; i<size; ++i){
			ss<<selected[i]<<'\t';
		}
		
		slotXML.set("PID", ss.str());
	}
	
		if(isPlotInExecution()){
			slotXML.set("PlotInExecution", std::to_string(1));
		}
		
		if(m_outputPlot != "."){
			slotXML.set("OutputPlot", m_outputPlot);
		}
	
	return;
};



//------------------------------------------------------------------------
//SELECTION BY BOX WITH SCALAR class    **********************************
//------------------------------------------------------------------------

/*!
 * Basic Constructor
 */
SelectionByBoxWithScalar::SelectionByBoxWithScalar(){
    m_name = "MiMMiNO.SelectionByBoxWithScalar";
};

/*!
 * Custom Constructor
 * \param[in] origin of the box->baricenter
 * \param[in] span   of the box, width/height/depth
 * \param[in] target    pointer to a target geometry
 */
SelectionByBoxWithScalar::SelectionByBoxWithScalar(darray3E origin, darray3E span, MimmoObject * target){
    m_name = "MiMMiNO.SelectionByBoxWithScalar";
    SelectionByBox(origin, span, target);
};

/*!
 * Destructor
 */
SelectionByBoxWithScalar::~SelectionByBoxWithScalar(){};

/*!
 * Copy Constructor
 */
SelectionByBoxWithScalar::SelectionByBoxWithScalar(const SelectionByBoxWithScalar & other){
    *this = other;
};

/*!
 * Copy operator
 */
SelectionByBoxWithScalar & SelectionByBoxWithScalar::operator=(const SelectionByBoxWithScalar & other){
    *(static_cast<SelectionByBox * >(this)) = *(static_cast<const SelectionByBox *>(&other));
    return *this;
};

/*!
 * It builds the input/output ports of the object
 */
void SelectionByBoxWithScalar::buildPorts(){

    bool built = true;

    //inheritance
    SelectionByBox::buildPorts();

    //input
    built = (built && createPortIn<dvector1D, SelectionByBoxWithScalar>(this, &SelectionByBoxWithScalar::setField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));

    //output
    built = (built && createPortOut<dvector1D, SelectionByBoxWithScalar>(this, &SelectionByBoxWithScalar::getField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));

    m_arePortsBuilt = built;
};

/*!
 * Clear content of the class
 */
void SelectionByBoxWithScalar::clear(){
    m_field.clear();
    SelectionByBox::clear();
};

/*! It sets the starting scalar field attached to the whole patch.
 * \param[in] field Scalar field.
 */
void
SelectionByBoxWithScalar::setField(dvector1D field){
    m_field = field;
}

/*! It gets the scalar field attached to the extracted patch.
 * \return Scalar field.
 */
dvector1D
SelectionByBoxWithScalar::getField(){
    return (m_field);
}


/*!
 * Execute your object. A selection is extracted and trasferred in
 * an indipendent MimmoObject structure pointed by m_subpatch member.
 * The extracted field attached to the selection is built starting from the
 * intial whole scalar field given as input and stored in member
 * m_field (modified after the execution).
 */
void SelectionByBoxWithScalar::execute(){

    SelectionByBox::execute();

    if (m_field.size() != 0){
        m_field.resize(getGeometry()->getNVertex(), 0.0);
        bitpit::PiercedVector<bitpit::Vertex> vertices = m_subpatch->getVertices();
        dvector1D field_tmp(vertices.size());
        for (auto vertex : vertices){
            field_tmp[m_subpatch->getMapDataInv(vertex.getId())] = m_field[getGeometry()->getMapDataInv(vertex.getId())];
        }
        m_field = field_tmp;
    }
}


/*!
 * Plot optional result of the class in execution, that is the selected patch
 * as standard vtk unstructured grid and the related scalar field.
 */
void SelectionByBoxWithScalar::plotOptionalResults(){
    if(getPatch()->isEmpty()) return;
    std::string name = m_name + "_" + std::to_string(getClassCounter()) +  "_Patch";

    getPatch()->getPatch()->getVTK().addData("field", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, m_field);

    getPatch()->getPatch()->write(name);
}

