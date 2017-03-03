/*---------------------------------------------------------------------------*\
 * 
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of MiMMO.
 *
 *  MiMMO is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  MiMMO is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with MiMMO. If not, see <http://www.gnu.org/licenses/>.
 *
 \ *---------------------------------------------------------------------------*/

 #include "ReconstructFields.hpp"
 using namespace mimmo;
 
REGISTER(BaseManipulation, ReconstructScalar,"MiMMO.ReconstructScalar");

 /*!
  * Constructor
  */
ReconstructScalar::ReconstructScalar(){
	m_name = "MiMMO.ReconstructScalar";
	m_overlapCriterium = OverlapMethod::MAX;
	buildPorts();
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ReconstructScalar::ReconstructScalar(const bitpit::Config::Section & rootXML){
	
	m_name = "MiMMO.ReconstructScalar";
	m_overlapCriterium = OverlapMethod::MAX;
	buildPorts();
	
	std::string fallback_name = "ClassNONE";	
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::trim(input);
	if(input == "MiMMO.ReconstructScalar"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml MiMMO::ReconstructScalar constructor. No valid xml data found"<<std::endl;
	};
}

/*!
 * Destructor
 */
ReconstructScalar::~ReconstructScalar(){
	clear();
}

/*!
 * Copy Constructor
 */
ReconstructScalar::ReconstructScalar(const ReconstructScalar & other){
	*this = other;
}

/*!
 * Copy Operator
 */
ReconstructScalar & ReconstructScalar::operator=(const ReconstructScalar & other){
	*(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(&other));
	m_overlapCriterium = other.m_overlapCriterium;
	m_subpatch = other.m_subpatch;
	m_result = other.m_result;
	return *this;
}

/*!
 * Gets actually used criterium for overlap regions of given fields
 */
OverlapMethod	ReconstructScalar::getOverlapCriteriumENUM(){
	return m_overlapCriterium;
}

/*!
 * Gets actually used criterium for overlap regions of given fields
 */
int	ReconstructScalar::getOverlapCriterium(){
	return static_cast<int>(m_overlapCriterium);
}

/*!
 * Return data pointed for a given subpatch mesh
 * \param[in] patch	pointer to a subpatch
 * \return data of scalar field associated to the patch, if any. 
 */
dvector1D ReconstructScalar::getData(MimmoObject * patch ){

	std::unordered_map<MimmoObject *, dvector1D *>::iterator it = m_subpatch.find(patch);

	dvector1D result;
	if(it == m_subpatch.end())	return result;

	return *(it->second);
}

/*!
 * Return number of fields data actually set in your class
 */
int ReconstructScalar::getNData(){

	return m_subpatch.size();
}

/*!
 * Return your result field
 */
dvector1D	ReconstructScalar::getResultField(){
	return(m_result);
}; 

/*!
 * Return list of sub-patch meshes actually stored in the class as a vector of copied pointers.
 */
std::vector<MimmoObject	*>	ReconstructScalar::whichSubMeshes(){
	std::vector<MimmoObject	*> result(getNData());
	int counter=0;
	for (auto && pairInd : m_subpatch){
		result[counter] = pairInd.first;
		++counter;
	}
	return result;
};

/*!
 * Set overlap criterium for multi-fields reconstruction. See OverlapMethod enum
 * Class default is OverlapMethod::MAX. Enum overloading
 * \param[in] funct	type of overlap method
 */
void 		ReconstructScalar::setOverlapCriteriumENUM( OverlapMethod funct){
	setOverlapCriterium(static_cast<int>(funct));
};

/*!
 * Set overlap criterium for multi-fields reconstruction. See OverlapMethod enum
 * Class default is OverlapMethod::MAX. 
 * \param[in] funct	type of overlap method
 */
void 		ReconstructScalar::setOverlapCriterium( int funct){
	if(funct <1 ||funct > 4)	return;
	m_overlapCriterium = static_cast<OverlapMethod>(funct);
};


/*!
 * Insert in the list data field of a sub-patch, as typedef pField
 * (pointer to the sub-patch mesh, pointer to the sub-patch field)
 * \param[in] field	sub-patch to be inserted
 */
void		ReconstructScalar::setData( pField  field){
	m_subpatch.insert(field);
};

/*!
 * Insert in the list data field of a sub-patch, as typedef pField
 * (pointer to the sub-patch mesh, pointer to the sub-patch field)
 * \param[in] fieldList	vector of sub-patch to be inserted
 */
void 		ReconstructScalar::setData(std::vector<pField>  fieldList){
	for(auto && data : fieldList){
		setData(data);
	}
};
/*!
 * Remove a data field on the list by passing as key its pointer to sub-patch mesh
 */
void		ReconstructScalar::removeData(MimmoObject * patch){
	std::unordered_map<MimmoObject *, dvector1D *>::iterator it = m_subpatch.find(patch);
	if(it != m_subpatch.end()){
		m_subpatch.erase(it);
	}
};

/*!
 * Remove all data present in the list
 */
void		ReconstructScalar::removeAllData(){
	m_subpatch.clear();
	m_result.clear();
};

/*!
 * Clear your class data and restore defaults settings
 */
void 		ReconstructScalar::clear(){
	BaseManipulation::clear();
	removeAllData();
	m_overlapCriterium = OverlapMethod::MAX;
}

/*!
 * Plot data (resulting field data) on vtu unstructured grid file
 * \param[in]	dir		output directory
 * \param[in]	name	output filename
 * \param[in]	flag    writing codex flag, false ascii, binary true
 * \param[in]	counter int counter, identifying your output name
 */
void	ReconstructScalar::plotData(std::string dir, std::string name, bool flag, int counter){
	
	if(getGeometry() == NULL || getGeometry()->isEmpty())	return;
	dvecarr3E points = getGeometry()->getVertexCoords();
    ivector2D connectivity;
    bitpit::VTKElementType cellType;
    if (getGeometry()->getType() != 3){
        connectivity = getGeometry()->getCompactConnectivity();
        cellType = desumeElement(getGeometry()->getType(), connectivity);
    }
    else{
        int np = points.size();
        connectivity.resize(np);
        for (int i=0; i<np; i++){
            connectivity[i].resize(1);
            connectivity[i][0] = i;
            cellType = bitpit::VTKElementType::VERTEX;
        }
    }
	bitpit::VTKUnstructuredGrid output(dir,name,cellType);
	output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points);
	output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
	output.setDimensions(connectivity.size(), points.size());

    output.addData("field", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, m_result);

    std::vector<long> ids(points.size());
    long ID;
    for (auto vertex : getGeometry()->getVertices()){
        ID = vertex.getId();
        ids[getGeometry()->getMapDataInv(ID)] = ID;
    }

    output.addData("ID", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, ids);

	output.setCounter(counter);
	output.setCodex(bitpit::VTKFormat::APPENDED);
	if(!flag) output.setCodex(bitpit::VTKFormat::ASCII);

	output.write();
};


/*! Reconstruct your field
 * \return reconstructed field in m_result member;
 */
void 	ReconstructScalar::execute(){

	if(getGeometry() == NULL)	return;
	
	liimap & vMotherMap = getGeometry()->getMapDataInv();

	m_result.resize(getGeometry()->getPatch()->getVertexCount(), 0.0);
	
	if(m_subpatch.empty()) return;
	
	std::unordered_map<long, dvector1D> map = checkOverlapping();

	for(auto && obj : map){
		m_result[vMotherMap[obj.first]] = overlapFields(obj.second);
		obj.second.clear();
	}
	
}

/*!
 * Plot Optional results in execution, that is the reconstructed scalar field .
 */
void 	ReconstructScalar::plotOptionalResults(){
	std::string dir = m_outputPlot;
	std::string name = m_name;
	plotData(dir, name, true, getClassCounter());
}

//PRIVATE MEMBER OF THE CLASS

/*!
 * Overlap concurrent value of different fields in the same node. Overlap Method is specified
 * in the class set
 *\param[in] locField list of value of concurrent field. If value is unique, simply assigns it 
 *\return	assigned value
 */
//DEVELOPERS REMIND if more overlap methods are added refer to this method to implement them
double 	ReconstructScalar::overlapFields(dvector1D & locField){
	int size  = locField.size();
	if(size < 1) return 0.0;

	double value = locField[0];
	if(size ==1 )return value;

	switch(m_overlapCriterium){
		case OverlapMethod::MAX :
			for(auto && loc : locField){
				value = std::fmax(loc,value);
			}
			break;
		case OverlapMethod::MIN :
			for(auto && loc : locField){
				value = std::fmin(loc,value);
			}
			break;
		case OverlapMethod::AVERAGE :
			value = 0.0;
			for(auto && loc : locField){
				value += loc/double(size);
			}
			break;
		case OverlapMethod::SUM :
			value = 0.0;
			for(auto && loc : locField){
				value += loc;
			}
			break;
		default : //never been reached
			break;
	}

	return value;
};

/*!
 * Check your sub-patch fields and reorder them in an unordered map carrying as
 * key the bitpit::PatchKernel ID of the mother mesh vertex and as value a vector
 * of doubles with all different field values concurring in that vertex
 *\return reordered map
 */
std::unordered_map<long, dvector1D>		ReconstructScalar::checkOverlapping(){

	std::unordered_map<long, dvector1D> result;
	int counter;

	for(auto && pairInd : m_subpatch){

		livector1D & vMap = pairInd.first->getMapData();
		counter = 0;
		dvector1D field = *(pairInd.second);
		for(int i=0; i<field.size(); ++i){
			result[vMap[counter]].push_back(field[i]);
			++counter;
		}
	}

	return(result);
};

/*!
 * Desume Element type from passed typeGeom and connectivity. Return undefined type for unexistent 
 * or unsupported element, or mixed element type connectivity. NEED TO BE MOVED IN MimmoObject
 */
bitpit::VTKElementType	ReconstructScalar::desumeElement(int typeGeom, ivector2D & conn){
	bitpit::VTKElementType result = bitpit::VTKElementType::UNDEFINED;
	if(conn.empty())	return	result;

	switch(typeGeom){
		case	1:
				if(conn[0].size() == 3)		result = bitpit::VTKElementType::TRIANGLE;
				if(conn[0].size() == 4)		result = bitpit::VTKElementType::QUAD;
			break;
		case	2:
				if(conn[0].size() == 4)		result = bitpit::VTKElementType::TETRA;
				if(conn[0].size() == 8)		result = bitpit::VTKElementType::HEXAHEDRON;
			break;
		default : //never been reached
			break;
	}

	return result;
};


/*! 
 * It builds the input/output ports of the object
 */
void ReconstructScalar::buildPorts(){

	bool built = true;

	//input
	built = (built && createPortIn<MimmoObject *, ReconstructScalar>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	built = (built && createPortIn<std::pair<MimmoObject *, dvector1D *>,ReconstructScalar>(this, &mimmo::ReconstructScalar::setData, PortType::M_PAIRSCAFIELD, mimmo::pin::containerTAG::PAIR, mimmo::pin::dataTAG::MIMMO_VECFLOAT_));
	built = (built && createPortIn<std::vector<std::pair<MimmoObject *, dvector1D *>>,ReconstructScalar>(this, &mimmo::ReconstructScalar::setData, PortType::M_VECPAIRSF, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::PAIRMIMMO_VECFLOAT_));

	//output
	built = (built && createPortOut<dvector1D, ReconstructScalar>(this, &ReconstructScalar::getResultField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<MimmoObject *, ReconstructScalar>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	m_arePortsBuilt = built;
};



/*!
 * Get infos from a XML bitpit::Config::section. The parameters available are
 * 
 *  --> Absorbing data:
 * Priority  : uint marking priority in multi-chain execution;
 * OverlapCriterium  : set how to treat fields in the overlapped region 1-MaxVal, 2-MinVal, 3-AverageVal, 4-Summing
 * PlotInExecution : plot optional results in execution
 * OutputPlot : path to optional results
 * PlotInExecution : boolean 0/1 print optional results of the class.
 * OutputPlot : target directory for optional results writing. 
 * 
 * Fields and Geometry are mandatorily passed through ports. 
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ReconstructScalar::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
	
	//start absorbing
	if(slotXML.hasOption("Priority")){
		std::string input = slotXML.get("Priority");
		int value =0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss>>value;
		}
		setPriority(value);
	}; 
	
	if(slotXML.hasOption("OverlapCriterium")){
		std::string input = slotXML.get("OverlapCriterium");
		input = bitpit::utils::trim(input);
		int value = 1;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
			value = std::min(std::max(1, value),4);
		}
		setOverlapCriterium(value);
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
 * 
 * --> Flushing data// how to write it on XML:
 * ClassName : name of the class as "MiMMO.ReconstructScalar"
 * Priority  : uint marking priority in multi-chain execution;
 * OverlapCriterium  : set how to treat fields in the overlapped region 1-MaxVal, 2-MinVal, 3-AverageVal, 4-Summing
 * PlotInExecution : plot optional results in execution
 * OutputPlot : path to optional results
 * PlotInExecution : boolean 0/1 print optional results of the class.
 * OutputPlot : target directory for optional results writing.
 * 
 * Fields and Geometry are mandatorily passed through ports. 
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ReconstructScalar::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("Priority", std::to_string(getPriority()));
	
	int value = static_cast<int>(m_overlapCriterium);
	slotXML.set("OverlapCriterium", std::to_string(value));
	
	if(isPlotInExecution()){
		slotXML.set("PlotInExecution", std::to_string(1));
	}
	
	if(m_outputPlot != "."){
		slotXML.set("OutputPlot", m_outputPlot);
	}
	
	return;
};


