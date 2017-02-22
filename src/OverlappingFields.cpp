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

 #include "OverlappingFields.hpp"
 using namespace mimmo;
 using namespace mimmino;

 /*!
  * Constructor
  */
OverlapScalarFields::OverlapScalarFields(){
	m_name = "MiMMiNO.OverlapScalarFields";
	m_overlapCriterium = OverlapMethod::SUM;
}

/*!
 * Destructor
 */
OverlapScalarFields::~OverlapScalarFields(){
	clear();
}

/*!
 * Copy Constructor
 */
OverlapScalarFields::OverlapScalarFields(const OverlapScalarFields & other){
	*this = other;
}

/*!
 * Copy Operator
 */
OverlapScalarFields & OverlapScalarFields::operator=(const OverlapScalarFields & other){
	*(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(&other));
	m_overlapCriterium = other.m_overlapCriterium;
	m_originals = other.m_originals;
	return *this;
}

/*!
 * Gets actually used criterium for overlap regions of given fields
 */
OverlapMethod	OverlapScalarFields::getOverlapCriteriumENUM(){
	return m_overlapCriterium;
}

/*!
 * Gets actually used criterium for overlap regions of given fields
 */
int	OverlapScalarFields::getOverlapCriterium(){
	return static_cast<int>(m_overlapCriterium);
}

/*!
 * Return overlapped data pointed for a given subpatch mesh
 * \param[in] patch	pointer to a geometry subpatch
 * \return data of processed scalar field associated to the patch, if any. 
 */
dvector1D OverlapScalarFields::getResultData(MimmoObject * patch ){
	dvector1D result;
	if(m_results.count(patch) == 0) return result;
	return m_results[patch];
}

/*!
 * Return effective number of fields that would count after overlapping.
 */
int OverlapScalarFields::getNEffectiveFields(){
	return m_originals.size();
}

/*!
 * Return all linked fields in your class. Multiple fields associated to a unique
 * geometry that need to be overlapped are counted too.
 */
int OverlapScalarFields::getNLinkedFields(){
	
	int sum = 0;
	
	for(auto & val : m_originals){
		sum +=val.second.size();
	}
	return sum;
}

/*!
 * Return list of geometry pointers actually linked
 */
std::vector<MimmoObject * > 
OverlapScalarFields::whichGeometriesLinked(){
	
	std::vector<MimmoObject * > res(m_originals.size()); 
	int counter = 0;
	for(auto & val : m_originals){
		res[counter] = val.first;
		++counter;
	}
	return res;
}

/*!
 * Return overlapped fields after class execution as map of geometry pointers 
 * and data field pointers which refers to
 */
std::unordered_map<MimmoObject*, dvector1D* > 
OverlapScalarFields::getDataFieldMap(){
	
	std::unordered_map<MimmoObject*, dvector1D* > res; 
	for(auto & val : m_results){
		res[val.first] = &(val.second);
	}
	return res;
}

/*!
 * Return overlapped fields after class execution as a list of pairs having geometry pointers 
 * as key and data field pointers which refers to as argument
 */
std::vector<std::pair<MimmoObject*, dvector1D* > >
OverlapScalarFields::getDataFieldList(){
	
	std::vector<std::pair<MimmoObject*, dvector1D* > > res(m_results.size()); 
	int counter = 0;
	for(auto & val : m_results){
		res[counter] = std::make_pair(val.first, &(val.second));
		++counter;
	}
	return res;
}

/*!
 * Set overlap criterium for multi-fields overlapping. See OverlapMethod enum
 * Class default is OverlapMethod::SUM. Enum overloading
 * \param[in] funct	type of overlap method
 */
void 		OverlapScalarFields::setOverlapCriteriumENUM( OverlapMethod funct){
	setOverlapCriterium(static_cast<int>(funct));
	m_results.clear();
};

/*!
 * Set overlap criterium for multi-fields overlapping. See OverlapMethod enum
 * Class default is OverlapMethod::SUM. 
 * \param[in] funct	type of overlap method
 */
void 		OverlapScalarFields::setOverlapCriterium( int funct){
	if(funct <1 ||funct > 4)	return;
	m_overlapCriterium = static_cast<OverlapMethod>(funct);
	m_results.clear();
};


/*!
 * Append data fields of a geometry in a list , as a pair
 * (pointer to the geometry, pointer to the field it refers to).
 * \param[in] field	field to be inserted
 */
void		OverlapScalarFields::setAddDataField( std::pair<MimmoObject*, dvector1D*> field){
	
	if(field.first == NULL || field.second == NULL) return;
	if(field.first->isEmpty() || field.second->empty()) return;
	
	m_originals[field.first].push_back(field.second);
	m_results.clear();
};

/*!
 * Append data fields of a geometry in a list , as a map of
 * (pointer to the geometry, pointer to the field it refers to).
 * \param[in] fieldMap	map of fields to be inserted
 */
void		OverlapScalarFields::setDataFieldMap( std::unordered_map<MimmoObject*, dvector1D*> fieldMap){
	
	for(auto & val : fieldMap){
		setAddDataField(val);
	}
};

/*!
 * Append data fields of a geometry in a list , as a vector of pairs of
 * (pointer to the geometry, pointer to the field it refers to).
 * \param[in] fieldMap	list of fields to be inserted
 */
void		OverlapScalarFields::setDataFieldList(std::vector<std::pair<MimmoObject*, dvector1D*> > fieldList){
	
	for(auto & val : fieldList){
		setAddDataField(val);
	}
};

/*!
 * Remove a data field on the list by passing as key its pointer to geometry mesh
 */
void		OverlapScalarFields::removeData(MimmoObject * patch){
	std::unordered_map<MimmoObject *, std::vector<dvector1D *> >::iterator it = m_originals.find(patch);
	if(it != m_originals.end()){
		m_originals.erase(it);
	}
	m_results.clear();
};

/*!
 * Remove all data present in the list
 */
void		OverlapScalarFields::removeAllData(){
	m_originals.clear();
	m_results.clear();
};

/*!
 * Clear your class data and restore defaults settings
 */
void 		OverlapScalarFields::clear(){
	BaseManipulation::clear();
	removeAllData();
	m_overlapCriterium = OverlapMethod::SUM;
}

/*! 
 * It builds the input/output ports of the object
 */
void OverlapScalarFields::buildPorts(){
	
	bool built = true;
	
	//input
	built = (built && createPortIn<std::pair<MimmoObject *, dvector1D *>, OverlapScalarFields>(this, &mimmino::OverlapScalarFields::setAddDataField, PortType::M_PAIRSCAFIELD, mimmo::pin::containerTAG::PAIR, mimmo::pin::dataTAG::MIMMO_VECFLOAT_));
	built = (built && createPortIn<std::unordered_map<MimmoObject *, dvector1D *>,OverlapScalarFields>(this, &mimmino::OverlapScalarFields::setDataFieldMap, PortType::M_UMGEOSFD, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::MIMMO_VECFLOAT_));
	built = (built && createPortIn<int, OverlapScalarFields>(this, &OverlapScalarFields::setOverlapCriterium, CAMILOPortType::C_OVERLAPMTH, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::INT));
	built = (built && createPortIn<std::vector<std::pair<MimmoObject *, dvector1D *>>,OverlapScalarFields>(this, &mimmino::OverlapScalarFields::setDataFieldList, CAMILOPortType::C_VECPAIRSF, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::PAIRMIMMO_VECFLOAT_));
	
	//output
	built = (built && createPortOut<std::unordered_map<MimmoObject *, dvector1D *>,OverlapScalarFields>(this, &mimmino::OverlapScalarFields::getDataFieldMap, PortType::M_UMGEOSFD, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::MIMMO_VECFLOAT_));
	built = (built && createPortOut<int, OverlapScalarFields>(this,  &OverlapScalarFields::getOverlapCriterium, CAMILOPortType::C_OVERLAPMTH, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::INT));
	built = (built && createPortOut<std::vector<std::pair<MimmoObject *, dvector1D *>>,OverlapScalarFields>(this, &mimmino::OverlapScalarFields::getDataFieldList, CAMILOPortType::C_VECPAIRSF, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::PAIRMIMMO_VECFLOAT_));
	m_arePortsBuilt = built;
};




/*!
 * Plot overlapped data field of a target geometry. If geometry is not listed into the class, plot nothing.
 * If overlapped data field is not available into results (because of plotting before proper execution), returns
 * a vtu file with 0 field.
 * \param[in]	dir		output directory
 * \param[in]	name	output filename
 * \param[in]	flag    writing codex flag, false ascii, binary true
 * \param[in]	counter int counter, identifying your output name
 * \param[in]	geo     pointer to a geometry listed to the class.
 */
void	OverlapScalarFields::plotData(std::string dir, std::string name, bool flag, int counter, mimmo::MimmoObject * geo){
	
	if(m_originals.count(geo) == 0) return;
	dvector1D field;
	bool flagCalc = (m_results.count(geo) > 0);


	dvecarr3E points = geo->getVertexCoords();
    ivector2D connectivity;
    bitpit::VTKElementType cellType;
    if (geo->getType() != 3){
        connectivity = geo->getCompactConnectivity();
        cellType = desumeElement(geo->getType(), connectivity);
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

	if(flagCalc)	field = m_results[geo];
	else			field.resize(points.size(), 0.0);
	
    output.addData("field", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, field);

    std::vector<long> ids(points.size());
    long ID;
    for (auto vertex : geo->getVertices()){
        ID = vertex.getId();
        ids[geo->getMapDataInv(ID)] = ID;
    }

    output.addData("ID", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, ids);

	output.setCounter(counter);
	output.setCodex(bitpit::VTKFormat::APPENDED);
	if(!flag) output.setCodex(bitpit::VTKFormat::ASCII);

	output.write();
};

/*!
 * Plot all overlapped data field of all target geometries. Results will returned with the target filename and 
 * a consecutive numeration. If class is not executed and overlapped fields are not available plot nothing
 * \param[in]	dir		output directory
 * \param[in]	name	output unique filename
 * \param[in]	flag    writing codex flag, false ascii, binary true
 */
void	OverlapScalarFields::plotAllData(std::string dir, std::string name, bool flag){
	
	int counter = 0;
	for(auto & val : m_results){
		plotData(dir, name, flag, counter, val.first);
		++counter;
	}
};


/*! Overlap all your fields listed
 * \return overlapped fields in m_results member;
 */
void 	OverlapScalarFields::execute(){

	if(m_originals.empty())	return;
	m_results.clear();
	int size, listsize;
	dvector1D temp;
	dvector1D results;
	for(auto & obj : m_originals){
	
		size = obj.first->getPatch()->getVertexCount();
		//careful resizing
		results.resize(size, 0.0);
		for(auto &vec : obj.second)	vec->resize(size,0.0);
		
		listsize = obj.second.size();
		temp.resize(listsize);

		for(int i=0; i<size; ++i){
			for(int j=0; j<listsize; j++)	temp[j] = (*(obj.second[j]))[i];
				results[i] = overlapFields(temp);
		}
		
		m_results[obj.first] = results;
	}
	return;	
}

/*!
 * Plot Optional results in execution, that is the series of overlapped field .
 */
void 	OverlapScalarFields::plotOptionalResults(){
	std::string dir = m_outputPlot;
	std::string name = m_name;
	plotAllData(dir, name, true);
}

//PRIVATE MEMBER OF THE CLASS

/*!
 * Overlap concurrent value of different fields in the same node. Overlap Method is specified
 * in the class set
 *\param[in] locField list of value of concurrent field. If value is unique, simply assigns it 
 *\return	assigned value
 */
//DEVELOPERS REMIND if more overlap methods are added refer to this method to implement them
double 	OverlapScalarFields::overlapFields(dvector1D & locField){
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
 * Desume Element type from passed typeGeom and connectivity. Return undefined type for unexistent 
 * or unsupported element, or mixed element type connectivity. NEED TO BE MOVED IN MimmoObject
 */
bitpit::VTKElementType	OverlapScalarFields::desumeElement(int typeGeom, ivector2D & conn){
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
 * Get infos from a XML bitpit::Config::section. The parameters available are
 * 
 * 1) OverlapCriterium  - set how to treat fields in the overlapped region 1-MaxVal, 2-MinVal, 3-AverageVal, 4-Summing
 * 2) PlotInExecution - plot optional results in execution
 * 3) OutputPlot - path to store 3)
 * 
 * Fields and Geometry are mandatorily passed through ports. 
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void OverlapScalarFields::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	//start absorbing
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
 * 1) OverlapCriterium  - set how to treat fields in the overlapped region 1-MaxVal, 2-MinVal, 3-AverageVal, 4-Summing
 * 2) PlotInExecution - plot optional results in execution
 * 3) OutputPlot - path to store 3)
 * 
 * Fields and Geometry are mandatorily passed through ports. 
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void OverlapScalarFields::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
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



/*
 * ***************************************************
 * OVERLAPVECTORFIELDS CLASS IMPLEMENTATION
 *****************************************************
 */

/*!
 * Constructor
 */
OverlapVectorFields::OverlapVectorFields(){
	m_name = "MiMMiNO.OverlapVectorFields";
	m_overlapCriterium = OverlapMethod::SUM;
}

/*!
 * Destructor
 */
OverlapVectorFields::~OverlapVectorFields(){
	clear();
}

/*!
 * Copy Constructor
 */
OverlapVectorFields::OverlapVectorFields(const OverlapVectorFields & other){
	*this = other;
}

/*!
 * Copy Operator
 */
OverlapVectorFields & OverlapVectorFields::operator=(const OverlapVectorFields & other){
	*(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(&other));
	m_overlapCriterium = other.m_overlapCriterium;
	m_originals = other.m_originals;
	return *this;
}

/*!
 * Gets actually used criterium for overlap regions of given fields
 */
OverlapMethod	OverlapVectorFields::getOverlapCriteriumENUM(){
	return m_overlapCriterium;
}

/*!
 * Gets actually used criterium for overlap regions of given fields
 */
int	OverlapVectorFields::getOverlapCriterium(){
	return static_cast<int>(m_overlapCriterium);
}

/*!
 * Return overlapped data pointed for a given subpatch mesh
 * \param[in] patch	pointer to a geometry subpatch
 * \return data of processed vector field associated to the patch, if any. 
 */
dvecarr3E OverlapVectorFields::getResultData(MimmoObject * patch ){
	dvecarr3E result;
	if(m_results.count(patch) == 0) return result;
	return m_results[patch];
}

/*!
 * Return effective number of fields that would count after overlapping.
 */
int OverlapVectorFields::getNEffectiveFields(){
	return m_originals.size();
}

/*!
 * Return all linked fields in your class. Multiple fields associated to a unique
 * geometry that need to be overlapped are counted too.
 */
int OverlapVectorFields::getNLinkedFields(){
	
	int sum = 0;
	
	for(auto & val : m_originals){
		sum +=val.second.size();
	}
	return sum;
}

/*!
 * Return list of geometry pointers actually linked
 */
std::vector<MimmoObject * > 
OverlapVectorFields::whichGeometriesLinked(){
	
	std::vector<MimmoObject * > res(m_originals.size()); 
	int counter = 0;
	for(auto & val : m_originals){
		res[counter] = val.first;
		++counter;
	}
	return res;
}

/*!
 * Return overlapped fields after class execution as map of geometry pointers 
 * and data field pointers which refers to
 */
std::unordered_map<MimmoObject*, dvecarr3E* > 
OverlapVectorFields::getDataFieldMap(){
	
	std::unordered_map<MimmoObject*, dvecarr3E* > res; 
	for(auto & val : m_results){
		res[val.first] = &(val.second);
	}
	return res;
}

/*!
 * Return overlapped fields after class execution as a list of pairs having geometry pointers 
 * as key and data field pointers which refers to as argument
 */
std::vector<std::pair<MimmoObject*, dvecarr3E* > >
OverlapVectorFields::getDataFieldList(){
	
	std::vector<std::pair<MimmoObject*, dvecarr3E* > > res(m_results.size()); 
	int counter = 0;
	for(auto & val : m_results){
		res[counter] = std::make_pair(val.first, &(val.second));
		++counter;
	}
	return res;
}

/*!
 * Set overlap criterium for multi-fields overlapping. See OverlapMethod enum
 * Class default is OverlapMethod::SUM. Enum overloading
 * \param[in] funct	type of overlap method
 */
void 		OverlapVectorFields::setOverlapCriteriumENUM( OverlapMethod funct){
	setOverlapCriterium(static_cast<int>(funct));
	m_results.clear();
};

/*!
 * Set overlap criterium for multi-fields overlapping. See OverlapMethod enum
 * Class default is OverlapMethod::SUM. 
 * \param[in] funct	type of overlap method
 */
void 		OverlapVectorFields::setOverlapCriterium( int funct){
	if(funct <1 ||funct > 4)	return;
	m_overlapCriterium = static_cast<OverlapMethod>(funct);
	m_results.clear();
};


/*!
 * Append data fields of a geometry in a list , as a pair
 * (pointer to the geometry, pointer to the field it refers to).
 * \param[in] field	field to be inserted
 */
void		OverlapVectorFields::setAddDataField( std::pair<MimmoObject*, dvecarr3E*> field){
	
	if(field.first == NULL || field.second == NULL) return;
	if(field.first->isEmpty() || field.second->empty()) return;
	
	m_originals[field.first].push_back(field.second);
	m_results.clear();
};

/*!
 * Append data fields of a geometry in a list , as a map of
 * (pointer to the geometry, pointer to the field it refers to).
 * \param[in] fieldMap	map of fields to be inserted
 */
void		OverlapVectorFields::setDataFieldMap( std::unordered_map<MimmoObject*, dvecarr3E*> fieldMap){
	
	for(auto & val : fieldMap){
		setAddDataField(val);
	}
};

/*!
 * Append data fields of a geometry in a list , as a vector of pairs of
 * (pointer to the geometry, pointer to the field it refers to).
 * \param[in] fieldMap	list of fields to be inserted
 */
void		OverlapVectorFields::setDataFieldList(std::vector<std::pair<MimmoObject*, dvecarr3E*> > fieldList){
	
	for(auto & val : fieldList){
		setAddDataField(val);
	}
};

/*!
 * Remove a data field on the list by passing as key its pointer to geometry mesh
 */
void		OverlapVectorFields::removeData(MimmoObject * patch){
	std::unordered_map<MimmoObject *, std::vector<dvecarr3E *> >::iterator it = m_originals.find(patch);
	if(it != m_originals.end()){
		m_originals.erase(it);
	}
	m_results.clear();
};

/*!
 * Remove all data present in the list
 */
void		OverlapVectorFields::removeAllData(){
	m_originals.clear();
	m_results.clear();
};

/*!
 * Clear your class data and restore defaults settings
 */
void 		OverlapVectorFields::clear(){
	BaseManipulation::clear();
	removeAllData();
	m_overlapCriterium = OverlapMethod::SUM;
}

/*! 
 * It builds the input/output ports of the object
 */
void OverlapVectorFields::buildPorts(){
	
	bool built = true;
	
	//input
	built = (built && createPortIn<std::pair<MimmoObject *, dvecarr3E *>, OverlapVectorFields>(this, &mimmino::OverlapVectorFields::setAddDataField, PortType::M_PAIRVECFIELD, mimmo::pin::containerTAG::PAIR, mimmo::pin::dataTAG::MIMMO_VECARR3FLOAT_));
	built = (built && createPortIn<std::unordered_map<MimmoObject *, dvecarr3E *>,OverlapVectorFields>(this, &mimmino::OverlapVectorFields::setDataFieldMap, PortType::M_UMGEOVFD, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::MIMMO_VECARR3FLOAT_));
	built = (built && createPortIn<int, OverlapVectorFields>(this, &OverlapVectorFields::setOverlapCriterium, CAMILOPortType::C_OVERLAPMTH, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::INT));
	built = (built && createPortIn<std::vector<std::pair<MimmoObject *, dvecarr3E *>>,OverlapVectorFields>(this, &mimmino::OverlapVectorFields::setDataFieldList, CAMILOPortType::C_VECPAIRVF, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::PAIRMIMMO_VECARR3FLOAT_));
	
	//output
	built = (built && createPortOut<std::unordered_map<MimmoObject *, dvecarr3E *>,OverlapVectorFields>(this, &mimmino::OverlapVectorFields::getDataFieldMap, PortType::M_UMGEOVFD, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::MIMMO_VECARR3FLOAT_));
	built = (built && createPortOut<int, OverlapVectorFields>(this,  &OverlapVectorFields::getOverlapCriterium, CAMILOPortType::C_OVERLAPMTH, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::INT));
	built = (built && createPortOut<std::vector<std::pair<MimmoObject *, dvecarr3E *>>,OverlapVectorFields>(this, &mimmino::OverlapVectorFields::getDataFieldList, CAMILOPortType::C_VECPAIRVF, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::PAIRMIMMO_VECARR3FLOAT_));
	m_arePortsBuilt = built;
};




/*!
 * Plot overlapped data field of a target geometry. If geometry is not listed into the class, plot nothing.
 * If overlapped data field is not available into results (because of plotting before proper execution), returns
 * a vtu file with 0 field.
 * \param[in]	dir		output directory
 * \param[in]	name	output filename
 * \param[in]	flag    writing codex flag, false ascii, binary true
 * \param[in]	counter int counter, identifying your output name
 * \param[in]	geo     pointer to a geometry listed to the class.
 */
void	OverlapVectorFields::plotData(std::string dir, std::string name, bool flag, int counter, mimmo::MimmoObject * geo){
	
	if(m_originals.count(geo) == 0) return;
	dvecarr3E field;
	bool flagCalc = (m_results.count(geo) > 0);
	
	
	dvecarr3E points = geo->getVertexCoords();
	ivector2D connectivity;
	bitpit::VTKElementType cellType;
	if (geo->getType() != 3){
		connectivity = geo->getCompactConnectivity();
		cellType = desumeElement(geo->getType(), connectivity);
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
	
	if(flagCalc)	field = m_results[geo];
	else			field.resize(points.size(), {{0.0,0.0,0.0}});
	
	output.addData("vectorfield", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, field);
	
	std::vector<long> ids(points.size());
	long ID;
	for (auto vertex : geo->getVertices()){
		ID = vertex.getId();
		ids[geo->getMapDataInv(ID)] = ID;
	}
	
	output.addData("ID", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, ids);
	
	output.setCounter(counter);
	output.setCodex(bitpit::VTKFormat::APPENDED);
	if(!flag) output.setCodex(bitpit::VTKFormat::ASCII);
	
	output.write();
};

/*!
 * Plot all overlapped data field of all target geometries. Results will returned with the target filename and 
 * a consecutive numeration. If class is not executed and overlapped fields are not available plot nothing
 * \param[in]	dir		output directory
 * \param[in]	name	output unique filename
 * \param[in]	flag    writing codex flag, false ascii, binary true
 */
void	OverlapVectorFields::plotAllData(std::string dir, std::string name, bool flag){
	
	int counter = 0;
	for(auto & val : m_results){
		plotData(dir, name, flag, counter, val.first);
		++counter;
	}
};


/*! Overlap all your fields listed
 * \return overlapped fields in m_results member;
 */
void 	OverlapVectorFields::execute(){
	
	if(m_originals.empty())	return;
	m_results.clear();
	int size, listsize;
	dvecarr3E temp;
	dvecarr3E results;
	for(auto & obj : m_originals){
		
		size = obj.first->getPatch()->getVertexCount();
		//careful resizing
		results.resize(size, {{0.0,0.0,0.0}});
		for(auto &vec : obj.second)	vec->resize(size,{{0.0,0.0,0.0}});
		
		listsize = obj.second.size();
		temp.resize(listsize);
		
		for(int i=0; i<size; ++i){
			for(int j=0; j<listsize; j++)	temp[j] = (*(obj.second[j]))[i];
			results[i] = overlapFields(temp);
		}
		
		m_results[obj.first] = results;
	}
	return;	
}

/*!
 * Plot Optional results in execution, that is the series of overlapped field .
 */
void 	OverlapVectorFields::plotOptionalResults(){
	std::string dir = m_outputPlot;
	std::string name = m_name;
	plotAllData(dir, name, true);
}

//PRIVATE MEMBER OF THE CLASS

/*!
 * Overlap concurrent value of different fields in the same node. Overlap Method is specified
 * in the class set
 *\param[in] locField list of value of concurrent field. If value is unique, simply assigns it 
 *\return	assigned value
 */
//DEVELOPERS REMIND if more overlap methods are added refer to this method to implement them
darray3E 	OverlapVectorFields::overlapFields(dvecarr3E & locField){
	int size  = locField.size();
	if(size < 1) return {{0.0,0.0,0.0}};
	
	if(size ==1 )return locField[0];
	darray3E value;
	darray3E dir; dir.fill(0.0);
	double match;
	
	
	
	
	switch(m_overlapCriterium){
		case OverlapMethod::MAX :
			value = {{0.0,0.0,0.0}};
			for(auto && loc : locField){
				dir += loc;
			}
			
			dir /= norm2(dir);
			
			match = 1.e-18;
			for(auto && loc : locField){
				double dummy = dotProduct(loc, dir);
				match = std::fmax(match,dummy);
			}
			
			value = match*dir;	
			break;
			
		case OverlapMethod::MIN :
			
			value = {{0.0,0.0,0.0}};
			for(auto && loc : locField){
				dir += loc;
			}
			
			dir /= norm2(dir);
			
			match = 1.e18;
			for(auto && loc : locField){
				double dummy = dotProduct(loc, dir);
				match = std::fmin(match,dummy);
			}
			
			value = match*dir;	
			break;
			
		case OverlapMethod::AVERAGE :
			value = {{0.0,0.0,0.0}};
			for(auto && loc : locField){
				value += loc/double(size);
			}
			break;
			
		case OverlapMethod::SUM :
			value = {{0.0,0.0,0.0}};
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
 * Desume Element type from passed typeGeom and connectivity. Return undefined type for unexistent 
 * or unsupported element, or mixed element type connectivity. NEED TO BE MOVED IN MimmoObject
 */
bitpit::VTKElementType	OverlapVectorFields::desumeElement(int typeGeom, ivector2D & conn){
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
 * Get infos from a XML bitpit::Config::section. The parameters available are
 * 
 * 1) OverlapCriterium  - set how to treat fields in the overlapped region 1-MaxVal, 2-MinVal, 3-AverageVal, 4-Summing
 * 2) PlotInExecution - plot optional results in execution
 * 3) OutputPlot - path to store 3)
 * 
 * Fields and Geometry are mandatorily passed through ports. 
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void OverlapVectorFields::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	//start absorbing
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
 * 1) OverlapCriterium  - set how to treat fields in the overlapped region 1-MaxVal, 2-MinVal, 3-AverageVal, 4-Summing
 * 2) PlotInExecution - plot optional results in execution
 * 3) OutputPlot - path to store 3)
 * 
 * Fields and Geometry are mandatorily passed through ports. 
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void OverlapVectorFields::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
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
