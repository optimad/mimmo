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

#include "MeshSelection.hpp"
#include "levelSet.hpp"
#include <cstddef>
using namespace mimmo;

//------------------------------------------------------------------------
//SELECTION	BY PID class 	******************************************
//------------------------------------------------------------------------
// REGISTER_MANIPULATOR("MiMMO.SelectionByPID", "selectionbypid");

/*!
 * Basic Constructor
 */
SelectionByPID::SelectionByPID(){
	m_name = "MiMMO.SelectionByPID";
};

/*!
 * Custom Constructor
 * \param[in] pidlist list of pid to be included into selection
 * \param[in] target	pointer to taret geometry
 */
SelectionByPID::SelectionByPID(shivector1D & pidlist, MimmoObject * target){
	m_name = "MiMMO.SelectionByPID";
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
	int counter = 0;
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
* -->Absorbing Data
* Dual       : boolean to get straight what given by selection method or its exact dual
* nPID		: number of PID to be selected
* PID   		: set PID to select, separate by blank space. -1 select all available PID
* PlotInExecution : boolean 0/1 print optional results of the class.
* OutputPlot : target directory for optional results writing. 
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
 * --> Flushing data// how to write it on XML:
 * ClassName : name of the class as "MiMMO.Lattice"
 * ClassID	  : integer identifier of the class	
 * Dual       : boolean to get straight what given by selection method or its exact dual
 * nPID		: number of PID to be selected
 * PID   		: set PID to select, separate by blank space. -1 select all available PID
 * PlotInExecution : boolean 0/1 print optional results of the class.
 * OutputPlot : target directory for optional results writing.
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



