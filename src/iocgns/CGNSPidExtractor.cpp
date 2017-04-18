/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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
\*---------------------------------------------------------------------------*/

#include "CGNSPidExtractor.hpp"
#include <unordered_map>

using namespace std;
using namespace bitpit;
namespace mimmo{


/*!Default constructor of CGNSPidExtractor.
 */
CGNSPidExtractor::CGNSPidExtractor(){
	m_name 		= "MiMMINO.CGNSPidExtractor";
	m_force = false;
	m_targetpid.insert(0);
}

/*!Default destructor of CGNSPidExtractor.
 */
CGNSPidExtractor::~CGNSPidExtractor(){};

/*!Copy constructor of CGNSPidExtractor.
 */
CGNSPidExtractor::CGNSPidExtractor(const CGNSPidExtractor & other){
	*this = other;
};

/*!Assignement operator of CGNSPidExtractor.
 */
CGNSPidExtractor & CGNSPidExtractor::operator=(const CGNSPidExtractor & other){
	
	*(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(&other));
	m_force = other.m_force;
	m_targetpid = other.m_targetpid;
	return *this;
};

/*! It builds the input/output ports of the object
 */
void
CGNSPidExtractor::buildPorts(){
	bool built = true;
	built = (built && createPortIn<MimmoObject*, CGNSPidExtractor>(this, &CGNSPidExtractor::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	
	built = (built && createPortOut<MimmoObject*, CGNSPidExtractor>(this, &CGNSPidExtractor::getPatch, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	m_arePortsBuilt = built;
};



/*!
 * Return the extracted surface by PID. If not valid PID is found, return an exact copy of the original target geometry.
 */

MimmoObject*
CGNSPidExtractor::getPatch(){
	return m_patch.get();
}

/*!
 * Return current active PID for selectioning in the class.
 */
std::set<short>
CGNSPidExtractor::whatPIDActive(){
	return m_targetpid;
}

/*!
 * Return if the class is set to reduce (in execution) the resulting PID-selected patch in a homogeneous straight triangulation, 
 * if eventually not.
 */
bool
CGNSPidExtractor::isForcedToTriangulate(){
	return m_force;
}

/*!
 * Set current active PID for selectioning in the class. If PID is not
 * supported by the linked target geoemetry, return a stand-alone copy of the target geometry itself in execution.
 * \param[in]	val PID holded by your target geometry
 */
void
CGNSPidExtractor::addPID( short val){
	m_targetpid.insert(val);
}

/*!
 * Set list of current active PID for selectioning in the class. If no PID in list is not
 * supported by the linked target geoemetry, return a stand-alone copy of the target geometry itself in execution.
 * \param[in]	vval PID list holded by your target geometry
 */
void
CGNSPidExtractor::setPID( std::vector<short int> vval){
	for(auto & val : vval) addPID(val);
}

/*!
 * Set the class to reduce (in execution) the resulting PID-selected patch in a homogeneous straight triangulation, 
 * if eventually not.
 * \param[in]	flag if true, try to triangulate homogeneously the selected patch.
 */
void
CGNSPidExtractor::setForcedToTriangulate( bool flag){
	m_force = flag;
}

/*!
 * Set current geometry to an external surface boundary mesh, from a mimmo CGNS reader.
 */
void
CGNSPidExtractor::setGeometry(MimmoObject * geo){
	if(geo == NULL)	return;
	if(geo->isEmpty())	return;
	if(geo->getType() != 1)	return;
	BaseManipulation::setGeometry(geo);
}


/*!Execution command.
 * Extract your PIDed patch.
 */
void
CGNSPidExtractor::execute(){
	
	livector1D extracted;
	
	for(auto & val: m_targetpid){
		livector1D temp= getGeometry()->extractPIDCells(val);
		extracted.insert(extracted.end(), temp.begin(), temp.end());
	}
	std::unique_ptr<MimmoObject> patchTemp(new MimmoObject(1));
	
	
	MimmoObject * mother = getGeometry();
	
	if(extracted.empty()){ 
		patchTemp->setHARDCopy(getGeometry());
		
	}else{
		
		for(auto &val: extracted){
			
			livector1D conn = mother->getCellConnectivity(val);
			bitpit::ElementInfo::Type eletype = mother->getPatch()->getCell(val).getType();
			patchTemp->addConnectedCell(conn, eletype, val);
		}
			
		{
			darray3E temp;
			livector2D tempConn = patchTemp->getConnectivity();
            //std::set<long> ordIndex;
            std::map<long,long> ordIndex;
				
			for(auto & val2: tempConn){
			    for (auto val22 : val2){
//                    ordIndex.insert(val22);
                    ordIndex[val22] = val22;
			    }
            //ordIndex.insert(val2.begin(), val2.end());
			}
				
			for(auto & val3: ordIndex){
				temp = mother->getVertexCoords(val3.second);
				patchTemp->addVertex(temp,val3.second);
			}
		}
		
	}	

	
	//check if patch is forced to be triangulated 
	if(m_force){
		
		long maxID, newID;
		
		auto orderedCellID = patchTemp->getCells().getIds(true);
		maxID = orderedCellID[(int)orderedCellID.size()-1];
		newID = maxID+1;
		bitpit::ElementInfo::Type eletype, eletri = bitpit::ElementInfo::Type::TRIANGLE;
		
		for(auto &idcell : orderedCellID){
			
			livector1D conn = patchTemp->getCellConnectivity(idcell);
			eletype = mother->getPatch()->getCell(idcell).getType();
			
			if(eletype == bitpit::ElementInfo::Type::QUAD){
				//create new triangle connectivity
				livector1D conn1(3), conn2(3);
				conn1[0] = conn[0];
				conn1[1] = conn[1];
				conn1[2] = conn[2];
				
				conn2[0] = conn[0];
				conn2[1] = conn[2];
				conn2[2] = conn[3];
				
				patchTemp->getPatch()->deleteCell(idcell);
				patchTemp->addConnectedCell(conn1, eletri, idcell);
				patchTemp->addConnectedCell(conn2, eletri, newID);
				++newID;
			}
		}
	}

	m_patch = std::move(patchTemp);
	mother=NULL;
};


/*!
 * Plot optional result of the class in execution, that is the selected patch
 * as standard vtk unstructured grid.
 */
void CGNSPidExtractor::plotOptionalResults(){
	if(getPatch()->isEmpty()) return;
	std::string name = m_name + "_" + std::to_string(getClassCounter()) +  "_Patch";
	getPatch()->getPatch()->write(name);
}

}
