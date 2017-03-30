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
namespace mimmo {

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
	if(getPatch() == NULL) return;
	if(getPatch()->isEmpty()) return;
	
	dvecarr3E points = getPatch()->getVertexCoords();
	ivector2D connectivity;
	bitpit::VTKElementType cellType;
	
	std::string dir = m_outputPlot;
	std::string name = m_name + "_Patch";
		
	
	if (getPatch()->getType() != 3){
		connectivity = getPatch()->getCompactConnectivity();
	}
	else{
		int np = points.size();
		connectivity.resize(np);
		for (int i=0; i<np; i++){
			connectivity[i].resize(1);
			connectivity[i][0] = i;
			
		}
	}
	cellType = desumeElement(getPatch()->getType(), connectivity); 
	
	
	bitpit::VTKUnstructuredGrid output(dir,name,cellType);
	output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points);
	output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
	output.setDimensions(connectivity.size(), points.size());
	
	auto pids = getPatch()->getCompactPID();
	if(pids.size() > 0) output.addData("PID", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, pids);
	
	output.setCounter(getClassCounter());
	output.setCodex(bitpit::VTKFormat::APPENDED);
	
	output.write();
}

/*!
 * Desume Element type from passed typeGeom and connectivity. Return undefined type for unexistent 
 * or unsupported element, or mixed element type connectivity. NEED TO BE MOVED IN MimmoObject
 */
bitpit::VTKElementType	GenericSelection::desumeElement(int typeGeom, ivector2D & conn){
	bitpit::VTKElementType resultUND = bitpit::VTKElementType::UNDEFINED;
	bitpit::VTKElementType result;
	
	switch(typeGeom){
		case	1:
			if(conn.empty()) 			return resultUND;
			if(conn[0].size() == 3)		result = bitpit::VTKElementType::TRIANGLE;
			if(conn[0].size() == 4)		result = bitpit::VTKElementType::QUAD;
			break;
		case	2:
			if(conn.empty()) 			return resultUND;
			if(conn[0].size() == 4)		result = bitpit::VTKElementType::TETRA;
			if(conn[0].size() == 8)		result = bitpit::VTKElementType::HEXAHEDRON;
		case	3:
			result = bitpit::VTKElementType::VERTEX;
			break;
		case	4:
			result = bitpit::VTKElementType::LINE;
			break;
		default : 
			result =resultUND;
			break;
	}
	
	return result;
};


}