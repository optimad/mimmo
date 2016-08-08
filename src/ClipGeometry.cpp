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
\*---------------------------------------------------------------------------*/
#include "ClipGeometry.hpp"

using namespace mimmo;

/*!Default constructor of ClipGeometry
*/
ClipGeometry::ClipGeometry(){
	m_name = "MiMMO.ClipGeometry";
	m_plane.fill(0.0);
	m_insideout = false;
	m_patch.reset(nullptr);
};

/*!Default destructor of ClipGeometry
 */
ClipGeometry::~ClipGeometry(){};

/*!Copy constructor of ClipGeometry.
 */
ClipGeometry::ClipGeometry(const ClipGeometry & other){
	*this = other;
};

/*!Assignement operator of ClipGeometry.
 */
ClipGeometry & ClipGeometry::operator=(const ClipGeometry & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_plane = other.m_plane;
	m_insideout = other.m_insideout;
	return(*this);
};

/*! It builds the input/output ports of the object
 */
void
ClipGeometry::buildPorts(){
	bool built = true;
	
	built = (built && createPortIn<darray4E, ClipGeometry>(this, &mimmo::ClipGeometry::setClipPlane, PortType::M_PLANE, mimmo::pin::containerTAG::ARRAY4, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<bool, ClipGeometry>(this, &mimmo::ClipGeometry::setInsideOut, PortType::M_VALUEB, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::BOOL));
	built = (built && createPortIn<MimmoObject*, ClipGeometry>(this, &mimmo::ClipGeometry::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	
	built = (built && createPortOut<MimmoObject*, ClipGeometry>(this, &mimmo::ClipGeometry::getClippedPatch, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	m_arePortsBuilt = built;
};

/*! Return direction for clipping. If false take all parts of target geometry lying on the half positive space 
 *  delimited by the plane (where plane normal pointing), true the exact opposite.
 * \return boolean flag 
 */
bool 
ClipGeometry::isInsideOut(){
	return(m_insideout);
};

/*!
 * Return copy of a pointer to the clipped geometry treated as an indipendent MimmoObject (owned by the class).
 */
MimmoObject * 
ClipGeometry::getClippedPatch(){
	return(m_patch.get());
};

/*!
 * Return coefficients of the clipping plane in its implicit form a*x+b*y+c*z+d=0
 */
darray4E
ClipGeometry::getClipPlane(){
	return m_plane;
};

/*!
 * Set coefficients of the clipping plane in its implicit form a*x+b*y+c*z+d=0
 * \param[in]	plane array of 4 coefficients a,b,c,d.
 */
void
ClipGeometry::setClipPlane(darray4E plane){
	m_plane = plane;
};

/*! Set direction for clipping. If false take all parts of target geometry lying on the half positive space 
 *  delimited by the plane (where plane normal pointing), true the exact opposite.
 * \param[in] flag boolean
 */
void 
ClipGeometry::setInsideOut(bool flag){
	m_insideout = flag;
};

/*!Execution command. Clip geometry and save result in m_patch member.
 */
void
ClipGeometry::execute(){

	if(getGeometry() == NULL || getGeometry()->isEmpty()) return;
	
	m_patch.reset(nullptr);
	
	livector1D extracted = clipPlane();
	if(extracted.empty()) return;
	
	//create your subpatch.
	std::unique_ptr<MimmoObject> temp(new MimmoObject(getGeometry()->getType()));
	bitpit::PatchKernel * tri = getGeometry()->getPatch();
	bitpit::PiercedVector<bitpit::Vertex> & mapV = temp->getPatch()->getVertices();
	
	livector1D TT;
	
	long idV;
	int sizeCC;
	int counter=0;
	bitpit::ElementInfo::Type eltype;
	
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
		counter++;
	}
	
	m_patch = std::move(temp);
	tri = NULL;
	return;
};

/*!
 * Return ID of elements composing geometry after clipping. Can be vertex IDs if the geometry is a PointCloud
 * or cell IDs if the geometry is a superficial or volumetric tessellation
 */
livector1D ClipGeometry::clipPlane(){
	
	livector1D result;
	darray3E norm;
	double offset;
	int counter;
	double sig = 1.0 - 2.0*(int)isInsideOut();
	long iD; 
	
	
	for(int i=0; i<3; ++i)	norm[i] = m_plane[i];
	offset = m_plane[3];
	
	double normPlane = norm2(norm);
	if(normPlane < 1.E-18)	return result;
	norm /= normPlane;
	
	bitpit::PatchKernel * tri = getGeometry()->getPatch();
	
	if(getGeometry()->getType() != 1 && getGeometry()->getType() != 2 ){
		counter = 0;
		result.resize(tri->getVertexCount());
		for(auto vert : tri->getVertices()){
			iD = vert.getId();
			if(sig*(dotProduct(norm, vert.getCoords()) + offset) >0)	{
				result[counter] = iD;
				++counter;
			}	
		}
		
	}else{
		counter = 0;
		result.resize(tri->getCellCount());
		for(auto cell : tri->getCells()){
			iD = cell.getId();
			if(sig*(dotProduct(norm, tri->evalCellCentroid(iD)) + offset) >0)	{
				result[counter] = iD;
				++counter;
			}	
		}
	}
	result.resize(counter);
	return result;
};


/*!
 * Plot optional result of the class in execution, that is the clipped geometry
 * as standard vtk unstructured grid.
 */
void ClipGeometry::plotOptionalResults(){
	if(getClippedPatch()->isEmpty() ) return;
	std::string name = m_name + "_" + std::to_string(getClassCounter()) +  "_Patch";
	getClippedPatch()->getPatch()->write(name);
}


/*!
 * Get infos from a XML bitpit::Config::section. The parameters available are
 * 
 * 1) InsideOut	- boolean to get direction of clipping according to given plane
 * 2) ClipPlane	- array of 4 doubles identifying coefficients of implicit plane equation a*x+b*y+c*z+d=0
 * 
 * Geometry is mandatorily passed through ports. 
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ClipGeometry::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	//start absorbing
	if(slotXML.hasOption("InsideOut")){
		std::string input = slotXML.get("InsideOut");
		input = bitpit::utils::trim(input);
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setInsideOut(value);
	}
	
	if(slotXML.hasOption("ClipPlane")){
		std::string input = slotXML.get("ClipPlane");
		input = bitpit::utils::trim(input);
		darray4E temp = {{0.0,0.0,0.0,0.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp[0]>>temp[1]>>temp[2]>>temp[3];
			setClipPlane(temp);
		}else{
			setClipPlane(temp);
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
 * 1) InsideOut	- boolean to get direction of clipping according to given plane
 * 2) ClipPlane	- array of 4 doubles identifying coefficients of implicit plane equation a*x+b*y+c*z+d=0
 * 
 * Geometry is mandatorily passed through ports.  
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ClipGeometry::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
	int value = m_insideout;
	slotXML.set("InsideOut", std::to_string(value));
	
	{
		darray4E org = getClipPlane();
		std::stringstream ss;
		ss<<std::scientific<<org[0]<<'\t'<<org[1]<<'\t'<<org[2]<<'\t'<<org[3];
		slotXML.set("ClipPlane",ss.str());
	}
	
	if(isPlotInExecution()){
		slotXML.set("PlotInExecution", std::to_string(1));
	}
	
	if(m_outputPlot != "."){
		slotXML.set("OutputPlot", m_outputPlot);
	}
	
	return;
};




