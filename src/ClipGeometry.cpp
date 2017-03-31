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

namespace mimmo{

/*!Default constructor of ClipGeometry
*/
ClipGeometry::ClipGeometry(){
	m_name = "MiMMO.ClipGeometry";
	m_plane.fill(0.0);
	m_insideout = false;
	m_patch.reset(nullptr);
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ClipGeometry::ClipGeometry(const bitpit::Config::Section & rootXML){
	
	m_name = "MiMMO.ClipGeometry";
	m_plane.fill(0.0);
	m_insideout = false;
	m_patch.reset(nullptr);

	std::string fallback_name = "ClassNONE";	
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::trim(input);
	if(input == "MiMMO.ClipGeometry"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml MiMMO::ClipGeometry constructor. No valid xml data found"<<std::endl;
	};
}

/*!Default destructor of ClipGeometry
 */
ClipGeometry::~ClipGeometry(){};

/*!Copy constructor of ClipGeometry.
 */
ClipGeometry::ClipGeometry(const ClipGeometry & other):BaseManipulation(){
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

/*!
 * Set Plane for clipping by point and normal 
 * \param[in] origin point belonging to plane	
 * \param[in] normal plane normal
 *  
 */
void
ClipGeometry::setClipPlane(darray3E origin, darray3E normal){
	
	normal /= norm2(normal);
	double b = -1.0*dotProduct(origin, normal);
	
	m_plane[0] = normal[0];
	m_plane[1] = normal[1];
	m_plane[2] = normal[2];
	m_plane[3] = b;
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
	bitpit::ElementInfo::Type eltype;
	int PID;
	
	if (getGeometry()->getType() != 3){
		for(auto && idCell : extracted){
			
			bitpit::Cell & cell = tri->getCell(idCell);
			eltype = cell.getType();
			sizeCC = cell.getVertexCount();
			PID = cell.getPID();
			TT.resize(sizeCC);
			
			for(int i=0; i<sizeCC; ++i){
				idV = cell.getVertex(i);
				TT[i] = idV;
				
				if(!mapV.exists(idV))	temp->addVertex(tri->getVertexCoords(idV),idV);
			}
			temp->addConnectedCell(TT,eltype,short(PID), idCell);
			TT.clear();
		}
	}
	else{
		for(auto && idV : extracted){
			temp->addVertex(tri->getVertexCoords(idV),idV);
		}
	}
	
	m_patch = std::move(temp);
	tri = NULL;
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
	
	if(getGeometry()->getType() == 3){
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
	if(getClippedPatch() == NULL) return;
	if(getClippedPatch()->isEmpty()) return;
	
	dvecarr3E points = getClippedPatch()->getVertexCoords();
	ivector2D connectivity;
	bitpit::VTKElementType cellType;
	
	std::string dir = m_outputPlot;
	std::string name = m_name + "_Patch";
	
	
	if (getClippedPatch()->getType() != 3){
		connectivity = getClippedPatch()->getCompactConnectivity();
	}
	else{
		int np = points.size();
		connectivity.resize(np);
		for (int i=0; i<np; i++){
			connectivity[i].resize(1);
			connectivity[i][0] = i;
			
		}
	}
	cellType = getClippedPatch()->desumeElement(); 
	
	
	bitpit::VTKUnstructuredGrid output(dir,name,cellType);
	output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points);
	output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
	output.setDimensions(connectivity.size(), points.size());
	
	auto pids = getClippedPatch()->getCompactPID();
	if(pids.size() > 0) output.addData("PID", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, pids);
	
	output.setCounter(getClassCounter());
	output.setCodex(bitpit::VTKFormat::APPENDED);
	
	output.write();
}


/*!
 * Get infos from a XML bitpit::Config::section. The parameters available are
 * 
 *  --> Absorbing data:
 * 
 *   - <B>Priority</B>: uint marking priority in multi-chain execution; 
 *   - <B>InsideOut</B>: boolean to get direction of clipping according to given plane
 *   - <B>ClipPlane</B>: section defining the plane's normal and a point belonging to it
 *   - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 *   - <B>OutputPlot</B>: target directory for optional results writing.
 * 
 * Geometry is mandatorily passed through ports. 
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ClipGeometry::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
	
	//start absorbing
	BITPIT_UNUSED(name);
	
	if(slotXML.hasOption("Priority")){
		std::string input = slotXML.get("Priority");
		int value =0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss>>value;
		}
		setPriority(value);
	}; 
	
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
	
	if(slotXML.hasSection("ClipPlane")){
		const bitpit::Config::Section & planeXML = slotXML.getSection("ClipPlane");
		
		std::string input1 = planeXML.get("Point");
		std::string input2 = planeXML.get("Normal");
		input1 = bitpit::utils::trim(input1);
		input2 = bitpit::utils::trim(input2);
		
		darray3E temp1 = {{0.0,0.0,0.0}};
		darray3E temp2 = {{0.0,0.0,0.0}};
		
		if(!input1.empty()){
			std::stringstream ss(input1);
			ss>>temp1[0]>>temp1[1]>>temp1[2];
		}
		if(!input2.empty()){
			std::stringstream ss(input2);
			ss>>temp2[0]>>temp2[1]>>temp2[2];
		}
		
		setClipPlane(temp1, temp2);
	};	
	
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
 * * --> Flushing data// how to write it on XML:
 * - <B>ClassName</B>: name of the class as "MiMMO.ClipGeometry"
 * - <B>Priority</B>: uint marking priority in multi-chain execution; 
 * - <B>InsideOut</B>: boolean 0/1 to get direction of clipping according to given plane
 * - <B>ClipPlane</B>: section defining the plane's normal and a point belonging to it
 * 				<ClipPlane>
 * 					<Point>	0.0 0.0 0.0 </Point>
 * 					<Normal> 0.0 1.0 0.0 </Normal>
 * 				</ClipPlane>
 * - <B>PlotInExecution</B> : boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 * 
 * Geometry is mandatorily passed through ports.  
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ClipGeometry::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	BITPIT_UNUSED(name);
	
	slotXML.set("ClassName", m_name);
	slotXML.set("Priority", std::to_string(getPriority()));
	
	int value = m_insideout;
	slotXML.set("InsideOut", std::to_string(value));
	
	{
		darray4E org = getClipPlane();
		darray3E normal;
		darray3E point = {{0.0,0.0,0.0}};
		int imax = -1;
		double dum = 0.0;
		for(int i=0; i<3; ++i)	{
			normal[i] =org[i];
			if(abs(normal[i]) > dum) {
				imax = i;
				dum = abs(normal[i]);
			}
		}	
		if(imax != -1)	point[imax] = -1.0*org[3]/normal[imax];
		
		std::stringstream ss1, ss2;
		ss1<<std::scientific<<point[0]<<'\t'<<point[1]<<'\t'<<point[2];
		ss2<<std::scientific<<normal[0]<<'\t'<<normal[1]<<'\t'<<normal[2];
		
		bitpit::Config::Section & planeXML = slotXML.addSection("ClipPlane");
		planeXML.set("Point",ss1.str());
		planeXML.set("Normal",ss2.str());
	}
	
	if(isPlotInExecution()){
		slotXML.set("PlotInExecution", std::to_string(1));
	}
	
	if(m_outputPlot != "."){
		slotXML.set("OutputPlot", m_outputPlot);
	}
	
	return;
};

}


