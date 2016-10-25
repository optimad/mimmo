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
#include "SpecularPoints.hpp"

using namespace mimmo;

/*!Default constructor of SpecularPoints
*/
SpecularPoints::SpecularPoints(){
	m_name = "MiMMO.SpecularPoints";
	m_plane.fill(0.0);
    m_insideout = false;
    m_force = true;
};

/*!Default destructor of SpecularPoints
 */
SpecularPoints::~SpecularPoints(){};

/*!Copy constructor of SpecularPoints.
 */
SpecularPoints::SpecularPoints(const SpecularPoints & other){
	*this = other;
};

/*!Assignement operator of SpecularPoints.
 */
SpecularPoints & SpecularPoints::operator=(const SpecularPoints & other){
	*(static_cast<ProjectCloud*> (this)) = *(static_cast<const ProjectCloud*> (&other));
	m_plane = other.m_plane;
	m_scalar = other.m_scalar;
	m_vector = other.m_vector;
    m_insideout = other.m_insideout;
    m_force = other.m_force;
	return(*this);
};

/*! It builds the input/output ports of the object
 */
void
SpecularPoints::buildPorts(){
	bool built = true;

	built = (built && createPortIn<dvecarr3E, SpecularPoints>(this,&mimmo::SpecularPoints::setCoords, PortType::M_COORDS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dvecarr3E, SpecularPoints>(this, &mimmo::SpecularPoints::setVectorData, PortType::M_DISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dvector1D, SpecularPoints>(this, &mimmo::SpecularPoints::setScalarData, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<darray4E, SpecularPoints>(this, &mimmo::SpecularPoints::setPlane, PortType::M_PLANE, mimmo::pin::containerTAG::ARRAY4, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<MimmoObject*, SpecularPoints>(this, &mimmo::SpecularPoints::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortIn<bool, SpecularPoints>(this, &mimmo::SpecularPoints::setInsideOut, PortType::M_VALUEB, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::BOOL));
    built = (built && createPortIn<bool, SpecularPoints>(this, &mimmo::SpecularPoints::setForce, PortType::M_VALUEB2, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::BOOL));
	
	built = (built && createPortOut<dvecarr3E, SpecularPoints>(this, &mimmo::SpecularPoints::getCloudResult, PortType::M_COORDS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<dvecarr3E, SpecularPoints>(this, &mimmo::SpecularPoints::getCloudVectorData, PortType::M_DISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<dvector1D, SpecularPoints>(this, &mimmo::SpecularPoints::getCloudScalarData, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
	m_arePortsBuilt = built;
};

/*!
 * Returns the original scalar data field attached to cloud points
 */
dvector1D
SpecularPoints::getOriginalScalarData(){
	return(m_scalar);
};

/*!
 * Returns the original vector data field attached to cloud points
 */
dvecarr3E
SpecularPoints::getOriginalVectorData(){
	return(m_vector);
};

/*!
 * Returns the resulting scalar data field attached to mirrored cloud points
 */
dvector1D
SpecularPoints::getCloudScalarData(){
	return(m_scalarMirrored);
};

/*!
 * Returns the resulting vector data field attached to mirrored cloud points
 */
dvecarr3E
SpecularPoints::getCloudVectorData(){
	return(m_vectorMirrored);
};


/*!
 * Returns plane set up in the class, for mirroring 
 */
darray4E
SpecularPoints::getPlane(){
	return(m_plane);
};

/*!
 * Returns which half-space intercepeted by the plane is interested by mirroring. 
 * False represents the half-space where plane normal is directed, true the other one.
 */
bool
SpecularPoints::isInsideOut(){
    return m_insideout;
}

/*!
 * Returns if even the points belonging to the symmetry
 * plane are mirrored (i.e. duplicated).
 */
bool
SpecularPoints::isForce(){
    return m_force;
}

/*!
 * Set the original scalar data field attached to points that needs to be mirrored
 * \param[in] data scalar field;
 */
void
SpecularPoints::setScalarData(dvector1D data){
	m_scalar = data;
};

/*!
 * Set the original field data field attached to points that needs to be mirrored
 * \param[in] data scalar field;
 */

void 
SpecularPoints::setVectorData(dvecarr3E data){
	m_vector = data;
};


/*!
 * Set Plane for mirroring cloud points. All points not belonging to plane will be mirrored
 * \param[in] plane coefficients a,b,c,d of plane in its implicit form a*x+b*y+c*z+d = 0
 */
void
SpecularPoints::setPlane(darray4E plane){
	m_plane = plane;
};

/*!
 * Set Plane for mirroring cloud points. All points not belonging to plane will be mirrored
 * \param[in] origin points belonging to plane	
 * \param[in] normal plane normal
 *  
 */
void
SpecularPoints::setPlane(darray3E origin, darray3E normal){
	
	normal /= norm2(normal);
	double b = -1.0*dotProduct(origin, normal);
		
	m_plane[0] = normal[0];
	m_plane[1] = normal[1];
	m_plane[2] = normal[2];
	m_plane[3] = b;
};

/*!
 * Returns which half-space intercepeted by the plane is interested by mirroring. 
 * \param[in] flag false to select the half-space where plane normal is directed, true to select the other one.
 */
void
SpecularPoints::setInsideOut(bool flag){
	m_insideout = flag;
};

/*!
 * Set if even the points belonging to the symmetry
 * plane have to be mirrored (i.e. duplicated).
 * \param[in] flag true if the points on the symmetry plane have to be duplicated during mirroring.
 */
void
SpecularPoints::setForce(bool flag){
    m_force = flag;
}


/*!Execution command.Mirror the list of points linked, with data attached if any.
 * If a geometry is linked, project all resulting points on it.
 */
void
SpecularPoints::execute(){
	
	darray3E norm;
	double offset;
	for(int i=0; i<3; ++i)	norm[i] = m_plane[i];
	offset = m_plane[3];
	double normPlane = norm2(norm);
	if(normPlane < 1.E-18 || m_points.empty())	return;
	norm /= normPlane;
	bool project =!(getGeometry() == NULL || getGeometry()->isEmpty());

	
	//choosing margin for plane offset
    double margin = 1.e-12;
	if(project && getGeometry()->getType() == 1){
		double aTot = 0.0;
		int cellSize = getGeometry()->getNCells();
		bitpit::SurfaceKernel * tri = static_cast<bitpit::SurfaceKernel * >(getGeometry()->getPatch());
		for(auto &cell: tri->getCells()){
			aTot += tri->evalCellArea(cell.getId());
		}
		
		if(aTot > 0.0) {margin =  1.2*pow(aTot/((double)cellSize),0.5);}
		
	}else if(project && getGeometry()->getType() == 2){
		double vTot = 0.0;
		int cellSize = getGeometry()->getNCells();
		bitpit::VolumeKernel * tetra = static_cast<bitpit::VolumeKernel * >(getGeometry()->getPatch());
		for(auto &cell: tetra->getCells()){
			vTot += tetra->evalCellVolume(cell.getId());
		}
		
		if(vTot > 0.0) {margin =  1.2*pow(vTot/((double)cellSize),0.5);}
	}

	
	double sig = (1.0  - 2.0*((int)m_insideout));
	double distance;

	//mirroring.
	m_scalar.resize(m_points.size());
	m_vector.resize(m_points.size());
	
	int counterProj = m_points.size();
	int counterData = 0;
	
	m_proj = m_points;
	m_scalarMirrored = m_scalar;
	m_vectorMirrored = m_vector;
	
	m_proj.resize(2*counterProj);
	m_scalarMirrored.resize(2*counterProj);
	m_vectorMirrored.resize(2*counterProj);
	
	for(auto &val: m_points){
		distance = sig*(dotProduct(norm, val) + offset);
		if(distance > margin || m_force){
			m_proj[counterProj] = val - 2.0*distance*sig*norm;
			m_scalarMirrored[counterProj] = m_scalar[counterData];
			m_vectorMirrored[counterProj] = m_vector[counterData] -2.0*dotProduct(m_vector[counterData], sig*norm)*sig*norm;
			counterProj++;
		}
		counterData++;
	}	
	
	m_proj.resize(counterProj);
	m_scalarMirrored.resize(counterProj);
	m_vectorMirrored.resize(counterProj);
	
	if(project){
		if(!getGeometry()->isBvTreeBuilt())	getGeometry()->buildBvTree();
	
		//project points on surface.
		int counter = 0;
		for(auto &val : m_proj){
			m_proj[counter]= bvTreeUtils::projectPoint(&val, getGeometry()->getBvTree()); 
			++counter;
		}
	}	
	return;
};

/*!
 * Clear all content of the class
 */
void SpecularPoints::clear(){
	ProjectCloud::clear();
	m_scalar.clear();
	m_vector.clear();
	m_scalarMirrored.clear();
	m_vectorMirrored.clear();
	m_plane.fill(0.0);
	m_insideout = false;
};

/*!
 * Get infos from a XML bitpit::Config::section. The parameters available are
 * 
 *  
 * 1) Plane -> give mirror plane coefficients in inplicit form
 * 2) InsideOut-> half-space direction for mirroring
 * 
 * Coordinates and data attached are mandatorely passed through ports
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SpecularPoints::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){

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
	
	if(slotXML.hasOption("Plane")){
		std::string input = slotXML.get("Plane");
		input = bitpit::utils::trim(input);
		darray4E temp = {{0.0,0.0,0.0,0.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp[0]>>temp[1]>>temp[2]>>temp[3];
			setPlane(temp);
		}else{
			setPlane(temp);
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
 * Plot infos to a XML bitpit::Config::section. The parameters available are
 * 
 *  
 * 1) Plane -> give mirror plane coefficients in inplicit form
 * 2) InsideOut-> half-space direction for mirroring
 * 
 * Coordinates and data attached are mandatorely passed through ports
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SpecularPoints::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
	int value = m_insideout;
	slotXML.set("InsideOut", std::to_string(value));
	
	{
		darray4E org = getPlane();
		std::stringstream ss;
		ss<<std::scientific<<org[0]<<'\t'<<org[1]<<'\t'<<org[2]<<'\t'<<org[3];
		slotXML.set("Plane",ss.str());
	}
	
	if(isPlotInExecution()){
		slotXML.set("PlotInExecution", std::to_string(1));
	}
	
	if(m_outputPlot != "."){
		slotXML.set("OutputPlot", m_outputPlot);
	}
};

/*!
 * Plot as optional results the mirrored list of points with the updated 
 * data field associated to it
 */
void SpecularPoints::plotOptionalResults(){
		bitpit::VTKFormat codex = bitpit::VTKFormat::APPENDED;
		
		int size = m_proj.size();
		ivector1D conn(size);
		for(int i=0; i<size; i++){
			conn[i] = i;
		}
		std::string dir = "./";
		std::string file = m_name + "_" + std::to_string(getClassCounter());
		
		bitpit::VTKUnstructuredGrid vtk(dir, file, bitpit::VTKElementType::VERTEX);
		vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, m_proj) ;
		vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, conn) ;
		vtk.setDimensions(size, size);
		vtk.setCodex(codex);
	
		std::string scalarfield = "ScalarData";
		std::string vectorfield = "VectorData";
	
		vtk.addData( scalarfield, bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, m_scalarMirrored ) ;
		vtk.addData( vectorfield, bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, m_vectorMirrored ) ;

		vtk.write();
};



