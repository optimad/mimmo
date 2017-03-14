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

namespace mimmo{

/*!
 * Default constructor of SpecularPoints
*/
SpecularPoints::SpecularPoints(){
	m_name = "MiMMO.SpecularPoints";
	m_plane.fill(0.0);
    m_insideout = false;
    m_force = true;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
SpecularPoints::SpecularPoints(const bitpit::Config::Section & rootXML){
	
	m_name = "MiMMO.SpecularPoints";
	m_plane.fill(0.0);
	m_insideout = false;
	m_force = true;
	
	std::string fallback_name = "ClassNONE";	
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::trim(input);
	if(input == "MiMMO.SpecularPoints"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml MiMMO::SpecularPoints constructor. No valid xml data found"<<std::endl;
	};
}

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
 *  --> Absorbing data:
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>InsideOut</B>: boolean to get direction of clipping according to given plane
 * - <B>Plane</B>: section defining the plane's normal and a point belonging to it
 * - <B>Force</B>: boolean 0/1. If 1, force mirroring of points that lies on the plane.
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 * 
 * Coordinates and data attached are mandatorely passed through ports
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SpecularPoints::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

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
	
	if(slotXML.hasOption("Force")){
		std::string input = slotXML.get("Force");
		input = bitpit::utils::trim(input);
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setForce(value);
	}

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
	
	if(slotXML.hasSection("Plane")){
		const bitpit::Config::Section & planeXML = slotXML.getSection("Plane");
		
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
		
		setPlane(temp1, temp2);
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
 * Plot infos to a XML bitpit::Config::section. The parameters available are
 * 
 *  
 * --> Flushing data// how to write it on XML:
 * - <B>ClassName</B>: name of the class as "MiMMO.SpecularPoints"
 * - <B>Priority</B>: uint marking priority in multi-chain execution; 
 * - <B>Force</B>: boolean 0/1. If 1, force mirroring of points that lies on the plane.
 * - <B>InsideOut</B>: boolean 0/1 to get direction of clipping according to given plane
 * - <B>Plane</B>: section defining the plane's normal and a point belonging to it
 * 				<Plane>
 * 					<Point>	0.0 0.0 0.0 </Point>
 * 					<Normal> 0.0 1.0 0.0 </Normal>
 * 				</Plane>
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 * 
 * Coordinates and data attached are mandatorely passed through ports
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SpecularPoints::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	slotXML.set("ClassName", m_name);
	slotXML.set("Priority", std::to_string(getPriority()));
	
	int value = m_force;
	slotXML.set("Force", std::to_string(value));

	value = m_insideout;
	slotXML.set("InsideOut", std::to_string(value));

	
	{
		darray4E org = getPlane();
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
		
		bitpit::Config::Section & planeXML = slotXML.addSection("Plane");
		planeXML.set("Point",ss1.str());
		planeXML.set("Normal",ss2.str());
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

}

