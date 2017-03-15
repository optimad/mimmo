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

#include "IOCloudPoints.hpp"
#include "Operators.hpp"
#include <fstream>

using namespace std;
namespace mimmo {

/*!
 * Default constructor of IOCloudPoints.
 * \param[in] readMode True if the object is in read mode, false if in Write mode.
 */
IOCloudPoints::IOCloudPoints(bool readMode){
	m_name 		= "MiMMO.IOCloudPoints";
	m_read 		= readMode;
	m_template 	= false;
	m_filename 	= m_name+"_source.dat";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
IOCloudPoints::IOCloudPoints(const bitpit::Config::Section & rootXML){
	
	m_name 		= "MiMMO.IOCloudPoints";
	m_template 	= false;
	m_filename 	= m_name+"_source.dat";
	m_read 		= true;
	
	std::string fallback_name = "ClassNONE";	
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::trim(input);

	std::string fallback_name2 = "1";	
	std::string input2 = rootXML.get("IOmode", fallback_name2);
	input2 = bitpit::utils::trim(input2);
	
	m_read = bool(std::atoi(input2.c_str()));

	if(input == "MiMMO.IOCloudPoints"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml MiMMO::IOCloudPoints constructor. No valid xml data found"<<std::endl;
	};
}


IOCloudPoints::~IOCloudPoints(){};

/*!
 * Copy constructor of IOCloudPoints.
 */
IOCloudPoints::IOCloudPoints(const IOCloudPoints & other){
	*this = other;
};

/*!Assignement operator of IOCloudPoints.
 */
IOCloudPoints & IOCloudPoints::operator=(const IOCloudPoints & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_read 		= other.m_read;
	m_filename 	= other.m_filename;
	m_template = other.m_template;
	
	//data structure is not copied
	return *this;
};


/*! 
 * It builds the input/output ports of the object
 */
void
IOCloudPoints::buildPorts(){
	bool built = true;
	
	built = (built && createPortIn<dvecarr3E, IOCloudPoints>(this, &mimmo::IOCloudPoints::setPoints, PortType::M_COORDS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dvecarr3E, IOCloudPoints>(this, &mimmo::IOCloudPoints::setVectorField, PortType::M_DISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dvector1D, IOCloudPoints>(this, &mimmo::IOCloudPoints::setScalarField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<livector1D, IOCloudPoints>(this, &mimmo::IOCloudPoints::setLabels, PortType::M_VECTORLI, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::LONG));
		
	built = (built && createPortOut<dvecarr3E, IOCloudPoints>(this, &mimmo::IOCloudPoints::getPoints, PortType::M_COORDS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<dvecarr3E, IOCloudPoints>(this, &mimmo::IOCloudPoints::getVectorField, PortType::M_DISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<dvector1D, IOCloudPoints>(this, &mimmo::IOCloudPoints::getScalarField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<livector1D, IOCloudPoints>(this, &mimmo::IOCloudPoints::getLabels, PortType::M_VECTORLI, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::LONG));
	
	m_arePortsBuilt = built;
}

/*!
 * Return the points coordinate stored in the class
 */
dvecarr3E
IOCloudPoints::getPoints(){
	return m_points; 
};

/*!
 * Return the vector field stored in the class. Field size is always checked and automatically 
 * fit to point list during the class execution 
 * 
 */
dvecarr3E
IOCloudPoints::getVectorField(){
	return m_vectorfield; 
};

/*!
 * Return the scalar field stored in the class. Field size is always checked and automatically 
 * fit to point list during the class execution 
 * 
 */
dvector1D
IOCloudPoints::getScalarField(){
	return m_scalarfield; 
};


/*!
 * Return the labels attached to displacements actually stored into the class.
 * Labels are always checked and automatically fit to point list during the class execution 
 */
livector1D
IOCloudPoints::getLabels(){
	return m_labels; 
};

/*!
 * Return if template option is active. This method is only meant in class working in Write mode
 */
bool
IOCloudPoints::isTemplate(){
	return m_template; 
};

/*!
 * It sets the name of the input/output file.
 * \param[in] filename absolute path to the I/O file.
 */
void
IOCloudPoints::setFilename(std::string filename){
	m_filename = filename;
};

/*!
 * Set the point list into the class 
 * The method is not active in Read mode.
 * \param[in] points list of 3D point
 */
void
IOCloudPoints::setPoints(dvecarr3E points){
	if(m_read) return;
	m_points = points;
};

/*!
 * It sets the labels attached to each point.
 * The method is not active in Read mode.
 * \param[in] labels list of label ids 
 */
void
IOCloudPoints::setLabels(livector1D labels){
	if(m_read) return;
	m_labels = labels;
};

/*!
 * It sets the vector field associated to point.
 * The method is not active in Read mode.
 * \param[in] vectorfield vector field 
 */
void
IOCloudPoints::setVectorField(dvecarr3E vectorfield){
	if(m_read) return;
	m_vectorfield = vectorfield;
};

/*!
 * It sets the scalar field associated to point.
 * The method is not active in Read mode.
 * \param[in] scalarfield scalar field 
 */
void
IOCloudPoints::setScalarField(dvector1D scalarfield){
	if(m_read) return;
	m_scalarfield = scalarfield;
};

/*!
 * Enables the template writing mode. The method is not active in Read mode.
 * \param[in] flag true to enable the template writing
 */
void
IOCloudPoints::setTemplate(bool flag){
	if(m_read) return;
	m_template = flag;
};

/*!
 * Clear all data stored in the class
 */
void
IOCloudPoints::clear(){
	m_vectorfield.clear();
	m_scalarfield.clear();
	m_points.clear();
	m_labels.clear();
	m_template 	= false;
	m_filename 	= m_name+"_source.dat";
}

/*!
 * Execution command. Read data from or Write data on linked filename  
 */
void
IOCloudPoints::execute(){
	
	if(m_read){
		read();
	}else{
		write();
	}
};

/*!
 * Get settings of the class from bitpit::Config::Section slot.
 * Displacements and labels are eventually passed only through ports. 
 * 
 * --> Absorbing data:
 * - <B>IOmode</B>: 1/0 enable Read and Write mode,respectively	
 * - <B>Priority</B>: uint marking priority in multi-chain execution; 
 * - <B>Filename</B>: path to your current file data
 * - <B>Template</B>: path to your current file data
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.  
 * 
 * \param[in] slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
void IOCloudPoints::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
	
	std::string input; 

	//checking IOmode
	
	if(slotXML.hasOption("IOmode")){
		input = slotXML.get("IOmode");
		bool value = true;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss>>value;
		}
		if (m_read != value){
			std::cout<< "Warning in class "<<m_name<<": cannot absorb xml parameters for class IOmode mismatching"<<std::endl;
			return;
		}
	}; 
	
	if(slotXML.hasOption("Priority")){
		input = slotXML.get("Priority");
		int value =0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss>>value;
		}
		setPriority(value);
	}; 

	if(slotXML.hasOption("Filename")){
		std::string input = slotXML.get("Filename");
		input = bitpit::utils::trim(input);
		setFilename(input);
	}; 
	
	
	if(slotXML.hasOption("Template")){
		std::string input = slotXML.get("Template");
		input = bitpit::utils::trim(input);
		bool temp = false;
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp;
		}
		setTemplate(temp);
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
	
}

/*!
 * Write settings of the class from bitpit::Config::Section slot.
 * Displacements and labels are eventually passed only through ports.
 * 
 * --> Flushing data// how to write it on XML:
 * - <B>ClassName</B>: name of the class as "MiMMO.IOCloudPoints"
 * - <B>IOmode</B>: 1/0 enable Read and Write mode,respectively
 * - <B>Priority</B>: uint marking priority in multi-chain execution; 
 * - <B>Filename</B>: path to your current file data
 * - <B>Template</B>: path to your current file data
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.  
 * 
 * \param[in] slotXML bitpit::Config::Section which writes to
 * \param[in] name   name associated to the slot
 */
void IOCloudPoints::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("IOmode", std::to_string(int(m_read)));
	slotXML.set("Priority", std::to_string(getPriority()));
	slotXML.set("Filename", m_filename);
	slotXML.set("Template", std::to_string(int(m_template)));
	if(isPlotInExecution()){
		slotXML.set("PlotInExecution", std::to_string(1));
	}
	
	if(m_outputPlot != "."){
		slotXML.set("OutputPlot", m_outputPlot);
	}
};	


/*!
 * Plot cloud point and store it in *.vtu file
 */
void IOCloudPoints::plotOptionalResults(){
	
	std::string path = m_outputPlot;
	std::string name = m_name +".Cloud_" + std::to_string(getClassCounter());
	bitpit::VTKUnstructuredGrid output( path, name, bitpit::VTKElementType::VERTEX);
	
	int size = m_points.size();
	ivector1D conn(size);
	for(int i=0; i<size; ++i)	conn[i] = i;
	output.setGeomData(bitpit::VTKUnstructuredField::POINTS, m_points);
	output.setGeomData(bitpit::VTKUnstructuredField::CONNECTIVITY, conn);
	output.setDimensions(size, size);
	
	std::string sfield = "scalarfield";
	dvector1D scafield = m_scalarfield;
	scafield.resize(size, 0.0);

	std::string vfield = "vectorfield";
	dvecarr3E vecfield = m_vectorfield;
	vecfield.resize(size, {{0.0,0.0,0.0}});
	
	output.addData( sfield, bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, scafield ) ;
	output.addData( vfield, bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, vecfield ) ;
	output.write() ;
}

/*!
 * Read displacement data from file
 */
void IOCloudPoints::read(){

	std::unordered_map<long, int> mapP;	
	std::ifstream reading;
	reading.open(m_filename.c_str());
	if(reading.is_open()){
		
		m_points.clear();
		m_labels.clear();
		m_scalarfield.clear();
		m_vectorfield.clear();
		
		std::string keyword, line;
		long label;
		darray3E dtrial;
		double temp;
		
		do{
			line.clear();
			keyword.clear();
			std::getline(reading,line);
			line = bitpit::utils::trim(line);
			std::stringstream ss(line);
			ss>>keyword; 
			keyword = bitpit::utils::trim(keyword);
			if(keyword == "$POINT")	{
				
				ss>>label;
				m_labels.push_back(label);
				
				dtrial.fill(0.0);
				ss>>dtrial[0]>>dtrial[1]>>dtrial[2];
				m_points.push_back(dtrial);
			}
			
		}while(!reading.eof());
		
		int counter = 0;
		for(auto &lab :m_labels){
			mapP[lab] = counter;
			++counter;
		}
		
		m_scalarfield.resize(m_points.size(),0.0);
		m_vectorfield.resize(m_points.size(),{{0.0,0.0,0.0}});
		
		reading.clear();
		reading.seekg(0, ios::beg);

		do{
			line.clear();
			keyword.clear();
			std::getline(reading,line);
			line = bitpit::utils::trim(line);
			std::stringstream ss(line);
			ss>>keyword; 
			keyword = bitpit::utils::trim(keyword);
			if(keyword == "$SCALARF")	{
				
				ss>>label>>temp;
				m_scalarfield[mapP[label]] = temp;
			}
			if(keyword == "$VECTORF")	{
				
				dtrial.fill(0.0);
				ss>>label>>dtrial[0]>>dtrial[1]>>dtrial[2];
				m_vectorfield[mapP[label]] = dtrial;
			}
			
		}while(!reading.eof());
		
	}else{
		std::cout<<"Error of "<<m_name<<" : cannot open "<<m_filename<< " requested. Exiting... "<<std::endl;
		exit(1);
	}
	
	reading.close();
};

/*!
 * Write displacement data to file
 */
void IOCloudPoints::write(){
		
	//assessing data;
	int sizelabels = m_labels.size();
	long maxlabel = 0;
	for(auto & val: m_labels){
		maxlabel = std::max(maxlabel,val);
	}
	
	int sizedispl = m_points.size();
	
	m_labels.resize(m_points.size(), -1);
	if(sizelabels < sizedispl){
		for(int i = sizelabels; i<sizedispl; ++i){
			m_labels[i] = maxlabel+1;
			++maxlabel;
		}
	}
	
	m_scalarfield.resize(m_points.size(), 0.0);
	m_vectorfield.resize(m_points.size(), {{0.0,0.0,0.0}});
	
	std::string filename;
	if(m_template){
		filename = "TEMPLATE_"+m_filename; 
	}else{
		filename = m_filename;
	}
	
	std::ofstream writing;
	writing.open(filename.c_str());
	std::string keyT1 = "{", keyT2 = "}";
	if(writing.is_open()){
		
		int counter = 0;
		for(auto & dd : m_points){
			writing<<"$POINT"<<'\t'<<m_labels[counter]<<'\t'<<dd[0]<<'\t'<<dd[1]<<'\t'<<dd[2]<<std::endl;
			++counter;
		}
		writing<<""<<std::endl;
		
		counter = 0;
		for(auto & dd : m_scalarfield){
			if(m_template){
				std::string str1 = keyT1+"s"+std::to_string(m_labels[counter])+keyT2;
				writing<<"$SCALARF"<<'\t'<<m_labels[counter]<<'\t'<<str1<<std::endl;
			}else{
				writing<<"$SCALARF"<<'\t'<<m_labels[counter]<<'\t'<<dd<<std::endl;
			}
			++counter;
		}
		writing<<""<<std::endl;
		
		counter = 0;
		for(auto & dd : m_vectorfield){
			if(m_template){
				std::string str1 = keyT1+"x"+std::to_string(m_labels[counter])+keyT2;
				std::string str2 = keyT1+"y"+std::to_string(m_labels[counter])+keyT2;
				std::string str3 = keyT1+"z"+std::to_string(m_labels[counter])+keyT2;
				
				writing<<"$VECTORF"<<'\t'<<m_labels[counter]<<'\t'<<str1<<'\t'<<str2<<'\t'<<str3<<std::endl;
			}else{
				writing<<"$VECTORF"<<'\t'<<m_labels[counter]<<'\t'<<dd[0]<<'\t'<<dd[1]<<'\t'<<dd[2]<<std::endl;
			}
			++counter;
		}
		writing<<""<<std::endl;
	}else{
		std::cout<<"Error of "<<m_name<<" : cannot open "<<m_filename<< " requested. Exiting... "<<std::endl;
		exit(1);
	}
	
	writing.close();
};


}

