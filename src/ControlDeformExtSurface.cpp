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
#include "ControlDeformExtSurface.hpp"
// #include <chrono>
// 
// using namespace std::chrono;

using namespace mimmo;

/*!Default constructor of ControlDeformExtSurface
*/
ControlDeformExtSurface::ControlDeformExtSurface(){
	m_name = "MiMMO.ControlDeformExtSurface";
	m_tolerance = 1.E-8;
	m_cellBackground = 20;
	m_allowedType.resize(2);
	m_allowedType[1].insert(FileType::STL);
	m_allowedType[1].insert(FileType::STVTU);
	m_allowedType[1].insert(FileType::SQVTU);
	m_allowedType[1].insert(FileType::NAS);
	
};

/*!Default destructor of ControlDeformExtSurface
 */
ControlDeformExtSurface::~ControlDeformExtSurface(){};

/*!Copy constructor of ControlDeformExtSurface.
 */
ControlDeformExtSurface::ControlDeformExtSurface(const ControlDeformExtSurface & other){
	*this = other;
};

/*!
 * Assignement operator of ControlDeformExtSurface. Create an exact copy of the class,
 * except for the deformation field referred to the target geometry.
 */
ControlDeformExtSurface & ControlDeformExtSurface::operator=(const ControlDeformExtSurface & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_geolist = other.m_geolist;
	m_allowedType = other.m_allowedType;
	m_cellBackground = other.m_cellBackground;
	m_tolerance = other.m_tolerance;
	//deformation and violation field are not copied
	return(*this);
};

/*! It builds the input/output ports of the object
 */
void
ControlDeformExtSurface::buildPorts(){
	bool built = true;
	
	built = (built && createPortIn<dvecarr3E, ControlDeformExtSurface>(this, &mimmo::ControlDeformExtSurface::setDefField, PortType::M_GDISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<MimmoObject*, ControlDeformExtSurface>(this, &mimmo::ControlDeformExtSurface::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	
	built = (built && createPortOut<double, ControlDeformExtSurface>(this, &mimmo::ControlDeformExtSurface::getViolation, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<dvector1D, ControlDeformExtSurface>(this, &mimmo::ControlDeformExtSurface::getViolationField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));	
	built = (built && createPortOut<std::pair<BaseManipulation*, double>, ControlDeformExtSurface>(this, &mimmo::ControlDeformExtSurface::getViolationPair, PortType::M_VIOLATION, mimmo::pin::containerTAG::PAIR, mimmo::pin::dataTAG::PAIRMIMMO_OBJFLOAT_));	
	m_arePortsBuilt = built;
};

/*! 
 * Return the value of violation of deformed geometry, after class execution. 
 *  If value is positive or at least zero, constraint violation occurs, penetrating or touching at least in one point the 
 *  constraint sourface. Otherwise, returning negative values means that no violation occurs 
 *  delimited by the plane (where plane normal pointing), true the exact opposite.
 * \return violation value
 */
double 
ControlDeformExtSurface::getViolation(){
	
	double result = -1.0E+18;
	for(auto & val : m_violationField){
		result = std::fmax(result, val);
	}
	
	return	result;
};


/*! 
 *  Return the value of violation of deformed geometry, after class execution. 
 *  A BaseManipulation object pointer, representing the sender of geometry to which violation is referred, 
 *  is attached to violation value; 
 *  See getViolation method doc for further details. 
 * \return std::pair structure reporting geometry sender pointer and violation value.
 */
std::pair<BaseManipulation*, double> 
ControlDeformExtSurface::getViolationPair(){
	
	//get map of Input ports of the class.
	std::map<short int, mimmo::PortIn*> mapPorts = getPortsIn();
	
	//get class who send geometry here - portID = 99 -> M_GEOM
	
	std::vector<BaseManipulation*> senders = mapPorts[99]->getLink();
	
	std::string etiq;
	if(senders.size() == 0){
		return	std::make_pair(this, getViolation());	
	}else{
		return	std::make_pair(senders[0], getViolation());
	}	
	
};



/*! 
 * Return the violation distances of each point of deformed geometry, on the geometry itself. The info is available after class execution. 
 *  If value is positive or at least zero, constraint violation occurs, penetrating or touching at least in one point the 
 *  constraint sourface. Otherwise, returning negative values means that no violation occurs 
 *  delimited by the plane (where plane normal pointing), true the exact opposite.
 * \return violation field values
 */
dvector1D 
ControlDeformExtSurface::getViolationField(){
	return(m_violationField);
};


/*!
 * Get Tolerance fixed in the class to consider a deformed body close 
 * to target constraints not touching/penetrating them.
 * \return 	fixed tolerance. Default is 1.E-12
 */
double
ControlDeformExtSurface::getToleranceWithinViolation(){
	return m_tolerance;
}

/*!
 * Get number of cells used for determining background grid spacing.
 * See ControlDeformExtSurface::setBackgroundDetails() method documentation;
 * \return 	numebr of cells. Default is 20; 
 */
int
ControlDeformExtSurface::getBackgroundDetails(){
	return m_cellBackground;
}

/*!
 * Set the deformative field associated to each point of the target geometry. 
 * Field resize occurs in execution, if point dimension between field and geoemetry does not match.
 * \param[in]	field of deformation
 */
void
ControlDeformExtSurface::setDefField(dvecarr3E field){
	m_defField.clear();
	m_violationField.clear();
	m_defField = field;
	m_violationField.resize(field.size(),-1.E+18);
};

/*!
 * Set link to target geometry for your selection. Reimplementation of 
 * GenericSelection::setGeometry();
 */
void	ControlDeformExtSurface::setGeometry( MimmoObject * target){
	
	if(target->getType() != 1){
		std::cout<<"ControlDeformExtSurface cannot support current geometry. It works only w/ 3D surface."<<std::endl;
		return;
	}
	m_geometry = target;
	
};


/*!
 * Set the tolerance to consider a deformed body penetrating target constraints
 * still within the non-violation regime.
 */
void 	ControlDeformExtSurface::setToleranceWithinViolation(double tol){
	m_tolerance = std::fmax(1.0E-18, tol);
}

/*!
 * Set the number of Cells NC necessary to determine the spacing of a background grids.
 * A background axis-aligned cartesian volume grid wrapping each constraint + deformed body is created to evaluate 
 * distances of each deformed point from constraint surface. Spacing of such grid is
 * calculated as dh = D/NC where dh is the spacing, D is the diagonal of grid bounding-box.  
 */
void 	ControlDeformExtSurface::setBackgroundDetails(int nCell){
	m_cellBackground = std::fmax(2, nCell);
}


/*!
 * Return the actual list of external geometry files selected as constraint to check your deformation
 */
const std::unordered_map<std::string, int> & ControlDeformExtSurface::getFiles() const{
	return	m_geolist;
}

/*!
 * Set a list of external geometry files as constraint to check your deformation violation
 * \param[in] files list external geometries to be read
 */
void	ControlDeformExtSurface::setFiles(std::unordered_map<std::string, int>  files){
	for(auto && val : files){
		addFile(val);
	}	
};

/*!
 * Add a new file to a list of external constraint geometry files 
 *\param[in] file of external geometry to be read
 */
void 	ControlDeformExtSurface::addFile(std::pair<std::string, int> file){
	int type = 1;
	if(m_allowedType[type].find(file.second) != m_allowedType[type].end()){									
		m_geolist.insert(file);
	}	
};

/*!
 * Remove an existent file to a list of external geometry files. If not in the list, do nothing 
 * \param[in] file to be removed from the list 
 */
void 	ControlDeformExtSurface::removeFile(std::string file){
	if(m_geolist.find(file) != m_geolist.end())	m_geolist.erase(file);
};

/*!
 * Empty your list of external constraint geometry files 
 */
void 	ControlDeformExtSurface::removeFiles(){
	m_geolist.clear();
};

/*!
 * Clear your class
 */
void ControlDeformExtSurface::clear(){
	removeFiles();
	m_defField.clear();
	m_violationField.clear();
	BaseManipulation::clear();
	m_tolerance = 1.E-12;
	m_cellBackground = 20;
};

/*!Execution command. Calculate violation value and store it in the class member m_violation
 */
void
ControlDeformExtSurface::execute(){
	
	MimmoObject * geo = getGeometry();
	if(geo == NULL || geo->isEmpty()) return;

	int nDFS = m_defField.size();
	m_defField.resize(geo->getNVertex(),darray3E{{0.0,0.0,0.0}});
	m_violationField.resize(nDFS,-1.E+18);
	
	dvector1D violationField;
	violationField.resize(nDFS);

	dvecarr3E points = geo->getVertexCoords();
	//evaluate deformable geometry barycenter
	darray3E geoBary = {{0.0,0.0,0.0}};
	double distBary = 0.0;
	
	for(auto & p: points)	geoBary += p;
	
	geoBary /= (double)geo->getNVertex();
	
	for(auto & p: points){
		distBary= std::fmax(distBary,norm2(p-geoBary));
	}	
	//adding deformation to points.
	points+= m_defField;

	//read external surfaces.
	std::vector<std::unique_ptr<MimmoGeometry> > extgeo;
	readGeometries(extgeo);
	
	if(extgeo.size() < 1)	return;

	//evaluating bounding box of deformation;
	darray3E bbMinDef, bbMaxDef, bbMin, bbMax;
	{
		int count2 = 0;
		for(auto & p : points){
			if(count2 == 0)	{
				bbMinDef = p;
				bbMaxDef = p;
			}else{
				for(int i=0; i<3; ++i ){
					bbMinDef[i] = std::fmin(bbMinDef[i], p[i]);
					bbMaxDef[i] = std::fmax(bbMaxDef[i], p[i]);
				}	
			}
			count2++;
		}
	}

	
	for(auto &gg : extgeo){

		MimmoObject * local = gg->getGeometry();
		if(!(local->isBvTreeBuilt()))	local->buildBvTree();
		
		double dist;
		double radius, radius_old;
		long id;
		darray3E normal;
		int count;
		
		double psign = 1.0;
		
		//get the whole bounding box of deformation + constraints
		bbMin = bbMinDef;
		bbMax = bbMaxDef;
		
		for(auto & p : local->getVertexCoords()){
			for(int i=0; i<3; ++i ){
				bbMin[i] = std::fmin(bbMin[i], p[i]);
				bbMax[i] = std::fmax(bbMax[i], p[i]);
			}	
		}
		//calculate dh as a tot part of bb diagonal	
		
		darray3E span = bbMax - bbMin;
		bbMin += -1.0*(0.05*span);
		span *= 1.1;
		double dh = norm2(span)/(double)m_cellBackground;
		iarray3E dim;
		for(int i=0; i<3; ++i){
			dim[i] = (int)(span[i]/dh + 0.5);
			dim[i] = std::max(2, dim[i]);
		}
		
		if((dim[0]+1)*(dim[1]+1)*(dim[2]+1) > 0.8*nDFS){
			
			radius = distBary;
			radius_old = radius;
				
			dist = evaluateSignedDistance(geoBary, local->getBvTree(), id, normal, radius);	
			
			if(dist > 0)	psign = -1.0;
			if(radius <1.0E-12)	radius = 1.0;
			
			count = 0; 
			for(auto &p : points){
				
				dist = evaluateSignedDistance(p, local->getBvTree(), id, normal, radius);
				if(radius < 1.0E-12)	radius = radius_old;
				else 					radius_old = radius;
				violationField[count] = psign*(dist);
				count++;
			}	
		}else{
			
			//instantiate a VolCartesian;
			bitpit::VolCartesian * mesh = new bitpit::VolCartesian(0,3,bbMin,span, dim);
			mesh->update();
			
			//calculate distance on point of background grid.
			dvecarr3E background_points(mesh->getVertexCount());
			int count = 0;
			for(auto & v: mesh->getVertices()){
				background_points[count] = v.getCoords(); 
				count++;
			}	
			dvector1D background_LS(background_points.size());
			
			radius = 0.5*norm2(span);
			radius_old = radius;
			
			dist = evaluateSignedDistance(geoBary, local->getBvTree(), id, normal, radius);	
			
			if(dist > 0)	psign = -1.0;
			if(radius < 1.0E-12)	radius = radius_old;
			
			count = 0; 
			for(auto &p : background_points){
				
				dist = evaluateSignedDistance(p, local->getBvTree(), id, normal, radius);
				if(radius < 1.0E-12)	radius = radius_old;
				else 					radius_old = radius;
				background_LS[count] = psign*(dist);
				count++;
			}	
			
			//interpolate on background to obtain distance extimated for deformed points	
			
			count = 0;
			for(auto & p : points){
				ivector1D stencil;
				dvector1D weights;
				
				int res = mesh->linearVertexInterpolation(p, stencil,weights);
				
				violationField[count] = 0.0;
				for(int j=0; j<res; ++j){
					violationField[count] += background_LS[stencil[j]]*weights[j]; 
				}
			++count;
			}
			
			delete mesh;
			mesh = NULL;
		}
		
		for(int i=0; i< nDFS; ++i){
			m_violationField[i] = std::fmax(m_violationField[i], (violationField[i]-m_tolerance));
		}
	}
	return;
};

/*!
 * Get infos from a XML bitpit::Config::section. The parameters available are
 * 
 * 1) Files	- external constraint surfaces list of file 
 * 2) Tolerance - within violation tolerance
 * 3) BGDetails - define spacing of background grid, dividing diagonal of box containing geometries by this int factor.
 * 
 * Geometry and its deformation fields are mandatorily passed through ports. 
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ControlDeformExtSurface::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	std::unordered_map<std::string, int> mapp;
	if(slotXML.hasSection("Files")){
		
		bitpit::Config::Section & filesXML = slotXML.getSection("Files");
		
		for(auto & subfile : filesXML.getSections()){
			std::string path;
			std::string tag;
			
			if(subfile.second->hasOption("fullpath"))	{
				path = subfile.second->get("fullpath");
				path = bitpit::utils::trim(path);
			}
			if(subfile.second->hasOption("tag")){
				tag = subfile.second->get("tag");
				tag = bitpit::utils::trim(tag);
				//check tag;
				auto maybe_tag = FileType::_from_string_nothrow(tag.c_str());
				if(!maybe_tag)	tag.clear();
										 else	tag = maybe_tag->_to_string();
			}	
		
		if(!path.empty() && !tag.empty()){
			mapp[path] = (int) FileType::_from_string(tag.c_str());
			}	
		}
	
		setFiles(mapp);
	
	}
	
	if(slotXML.hasOption("Tolerance")){
		std::string input = slotXML.get("Tolerance");
		input = bitpit::utils::trim(input);
		double value = 1.E-8;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setToleranceWithinViolation(value);
	}
	
	if(slotXML.hasOption("BGDetails")){
		std::string input = slotXML.get("BGDetails");
		input = bitpit::utils::trim(input);
		int value = 20;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setBackgroundDetails(value);
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
 * 1) Files	- external constraint surfaces list of file 
 * 2) Tolerance - within violation tolerance
 * 3) BGDetails - define spacing of background grid, dividing diagonal of box containing geometries by this int factor.
 * 
 * Geometry and its deformation fields are mandatorily passed through ports. 
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ControlDeformExtSurface::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
	bitpit::Config::Section & filesXML = slotXML.addSection("Files");
	
	int counter = 0;
	for(auto & file : m_geolist){
		std::string name = "file"+std::to_string(counter);
		bitpit::Config::Section & local = filesXML.addSection(name);
		local.set("fullpath", file.first);
		std::string typetag = (FileType::_from_integral(file.second))._to_string(); 
		local.set("tag", typetag);
		++counter;
	}
	
	std::stringstream ss;
	ss<<std::scientific<<m_tolerance;
	slotXML.set("Tolerance", ss.str());
	
	slotXML.set("BGDetails", std::to_string(m_cellBackground));
	
	if(isPlotInExecution()){
		slotXML.set("PlotInExecution", std::to_string(1));
	}
	
	if(m_outputPlot != "."){
		slotXML.set("OutputPlot", m_outputPlot);
	}
	return;
};

/*!
 * Read all external geoemetries from files (whose name is stored in m_geolist) and return it 
 * in a list of unique pointers pointing to MimmoGeometry objects.
 * \param[out] extGeo list of read external constraint geoemetries.
 */
void ControlDeformExtSurface::readGeometries(std::vector<std::unique_ptr<MimmoGeometry> > & extGeo){
	
	extGeo.resize(m_geolist.size());
	int counter = 0;
	for(auto & geoinfo : m_geolist){
		
		svector1D info = extractInfo(geoinfo.first);
		std::unique_ptr<MimmoGeometry> geo (new MimmoGeometry());
		geo->setRead(true);
		geo->setWrite(false);
		geo->setReadDir(info[0]);
		geo->setReadFilename(info[1]);
		geo->setReadFileType(geoinfo.second);
		geo->setBuildBvTree(true);
		geo->execute();
		
		if(geo->getGeometry()->getNVertex() == 0 || geo->getGeometry()->getNCells() == 0 || !geo->getGeometry()->isBvTreeSupported()){ 
			std::cout<<"failed to read geometry in ControlDeformExtSurface::readGeometries. Skipping file..."<<std::endl;
		}else{
			extGeo[counter] = std::move(geo);
			++counter;
		}
	}	
	
	extGeo.resize(counter);
	return;
};

/*!
 * Extract root dir/filename/tag from an absolute file pattern
 */
svector1D ControlDeformExtSurface::extractInfo(std::string file){
	
	std::string root, name, tag,temp;
	std::string key1=".", key2="/\\";
	
	std::size_t found = file.find_last_of(key2); 
	root = file.substr(0, found);
	temp = file.substr(found+1);
	
	found = temp.find_last_of(key1);
	name = temp.substr(0,found);
	tag = temp.substr(found+1);
	
	svector1D result(3);
	result[0] = root;
	result[1] = name;
	result[2] = tag;
	
	return 	result;
}

/*!
 * Evaluate Signed Distance for a point from given BvTree of a open/closed geometry 3D surface. 
 * Return distance from target geometry with sign. Positive distance is returned, 
 * if the point is placed on the same side of simplex support normal. Negative otherwise.
 * Wrapper to bvTreeUtils::signedDistance
 *
 * \param[in] point 3D target point
 * \param[in] tree	BvTree of target geometry
 * \param[in,out] id returns id of support geometry simplex from which distance is evaluated 
 * \param[in,out] normal returns normal of support geometry simplex from which distance is evaluated 
 * \param[in,out] initRadius guess initial distance.
 * 
 * \return signed distance from target surface.  
 */
double ControlDeformExtSurface::evaluateSignedDistance(darray3E &point, mimmo::BvTree * tree, long & id, darray3E & normal, double &initRadius){
	
	double dist = 1.0E+18;
	double rate = 0.02;
	int kmax = 500;
	int kiter = 0;
	bool flag = true;
	
	while(flag && kiter < kmax){
		dist = bvTreeUtils::signedDistance(&point, tree, id, normal, initRadius);
		flag = (dist == 1.0E+18);
		if(flag)	initRadius *= (1.0+ rate*((double)flag));
		kiter++;
	}
	
	return dist;
}

/*!
 * Plot optional results in execution, that is the violation distance field on current deformed geometry.
 * Reimeplemented from BaseManipulation::plotOptionalResults;
 */
void ControlDeformExtSurface::plotOptionalResults(){
	
	dvecarr3E  points = getGeometry()->getVertexCoords();
	m_defField.resize(points.size());
	points+=m_defField;
	ivector2D connectivity = getGeometry()->getCompactConnectivity();
	
	bitpit::VTKFormat codex = bitpit::VTKFormat::APPENDED;
	bitpit::VTKElementType  elDM = bitpit::VTKElementType::TRIANGLE;
	
	std::string name = m_name +std::to_string(getClassCounter())+ "_ViolationField";
	bitpit::VTKUnstructuredGrid output(m_outputPlot, name, elDM);
	output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points) ;
	output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity) ;
	output.setDimensions(connectivity.size(), points.size());
	//output.setCodex(codex);
	
	std::string sdfstr = "Violation Distance Field";
	output.addData(sdfstr, bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, m_violationField);
	output.write();
}
