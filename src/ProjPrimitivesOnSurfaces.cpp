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

#include "ProjPrimitivesOnSurfaces.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;

/*!Default constructor of ProjPrimitivesOnSurfaces.
 */
ProjPrimitivesOnSurfaces::ProjPrimitivesOnSurfaces(){
	m_name 		= "";
	m_nC = 1000;
	
}

/*!
 * Default destructor of ProjPrimitivesOnSurfaces.
 */
ProjPrimitivesOnSurfaces::~ProjPrimitivesOnSurfaces(){
	clear();
};

/*!Copy constructor of ProjPrimitivesOnSurfaces.Soft Copy of MimmoObject;
 */
ProjPrimitivesOnSurfaces::ProjPrimitivesOnSurfaces(const ProjPrimitivesOnSurfaces & other){
	*this = other;
};	

/*!
 * Assignement operator of ProjPrimitivesOnSurfaces. Soft copy of the target class (no resulting projected 
 * elements will be copied)
 */
ProjPrimitivesOnSurfaces & ProjPrimitivesOnSurfaces::operator=(const ProjPrimitivesOnSurfaces & other){
	clear();
	*(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(&other));
	
	m_topo = other.m_topo;
	m_nC = other.m_nC;
	m_buildBvTree = other.m_buildBvTree;
	m_buildKdTree = other.m_buildKdTree;
	
	//warning the internal data structure of projected element is not copied. Relaunch the execution eventually to fill it.
	return *this;
};

void
ProjPrimitivesOnSurfaces::buildPorts(){
	bool built = true;
	built = (built && createPortIn<MimmoObject*, ProjPrimitivesOnSurfaces>(this, &mimmo::ProjPrimitivesOnSurfaces::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	
	built = (built && createPortOut<MimmoObject*, ProjPrimitivesOnSurfaces>(this, &mimmo::ProjPrimitivesOnSurfaces::getProjectedElement, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	m_arePortsBuilt = built;
}


/*!
 * Return kind of topology the primitive element support. 1-1D, 2-2D
 */
int 
ProjPrimitivesOnSurfaces::getTopology(){
	return m_topo;
}

/*!
 * Get the actual number of the cells that will be used to represent a discrete mesh of the primitive projection.
 */
int
ProjPrimitivesOnSurfaces::getProjElementTargetNCells(){
	return m_nC;
}


/*!
 * Get your current projected primitive as a 3D mesh
 */
MimmoObject *
ProjPrimitivesOnSurfaces::getProjectedElement(){
	return m_patch.get();
}



/*!
 * Set a target external surface mesh where primitive projection need to be performed.
 * Topology of the geometry must be of superficial type, so that MimmoObject::getType() must return 1, otherwise 
 * nothing will be set.
 * \param[in] geo  pointer to MimmoObject
 */
void
ProjPrimitivesOnSurfaces::setGeometry(MimmoObject* geo){
		if(geo->isEmpty()) return;
		if(geo->getType() != 1)	return;
		m_geometry = geo;
};

/*!It sets if the BvTree of the projected primitive element data structure has to be built during execution.
 * \param[in] build If true the BvTree is built in execution and stored in
 * the related MimmoObject member.
 */
void
ProjPrimitivesOnSurfaces::setBuildBvTree(bool build){
	m_buildBvTree = build;
}

/*!It sets if the KdTree the projected primitive element data structure has to be built during execution.
* \param[in] build If true the KdTree is built in execution and stored in
* the related MimmoObject member.
*/
void
ProjPrimitivesOnSurfaces::setBuildKdTree(bool build){
	m_buildKdTree = build;
}

/*!It sets the number of the cells to be used to represent in a discrete mesh the primitive projection.
 * \param[in] nC number of cells
 */
void
ProjPrimitivesOnSurfaces::setProjElementTargetNCells(int nC){
	m_nC = nC;
}

/*!
 * Check if resulting primitive element is present or not.
 * True - no geometry evaluated, False otherwise.
 */
bool 
ProjPrimitivesOnSurfaces::isEmpty(){
	return (m_patch.get() == NULL);
}

/*!
 * Clear all stuffs in your class
 */
void
ProjPrimitivesOnSurfaces::clear(){
	m_patch.reset(nullptr);
	m_nC = 1000;
	BaseManipulation::clear();
};

/*!Execution command.
 * stitch together multiple geoemetry in the same object. 
 */
void
ProjPrimitivesOnSurfaces::execute(){
	if(m_geometry->isEmpty())	return;
	projection();
	
	if(!isEmpty()){
		if(m_buildBvTree)	m_patch->buildBvTree();
		if(m_buildKdTree)	m_patch->buildKdTree();
	}
};	

/*!
 * Plot resulting projected element in a vtu mesh file;
 */
void 
ProjPrimitivesOnSurfaces::plotOptionalResults(){
	if(isEmpty()) return;
	if(m_patch->getNCells() < 1)	return;
	std::string name = m_name + "_" + std::to_string(getClassCounter()) +  "_Patch";
	m_patch->getPatch()->write(name);
}


//***********************************************************//
// ProjSegmentOnSurface IMPLEMENTATION **********************//
//***********************************************************//

/*!
 * Default constructor of ProjPrimitivesOnSurfaces.
 */
ProjSegmentOnSurface::ProjSegmentOnSurface(){
	m_name 		= "MiMMO.ProjSegmentOnSurface";
	m_topo     = 1;
	m_pointA.fill(0.0);
	m_pointB.fill(1.0);
};

/*!
 * Default destructor of ProjSegmentOnSurface.
 */
ProjSegmentOnSurface::~ProjSegmentOnSurface(){
	clear();
};

/*!Copy constructor of ProjSegmentOnSurface.Soft Copy of MimmoObject;
 */
ProjSegmentOnSurface::ProjSegmentOnSurface(const ProjSegmentOnSurface & other){
	*this = other;
};	

/*!
 * Assignement operator of ProjSegmentOnSurface. Soft copy of the target class (no resulting projected 
 * elements will be copied)
 */
ProjSegmentOnSurface & ProjSegmentOnSurface::operator=(const ProjSegmentOnSurface & other){
	clear();
	*(static_cast<ProjPrimitivesOnSurfaces * >(this)) = *(static_cast<const ProjPrimitivesOnSurfaces * >(&other));
	
	m_pointA = other.m_pointA;
	m_pointB = other.m_pointB;
	return *this;
};

/*!
 * Clear all elements on your current class
 */
void
ProjSegmentOnSurface::clear(){
	ProjPrimitivesOnSurfaces::clear();
	m_pointA.fill(0.0);
	m_pointB.fill(1.0);
}

/*!
 * Set the primitive segment, by passing its extremal points
 * \param[in] pointA first extremal point
 * \param[in] pointB second extremal point
 */
void	
ProjSegmentOnSurface::setSegment(darray3E pointA, darray3E pointB){

	if(norm2(pointB - pointA) < 1.E-18)	return;

	m_pointA = pointA;
	m_pointB = pointB;
}

/*!
 * Set the primitive segment, by passing an extremal point, a direction vector and the segment length
 * \param[in] origin extremal point
 * \param[in] dir segment direction
 * \param[in] length lenght of te segment
 */
void	
ProjSegmentOnSurface::setSegment(darray3E origin, darray3E dir, double length){
	setSegment(origin, origin+length*dir);
}

/*!
 * Get settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML. Except of geometry parameter (which is instantiated internally
 * or passed by port linking), the class reads the following parameters:
 * 
 * 0) Segment - SegmentPoints
 * 1) nCells - number of discrete cells of 3D curve
 * 2) BvTree - evaluate bvTree true/false
 * 3) KdTree - evaluate kdTree ture/false
 * 
 * \param[in]	slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
void ProjSegmentOnSurface::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	std::string input; 
	int counter;
	
	if(slotXML.hasOption("Segment")){
		input = slotXML.get("Segment");
		darray3E p1,p2;
		p1.fill(0.0); 
		p2.fill(1.0);
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> p1[0]>>p1[1]>>p1[2]>>p2[0]>>p2[1]>>p2[2];
		}
		setSegment(p1,p2);
	}; 
	if(slotXML.hasOption("nCells")){
		input = slotXML.get("nCells");
		int value = 1000;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> value;
		}
		setProjElementTargetNCells(value);
	}; 
	
	
	if(slotXML.hasOption("BvTree")){
		input = slotXML.get("BvTree");
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> value;
		}
		setBuildBvTree(value);
	}; 
	
	if(slotXML.hasOption("KdTree")){
		input = slotXML.get("KdTree");
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> value;
		}
		setBuildKdTree(value);
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
 * Write settings of the class to bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::flushSectionXML. Except of geometry parameter (which is instantiated internally
 * or passed by port linking), the class writes the following parameters(if different from default):
 * 
 * 0) Segment - SegmentPoints
 * 1) nCells - number of discrete cells of 3D curve
 * 2) BvTree - evaluate bvTree true/false
 * 3) KdTree - evaluate kdTree ture/false
 * 
 * \param[in]	slotXML bitpit::Config::Section which writes to
 * \param[in] name   name associated to the slot 
 */
void ProjSegmentOnSurface::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	slotXML.set("Topology", m_topo);
	
	std::string output;
	
	{
		stringstream ss;
		ss<<m_pointA[0]<<m_pointA[1]<<m_pointA[2]<<m_pointB[0]<<m_pointB[1]<<m_pointB[2];
		slotXML.set("Segment", ss.str());
	}
	{
		output = std::to_string(m_nC);
		slotXML.set("nCells", output);
	}
	
	if(m_buildBvTree){
		output = std::to_string(m_buildBvTree);
		slotXML.set("BvTree", output);
	}
	
	if(m_buildKdTree){
		output = std::to_string(m_buildKdTree);
		slotXML.set("KdTree", output);
	}
	
	if(isPlotInExecution()){
		slotXML.set("PlotInExecution", std::to_string(1));
	}
	
	if(m_outputPlot != "."){
		slotXML.set("OutputPlot", m_outputPlot);
	}
	return;	
};	


/*!
 * core engine for projection.
 */
void
ProjSegmentOnSurface::projection(){	
	
	if(getGeometry()->isEmpty())	return;
	int counter = 0;
	std::unique_ptr<MimmoObject> dum(new MimmoObject(4));
	//reserving memory
	dum->getPatch()->reserveVertices(m_nC+1);
	dum->getPatch()->reserveCells(m_nC);

	//start filling connectivity of your object.
	bitpit::ElementInfo::Type eltype = bitpit::ElementInfo::Type::LINE;
	
	for(int i=0; i<m_nC; ++i){
		livector1D conn(2);
		conn[0] = i;
		conn[1] = i+1;
		dum->addConnectedCell(conn, eltype);
	}
	
	//create points ...
	dvecarr3E verts((m_nC + 1)), projs;
	{
		//create the vertices array, ordered from pointA to pointB
		double dx = norm2(m_pointB - m_pointA)/m_nC;
		counter = 0;
		for( auto && ele : verts){
			ele  = m_pointA + counter*dx;
			++counter;
		}
	}
	
	//...and projecting them onto target surface 
	if(!getGeometry()->isBvTreeBuilt())	getGeometry()->buildBvTree();
	counter = 0;
	projs.resize(verts.size());
	for(auto &val : verts){
		projs[counter]= bvTreeUtils::projectPoint(&val, getGeometry()->getBvTree()); 
		++counter;
	}
	
	//storing the projected points in the MImmoObject:
	long idS = 0;
	for(auto vv : verts){
		dum->addVertex(vv, idS);
		++idS;
	} 

	dum->cleanGeometry();
	m_patch = std::move(dum);
	return;
};

