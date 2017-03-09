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

#include "MRBF.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;


/*! Default Constructor.*/
MRBF::MRBF(){
	m_name = "MiMMO.MRBF";
	m_maxFields=-1;
	m_tol = 0.00001;
	setMode(MRBFSol::NONE);
	m_bfilter = false;
	m_SRRatio = -1.0;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
MRBF::MRBF(const bitpit::Config::Section & rootXML){
	
	m_name = "MiMMO.MRBF";
	m_maxFields=-1;
	m_tol = 0.00001;
	setMode(MRBFSol::NONE);
	m_bfilter = false;
	m_SRRatio = -1.0;

	std::string fallback_name = "ClassNONE";	
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::trim(input);
	if(input == "MiMMO.MRBF"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml MiMMO::MRBF constructor. No valid xml data found"<<std::endl;
	};
}

/*! Default Destructor */
MRBF::~MRBF(){};

/*! Copy Constructor
 *@param[in] other MRBF where copy from
 */
MRBF::MRBF(const MRBF & other){
	*this = other;
};

/*! It builds the input/output ports of the object
 */
void MRBF::buildPorts(){
	bool built = true;
	built = (built && createPortIn<dvecarr3E, MRBF>(this, &mimmo::MRBF::setDisplacements, PortType::M_DISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dvecarr3E, MRBF>(this, &mimmo::MRBF::setNode, PortType::M_COORDS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dvector1D, MRBF>(this, &mimmo::MRBF::setFilter, PortType::M_FILTER, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<double, MRBF>(this, &mimmo::MRBF::setSupportRadius, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<double, MRBF>(this, &mimmo::MRBF::setSupportRadiusValue, PortType::M_VALUED2, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<MimmoObject*, MRBF>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	
	built = (built && createPortOut<dvecarr3E, MRBF>(this, &mimmo::MRBF::getDisplacements, PortType::M_GDISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<std::pair<MimmoObject*, dvecarr3E*> , MRBF>(this, &mimmo::MRBF::getDeformedField, PortType::M_PAIRVECFIELD, mimmo::pin::containerTAG::PAIR, mimmo::pin::dataTAG::MIMMO_VECARR3FLOAT_));
	m_arePortsBuilt = built;
};

/*! Copy Operator
 * @param[in] other MRBF where copy from
 */
MRBF & MRBF::operator=(const MRBF & other){
	*(static_cast<RBF * > (this)) = *(static_cast <const RBF*>(&other));
	*(static_cast<BaseManipulation * > (this)) = *(static_cast <const BaseManipulation * >(&other));
	m_tol = other.m_tol;
	m_solver = other.m_solver;
	m_SRRatio  = other.m_SRRatio;
	m_supRIsValue = other.m_supRIsValue;
	m_bfilter = other.m_bfilter;
	if(m_bfilter)	m_filter = other.m_filter;
	
	return(*this);
};

/*!It sets the geometry linked by the manipulator object (overloading of base class method).
 * @param[in] geometry Pointer to geometry to be deformed by the manipulator object.
 */
void
MRBF::setGeometry(MimmoObject* geometry){
	m_geometry = geometry;
};

/*!It returns a pointer to the RBF node stored in the object.
 */
dvecarr3E*
MRBF::getNodes(){
	return(&m_node);
}

/*! 
 * Return actual solver set for RBF data fields interpolation in MRBF::execute
 * Reimplemented from RBF::getMode() of bitpit;
 */
MRBFSol
MRBF::getMode(){
	return m_solver;
};

/*!
 * Set type of solver set for RBF data fields interpolation/parameterization in MRBF::execute.
 * Reimplemented from RBF::setMode() of bitpit;
 * @param[in] solver type of MRBFSol enum; 
 */
void
MRBF::setMode(MRBFSol solver){
	m_solver = solver;
	if (m_solver == MRBFSol::NONE)	RBF::setMode(RBFMode::PARAM);
	else							RBF::setMode(RBFMode::INTERP);
};
/*!
 * Overloading of MRBF::setSolver(MRBFSol solver) with int input parameter
 * Reimplemented from RBF::setMode() of bitpit;
 * @param[in] int type of solver 1-WHOLE, 2-GREEDY, see MRBFSol enum; 
 */
void 
MRBF::setMode(int type){
	switch(type){
		case 1 : setMode(MRBFSol::WHOLE);
		   break;
		case 2 : setMode(MRBFSol::GREEDY);
		   break;
		default: setMode(MRBFSol::NONE);	
			break;
	}
};

/*! It gets current set filter field. See MRBF::setFilter
 * @return filter field.
 */
dvector1D	MRBF::getFilter(){
	return(m_filter);
};

/*! 
 * It gets current support radius ratio (or value if defined as absolute value) as set up in the class.
 * See MRBF::setSupportRadius and MRBF::setSupportRadiusValue method documentation.
 * Reimplemented from bitpit::RBF::getSupportRadius.
 * @return support radius ratio
 */
double	MRBF::getSupportRadius(){
    if (m_supRIsValue)
    {
        return(getSupportRadiusValue());
    }
    return(m_SRRatio);
};

/*! 
 * It gets the current real value of support radius for RBF kernels set up in the class. 
 * @return support radius value
 */
double	MRBF::getSupportRadiusValue(){
	return(RBF::getSupportRadius());
};

/*!
 * It gets if the support radius for RBF kernels is set up as absolute value (true) or
 * ratio of diagonal of bounding box of the geometry (false).
 * @return support radius is set as value?
 */

bool    MRBF::getIsSupportRadiusValue(){
    return(m_supRIsValue);
}

/*!
 * Return actual computed deformation field (if any) for the geometry linked.
 * If no field is actually present, return null pointers;
 * @return 	std::pair of pointers linking to actual geometry pointed by the class, and the computed deformation field on its vertices
 */
std::pair<MimmoObject * , dvecarr3E * >	MRBF::getDeformedField(){

	std::pair<MimmoObject *, dvecarr3E * > pairField;
	pairField.first = getGeometry();
	pairField.second = &m_displ;
	return pairField;
};

/*!
 * Return actual computed displacements field (if any) for the geometry linked.
 * @return 	The computed deformation field on the vertices of the linked geometry
 */
dvecarr3E
MRBF::getDisplacements(){
	return m_displ;
};


/*!Adds a RBF point to the total control node list and activate it.
 * Reimplemented from RBF::addNode of bitpit;
 * @param[in] node coordinates of control point.
 * @return RBF id.
 */
int MRBF::addNode(darray3E node){
	return(RBF::addNode(node));
};

/*!Adds a list of RBF points to the total control node list and activate them.
 * Reimplemented from RBF::addNode of bitpit;
 * @param[in] nodes coordinates of control points.
 * @return Vector of RBF ids.
 */
std::vector<int> MRBF::addNode(dvecarr3E nodes){
	return(RBF::addNode(nodes));
};

/*!Adds a set of RBF points to the total control node list extracting
 * the vertices stored in a MimmoObject container. Return a vector containing 
 * the RBF node int id.
 * Reimplemented from RBF::addNode of bitpit;
 * @param[in] geometry Pointer to MimmoObject that contains the geometry.
 * @return Vector of RBF ids.
 */
ivector1D MRBF::addNode(MimmoObject* geometry){
	if(geometry == NULL)	return	ivector1D(0);
	dvecarr3E vertex = geometry->getVertexCoords();
	return(RBF::addNode(vertex));
};


/*!Set a RBF point as unique control node and activate it.
 * @param[in] node coordinates of control point.
 */
void MRBF::setNode(darray3E node){
	removeAllNodes();
	RBF::addNode(node);
	return;
};

/*!Set a list of RBF points as control nodes and activate it.
 * @param[in] node coordinates of control points.
 */
void MRBF::setNode(dvecarr3E nodes){
	removeAllNodes();
	RBF::addNode(nodes);
	return;
};

/*!Set the RBF points as control nodes extracting
 * the vertices stored in a MimmoObject container.
 * @param[in] geometry Pointer to MimmoObject that contains the geometry.
 */
void MRBF::setNode(MimmoObject* geometry){
	if(geometry == NULL)	return ;
	removeAllNodes();
	dvecarr3E vertex = geometry->getVertexCoords();
	RBF::addNode(vertex);
	return;
};

/*! Sets filter field. Note: filter field is defined on nodes of the current linked geometry.
 * coherent size between field size and number of geometry vertices is expected.
 * @param[in] filter fields.
 */
void	MRBF::setFilter(dvector1D filter){
	m_filter.clear();
	m_bfilter = !(filter.empty());
	m_filter = filter;
};


/*! Find all possible duplicated nodes within a prescribed distance tolerance.
 * Default tolerance value is 1.0E-12;
 * @param[in] tol distance tolerance
 * @return	list of duplicated nodes.
 */
ivector1D MRBF::checkDuplicatedNodes(double tol){
	ivector1D marked;
	int sizeEff = getTotalNodesCount();
	if( sizeEff == 0 ) return marked;
	
	bvector1D check(sizeEff, false);
	
	for(int i=0; i<sizeEff; ++i){
		for(int j=i+1; j<sizeEff; ++j){
			double dist = norm2(m_node[j] - m_node[i]);
			if(!check[j] && dist <= tol){
				marked.push_back(i);
				check[j] = true;
			}
		}
	}	
	return(marked);	
}

/*! Erase all nodes passed by their RBF id list. If no list is provided, the method find all 
 * possible duplicated nodes within a default tolerance of 1.0E-12 and erase them, if any.
 * @param[in] list pointer to a list of id's of RBF candidate nodes
 * @return	boolean, true if all duplicated nodes are erased, false if one or more of them are not.
 */
bool MRBF::removeDuplicatedNodes(ivector1D * list){
	ivector1D marked;
	if(list==NULL){
		marked = checkDuplicatedNodes();
		list = &marked;
	}
	return(removeNode(*list));
}


/*! Set ratio a of support radius R of RBF kernel functions, according to the formula
 * R = a*D, where D is the diagonal of the Axis Aligned Bounding Box referred to the target
 * geometry. During the execution the correct value of R is applied.
 * The ratio a can have value between 0 and +inf (with 0 excluded), which corresponding to minimum locally narrowed
 * function, and almost flat functions (as sphere of infinite radius), respectively. 
 * Negative or zero values, bind the evaluation of R to the maximum displacement applied to RBF node, that is 
 * R is set proportional to the maximum displacement value.
 * @param[in] suppR new value of suppR
 */
void
MRBF::setSupportRadius(double suppR_){
    suppR_ = std::fmax(-1.0,suppR_);
    m_SRRatio = suppR_;
    m_supRIsValue = false;
}


/*! Set the value of the support radius R of RBF kernel functions.
 * During the execution the correct value of R is applied.
 * The support radius a can have value between 0 and +inf (with 0 excluded), which corresponding to minimum locally narrowed
 * function, and almost flat functions (as sphere of infinite radius), respectively.
 * Negative or zero values, bind the evaluation of R to the maximum displacement applied to RBF node, that is
 * R is set proportional to the maximum displacement value.
 * @param[in] suppR_ new value of support radius.
 */
void
MRBF::setSupportRadiusValue(double suppR_){
    suppR_ = std::fmax(-1.0,suppR_);
    m_SRRatio = suppR_;
    m_supRIsValue = true;
}

/*!It sets the tolerance for greedy - interpolation algorithm.
 * Tolerance infos are not used in MRBFSol::NONE mode. 
 * @param[in] tol Target tolerance.
 */
void MRBF::setTol(double tol){
	m_tol = tol;
}

/*!
 * Set a field  of 3D displacements on your RBF Nodes. According to MRBFSol mode
 * active in the class set: displacements as direct RBF weights coefficients in MRBFSol::NONE mode,
 * or interpolate displacements to get the best fit weights in other modes MRBFSol::GREEDY/WHOLE
 * Displacements size may not match the actual number of RBF nodes stored in the class.
 * To ensure consistency call fitDataToNodes() method inherited from RBF class.
 * 
 * @param[in] displ list of nodal displacements
 */
void MRBF::setDisplacements(dvecarr3E displ){
	int size = displ.size();
	if(size != getTotalNodesCount()){
		std::cout << "MiMMO : WARNING : " << getName() << " sets displacements with size (" << size << ") that does not fit number of RBF nodes ("<< getTotalNodesCount() << ")" << std::endl;
	}
	
	removeAllData();
	
	dvector1D temp(size);
	for(int loc=0; loc<3; ++loc){
		for(int i=0; i<size; ++i){
			temp[i] = displ[i][loc];
		}
		addData(temp);
	}
}

/*!Clean all except nodal RBF and its displacements. Use apposite methods RemoveAll*** */
void MRBF::clear(){
	BaseManipulation::clear();
	clearFilter();
	m_tol = 0.00001;
	m_SRRatio = -1;
};

/*!Clean filter field */
void MRBF::clearFilter(){
	m_filter.clear();
	m_bfilter = false;
};

/*!
 * Set a field  of n-Dim weights on your RBF Nodes. Supported only in MRBFSol::NONE mode.
 * Weights total number may not match the actual number of RBF nodes stored in the class.
 * To ensure consistency call fitDataToNodes() method inherited from RBF class.
 * 
 * @param[in] displ list of nodal weights
 */
void
MRBF::setWeight(dvector2D value){
	if(m_solver != MRBFSol::NONE)	return;
		   
	int size = value.size();
	if(size != getTotalNodesCount()){
		std::cout << "MiMMO : WARNING : " << getName() << " sets weights with size (" << size << ") that does not fit number of RBF nodes ("<< getTotalNodesCount() << ")" << std::endl;
	}
	
	removeAllData();
	
	dvector1D temp(size);
	int sizeLoc = 0;
	if(!(value.empty()))	sizeLoc = value[0].size();
	for(int loc=0; loc<sizeLoc; ++loc){
		for(int i=0; i<size; ++i){
			temp[i] = value[i][loc];
		}
		addData(temp);
	}
}

/*!Execution of RBF object. It evaluates the displacements (values) over the point of the
 * linked geometry, given as result of RBF technique implemented in bitpit::RBF base class.
 * The result is stored in the m_displ member.
 *
 */
void MRBF::execute(){

	MimmoObject * container = getGeometry();
	if(container->isEmpty() ) return;

	int size, sizeF = getDataCount();
	for (int i=0; i<sizeF; i++){
		
		if(m_solver == MRBFSol::NONE)	size = m_weight[i].size();
		else							size = m_value[i].size();
		
		if(size != getTotalNodesCount()){
			std::cout << "MiMMO : WARNING : " << getName() << " has displacements of " << i << " field with size (" << size << ") that does not fit number of RBF nodes ("<< getTotalNodesCount() << ")" << std::endl;
			fitDataToNodes(i);
		}
	}

	double bboxDiag;
	{
		darray3E pmin, pmax;
		container->getPatch()->getBoundingBox(pmin, pmax);
		bboxDiag= norm2(pmax - pmin);
	}
	
	//Checking supportRadius.
	double distance = 0.0;
	if(m_SRRatio <=0.0){ //get maximum weight/value displ and assign support radius a 3 times this value.

		for(int i=0; i<size; ++i){
			
			dvector1D data(sizeF);

			for(int j=0; j<sizeF; ++j){
				if(m_solver == MRBFSol::NONE)	data[j] = m_weight[j][i];
				else							data[j] = m_value[j][i];
			}	

			distance = std::max(distance, norm2(data));
		}
		
		distance *=3.0;
	
	}else{
	    if (m_supRIsValue){
	        distance = m_SRRatio;
	    }
	    else{
	        distance = m_SRRatio * bboxDiag;
	    }
	}
	
	//TODO remove it (I can't use a support radius 0....why?)
	if(distance <=1.E-18){ //checkSupportRadius if too small, set it to the semidiagonal value of the geometry AABB
		distance = 0.5*bboxDiag;
	}
	
	const double radius = distance;
	RBF::setSupportRadius(radius);
	
   if (m_solver == MRBFSol::WHOLE)	solve();
   if (m_solver == MRBFSol::GREEDY)	greedy(m_tol);

	int nv = container->getNVertex();
	dvecarr3E vertex = container->getVertexCoords();

	m_displ.resize(nv, darray3E{0,0,0});
	dvector1D displ;
	for(int i=0; i<nv; ++i){
		displ = RBF::evalRBF(vertex[i]);
		for (int j=0; j<3; j++) m_displ[i][j] = displ[j];
	}
	
	//if m_filter is active;
	if(m_bfilter){
		m_filter.resize(nv,1.0);
		int counter = 0;
		for (auto && vec : m_displ){
			vec = vec * m_filter[counter];
			++counter;
		}
	}

	return;
};


/*!
 * Method to absorb parameter infos from an XML parser class of bitpit. 
 * The sensible parameters are:
 * 
 * --> Absorbing data:
 * - Priority  : uint marking priority in multi-chain execution;
 * - Mode : mode of usage of the class 0-parameterizator class, 1-regular interpolator class, 2- greedy interpolator class )
 * - SupportRadius : local radius of RBF function for each nodes
 * - RBFShape : shape of RBF function 1(wendlandc2), 2(linear), 3(gauss90), 4(gauss95),5(gauss99) 
 * - Tolerance : greedy engine tolerance (meant for mode 2);
 * 
 * RBF node list, Filter to deformation, Geometry and RBF nodal displacements are passed through port linking
 * 
 * \param[in] slotXML	reference to a Section slot of bitpit::Config class.
 * \param[in] name   name associated to the slot
 */
void  MRBF::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
	
	std::string input; 
	
	if(slotXML.hasOption("Priority")){
		input = slotXML.get("Priority");
		int value =0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss>>value;
		}
		setPriority(value);
	};
	
	if(slotXML.hasOption("RBFShape")){
		input = slotXML.get("RBFShape");
		int value =1;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss>>value;
		}
		value = std::max(1, value);
		if(value > 5)	value =1;
	   bitpit::RBFBasisFunction shapetype = static_cast<bitpit::RBFBasisFunction>(value);
	   setFunction(shapetype);
	}; 
	
	if(slotXML.hasOption("Mode")){
		input = slotXML.get("Mode");
		int value = 0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> value;
			value = std::max(value, 0);
			if(value > 2) value = 0;
		}
		setMode(value);
	}; 
	
	if(slotXML.hasOption("SupportRadius")){
		input = slotXML.get("SupportRadius");
		double value = -1.0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> value;
			setSupportRadius(value);
		}
	}; 
	
	m_tol = 1.0E-6;
	if(slotXML.hasOption("Tolerance")){
		input = slotXML.get("Tolerance");
		input = bitpit::utils::trim(input);
		double value = m_tol;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> value;
			if(value > 0.0)	setTol(value);
		}
	}; 
	
	return;
}

/*!
 * Method to flush parameter infos to an XML parser class of bitpit. 
 * The sensible parameters are:
 * 
 * --> Flushing data// how to write it on XML:
 * - ClassName : name of the class as "MiMMO.MRBF"
 * - Priority  : uint marking priority in multi-chain execution;
 * - Mode : mode of usage of the class 0-parameterizator class, 1-regular interpolator class, 2- greedy interpolator class )
 * - SupportRadius : local radius of RBF function for each nodes
 * - RBFShape : shape of RBF function 1(wendlandc2), 2(linear), 3(gauss90), 4(gauss95),5(gauss99) 
 * - Tolerance : greedy engine tolerance (meant for mode 2);
 * 
 * RBF node list, Filter to deformation, Geometry and RBF nodal displacements are passed through port linking
 * In any case, if different by default, such parameters are always written, even.
 * 
 * \param[in] slotXML	reference to a Section slot of bitpit::Config class.
 * \param[in] name   name associated to the slot
 */
void  MRBF::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("Priority", std::to_string(getPriority()));
	
	std::string input;
	input = std::to_string(static_cast<int>(m_solver));
	slotXML.set("Mode", input);
	
	{
		std::stringstream ss;
		ss<<std::scientific<<getSupportRadius();
		slotXML.set("SupportRadius", ss.str());
	}
	
	bitpit::RBFBasisFunction type = getTypeFunction();
	if(type != RBFBasisFunction::CUSTOM){
		int val = static_cast<int>(type);
		slotXML.set("RBFShape", std::to_string(val));
	}
	
	//checking if not default and if not connected to a port
	if(m_tol != 1.0E-6 ){
		std::stringstream ss;
		ss<<std::scientific<<m_tol;
		slotXML.set("Tolerance", ss.str());
	}

	return;
}

