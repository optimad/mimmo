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

#include "MRBFStyleObj.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;


/*! Default Constructor.*/
MRBFStyleObj::MRBFStyleObj(){
	m_name = "MiMMO.MRBFStyleObj";
	m_maxFields=-1;
	m_tol = 0.00001;
	setMode(MRBFSol::NONE);
	m_bfilter = false;
	m_SRRatio = -1.0;
};

/*! Default Destructor */
MRBFStyleObj::~MRBFStyleObj(){};

/*! Copy Constructor
 *@param[in] other MRBFStyleObj where copy from
 */
MRBFStyleObj::MRBFStyleObj(const MRBFStyleObj & other){
	*this = other;
};

/*! It builds the input/output ports of the object
 */
void MRBFStyleObj::buildPorts(){
	bool built = true;
	built = (built && createPortIn<dvecarr3E, MRBFStyleObj>(this, &mimmo::MRBFStyleObj::setDisplacements, PortType::M_DISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<MimmoObject*, MRBFStyleObj>(this, &mimmo::MRBFStyleObj::setAddNode, PortType::M_GEOM2, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	built = (built && createPortIn<std::vector<MimmoObject*>, MRBFStyleObj>(this, &mimmo::MRBFStyleObj::setAddNode, PortType::M_VECGEOM, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::MIMMO_));
	built = (built && createPortIn<dvector1D, MRBFStyleObj>(this, &mimmo::MRBFStyleObj::setFilter, PortType::M_FILTER, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<double, MRBFStyleObj>(this, &mimmo::MRBFStyleObj::setSupportRadius, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<double, MRBFStyleObj>(this, &mimmo::MRBFStyleObj::setSupportRadiusValue, PortType::M_VALUED2, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<MimmoObject*, MRBFStyleObj>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	
	built = (built && createPortOut<dvecarr3E, MRBFStyleObj>(this, &mimmo::MRBFStyleObj::getDisplacements, PortType::M_GDISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<std::pair<MimmoObject*, dvecarr3E*> , MRBFStyleObj>(this, &mimmo::MRBFStyleObj::getDeformedField, PortType::M_PAIRVECFIELD, mimmo::pin::containerTAG::PAIR, mimmo::pin::dataTAG::MIMMO_VECARR3FLOAT_));
	m_arePortsBuilt = built;
};

/*! Copy Operator
 * @param[in] other MRBFStyleObj where copy from
 */
MRBFStyleObj & MRBFStyleObj::operator=(const MRBFStyleObj & other){
	*(static_cast<RBFKernel * > (this)) = *(static_cast <const RBFKernel*>(&other));
	*(static_cast<BaseManipulation * > (this)) = *(static_cast <const BaseManipulation * >(&other));
	m_tol = other.m_tol;
	m_solver = other.m_solver;
	m_SRRatio  = other.m_SRRatio;
	m_node = other.m_node;
	m_bboxes = other.m_bboxes;
	m_supRIsValue = other.m_supRIsValue;
	m_bfilter = other.m_bfilter;
	if(m_bfilter)	m_filter = other.m_filter;
	return(*this);
};

/*!It sets the geometry linked by the manipulator object (overloading of base class method).
 * @param[in] geometry Pointer to geometry to be deformed by the manipulator object.
 */
void
MRBFStyleObj::setGeometry(MimmoObject* geometry){
	m_geometry = geometry;
};

/*!
 * It returns a pointer to the RBF nodes stored in the object.
 */
std::vector<MimmoObject*>
MRBFStyleObj::getNodes(){
	return(m_node);
}

/*!
 * Return the number of nodes set for your class
 */
int
MRBFStyleObj::getTotalNodesCount(){
	return m_node.size();
}

/*! 
 * Return actual solver set for RBF data fields interpolation in MRBFStyleObj::execute
 * Reimplemented from RBFKernel::getMode() of bitpit;
 */
MRBFSol
MRBFStyleObj::getMode(){
	return m_solver;
};

/*!
 * Set type of solver set for RBF data fields interpolation/parameterization in MRBFStyleObj::execute.
 * Reimplemented from RBFKernel::setMode() of bitpit;
 * @param[in] solver type of MRBFSol enum; 
 */
void
MRBFStyleObj::setMode(MRBFSol solver){
	m_solver = solver;
	if (m_solver == MRBFSol::NONE)	RBFKernel::setMode(RBFMode::PARAM);
	else							RBFKernel::setMode(RBFMode::INTERP);
};
/*!
 * Overloading of MRBFStyleObj::setSolver(MRBFSol solver) with int input parameter
 * Reimplemented from RBFKernel::setMode() of bitpit;
 * @param[in] int type of solver 1-WHOLE, 2-GREEDY, see MRBFSol enum; 
 */
void 
MRBFStyleObj::setMode(int type){
	switch(type){
		case 1 : setMode(MRBFSol::WHOLE);
		   break;
		case 2 : setMode(MRBFSol::GREEDY);
		   break;
		default: setMode(MRBFSol::NONE);	
			break;
	}
};

/*! It gets current set filter field. See MRBFStyleObj::setFilter
 * @return filter field.
 */
dvector1D	MRBFStyleObj::getFilter(){
	return(m_filter);
};

/*! 
 * It gets current support radius ratio (or value if defined as absolute value) as set up in the class.
 * See MRBFStyleObj::setSupportRadius and MRBFStyleObj::setSupportRadiusValue method documentation.
 * Reimplemented from bitpit::RBFKernel::getSupportRadius.
 * @return support radius ratio
 */
double	MRBFStyleObj::getSupportRadius(){
    if (m_supRIsValue){
        return(getSupportRadiusValue());
    }
    return(m_SRRatio);
};

/*! 
 * It gets the current real value of support radius for RBF kernels set up in the class. 
 * @return support radius value
 */
double	MRBFStyleObj::getSupportRadiusValue(){
	return(RBFKernel::getSupportRadius());
};

/*!
 * It gets if the support radius for RBF kernels is set up as absolute value (true) or
 * ratio of diagonal of bounding box of the geometry (false).
 * @return support radius is set as value?
 */

bool    MRBFStyleObj::getIsSupportRadiusValue(){
    return(m_supRIsValue);
}

/*!
 * Return actual computed deformation field (if any) for the geometry linked.
 * If no field is actually present, return null pointers;
 * @return 	std::pair of pointers linking to actual geometry pointed by the class, and the computed deformation field on its vertices
 */
std::pair<MimmoObject * , dvecarr3E * >	MRBFStyleObj::getDeformedField(){

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
MRBFStyleObj::getDisplacements(){
	return m_displ;
};


/*!Adds a RBF node to the total control node list and activate it. Appending mode.
 * @param[in] node pointer to a MimmoObject container holding the mesh/node  
 */
void MRBFStyleObj::setAddNode(MimmoObject * node){
	if(node->isEmpty())	return;
	//save bounding boxes adding nodes
	darray3E minP, maxP;
	node->getPatch()->getBoundingBox(minP, maxP);
	m_bboxes[node].resize(2);
	m_bboxes[node][0] = minP;
	m_bboxes[node][1] = maxP;
	//store node
	m_node.push_back(node);
	m_activeNodes.push_back(true);
	m_nodes++;
	return;
};

/*!Adds a list of RBF nodes to the total control node list and activate them. Appending mode.
 * @param[in] nodes list of pointers to MimmoObject containers holding the meshes/nodes
 */
void MRBFStyleObj::setAddNode(std::vector<MimmoObject*> nodelist){
	for(auto & val: nodelist){
		setAddNode(val);
	}
	return;
};

/*!Set a list of mesh objects as RBF control nodes and activate it. Any previous list will
 * be thrown away and substituted by the current list.
 * @param[in] node list of pointers to MimmoObject containers holding the meshes/nodes
 */
void MRBFStyleObj::setNode(std::vector<MimmoObject*> nodelist){
	removeAllNodes();
	setAddNode(nodelist);
	return;
};

/*! Sets filter field. Note: filter field is defined on nodes of the current linked geometry.
 * coherent size between field size and number of geometry vertices is expected.
 * @param[in] filter fields.
 */
void	MRBFStyleObj::setFilter(dvector1D filter){
	m_filter.clear();
	m_bfilter = !(filter.empty());
	m_filter = filter;
};



/*! Remove pre-existent node. MRBFStyleObj Node list is resized and renumbered after extraction.
 * Supported in both modes.
 * @param[in] id id of node
 * @return boolean, true if successfully extracted, false otherwise
 */
bool MRBFStyleObj::removeNode(int id){
	
	if(id < 0 || id >=m_nodes) return false;
	
	m_nodes--;
	auto key = *(m_node.begin()+id);
	m_node.erase(m_node.begin()+id);
	m_activeNodes.erase(m_activeNodes.begin()+id);
	m_bboxes.erase(key);
	return(true);
}

/*! Remove pre-existent set of nodes. MRBFStyleObj nodal list is resized and renumbered after extraction.
 *  Supported in both modes.
 * @param[in] list id list of candidates to extraction
 * @return boolean, true if all nodes are successfully extracted, false if any of them or none are extracted
 */
bool MRBFStyleObj::removeNode(std::vector<int> & list){
	
	std::set<int> setList;
	for(auto && id : list) setList.insert(id);
	
	int extracted = 0;
	for(auto && id : setList){
		if(id>=0 && id <m_nodes){;
			m_nodes--;
			int index = id-extracted;
			auto key = *(m_node.begin() + index);
			m_node.erase(m_node.begin() + index);
			m_activeNodes.erase(m_activeNodes.begin() + index);
			m_bboxes.erase(key);
			extracted++;
		}
	}
	return(extracted == (int)(list.size()));
}

/*!
 * Remove all nodes in MRBFStyleObj nodal list. Supported in both modes.
 */
void MRBFStyleObj::removeAllNodes(){
	m_nodes = 0;
	m_node.clear();
	m_activeNodes.clear();
	m_bboxes.clear();
}


/*! Find all possible duplicated nodes checking their address signature
 * @return	list of duplicated nodes.
 */
ivector1D MRBFStyleObj::checkDuplicatedNodes(){
	ivector1D marked;
	int sizeEff = getTotalNodesCount();
	if( sizeEff == 0 ) return marked;
	
	bvector1D check(sizeEff, false);
	
	for(int i=0; i<sizeEff; ++i){
		for(int j=i+1; j<sizeEff; ++j){
			if(!check[j] && m_node[j] == m_node[i]){
				marked.push_back(j);
				check[j] = true;
			}
		}
	}	
	return(marked);	
}

/*! Erase all nodes passed by their RBF id list. If no list is provided, the method find all 
 * possible duplicated nodes and erase them, if any.
 * @param[in] list pointer to a list of id's of RBF candidate nodes
 * @return	boolean, true if all duplicated nodes are erased, false if one or more of them are not.
 */
bool MRBFStyleObj::removeDuplicatedNodes(ivector1D * list){
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
MRBFStyleObj::setSupportRadius(double suppR_){
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
MRBFStyleObj::setSupportRadiusValue(double suppR_){
    suppR_ = std::fmax(-1.0,suppR_);
    m_SRRatio = suppR_;
    m_supRIsValue = true;
}

/*!It sets the tolerance for greedy - interpolation algorithm.
 * Tolerance infos are not used in MRBFSol::NONE mode. 
 * @param[in] tol Target tolerance.
 */
void MRBFStyleObj::setTol(double tol){
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
void MRBFStyleObj::setDisplacements(dvecarr3E displ){
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
void MRBFStyleObj::clear(){
	BaseManipulation::clear();
	clearFilter();
	m_tol = 0.00001;
	m_SRRatio = -1;
};

/*!Clean filter field */
void MRBFStyleObj::clearFilter(){
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
MRBFStyleObj::setWeight(dvector2D value){
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
 * linked geometry, given as result of RBF technique implemented in bitpit::RBFKernel base class.
 * The result is stored in the m_displ member.
 *
 */
void MRBFStyleObj::execute(){

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
	RBFKernel::setSupportRadius(radius);
	
   if (m_solver == MRBFSol::WHOLE)	solve();
   if (m_solver == MRBFSol::GREEDY)	greedy(m_tol);

	int nv = container->getNVertex();
	dvecarr3E vertex = container->getVertexCoords();

	m_displ.resize(nv, darray3E{0,0,0});
	dvector1D displ;
	for(int i=0; i<nv; ++i){
		displ = evalRBF(vertex[i]);
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
 * 1)	Mode - mode of usage of the class 0-parameterizator class, 1-regular interpolator class, 2- greedy interpolator class )
 * 2)	SupportRadius	- local radius of RBF function for each nodes
 * 3)	Tolerance - greedy engine tolerance (meant for mode 2);
 * 
 * RBF node list, Filter to deformation, Geometry and RBF nodal displacements are passed through port linking
 * Sometimes 2) and 3) are not parameters and can be passed through respective port.
 * 
 * \param[in] slotXML	reference to a Section slot of bitpit::Config class.
 * \param[in] name   name associated to the slot
 */
void  MRBFStyleObj::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	std::string input; 
	
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
	
	if(slotXML.hasOption("SupportRadiusValue")){
		input = slotXML.get("SupportRadiusValue");
		double value = -1.0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> value;
			setSupportRadiusValue(value);
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
 * 1)	Mode - mode of usage of the class 0-parameterizator class, 1-regular interpolator class, 2- greedy interpolator class )
 * 2)	SupportRadius	- local radius of RBF function for each nodes
 * 3)	Tolerance - greedy engine tolerance (meant for mode 2);
 * 
 * RBF node list, Filter to deformation, Geometry and RBF nodal displacements are passed through port linking
 * Sometimes 2) and 3) are not parameters and can be passed through respective port.
 * In any case, if different by default, such parameters are always written, even.
 * 
 * \param[in] slotXML	reference to a Section slot of bitpit::Config class.
 * \param[in] name   name associated to the slot
 */
void  MRBFStyleObj::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
	
	std::string input;
	input = std::to_string(static_cast<int>(m_solver));
	slotXML.set("Mode", input);
	
	{
		std::stringstream ss;
		ss<<std::scientific<<getSupportRadiusValue();
		slotXML.set("SupportRadiusValue", ss.str());
	}
	//checking if not default and if not connected to a port
	if(m_tol != 1.0E-6 ){
		std::stringstream ss;
		ss<<std::scientific<<m_tol;
		slotXML.set("Tolerance", ss.str());
	}

	return;
}

/*!
 * Evaluate the distance between two RBF nodes of the class  
 */
 double MRBFStyleObj::calcDist(int i, int j){
	 
	 double distMin = 1.0E+18;
	 if(!m_node[i]->isBvTreeBuilt())	m_node[i]->buildBvTree();
	 if(!m_node[j]->isBvTreeBuilt())	m_node[j]->buildBvTree();
	 
	 //find an indicative distance of search between geometries
	 double initRadius, workRadius;
	 double dist = 1.E+18;
	 long id;
	 
	{
		 darray3E cA,cB;
		 cA = 0.5*(m_bboxes[m_node[i]][0]+m_bboxes[m_node[i]][1]);
		 cB = 0.5*(m_bboxes[m_node[j]][0]+m_bboxes[m_node[j]][1]);
		 initRadius = norm2(cB - cA);
	}
	
	//search first set of points against second figure
	for(auto & val : m_node[i]->getVertexCoords()){
		 workRadius= initRadius;
		 dist = bvTreeUtils::distance(&val, m_node[j]->getBvTree(), id, workRadius);
		 if(dist < distMin)	distMin = dist;	
	}
	
	//search second set of points against first figure
	for(auto & val : m_node[j]->getVertexCoords()){
		workRadius= initRadius;
		dist = bvTreeUtils::distance(&val, m_node[i]->getBvTree(), id, workRadius);
		if(dist < distMin)	distMin = dist;	
	}
	
	return distMin;
	 
} 

/*!
 * Evaluate the distance between a 3D point and a RBF node of the class  
 */
double MRBFStyleObj::calcDist(const darray3E & point, int j){
	
	double dist = 1.0E+18;
	double rate = 0.02;
	int kmax = 1000;
	int kiter = 0;
	bool flag = true;
	long id;
	auto pp = point;
	
	if(!m_node[j]->isBvTreeBuilt())	m_node[j]->buildBvTree();
	double initRadius = norm2(m_bboxes[m_node[j]][1] - m_bboxes[m_node[j]][0]);
	
	while(flag && kiter < kmax){
		dist = bvTreeUtils::distance(&pp, m_node[j]->getBvTree(), id, initRadius);
		flag = (dist == 1.0E+18);
		if(flag)	initRadius *= (1.0+ rate*((double)flag));
		kiter++;
	}
	
	return dist;
} 
