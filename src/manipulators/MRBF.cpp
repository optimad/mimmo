/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
 *
 *  mimmo is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmo is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
 *
 \ *---------------------------------------------------------------------------*/

#include "MRBF.hpp"

namespace mimmo{


/*! Default Constructor.*/
MRBF::MRBF(MRBFSol mode){
    m_name = "mimmo.MRBF";
    setMode(mode);
    m_maxFields=-1;
    m_tol = 1.0E-06;
    m_bfilter = false;
    m_supportRadiusValue = -1.0;
    m_srIsReal = false;
    m_functype = -1;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
MRBF::MRBF(const bitpit::Config::Section & rootXML){

	m_name = "mimmo.MRBF";
	m_maxFields=-1;
	m_tol = 1.0E-06;
	m_bfilter = false;
    m_supportRadiusValue = -1.0;
    m_srIsReal = false;
    m_functype = -1;

    setMode(MRBFSol::NONE);

	std::string fallback_name = "ClassNONE";
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::string::trim(input);
	if(input == "mimmo.MRBF"){
        std::string fallback_mode = "0";
        std::string input2 = rootXML.get("Mode", fallback_mode);
        input2 = bitpit::utils::string::trim(input2);
        int mode_int = std::stoi(input2);
        mode_int = std::min(2, std::max(0, mode_int));
        setMode(static_cast<MRBFSol>(mode_int));
        absorbSectionXML(rootXML);
	}else{
		warningXML(m_log, m_name);
	};
}

/*! Default Destructor */
MRBF::~MRBF(){};

/*! Copy Constructor. Result geometry displacement are not copied.
 *\param[in] other MRBF where copy from
 */
MRBF::MRBF(const MRBF & other):BaseManipulation(other), bitpit::RBF(other){
	m_tol = other.m_tol;
	m_solver = other.m_solver;
    m_bfilter = other.m_bfilter;
	m_supportRadiusValue  = other.m_supportRadiusValue;
	m_srIsReal = other.m_srIsReal;
    m_supportRadii = other.m_supportRadii;
    m_functype = other.m_functype;
	if(m_bfilter)    m_filter = other.m_filter;
};

/*! Assignment operator. Result geometry displacement are not copied.
 * \param[in] other MRBF where copy from
 */
MRBF& MRBF::operator=(MRBF other){
	swap(other);
	return *this;
};

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void MRBF::swap(MRBF & x) noexcept
{
	std::swap(m_tol, x.m_tol);
	std::swap(m_solver, x.m_solver);
    m_filter.swap(x.m_filter);
    std::swap(m_bfilter, x.m_bfilter);
	std::swap(m_supportRadiusValue , x.m_supportRadiusValue);
	std::swap(m_srIsReal, x.m_srIsReal);
    std::swap(m_supportRadii, x.m_supportRadii);
    std::swap(m_functype, x.m_functype);
	m_displ.swap(x.m_displ);
    std::swap(m_effectiveSR, x.m_effectiveSR);

    RBF::swap(x);

	BaseManipulation::swap(x);
}

/*! It builds the input/output ports of the object
 */
void
MRBF::buildPorts(){
	bool built = true;
    built = (built && createPortIn<MimmoObject*, MRBF>(&m_geometry, M_GEOM, true));
    built = (built && createPortIn<dvecarr3E, MRBF>(this, &mimmo::MRBF::setNode, M_COORDS));
	built = (built && createPortIn<dvecarr3E, MRBF>(this, &mimmo::MRBF::setDisplacements, M_DISPLS));
	built = (built && createPortIn<dmpvector1D*, MRBF>(this, &mimmo::MRBF::setFilter, M_FILTER));
    built = (built && createPortIn<std::vector<double>, MRBF>(this, &mimmo::MRBF::setVariableSupportRadii, M_DATAFIELD));

	built = (built && createPortOut<dmpvecarr3E*, MRBF>(this, &mimmo::MRBF::getDisplacements, M_GDISPLS));
	built = (built && createPortOut<MimmoObject*, MRBF>(this, &BaseManipulation::getGeometry, M_GEOM));
	m_arePortsBuilt = built;
};

/*!It returns a pointer to the RBF node stored in the object.
 * \return pointer to nodes list
 */
dvecarr3E*
MRBF::getNodes(){
	return(&m_node);
}

/*!
 * Return actual solver set for RBF data fields interpolation in MRBF::execute
 * Reimplemented from RBF::getMode() of bitpit;
 * \return solver mode
 */
MRBFSol
MRBF::getMode(){
	return m_solver;
};



/*! It gets current set filter field. See MRBF::setFilter
 * \return filter field.
 */
dmpvector1D*
MRBF::getFilter(){
	return &m_filter;
};

/*!
 * Return true if support radius is set using setSupportRadiusReal method.
   Return false if support radius is set with setSupportRadiusLocal method,
   or provided as variable with setVariableSupportRadii method.
 * \return boolean true/false
 */
bool
MRBF::isSupportRadiusReal(){
	return m_srIsReal;
}

/*!
 * Return true if support radius is provided as variable with setVariableSupportRadii method.
   Return false if support radius is set homogeneous with setSupportRadiusLocal/Real method.
 * \return boolean true/false
 */
bool
MRBF::isVariableSupportRadiusSet(){
	return !m_supportRadii.empty();
}

/*!
 * Return the list of support radii effectively used for the current RBF set (one for each RBF node)
   BEWARE : returned results are meaningful only after class execution.
 * \return list of effective support radii used
 */
dvector1D &
MRBF::getEffectivelyUsedSupportRadii(){
	return m_effectiveSR;
}

/*!
    \return the type of shape function hold by the class.
 */
int
MRBF::getFunctionType(){
    if(m_functype < 0){
        return static_cast<int>(RBF::getFunctionType());
    }else{
        return m_functype;
    }
}

/*!
 * Return actual computed displacements field (if any) for the geometry linked.
 * \return     The computed deformation field on the vertices of the linked geometry
 */
dmpvecarr3E*
MRBF::getDisplacements(){
	return &m_displ;
};


/*!Adds a RBF point to the total control node list and activate it.
 * Reimplemented from RBF::addNode of bitpit;
 * \param[in] node coordinates of control point.
 * \return RBF id.
 */
int
MRBF::addNode(darray3E node){
	return(RBF::addNode(node));
};

/*!Adds a list of RBF points to the total control node list and activate them.
 * Reimplemented from RBF::addNode of bitpit;
 * \param[in] nodes coordinates of control points.
 * \return Vector of RBF ids.
 */
std::vector<int>
MRBF::addNode(dvecarr3E nodes){
	return(RBF::addNode(nodes));
};


/*!Set a RBF point as unique control node and activate it.
 * \param[in] node coordinates of control point.
 */
void
MRBF::setNode(darray3E node){
	removeAllNodes();
	RBF::addNode(node);

};

/*!Set a list of RBF points as control nodes and activate it.
 * \param[in] nodes coordinates of control points.
 */
void
MRBF::setNode(dvecarr3E nodes){
	removeAllNodes();
	RBF::addNode(nodes);
};

/*!Set the RBF points as control nodes extracting
 * the vertices stored in a MimmoObject container.
 * \param[in] geometry Pointer to MimmoObject that contains the geometry.
 */
void
MRBF::setNode(MimmoObject* geometry){
	if(!geometry)    return ;
	removeAllNodes();
	dvecarr3E vertex = geometry->getVerticesCoords();
	RBF::addNode(vertex);

};

/*! Sets filter field. Note: filter field is defined on nodes of the current linked geometry.
 * coherent size between field size and number of geometry vertices is expected.
 * \param[in] filter fields.
 */
void
MRBF::setFilter(dmpvector1D * filter){
    if(!filter) return;
    m_filter.clear();
	m_bfilter = !(filter->empty());
	m_filter = *filter;
};


/*! Find all possible duplicated nodes within a prescribed distance tolerance.
 * Default tolerance value is 1.0E-12;
 * \param[in] tol distance tolerance
 * \return    list of duplicated nodes.
 */
ivector1D
MRBF::checkDuplicatedNodes(double tol){
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
 * \param[in] list pointer to a list of id's of RBF candidate nodes
 * \return    boolean, true if all duplicated nodes are erased, false if one or more of them are not.
 */
bool
MRBF::removeDuplicatedNodes(ivector1D * list){
	ivector1D marked;
	if(list==NULL){
		marked = checkDuplicatedNodes();
		list = &marked;
	}
	return(removeNode(*list));
}


/*! Set ratio a of support radius R of RBF kernel functions (HOMOGENEOUS FOR ALL OF THEM),
   according to the formula R = a*D, where D is the diagonal of the Axis Aligned
   Bounding Box referred to the targetgeometry.
   During the execution the correct value of R is applied.
 * The ratio a can have value between ]0,+inf), which corresponding to minimum locally narrowed
 * function, and almost flat functions (as sphere of infinite radius), respectively.
 * Negative or zero values, bind the evaluation of R to the maximum displacement applied
   to RBF node, that is R is set proportional to the maximum displacement value.
 * \param[in] suppR_ new value of suppR
 */
void
MRBF::setSupportRadiusLocal(double suppR_){
    suppR_ = std::max(-1.0,suppR_);
    m_supportRadiusValue = suppR_;
    m_srIsReal = false;
    m_supportRadii.clear();
}


/*! Set the real physical value of the support radius R of RBF kernel functions (HOMOGENEOUS FOR ALL OF THEM).
 * During the execution the correct value of R is applied.
 * The support radius a can have value between ]0,+inf), which corresponds to minimum locally narrowed
 * function and almost flat functions (as sphere of infinite radius), respectively.
 * Negative or zero values, bind the evaluation of R to the maximum displacement
   applied to RBF node, that is  R is set proportional to the maximum displacement value.
 * \param[in] suppR_ new value of support radius.
 */
void
MRBF::setSupportRadiusReal(double suppR_){
	suppR_ = std::max(-1.0,suppR_);
	m_supportRadiusValue = suppR_;
	m_srIsReal = true;
    m_supportRadii.clear();
}

/*! Legacy method, it does exactly as setSupportRadiusReal method.
 * \param[in] suppR_ new value of support radius.
 */
void
MRBF::setSupportRadiusValue(double suppR_){
    setSupportRadiusReal(suppR_);
}

/*! Set a list of real physical values of the support radius R, one for each RBF
    kernel functions. Method is available only in MRBFSol::NONE mode of the class.
 * Support radii a can have value between ]0,+inf), which corresponds to minimum locally narrowed
 * function and almost flat functions (as sphere of infinite radius), respectively.
 * Negative values are set to minimum allowed support radius value.
   List size may not fit the number of RBF nodes: in that case automatic resize will be performed
   during execution.
 * \param[in] sradii non empty list of variable support radii (otherwise method does nothing)
 */
void
MRBF::setVariableSupportRadii(dvector1D sradii){
    if(sradii.empty() || m_solver != MRBFSol::NONE) return;
    m_supportRadii = sradii;
    m_supportRadiusValue = -1.0;
    m_srIsReal = false;
}

/*!It sets the tolerance for GREEDY mode - interpolation algorithm.
 * Tolerance infos are not used in MRBFSol::NONE/WHOLE mode.
 * \param[in] tol Target tolerance.
 */
void
MRBF::setTol(double tol){
	m_tol = tol;
}

/*!
 * Set a field  of 3D displacements on your RBF Nodes. According to MRBFSol mode
 * active in the class set: displacements as direct RBF weights coefficients in MRBFSol::NONE mode,
 * or interpolate displacements to get the best fit weights in other modes MRBFSol::GREEDY/WHOLE
 * Displacements size may not match the actual number of RBF nodes stored in the class.
 * To ensure consistency call fitDataToNodes() method inherited from RBF class.
 *
 * \param[in] displ list of nodal displacements
 */
void
MRBF::setDisplacements(dvecarr3E displ){
	int size = displ.size();

	removeAllData();

	dvector1D temp(size);
	for(int loc=0; loc<3; ++loc){
		for(int i=0; i<size; ++i){
			temp[i] = displ[i][loc];
		}
		addData(temp);
	}
}


/*!
 * Sets the rbf function to be used. Supported in both modes. (Overloading for mimmo rbf functions)
 * @param[in] bfunc basis function to be used
 */
void
MRBF::setFunction( const MRBFBasisFunction &bfunc )
{
	switch(bfunc){

	case( MRBFBasisFunction::HEAVISIDE10):
            RBF::setFunction(heaviside10);
            m_functype = static_cast<int>(MRBFBasisFunction::HEAVISIDE10);
            break;

	case( MRBFBasisFunction::HEAVISIDE50):
            RBF::setFunction(heaviside50);
            m_functype = static_cast<int>(MRBFBasisFunction::HEAVISIDE50);
            break;

	case( MRBFBasisFunction::HEAVISIDE100):
            RBF::setFunction(heaviside100);
            m_functype = static_cast<int>(MRBFBasisFunction::HEAVISIDE100);
            break;

	case( MRBFBasisFunction::HEAVISIDE1000):
            RBF::setFunction(heaviside1000);
            m_functype = static_cast<int>(MRBFBasisFunction::HEAVISIDE1000);
            break;

    case( MRBFBasisFunction::DSIGMOID):
            RBF::setFunction(dsigmoid);
            m_functype = static_cast<int>(MRBFBasisFunction::DSIGMOID);
            break;

	default:
            bitpit::RBF::setFunction(bitpit::RBFBasisFunction::WENDLANDC2);
            m_functype = -1;
            break;
	}
}

/*!
 * Sets the rbf function to be used. Supported in both modes. (Base method)
 * @param[in] bfunc basis function to be used
 */
void
MRBF::setFunction( const bitpit::RBFBasisFunction &bfunc )
{
	bitpit::RBF::setFunction(bfunc);
    m_functype = -1;
}



/*!Clean all except nodal RBF and its displacements. Use apposite methods RemoveAll*** */
void
MRBF::clear(){
	BaseManipulation::clear();
	clearFilter();
	m_tol = 0.00001;
	m_supportRadiusValue = -1.0;
    m_srIsReal = false;
    m_supportRadii.clear();
};

/*!Clean filter field only*/
void
MRBF::clearFilter(){
	m_filter.clear();
	m_bfilter = false;
};



/*!Execution of RBF object. It evaluates the displacements (values) over the point of the
 * linked geometry, given as result of RBF technique implemented in bitpit::RBF base class.
 * The result is stored in the m_displ member.
 *
 */
void
MRBF::execute(){

	MimmoObject * container = getGeometry();
    if(container == NULL){
        (*m_log)<<m_name + " : NULL pointer to linked geometry found"<<std::endl;
        throw std::runtime_error(m_name + "NULL pointer to linked geometry found");
    }

    if(container->isEmpty()){
        (*m_log)<<m_name + " : empty linked geometry found"<<std::endl;
    }

    //prepare m_displs;
    m_displ.clear();
	m_displ.setDataLocation(mimmo::MPVLocation::POINT);
	m_displ.reserve(getGeometry()->getNVertices());
	m_displ.setGeometry(getGeometry());


    //resize displacements.
	int size = 0;
	int sizeF = getDataCount();

	for (int i=0; i<sizeF; i++){

		if(m_solver == MRBFSol::NONE)    size = m_weight[i].size();
		else                            size = m_value[i].size();

		if(size != getTotalNodesCount()){
			(*m_log) << "warning: " << getName() << " has displacements of " << i << " field with size (" << size << ") that does not fit number of RBF nodes ("<< getTotalNodesCount() << ")" << std::endl;
			fitDataToNodes(i);
		}
	}

    //compute the support radius in m_effectiveSR
    //and push homogeneous support radius info to the base class,
    //in case of Mode WHOLE/GREEDY
    computeEffectiveSupportRadiusList();

    //calculate weights for interpolation modes. This is not required
    // in parameterization mode MRBFSol::NONE.
	if (m_solver == MRBFSol::WHOLE)    solve();
	if (m_solver == MRBFSol::GREEDY)    greedy(m_tol);

    // get deformation using own class evalRBF.
    for(const auto & vertex : container->getVertices()){
        m_displ.insert(vertex.getId(), evalRBF(vertex.getCoords()));
    }

    //apply m_filter if it's active;
    if(m_bfilter){
        checkFilter();
        for (auto it=m_displ.begin(); it!=m_displ.end(); ++it){
            (*it) *= m_filter.at(it.getId());
        }
    }

};

/*!
 * Directly apply deformation field to target geometry.
 */
void
MRBF::apply(){
    _apply(m_displ); //base manipulation utility method.
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
MRBF::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);

	std::string input;

	BaseManipulation::absorbSectionXML(slotXML, name);

    m_tol = 1.0E-6;
	if(slotXML.hasOption("Tolerance")){
		input = slotXML.get("Tolerance");
		input = bitpit::utils::string::trim(input);
		double value = m_tol;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::string::trim(input));
			ss >> value;
			if(value > 0.0)    setTol(value);
		}
	};

	if(slotXML.hasOption("RBFShape")){
		input = slotXML.get("RBFShape");
		int value =1;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::string::trim(input));
			ss>>value;
		}
		value = std::max(1, value);
		if(value > 100){
			if (value > 104){
				value = 104;
			}
			MRBFBasisFunction shapetype = static_cast<MRBFBasisFunction>(value);
			setFunction(shapetype);
		}
		else{
			if(value > 13){
				value =1;
			}
			bitpit::RBFBasisFunction shapetype = static_cast<bitpit::RBFBasisFunction>(value);
			setFunction(shapetype);
		}
	};

    if(slotXML.hasOption("SupportRadiusLocal")){
		input = slotXML.get("SupportRadiusLocal");
		double value = -1.0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::string::trim(input));
			ss >> value;
			setSupportRadiusLocal(value);
		}
	};
    //LEGACY ENTRY old version of SupportRadiusLocal
    if(slotXML.hasOption("SupportRadius")){
		input = slotXML.get("SupportRadius");
		double value = -1.0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::string::trim(input));
			ss >> value;
			setSupportRadiusLocal(value);
		}
	};

	if(slotXML.hasOption("SupportRadiusReal")){
		input = slotXML.get("SupportRadiusReal");
		double value = -1.0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::string::trim(input));
			ss >> value;
			setSupportRadiusReal(value);
		}
	};
}

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
MRBF::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);

	BaseManipulation::flushSectionXML(slotXML, name);

	slotXML.set("Mode", std::to_string(static_cast<int>(m_solver)));

    if(m_tol != 1.0E-6){
        std::stringstream ss;
        ss<<std::scientific<<m_tol;
        slotXML.set("Tolerance", ss.str());
    }
    int type = getFunctionType();
	if(type != static_cast<int>(bitpit::RBFBasisFunction::CUSTOM)){
		slotXML.set("RBFShape", std::to_string(type));
	}

    if(m_supportRadii.empty()){
        std::stringstream ss;
        ss<<std::scientific<<m_supportRadiusValue;
        if(m_srIsReal){
            slotXML.set("SupportRadiusReal", ss.str());
        }else{
            slotXML.set("SupportRadiusLocal", ss.str());
        }
    }
}

/*!
 * Check if the filter is related to the target geometry.
 * If not create a unitary filter field.
 */
void
MRBF::checkFilter(){
	bool check = m_filter.getDataLocation() == mimmo::MPVLocation::POINT;
	check = check && m_filter.completeMissingData(0.0);
	check = check && m_filter.getGeometry() == getGeometry();

	if (!check){
		m_log->setPriority(bitpit::log::Verbosity::DEBUG);
		(*m_log)<<"Not valid filter found in "<<m_name<<". Proceeding with default unitary field"<<std::endl;
		m_log->setPriority(bitpit::log::Verbosity::NORMAL);

		m_filter.clear();
		m_filter.setGeometry(m_geometry);
		m_filter.setDataLocation(mimmo::MPVLocation::POINT);
		m_filter.reserve(getGeometry()->getNVertices());
		for (const auto & vertex : getGeometry()->getVertices()){
			m_filter.insert(vertex.getId(), 1.0);
		}
	}
}

/*!
 * Set a field  of n-Dim weights on your RBF Nodes. Supported only in MRBFSol::NONE mode.
 * Weights total number may not match the actual number of RBF nodes stored in the class.
 * To ensure consistency call fitDataToNodes() method inherited from RBF class.
 *
 * \param[in] value list of nodal weights
 */
void
MRBF::setWeight(dvector2D value){
	if(m_solver != MRBFSol::NONE)    return;

	int size = value.size();

	removeAllData();

	dvector1D temp(size);
	int sizeLoc = 0;
	if(!(value.empty()))    sizeLoc = value[0].size();
	for(int loc=0; loc<sizeLoc; ++loc){
		for(int i=0; i<size; ++i){
			temp[i] = value[i][loc];
		}
		addData(temp);
	}
}

/*! Plot your current rbf nodes as a point cloud to *vtu file.
 * \param[in] directory output directory
 * \param[in] filename  output filename w/out tag
 * \param[in] counterFile   integer identifier of the file
 * \param[in] binary     boolean flag for 0-"ascii" or 1-"appended" writing
 * \param[in] deformed  boolean flag for plotting 0-"original points", 1-"moved points"
 */
void
MRBF::plotCloud(std::string directory, std::string filename, int counterFile, bool binary, bool deformed){

	int nnodes = getTotalNodesCount();
	nnodes = std::min(nnodes, int(m_displ.size()));
	dvecarr3E* nodes_ = getNodes();
	dvecarr3E nodes(nnodes);
	dvecarr3E data(nnodes);
	for(int i=0; i<nnodes; ++i){
		for(int j=0; j<3; ++j){
			if(m_solver == MRBFSol::NONE)   data[i][j] = m_weight[j][i];
			else                            data[i][j] = m_value[j][i];
		}
	}
	if(deformed){
		for(int i=0; i<nnodes; ++i){
			nodes[i] = (*nodes_)[i] + data[i];
		}
	}else{
		for(int i=0; i<nnodes; ++i){
			nodes[i] = (*nodes_)[i];
		}
	}

	bitpit::VTKFormat codex = bitpit::VTKFormat::ASCII;
	if(binary){codex=bitpit::VTKFormat::APPENDED;}

	ivector1D conn(nnodes);
	{
		int counter = 0;
		for(auto & val: conn){
			val = counter;
			++counter;
		}
	}
	bitpit::VTKUnstructuredGrid vtk(directory, filename, bitpit::VTKElementType::VERTEX);
	vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, nodes) ;
	vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, conn) ;
	vtk.setDimensions(conn.size(), nnodes);
    //use connectivity as node labels also;
    vtk.addData("labels",bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, conn) ;
	vtk.setCodex(codex);
	if(counterFile>=0){vtk.setCounter(counterFile);}
	vtk.write();
};

/*!
 * Plot Optional results of the class. It plots the RBF control nodes as a point cloud
 * in *.vtu format, for both original/moved control nodes.
 */
void
MRBF::plotOptionalResults(){

	std::string dir = m_outputPlot;
	std::string nameCloud = m_name;
	std::string nameCloudD = m_name+"_moved";

	plotCloud(dir, nameCloud, getId(), true, false );
	plotCloud(dir, nameCloudD, getId(), true, true );
}

/*!
 * Compute the support radius according on what preferences are collected by class interface methods.
   The m_effectiveSR vector will be filled (structure used by evalRBF method).
 */
void
MRBF::computeEffectiveSupportRadiusList(){

    m_effectiveSR.clear();
    int nodesCount = getTotalNodesCount();
    //First step check variable support radii.
    //if not empty, pass as it is to m_effectiveSR, resize to the total count of nodes and return.
    // if class is in MRBFSol::WHOLE/GREEDY m_supportRadii should already be empty
    // see control of setVariableSupportRadii
    if(!m_supportRadii.empty()){
        m_effectiveSR = m_supportRadii;
        m_effectiveSR.resize(nodesCount, std::numeric_limits<double>::min());
        return;
    }

    // now you are working with singular supportRadius
    double candidateRadius;
    //check if it's negative.
    if(m_supportRadiusValue < std::numeric_limits<double>::min()) {
    //get maximum weight/value displ and assign support radius a 3 times this value.
        int dataCount = getDataCount();
        dvector1D data(dataCount);
        for(int i=0; i<nodesCount; ++i){
            for(int j=0; j<dataCount; ++j){
                if(m_solver == MRBFSol::NONE)    data[j] = m_weight[j][i];
                else                             data[j] = m_value[j][i];
            }
            candidateRadius = std::max(candidateRadius, norm2(data));
        }
        candidateRadius *=3.0;
    }else{

        candidateRadius = m_supportRadiusValue;

        // check if its local(based on AABB)
        if (!m_srIsReal){
            double bboxDiag;
            darray3E pmin, pmax;
            getGeometry()->getPatch()->getBoundingBox(pmin, pmax);
            bboxDiag= norm2(pmax - pmin);
            candidateRadius *= bboxDiag;
        }
    }

    //create homogeneous vector of support radii with this candidateRadius;
    m_effectiveSR.resize(nodesCount, candidateRadius);

    //If you are in WHOLE/GREEDY mode, push to the class base
    // the info on the candidate radius
    if(m_solver != MRBFSol::NONE){
        setSupportRadius(candidateRadius);
    }
}


/*!
 * Evaluates the displacements value with RBF . Supported in all modes.
 * Use weights, RBF node positions and m_effectiveSR (support radius structure) of each RBF node
  to retrive the deformation field.
 *
 * \param[in] point point where to evaluate the basis
 * \return array containing interpolated/parameterized values of displacements.
 *
 */
std::array<double,3>
MRBF::evalRBF( const std::array<double,3> &point){

    std::array<double,3> values;
    values.fill(0.0);
    int                 i, j;
    double              dist, basis;

    for( i=0; i<m_nodes; ++i ){

        if( m_activeNodes[i] ) {

            dist = norm2(point - m_node[i]) / m_effectiveSR[i];
            basis = evalBasis( dist );

            for( j=0; j<3; ++j) {
                values[j] += basis * m_weight[j][i];
            }
        }
    }

    return values;
}
/*!
 * Set type of solver set for RBF data fields interpolation/parameterization in MRBF::execute.
 * Reimplemented from RBF::setMode() of bitpit;
 * \param[in] solver type of MRBFSol enum;
 */
void
MRBF::setMode(MRBFSol solver){
	m_solver = solver;
	if (m_solver == MRBFSol::NONE)    RBF::setMode(bitpit::RBFMode::PARAM);
	else                            RBF::setMode(bitpit::RBFMode::INTERP);
};

/*!
 * Non compact sharp heaviside 0.5*(1.+tanh(k*x)) with k = 10
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double
heaviside10( double dist )
{
	return 1./(1.+std::exp(-10.*(1.-dist)));
}

/*!
 * Non compact sharp heaviside 0.5*(1.+tanh(k*x)) with k = 50
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double
heaviside50( double dist )
{
	return 1./(1.+std::exp(-50.*(1.-dist)));
}

/*!
 * Non compact sharp heaviside 0.5*(1.+tanh(k*x)) with k = 100
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double
heaviside100( double dist )
{
	return 1./(1.+std::exp(-100.*(1.-dist)));
}

/*!
 * Non compact sharp heaviside 0.5*(1.+tanh(k*x)) with k = 1000
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double
heaviside1000( double dist )
{
	return 1./(1.+std::exp(-1000.*(1.-dist)));
}

/*!
 * Non compact sigmoid derivative e^(-10.*x) / (1.+ e^(-10.*x))^2
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double
dsigmoid( double dist )
{
    return std::exp(-10.*dist)/std::pow((1.+std::exp(-10.*dist)),int(2));
}

}
