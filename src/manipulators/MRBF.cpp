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

using namespace std;
using namespace bitpit;
namespace mimmo{


/*! Default Constructor.*/
MRBF::MRBF(){
	m_name = "mimmo.MRBF";
	m_maxFields=-1;
	m_tol = 0.00001;
	setMode(MRBFSol::NONE);
	m_bfilter = false;
	m_SRRatio = -1.0;
    m_functype = -1;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
MRBF::MRBF(const bitpit::Config::Section & rootXML){

	m_name = "mimmo.MRBF";
	m_maxFields=-1;
	m_tol = 0.00001;
	setMode(MRBFSol::NONE);
	m_bfilter = false;
	m_SRRatio = -1.0;
    m_functype = -1;

	std::string fallback_name = "ClassNONE";
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::string::trim(input);
	if(input == "mimmo.MRBF"){
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
	m_SRRatio  = other.m_SRRatio;
	m_supRIsValue = other.m_supRIsValue;
	m_bfilter = other.m_bfilter;
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
	std::swap(m_SRRatio , x.m_SRRatio);
	std::swap(m_supRIsValue, x.m_supRIsValue);
	std::swap(m_bfilter, x.m_bfilter);
    std::swap(m_functype, x.m_functype);
	//    std::swap(m_filter, x.m_filter);
	//    std::swap(m_displ, x.m_displ);
	m_filter.swap(x.m_filter);
	m_displ.swap(x.m_displ);
	RBF::swap(x);
	BaseManipulation::swap(x);
}

/*! It builds the input/output ports of the object
 */
void
MRBF::buildPorts(){
	bool built = true;
	built = (built && createPortIn<dmpvecarr3E, MRBF>(this, &mimmo::MRBF::setDisplacements, M_GDISPLS));
	built = (built && createPortIn<dvecarr3E, MRBF>(this, &mimmo::MRBF::setDisplacements, M_DISPLS));
	built = (built && createPortIn<dvecarr3E, MRBF>(this, &mimmo::MRBF::setNode, M_COORDS));
	built = (built && createPortIn<MimmoObject*, MRBF>(this, &mimmo::MRBF::setNode, M_GEOM2));
	built = (built && createPortIn<dmpvector1D, MRBF>(this, &mimmo::MRBF::setFilter, M_FILTER));
	built = (built && createPortIn<double, MRBF>(this, &mimmo::MRBF::setSupportRadius, M_VALUED));
	built = (built && createPortIn<double, MRBF>(this, &mimmo::MRBF::setSupportRadiusValue, M_VALUED2));
	built = (built && createPortIn<MimmoObject*, MRBF>(&m_geometry, M_GEOM, true));

	built = (built && createPortOut<dmpvecarr3E, MRBF>(this, &mimmo::MRBF::getDisplacements, M_GDISPLS));
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

/*!
 * Set type of solver set for RBF data fields interpolation/parameterization in MRBF::execute.
 * Reimplemented from RBF::setMode() of bitpit;
 * \param[in] solver type of MRBFSol enum;
 */
void
MRBF::setMode(MRBFSol solver){
	m_solver = solver;
	if (m_solver == MRBFSol::NONE)    RBF::setMode(RBFMode::PARAM);
	else                            RBF::setMode(RBFMode::INTERP);
};
/*!
 * Overloading of MRBF::setSolver(MRBFSol solver) with int input parameter
 * Reimplemented from RBF::setMode() of bitpit;
 * \param[in] type of solver 1-WHOLE, 2-GREEDY, see MRBFSol enum;
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
 * \return filter field.
 */
dmpvector1D
MRBF::getFilter(){
	return(m_filter);
};

/*!
 * It gets current support radius ratio (or value if defined as absolute value) as set up in the class.
 * See MRBF::setSupportRadius and MRBF::setSupportRadiusValue method documentation.
 * Reimplemented from bitpit::RBF::getSupportRadius.
 * \return support radius ratio
 */
double
MRBF::getSupportRadius(){
	if (m_supRIsValue)
	{
		return(getSupportRadiusValue());
	}
	return(m_SRRatio);
};

/*!
 * It gets the current real value of support radius for RBF kernels set up in the class.
 * \return
 * support radius value
 */
double
MRBF::getSupportRadiusValue(){
	return(m_SRRatio);
};

/*!
 * It gets if the support radius for RBF kernels is set up as absolute value (true) or
 * ratio of diagonal of bounding box of the geometry (false).
 * \return support radius is set as value?
 */

bool
MRBF::getIsSupportRadiusValue(){
	return(m_supRIsValue);
}

/*!
 * Return actual computed displacements field (if any) for the geometry linked.
 * \return     The computed deformation field on the vertices of the linked geometry
 */
dmpvecarr3E
MRBF::getDisplacements(){
	return m_displ;
};


/*!Adds a RBF point to the total control node list and activate it.
 * Reimplemented from RBF::addNode of bitpit;
 * \param[in] node coordinates of control point.
 * \return RBF id.
 */
long
MRBF::addNode(darray3E node){
	return(RBF::addNode(node));
};

/*!Adds a list of RBF points to the total control node list and activate them.
 * Reimplemented from RBF::addNode of bitpit;
 * \param[in] nodes coordinates of control points.
 * \return Vector of RBF ids.
 */
std::vector<long>
MRBF::addNode(dvecarr3E nodes){
	return(RBF::addNode(nodes));
};

/*!Adds a set of RBF points to the total control node list extracting
 * the vertices stored in a MimmoObject container. Return a vector containing
 * the RBF node int id.
 * Reimplemented from RBF::addNode of bitpit;
 * \param[in] geometry Pointer to MimmoObject that contains the geometry.
 * \return Vector of RBF ids.
 */
livector1D
MRBF::addNode(MimmoObject* geometry){
	if(geometry == NULL)    return    livector1D(0);
	dvecarr3E vertex = geometry->getVerticesCoords();
	return(RBF::addNode(vertex));
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
	if(geometry == NULL)    return ;
	removeAllNodes();
	dvecarr3E vertex = geometry->getVerticesCoords();
	RBF::addNode(vertex);

};

/*! Sets filter field. Note: filter field is defined on nodes of the current linked geometry.
 * coherent size between field size and number of geometry vertices is expected.
 * \param[in] filter fields.
 */
void
MRBF::setFilter(dmpvector1D filter){
	m_filter.clear();
	m_bfilter = !(filter.empty());
	m_filter = filter;
};


/*! Find all possible duplicated nodes within a prescribed distance tolerance.
 * Default tolerance value is 1.0E-12;
 * \param[in] tol distance tolerance
 * \return    list of duplicated nodes.
 */
livector1D
MRBF::checkDuplicatedNodes(double tol){
	livector1D marked;
	long sizeEff = getTotalNodesCount();
	if( sizeEff == 0 ) return marked;

	bvector1D check(sizeEff, false);

	for(long i=0; i<sizeEff; ++i){
		for(long j=i+1; j<sizeEff; ++j){
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
MRBF::removeDuplicatedNodes(livector1D * list){
	livector1D marked;
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
 * \param[in] suppR_ new value of suppR
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
 * \param[in] suppR_ new value of support radius.
 */
void
MRBF::setSupportRadiusValue(double suppR_){
	suppR_ = std::fmax(-1.0,suppR_);
	m_SRRatio = suppR_;
	m_supRIsValue = true;
}

/*!It sets the tolerance for greedy - interpolation algorithm.
 * Tolerance infos are not used in MRBFSol::NONE mode.
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
	long size = displ.size();
	//    if(size != getTotalNodesCount()){
	//        (*m_log) << "warning: " << getName() << " sets displacements with size (" << size << ") that does not fit number of RBF nodes ("<< getTotalNodesCount() << ")" << std::endl;
	//    }

	removeAllData();

	dvector1D temp(size);
	for(int loc=0; loc<3; ++loc){
		for(long i=0; i<size; ++i){
			temp[i] = displ[i][loc];
		}
		addData(temp);
	}
}

/*!
 * Set a field  of 3D displacements on your RBF Nodes defined as points of a MimmoObject. According to MRBFSol mode
 * active in the class set: displacements as direct RBF weights coefficients in MRBFSol::NONE mode,
 * or interpolate displacements to get the best fit weights in other modes MRBFSol::GREEDY/WHOLE
 * Displacements size may not match the actual number of RBF nodes stored in the class.
 * To ensure consistency call fitDataToNodes() method inherited from RBF class.
 *
 * \param[in] displ list of nodal displacements
 */
void
MRBF::setDisplacements(dmpvecarr3E displ){
	long size = displ.size();
	//    if(size != getTotalNodesCount()){
	//        (*m_log) << "warning: " << getName() << " sets displacements with size (" << size << ") that does not fit number of RBF nodes ("<< getTotalNodesCount() << ")" << std::endl;
	//    }

	removeAllData();

	dvector1D temp(size);
	for(int loc=0; loc<3; ++loc){
		long i = 0;
		for (auto val : displ){
			temp[i] = val[loc];
			i++;
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

	default:
            RBF::setFunction(RBFBasisFunction::WENDLANDC2);
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
	RBF::setFunction(bfunc);
    m_functype = -1;
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

/*!Clean all except nodal RBF and its displacements. Use apposite methods RemoveAll*** */
void
MRBF::clear(){
	BaseManipulation::clear();
	clearFilter();
	m_tol = 0.00001;
	m_SRRatio = -1;
};

/*!Clean filter field */
void
MRBF::clearFilter(){
	m_filter.clear();
	m_bfilter = false;
};

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
	//    if(size != getTotalNodesCount()){
	//        (*m_log) << "warning: " << getName() << " sets weights with size (" << size << ") that does not fit number of RBF nodes ("<< getTotalNodesCount() << ")" << std::endl;
	//    }

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

/*!Execution of RBF object. It evaluates the displacements (values) over the point of the
 * linked geometry, given as result of RBF technique implemented in bitpit::RBF base class.
 * The result is stored in the m_displ member.
 *
 */
void
MRBF::execute(){

	MimmoObject * container = getGeometry();
    if(container == NULL){
//        throw std::runtime_error(m_name + "NULL pointer to linked geometry found");
        (*m_log)<<m_name + " : NULL pointer to linked geometry found"<<std::endl;
        return;
    }

    if(container->isEmpty()){
        (*m_log)<<m_name + " : empty linked geometry found"<<std::endl;
    }

    m_displ.clear();
	m_displ.setDataLocation(mimmo::MPVLocation::POINT);
	m_displ.reserve(getGeometry()->getNVertices());
	m_displ.setGeometry(getGeometry());

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
				if(m_solver == MRBFSol::NONE)    data[j] = m_weight[j][i];
				else                            data[j] = m_value[j][i];
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

//	std::cout << " support radius " << radius << std::endl;
//	std::cout << " solving " << std::endl;

	if (m_solver == MRBFSol::WHOLE)    solve();
	if (m_solver == MRBFSol::GREEDY)    greedy(m_tol);

//	std::cout << " solved " << std::endl;

//	for (int i=0; i<sizeF; i++){
//		for (auto val : m_weight[i])
//			std::cout << " wegiht " <<  val << std::endl;
//		for (auto val : m_value[i])
//			std::cout << " value " <<  val << std::endl;
//	}

	dvector1D displ;
	darray3E adispl;
	for(const auto & vertex : container->getVertices()){
		displ = RBF::evalRBF(vertex.getCoords());
		for (int j=0; j<3; ++j)
			adispl[j] = displ[j];
		m_displ.insert(vertex.getId(), adispl);
//		if (norm2(adispl)>0.)
//			std::cout << " insert  id " <<  vertex.getId() << " displ " << adispl << std::endl;
	}

	//if m_filter is active;
	if(m_bfilter){

		checkFilter();

		long int ID;
		for (const auto & vertex : m_geometry->getVertices()){
			ID = vertex.getId();
			m_displ[ID] = m_displ[ID] * m_filter[ID];
		}
	}

//	std::cout << " end execution  " <<  std::endl;

};

/*!
 * Directly apply deformation field to target geometry.
 */
void
MRBF::apply(){

	if (getGeometry() == NULL) return;
	darray3E vertexcoords;
	long int ID;
	for (const auto & vertex : m_geometry->getVertices()){
		vertexcoords = vertex.getCoords();
		ID = vertex.getId();
		if(m_displ.exists(ID))    vertexcoords += m_displ[ID];
		getGeometry()->modifyVertex(vertexcoords, ID);
	}

#if MIMMO_ENABLE_MPI
	getGeometry()->updatePointGhostExchangeInfo();
#endif

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


/*! Plot your current rbf nodes as a point cloud to *vtu file.
 * \param[in] directory output directory
 * \param[in] filename  output filename w/out tag
 * \param[in] counterFile   integer identifier of the file
 * \param[in] binary     boolean flag for 0-"ascii" or 1-"appended" writing
 * \param[in] deformed  boolean flag for plotting 0-"original points", 1-"moved points"
 */
void
MRBF::plotCloud(std::string directory, std::string filename, int counterFile, bool binary, bool deformed){

	long nnodes = getTotalNodesCount();
	nnodes = min(nnodes, long(m_displ.size()));
	dvecarr3E* nodes_ = getNodes();
	dvecarr3E nodes(nnodes);
	dvecarr3E data(nnodes);
	for(long i=0; i<nnodes; ++i){
		for(int j=0; j<3; ++j){
			if(m_solver == MRBFSol::NONE)   data[i][j] = m_weight[j][i];
			else                            data[i][j] = m_value[j][i];
		}
	}
	if(deformed){
		for(long i=0; i<nnodes; ++i){
			nodes[i] = (*nodes_)[i] + data[i];
		}
	}else{
		for(long i=0; i<nnodes; ++i){
			nodes[i] = (*nodes_)[i];
		}
	}

	bitpit::VTKFormat codex = bitpit::VTKFormat::ASCII;
	if(binary){codex=bitpit::VTKFormat::APPENDED;}

	livector1D conn(nnodes);
	{
		long counter = 0;
		for(auto & val: conn){
			val = counter;
			++counter;
		}
	}
	bitpit::VTKUnstructuredGrid vtk(directory, filename, bitpit::VTKElementType::VERTEX);
	vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, nodes) ;
	vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, conn) ;
	vtk.setDimensions(conn.size(), nnodes);
	vtk.setCodex(codex);
	if(counterFile>=0){vtk.setCounter(counterFile);}
	vtk.write();
};

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

	if(slotXML.hasOption("Mode")){
		input = slotXML.get("Mode");
		int value = 0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::string::trim(input));
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
			std::stringstream ss(bitpit::utils::string::trim(input));
			ss >> value;
			setSupportRadius(value);
		}
	};

	if(slotXML.hasOption("SupportRadiusReal")){
		input = slotXML.get("SupportRadiusReal");
		double value = -1.0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::string::trim(input));
			ss >> value;
			setSupportRadiusValue(value);
		}
	};

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

	std::string input;
	input = std::to_string(static_cast<int>(m_solver));
	slotXML.set("Mode", input);


	if(!m_supRIsValue){
		std::stringstream ss;
		ss<<std::scientific<<getSupportRadius();
		slotXML.set("SupportRadius", ss.str());
	}else{
		std::stringstream ss;
		ss<<std::scientific<<getSupportRadiusValue();
		slotXML.set("SupportRadiusReal", ss.str());
	}

	int type = getFunctionType();
	if(type != static_cast<int>(RBFBasisFunction::CUSTOM)){
		slotXML.set("RBFShape", std::to_string(type));
	}

	//checking if not default and if not connected to a port
	if(m_tol != 1.0E-6 ){
		std::stringstream ss;
		ss<<std::scientific<<m_tol;
		slotXML.set("Tolerance", ss.str());
	}
}

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



}
