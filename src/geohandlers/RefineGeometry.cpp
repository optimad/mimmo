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
\*---------------------------------------------------------------------------*/

#include "RefineGeometry.hpp"

using namespace std;
using namespace bitpit;
namespace mimmo{

/*!
 * Default constructor of RefineGeometry.
 */
RefineGeometry::RefineGeometry(){
    m_name  = "mimmo.RefineGeometry";
    m_type	= RefineType(2);
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
RefineGeometry::RefineGeometry(const bitpit::Config::Section & rootXML){

    std::string fallback_name = "ClassNONE";
    std::string input_name = rootXML.get("ClassName", fallback_name);
    input_name = bitpit::utils::string::trim(input_name);

    m_name = "mimmo.RefineGeometry";
    m_type = RefineType(2);

    if(input_name == "mimmo.RefineGeometry"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Default destructor of RefineGeometry.
 */
RefineGeometry::~RefineGeometry(){
    clear();
};

/*!
 * Copy constructor of RefineGeometry.
 */
RefineGeometry::RefineGeometry(const RefineGeometry & other):BaseManipulation(other){
    m_type = other.m_type;
};

/*!
 * Assignement operator of RefineGeometry. Soft copy of MimmoObject
 */
RefineGeometry & RefineGeometry::operator=(RefineGeometry other){
    swap(other);
    return *this;
};


/*!
 * Swap function of RefineGeometry.
 * \param[in] x object to be swapped.
 */
void RefineGeometry::swap(RefineGeometry & x ) noexcept
{
    std::swap(m_type,x.m_type);
    BaseManipulation::swap(x);
};


/*!
 * Building ports of the class
 */
void
RefineGeometry::buildPorts(){
    bool built = true;

	built = (built && createPortIn<MimmoObject*, RefineGeometry>(this, &BaseManipulation::setGeometry, M_GEOM, true));

    built = (built && createPortOut<MimmoObject*, RefineGeometry>(this, &mimmo::RefineGeometry::getGeometry, M_GEOM));
    m_arePortsBuilt = built;
}

/*!
 * Return kind of refinement set for the object.
 * \return refine type
 */
RefineType
RefineGeometry::getRefineType(){
    return m_type;
}

/*!
 * It sets refinement method of refine block
 * \param[in] type refine type
 *
 */
void
RefineGeometry::setRefineType(RefineType type){
	if (type != RefineType::TERNARY)
		throw std::runtime_error(m_name + " : refinement method not allowed");
	m_type = type;
};

/*!
 * It sets partition method of partition block
 * \param[in] type partition method
 *
 */
void
RefineGeometry::setRefineType(int type){
	if (type != 0)
		throw std::runtime_error(m_name + " : refinement method not allowed");

	m_type = RefineType(type);
};

/*!
 * Clear all stuffs in your class
 */
void
RefineGeometry::clear(){
    BaseManipulation::clear();
};


/*!Execution command.
 * It refines the target surface geometry.
 */
void
RefineGeometry::execute(){

    if(getGeometry() == nullptr){
        (*m_log)<<m_name + " : no geometry to refine found"<<std::endl;
        return;
    }

	if (getGeometry()->getType() != 1){
		(*m_log)<<m_name + " : data structure type different from Surface "<<std::endl;
		throw std::runtime_error (m_name + " : data structure type different from Surface");
	}

    if (m_type == RefineType::TERNARY){
    	//TERNARY REFINEMENT

    	// Ternary refinement force the geometry to be a triangulation by placing a new vertex
    	// on the mean point of each cell (the cells have to be convex)

    	// For Now this method works only in SERIAL mode.
    #if MIMMO_ENABLE_MPI
    	if (m_nprocs > 1){
    		//TODO provide implementation to deal with insertion/deletion of vertices and cells in parallel
    		(*m_log)<< "WARNING " <<m_name <<" : is not available yet in parallel process."<<std::endl;
    	}
    #else

    	MimmoObject * geometry = getGeometry();
    	long maxID, newID, newVertID;
    	const auto orderedCellID = geometry->getCells().getIds(true);
    	maxID = orderedCellID[(int)orderedCellID.size()-1];
    	newID = maxID+1;
    	{
    		const auto orderedVertID = geometry->getVertices().getIds(true);
    		newVertID = orderedVertID[(int)orderedVertID.size()-1] +1;
    	}

    	bitpit::ElementType eletype;
    	bitpit::ElementType eletri = bitpit::ElementType::TRIANGLE;
    	livector1D connTriangle(3);

    	for(const auto &idcell : orderedCellID){

    		livector1D conn = geometry->getCellConnectivity(idcell);
    		eletype = geometry->getPatch()->getCell(idcell).getType();
    		long pid = geometry->getPatch()->getCell(idcell).getPID();

    		switch (eletype){
    		case bitpit::ElementType::TRIANGLE:
    		case bitpit::ElementType::PIXEL:
    		case bitpit::ElementType::QUAD:
    		{
    			std::size_t startIndex = 0;
    			std::size_t nnewTri = conn.size() - startIndex;
    			//calculate barycenter and add it as new vertex
    			darray3E barycenter = geometry->getPatch()->evalCellCentroid(idcell);
    			geometry->addVertex(barycenter, newVertID);
    			//delete current polygon
    			geometry->getPatch()->deleteCell(idcell);
    			//insert new triangles from polygon subdivision
    			for(std::size_t i=0; i<nnewTri; ++i){
    				connTriangle[0] = newVertID;
    				connTriangle[1] = conn[ startIndex + std::size_t( i % nnewTri) ];
    				connTriangle[2] = conn[ startIndex + std::size_t( (i+1) % nnewTri ) ];
    				geometry->addConnectedCell(connTriangle, eletri, pid, newID);
    				++newID;
    			}
    			//increment label of vertices
    			++newVertID;
    		}
    		break;
    		case bitpit::ElementType::POLYGON:
    		{
    			std::size_t startIndex = 1;
    			std::size_t nnewTri = conn.size() - startIndex;
    			//calculate barycenter and add it as new vertex
    			darray3E barycenter = geometry->getPatch()->evalCellCentroid(idcell);
    			geometry->addVertex(barycenter, newVertID);
    			//delete current polygon
    			geometry->getPatch()->deleteCell(idcell);
    			//insert new triangles from polygon subdivision
    			for(std::size_t i=0; i<nnewTri; ++i){
    				connTriangle[0] = newVertID;
    				connTriangle[1] = conn[ startIndex + std::size_t( i % nnewTri) ];
    				connTriangle[2] = conn[ startIndex + std::size_t( (i+1) % nnewTri ) ];
    				geometry->addConnectedCell(connTriangle, eletri, pid, newID);
    				++newID;
    			}
    			//increment label of vertices
    			++newVertID;
    		}
    		break;
    		default:
    			throw std::runtime_error("unrecognized cell type in 3D surface mesh of CGNSPidExtractor");
    			break;
    		}
    	}
    #endif
    }

}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void RefineGeometry::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;


	if(slotXML.hasOption("RefineType")){
		std::string input = slotXML.get("RefineType");
		input = bitpit::utils::string::trim(input);
		int value = 2;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setRefineType(value);
	}

    BaseManipulation::absorbSectionXML(slotXML, name);

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void RefineGeometry::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BaseManipulation::flushSectionXML(slotXML, name);
    slotXML.set("RefineType", int(m_type));

};

/*!
 * Plot stitched geometry in *.vtu file as optional result;
 */
void
RefineGeometry::plotOptionalResults(){
    if(getGeometry() == nullptr) return;
    std::string name = m_outputPlot +"/"+ m_name + "_" + std::to_string(getId()) +  "_Patch";
    getGeometry()->getPatch()->write(name);
}

}
