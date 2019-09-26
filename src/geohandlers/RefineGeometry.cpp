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
    m_mode	= RefineMode(2);
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
    m_mode = RefineMode(2);

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
    m_mode = other.m_mode;
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
    std::swap(m_mode,x.m_mode);
    BaseManipulation::swap(x);
};


/*!
 * Building ports of the class
 */
void
RefineGeometry::buildPorts(){
    bool built = true;

	built = (built && createPortIn<MimmoObject*, RefineGeometry>(this, &BaseManipulation::setGeometry, M_GEOM, true));
    built = (built && createPortIn<MimmoObject*, RefineGeometry>(this, &mimmo::RefineGeometry::addExternalGeometry, M_GEOM2));

    built = (built && createPortOut<MimmoObject*, RefineGeometry>(this, &mimmo::RefineGeometry::getGeometry, M_GEOM));
    m_arePortsBuilt = built;
}

/*!
 * Return kind of refinement set for the object.
 * \return refine mode
 */
RefineMode
RefineGeometry::getRefineMode(){
    return m_mode;
}

/*!
 * It sets refinement method of refine block
 * \param[in] mode refine mode
 *
 */
void
RefineGeometry::setRefineMode(RefineMode mode){
	if (mode != RefineMode::GLOBAL && mode != RefineMode::SELECTION && mode != RefineMode::LINE)
		throw std::runtime_error(m_name + " : refinement method not allowed");
	m_mode = mode;
};

/*!
 * It sets partition method of partition block
 * \param[in] mode partition method
 *
 */
void
RefineGeometry::setRefineMode(int mode){
	if (mode != 0 && mode != 1 && mode != 2)
		throw std::runtime_error(m_name + " : refinement method not allowed");

	m_mode = RefineMode(mode);
};

/*!
 * Add an external geometry to be used for refining. Topology of the geometry must be coeherent
 * with refinement mode of the class;
 * \param[in] geo  Pointer to MimmoObject
 */
void
RefineGeometry::addExternalGeometry(MimmoObject* geo){
    if(geo == nullptr) return;
    if(geo->getType() != 1 && geo->getType() != 4){
		(*m_log)<<m_name + " : warning, topology of provided external geometry different from surface or line."<<std::endl;
		return;
    }
    else if (geo->getType() == 1){
    	if (m_mode != RefineMode::SELECTION){
    		(*m_log)<<m_name + " : warning, topology of provided external geometry surface with mode different from selection."<<std::endl;
    		return;
    	}
    }
    else if (geo->getType() == 4){
    	if (m_mode != RefineMode::LINE){
    		(*m_log)<<m_name + " : warning, topology of provided external geometry 3D line with mode different from line."<<std::endl;
    		return;
    	}
    }

    //Don't insert a duplicated geometry
    if(std::count(m_extgeo.begin(), m_extgeo.end(), geo))
    	return;

    m_extgeo.push_back(geo);

};


/*!
 * Clear all stuffs in your class
 */
void
RefineGeometry::clear(){
    m_extgeo.clear();
    BaseManipulation::clear();
};


/*!Execution command.
 * It refines the target surface geometry.
 */
void
RefineGeometry::execute(){
    if(m_extgeo.empty()){
        (*m_log)<<m_name + " : no external geometries to refine were found"<<std::endl;
    }


    if (m_mode == RefineMode::LINE){

    	// Loop on external geometries
    	for (MimmoObject* line : m_extgeo){

    		//Compute projection of line segments on geometry
    		for (bitpit::Cell & segment : line->getCells()){

    		}// end segments loop

    	}// end external lines loop
    }






//    std::unique_ptr<MimmoObject> dum(new MimmoObject(m_topo));
//#if MIMMO_ENABLE_MPI
//    if(m_nprocs > 1){
//        //TODO you need a strategy to stitch together partioned mesh, keeping a unique id
//        //throughout cells and ids. For example, communicate the number of Global cells and Global verts
//        //(or min/max ids) of each patch to all communicators, and using them to organize offsets
//        //for safe inserting elements.
//        (*m_log)<<"WARNING "<< m_name<<" : stitching not available yet for MPI version with procs > 1"<<std::endl;
//        m_patch = std::move(dum);
//        return;
//    }
//#endif
//
//    long nCells = 0;
//    long nVerts = 0;
//
//    for(auto &obj : m_extgeo){
//        nCells += obj.first->getNCells();
//        nVerts += obj.first->getNVertices();
//    }
//
//    //reserving memory
//    dum->getPatch()->reserveVertices(nVerts);
//    dum->getPatch()->reserveCells(nCells);
//
//    long cV = 0, cC = 0;
//
//    std::unordered_map<MimmoObject*,long>    map_pidstart;
//    //initialize map
//    for(auto & obj : m_extgeo){
//        map_pidstart[obj.first] = 0;
//    }
//
//    if(m_repid){
//        long pidmax = -1;
//        for(auto & obj : m_extgeo){
//            auto pidlist = obj.first->getPIDTypeList();
//            std::vector<long> temp(pidlist.begin(), pidlist.end());
//            std::sort(temp.begin(), temp.end());
//            map_pidstart[obj.first] = pidmax + 1;
//            pidmax += (*(temp.rbegin())+1);
//        }
//    }
//
//
//    //start filling your stitched object.
//    {
//        //optional vars;
//        long vId, cId;
//        long PID;
//        bitpit::ElementType eltype;
//
//        for(auto &obj : m_extgeo){
//
//            std::unordered_map<long,long> mapVloc;
//            auto originalPidNames = obj.first->getPIDTypeListWNames();
//            std::unordered_map<long, long> newOldPid;
//            for(const auto & val : obj.first->getPIDTypeList()){
//                newOldPid[val + map_pidstart[obj.first]] = val;
//            }
//            //start extracting/reversing vertices of the current obj
//            for(const auto & vv : obj.first->getVertices()){
//                vId = vv.getId();
//                dum->addVertex(obj.first->getVertexCoords(vId), cV);
//                //update map;
//                mapVloc[vId] = cV;
//                cV++;
//            }
//
//            //start extracting/reversing cells of the current obj
//            for(const auto & cc : obj.first->getCells()){
//                cId = cc.getId();
//                PID = cc.getPID() + map_pidstart[obj.first];
//                eltype = cc.getType();
//                //get the local connectivity and update with new vertex numbering;
//                livector1D conn = obj.first->getCellConnectivity(cId);
//                livector1D connloc(conn.size());
//
//                if(eltype == bitpit::ElementType::POLYGON){
//                    std::size_t size = conn.size();
//                    connloc[0] = conn[0];
//                    for(std::size_t i = 1; i < size; ++i){
//                        connloc[i] = mapVloc[conn[i]];
//                    }
//
//                }else if(eltype == bitpit::ElementType::POLYHEDRON){
//                    connloc[0] = conn[0];
//                    for(int nF = 0; nF < conn[0]-1; ++nF){
//                        int facePos = cc.getFaceStreamPosition(nF);
//                        int beginVertexPos = facePos + 1;
//                        int endVertexPos   = facePos + 1 + conn[facePos];
//                        connloc[facePos] = conn[facePos];
//                        for (int i=beginVertexPos; i<endVertexPos; ++i){
//                            connloc[i] = mapVloc[conn[i]];
//                        }
//                    }
//                }else{
//                    int ic = 0;
//                    for (const auto & v : conn){
//                        connloc[ic] = mapVloc[v];
//                        ++ic;
//                    }
//                }
//                dum->addConnectedCell(connloc, eltype, PID, cC);
//                //update map;
//                cC++;
//            }
//
//            dum->resyncPID();
//            for(const auto & val: dum->getPIDTypeList()){
//                dum->setPIDName(val, originalPidNames[newOldPid[val]]);
//            }
//        }
//    }//scope for optional vars;
//
//    m_patch = std::move(dum);
//    m_patch->cleanGeometry();
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void RefineGeometry::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;


	if(slotXML.hasOption("RefineMode")){
		std::string input = slotXML.get("RefineMode");
		input = bitpit::utils::string::trim(input);
		int value = 2;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setRefineMode(value);
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
    slotXML.set("RefineMode", int(m_mode));

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
