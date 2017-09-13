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

#include "CGNSPidExtractor.hpp"
#include <unordered_map>

using namespace std;
using namespace bitpit;
namespace mimmo{

/*!
 * Default constructor of CGNSPidExtractor.
 */
CGNSPidExtractor::CGNSPidExtractor(){
    m_name         = "mimmo.CGNSPidExtractor";
    m_force = false;
    m_targetpid.insert(0);
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
CGNSPidExtractor::CGNSPidExtractor(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.CGNSPidExtractor";
    m_force = false;
    m_targetpid.insert(0);

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.CGNSPidExtractor"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Default destructor of CGNSPidExtractor.
 */
CGNSPidExtractor::~CGNSPidExtractor(){};

/*!
 * Copy constructor of CGNSPidExtractor.
 * No result patch is copied.
 */
CGNSPidExtractor::CGNSPidExtractor(const CGNSPidExtractor & other):BaseManipulation(other){
    m_force = other.m_force;
    m_targetpid = other.m_targetpid;
};

/*!
 * Assignment operator. No result patch is copied.
 */
CGNSPidExtractor & CGNSPidExtractor::operator=(CGNSPidExtractor other){
    swap(other);
    return *this;
}

/*!
 * Swap function. No result patch member is swapped.
 * \param[in] x object to be swapped.
 */
void CGNSPidExtractor::swap(CGNSPidExtractor & x) noexcept
{
    std::swap(m_force, x.m_force);
    std::swap(m_targetpid, x.m_targetpid);
    BaseManipulation::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
CGNSPidExtractor::buildPorts(){
    bool built = true;
    built = (built && createPortIn<MimmoObject*, CGNSPidExtractor>(this, &CGNSPidExtractor::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_, true));

    built = (built && createPortOut<MimmoObject*, CGNSPidExtractor>(this, &CGNSPidExtractor::getPatch, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    m_arePortsBuilt = built;
};



/*!
 * Return the extracted surface by PID. If not valid PID is found, return an exact copy of the original target geometry.
 * \return pointer to extracted surface
 */
MimmoObject*
CGNSPidExtractor::getPatch(){
    return m_patch.get();
}

/*!
 * Return current active PIDs for selection in the class.
 * \return list of target PIDs
 */
std::set<short>
CGNSPidExtractor::whatPIDActive(){
    return m_targetpid;
}

/*!
 * Return if the class is set to reduce (in execution)
 * the resulting PID-selected patch in a homogeneous straight triangulation,
 * if eventually not.
 * \return the surface has to be triangulated
 */
bool
CGNSPidExtractor::isForcedToTriangulate(){
    return m_force;
}

/*!
 * Set current active PID for selectioning in the class. If PID is not
 * supported by the linked target geoemetry, return a stand-alone copy of the target geometry itself in execution.
 * \param[in]    val target PID to insert in the extraction list
 */
void
CGNSPidExtractor::addPID( short val){
    if(val > 0) m_targetpid.insert(val);
}

/*!
 * Set list of current active PID for selectioning in the class. If no PID in list is not
 * supported by the linked target geoemetry, return a stand-alone copy of the target geometry itself in execution.
 * \param[in]   vval list of target PIDs to use during the extraction
 */
void
CGNSPidExtractor::setPID( std::vector<short int> vval){
    for(const auto & val : vval) addPID(val);
}

/*!
 * Set the class to reduce (in execution) the resulting PID-selected patch in a homogeneous straight triangulation, 
 * if eventually not.
 * \param[in]   flag if true, try to triangulate homogeneously the selected patch.
 */
void
CGNSPidExtractor::setForcedToTriangulate( bool flag){
    m_force = flag;
}

/*!
 * Set current geometry to an external surface boundary mesh, from a mimmo CGNS reader.
 * \param[in]   geo pointer to input geometry.
 */
void
CGNSPidExtractor::setGeometry(MimmoObject * geo){
    if(geo == NULL)    return;
    if(geo->isEmpty())    return;
    if(geo->getType() != 1)    return;
    BaseManipulation::setGeometry(geo);
}


/*!
 * Execution command.
 * Extract patch with target PIDs.
 */
void
CGNSPidExtractor::execute(){
    
    if (getGeometry() == NULL) {
        (*m_log)<<"Error in CGNSPidExtractor: found NULL linked geometry. Does not extract anything."<<std::endl;
        return;
    }
    livector1D extracted;

    for(const auto & val: m_targetpid){
        livector1D temp= getGeometry()->extractPIDCells(val);
        extracted.insert(extracted.end(), temp.begin(), temp.end());
    }
    std::unique_ptr<MimmoObject> patchTemp(new MimmoObject(1));


    MimmoObject * mother = getGeometry();

    if(extracted.empty()){
        patchTemp->setHARDCopy(getGeometry());

    }else{

        for(const auto &val: extracted){

            livector1D conn = mother->getCellConnectivity(val);
            bitpit::ElementInfo::Type eletype = mother->getPatch()->getCell(val).getType();
            patchTemp->addConnectedCell(conn, eletype, val);
        }

        {
            darray3E temp;
            livector2D tempConn = patchTemp->getConnectivity();
            std::map<long,long> ordIndex;

            for(const auto & val2: tempConn){
                for (auto val22 : val2){
                    ordIndex[val22] = val22;
                }
            }

            for(const auto & val3: ordIndex){
                temp = mother->getVertexCoords(val3.second);
                patchTemp->addVertex(temp,val3.second);
            }
        }

    }

    //check if patch is forced to be triangulated
    if(m_force){

        long maxID, newID;

        const auto orderedCellID = patchTemp->getCells().getIds(true);
        maxID = orderedCellID[(int)orderedCellID.size()-1];
        newID = maxID+1;
        bitpit::ElementInfo::Type eletype, eletri = bitpit::ElementInfo::Type::TRIANGLE;

        for(const auto &idcell : orderedCellID){

            livector1D conn = patchTemp->getCellConnectivity(idcell);
            eletype = mother->getPatch()->getCell(idcell).getType();

            if(eletype == bitpit::ElementInfo::Type::QUAD){
                //create new triangle connectivity
                livector1D conn1(3), conn2(3);
                conn1[0] = conn[0];
                conn1[1] = conn[1];
                conn1[2] = conn[2];

                conn2[0] = conn[0];
                conn2[1] = conn[2];
                conn2[2] = conn[3];

                patchTemp->getPatch()->deleteCell(idcell);
                patchTemp->addConnectedCell(conn1, eletri, idcell);
                patchTemp->addConnectedCell(conn2, eletri, newID);
                ++newID;
            }
        }
    }

    m_patch = std::move(patchTemp);
    mother=NULL;
};


/*!
 * Plot optional result of the class in execution, that is the selected patch
 * as standard vtk unstructured grid.
 */
void CGNSPidExtractor::plotOptionalResults(){
    if(getPatch() == NULL)    return;
    if(getPatch()->isEmpty()) return;
    std::string name = m_name + "_" + std::to_string(getClassCounter()) +  "_Patch";
    getPatch()->getPatch()->write(name);
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void 
CGNSPidExtractor::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    BaseManipulation::absorbSectionXML(slotXML, name);
    
    std::vector<short> temp;
    if(slotXML.hasOption("nPID")){
        input = slotXML.get("nPID");
        int value = 0;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
            value = std::max(0, value);
        }
        temp.resize(value, -1);
    };

    if(slotXML.hasOption("PID") && !temp.empty()){
        input = slotXML.get("PID");
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            for(auto & val: temp){
                ss >> val;
            }
        }
        setPID(temp);
    };

    if(slotXML.hasOption("ForcedToTriangulate")){
        input = slotXML.get("ForcedToTriangulate");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setForcedToTriangulate(value);
    };

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void 
CGNSPidExtractor::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    
    BaseManipulation::flushSectionXML(slotXML, name);
    
    slotXML.set("nPID", std::to_string(int(m_targetpid.size())));

    std::stringstream output;
    auto it = m_targetpid.begin();
    auto itE = m_targetpid.end();
    while(it != itE)    {
        output<<*it<<" ";
        ++it;
    }    
    slotXML.set("PID", output.str());
    slotXML.set("ForcedToTriangulate", std::to_string(int(m_force)));


};


}
