/*----------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Optimad Engineering S.r.l. ("COMPANY") CONFIDENTIAL
 *  Copyright (c) 2015-2021 Optimad Engineering S.r.l., All Rights Reserved.
 *
 *  --------------------------------------------------------------------------
 *
 *  NOTICE:  All information contained herein is, and remains the property
 *  of COMPANY. The intellectual and technical concepts contained herein are
 *  proprietary to COMPANY and may be covered by Italian and Foreign Patents,
 *  patents in process, and are protected by trade secret or copyright law.
 *  Dissemination of this information or reproduction of this material is
 *  strictly forbidden unless prior written permission is obtained from
 *  COMPANY. Access to the source code contained herein is hereby forbidden
 *  to anyone except current COMPANY employees, managers or contractors who
 *  have executed Confidentiality and Non-disclosure agreements explicitly
 *  covering such access.
 *
 *  The copyright notice above does not evidence any actual or intended
 *  publication or disclosure of this source code, which includes information
 *  that is confidential and/or proprietary, and is a trade secret, of
 *  COMPANY. ANY REPRODUCTION, MODIFICATION, DISTRIBUTION, PUBLIC PERFORMANCE,
 *  OR PUBLIC DISPLAY OF OR THROUGH USE  OF THIS  SOURCE CODE  WITHOUT THE
 *  EXPRESS WRITTEN CONSENT OF COMPANY IS STRICTLY PROHIBITED, AND IN
 *  VIOLATION OF APPLICABLE LAWS AND INTERNATIONAL TREATIES. THE RECEIPT OR
 *  POSSESSION OF THIS SOURCE CODE AND/OR RELATED INFORMATION DOES NOT CONVEY
 *  OR IMPLY ANY RIGHTS TO REPRODUCE, DISCLOSE OR DISTRIBUTE ITS CONTENTS, OR
 *  TO MANUFACTURE, USE, OR SELL ANYTHING THAT IT  MAY DESCRIBE, IN WHOLE OR
 *  IN PART.
 *
\*----------------------------------------------------------------------------*/

#include "ProjPatchOnSurface.hpp"
#include <SkdTreeUtils.hpp>

namespace mimmo{

/*!
 * Default constructor
 */
ProjPatchOnSurface::ProjPatchOnSurface(){
    m_name         = "mimmo.ProjPatchOnSurface";
    m_topo     = 2; //2D- surface tessellation
    m_cobj = nullptr;
    m_workingOnTarget = false;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ProjPatchOnSurface::ProjPatchOnSurface(const bitpit::Config::Section & rootXML){

    m_name         = "mimmo.ProjPatchOnSurface";
    m_topo     = 2; //2D- surface tessellation
    m_cobj = nullptr;
    m_workingOnTarget = false;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.ProjPatchOnSurface"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Default destructor
 */
ProjPatchOnSurface::~ProjPatchOnSurface(){
    clear();
};

/*!Copy constructor. No resulting projected patch is copied
 */
ProjPatchOnSurface::ProjPatchOnSurface(const ProjPatchOnSurface & other):ProjPrimitivesOnSurfaces(other){
    m_cobj = other.m_cobj;
    m_workingOnTarget = other.m_workingOnTarget;

};

/*!
 * Assignement operator. No resulting projected patch is copied
 */
ProjPatchOnSurface & ProjPatchOnSurface::operator=(ProjPatchOnSurface other){
    swap(other);
    return *this;
};

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void ProjPatchOnSurface::swap(ProjPatchOnSurface &x) noexcept
{
    std::swap(m_cobj, x.m_cobj);
    std::swap(m_workingOnTarget, x.m_workingOnTarget);
    ProjPrimitivesOnSurfaces::swap(x);
}

/*!
 * Clear all elements on your current class
 */
void
ProjPatchOnSurface::clear(){
    ProjPrimitivesOnSurfaces::clear();
    m_cobj = nullptr;
    m_workingOnTarget = false;
}

/*!
 * Link the target Patch to be projected on the reference surface.
 * Any patch type is allowed, except for volume meshes (MimmoObject of type 2).
   In that case nothing will be linked.
 * \param[in] geo pointer to patch that need to be projected
 */
void
ProjPatchOnSurface::setPatch(MimmoSharedPointer<MimmoObject> geo){

    if(geo == nullptr)    return;
    int type = geo->getType();
    if(type == 2 ) return;
    m_cobj = geo;
    m_topo = 2;
    if(type == 3) m_topo = 0;
    if(type == 4) m_topo = 1;

}

/*!
    Return true, if the class will apply the modifications directly to the target surface
    container linked with setPatch. False, if projection is returned in an indipendent container.
    \return true or false if is working on target.
*/
bool
ProjPatchOnSurface::isWorkingOnTarget(){
    return m_workingOnTarget;
}

/*!
    Set to true, if the class will apply the modifications directly to the target surface
    container linked with setPatch. To False instead, if projection is returned in an indipendent container.
    param[in] flag true or false to set class working on target.
*/
void
ProjPatchOnSurface::setWorkingOnTarget(bool flag){
    m_workingOnTarget = flag;
}


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ProjPatchOnSurface::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    BaseManipulation::absorbSectionXML(slotXML, name);

     if(slotXML.hasOption("SkdTree")){
        input = slotXML.get("SkdTree");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setBuildSkdTree(value);
    };

    if(slotXML.hasOption("KdTree")){
        input = slotXML.get("KdTree");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setBuildKdTree(value);
    };

    if(slotXML.hasOption("WorkingOnTarget")){
       input = slotXML.get("WorkingOnTarget");
       bool value = false;
       if(!input.empty()){
           std::stringstream ss(bitpit::utils::string::trim(input));
           ss >> value;
       }
       setWorkingOnTarget(value);
   };

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ProjPatchOnSurface::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    if(m_buildSkdTree){
        slotXML.set("SkdTree", std::to_string(m_buildSkdTree));
    }

    if(m_buildKdTree){
        slotXML.set("KdTree", std::to_string(m_buildKdTree));
    }
    slotXML.set("WorkingOnTarget", std::to_string(m_workingOnTarget));

};


/*!
 * Building ports of the class
 */
void
ProjPatchOnSurface::buildPorts(){
    ProjPrimitivesOnSurfaces::buildPorts();
    bool built = m_arePortsBuilt;

    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, ProjPatchOnSurface>(this, &mimmo::ProjPatchOnSurface::setPatch, M_GEOM2, true));
    m_arePortsBuilt = built;
}


/*!
 * Core engine for projection.
 */
void
ProjPatchOnSurface::projection(){
    if(m_workingOnTarget){
        m_patch = m_cobj;
    }else{
        m_patch = m_cobj->clone();
    }
    //...and projecting them onto target surface
    getGeometry()->buildSkdTree();

    dvecarr3E points;
    livector1D originalIds;
    points.reserve(m_patch->getNVertices());
    originalIds.reserve(points.size());
    for(bitpit::Vertex & vert : m_patch->getVertices()){
        points.push_back(vert.getCoords());
        originalIds.push_back(vert.getId());
     }

    std::size_t npoints = points.size();
    dvecarr3E projs(npoints);
    livector1D cell_ids(npoints);
#if MIMMO_ENABLE_MPI
    ivector1D ranks(npoints);
    double radius =  std::numeric_limits<double>::max();
    skdTreeUtils::projectPointGlobal(npoints, points.data(), getGeometry()->getSkdTree(), projs.data(), cell_ids.data(), ranks.data(), radius, false);
#else
    skdTreeUtils::projectPoint(npoints, points.data(), getGeometry()->getSkdTree(), projs.data(), cell_ids.data());
#endif

    std::size_t counter = 0;
    for(long id: originalIds){
        m_patch->modifyVertex(projs[counter], id);
        counter++;
    }

    m_patch->cleanPatchInfo();
#if MIMMO_ENABLE_MPI
    m_patch->resetPointGhostExchangeInfo();
#endif
    m_patch->update();

};


}
