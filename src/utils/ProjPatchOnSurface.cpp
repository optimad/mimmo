/*----------------------------------------------------------------------------*\
 *
 *  mimic
 *
 *  Optimad Engineering S.r.l. ("COMPANY") CONFIDENTIAL
 *  Copyright (c) 2015-2017 Optimad Engineering S.r.l., All Rights Reserved.
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
#include "mimic_common.hpp"

namespace mimmo{

/*!
 * Default constructor
 */
ProjPatchOnSurface::ProjPatchOnSurface(){
    m_name         = "mimmo.ProjPatchOnSurface";
    m_topo     = 2;
    m_cobj = NULL;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ProjPatchOnSurface::ProjPatchOnSurface(const bitpit::Config::Section & rootXML){

    m_name         = "mimmo.ProjPatchOnSurface";
    m_topo     = 2;
    m_cobj = NULL;

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
    ProjPrimitivesOnSurfaces::swap(x);
}

/*!
 * Clear all elements on your current class
 */
void
ProjPatchOnSurface::clear(){
    ProjPrimitivesOnSurfaces::clear();
    m_cobj = NULL;
}

/*!
 * Link as reference surface Patch an external surface mesh passed in a MimmoObject container.
 * If the MimmoObject is not a surface (type 1), nothing will be linked.
 * \param[in] geo pointer to geometry container
 */
void
ProjPatchOnSurface::setPatch(MimmoObject * geo){

    if(geo == NULL)    return;
    auto type = geo->getType();
    if(type != 1 )    return;
    m_cobj = geo;
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

    if(slotXML.hasOption("BvTree")){
        input = slotXML.get("BvTree");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setBuildSkdTree(value);
    };

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

};


/*!
 * Building ports of the class
 */
void
ProjPatchOnSurface::buildPorts(){
    bool built = true;

    built = (built && createPortIn<MimmoObject*, ProjPatchOnSurface>(this, &mimmo::ProjPatchOnSurface::setPatch, M_GEOM2, true));
    m_arePortsBuilt = built;
    ProjPrimitivesOnSurfaces::buildPorts();
}


/*!
 * Core engine for projection.
 */
void
ProjPatchOnSurface::projection(){
    std::unique_ptr<MimmoObject> dum = m_cobj->clone();
    //...and projecting them onto target surface
    if(!getGeometry()->isSkdTreeSync())    getGeometry()->buildSkdTree();
    bitpit::PiercedVector<bitpit::Vertex> & verts = dum->getVertices();
    darray3E point, proj;
    for(auto it = verts.begin(); it != verts.end(); ++it){
        point = it->getCoords();
        proj = skdTreeUtils::projectPoint(&point, getGeometry()->getSkdTree());
        it->setCoords(proj);
    }
    m_patch = std::move(dum);
};




}
