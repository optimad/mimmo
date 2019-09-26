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

#include "ProjSegmentOnSurface.hpp"
#include "mimmo_core.hpp"

using namespace std;
using namespace bitpit;
namespace mimmo{

/*!
 * Default constructor
 */
ProjSegmentOnSurface::ProjSegmentOnSurface(){
    m_name         = "mimmo.ProjSegmentOnSurface";
    m_topo     = 1;
    m_pointA.fill(0.0);
    m_pointB.fill(1.0);
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ProjSegmentOnSurface::ProjSegmentOnSurface(const bitpit::Config::Section & rootXML){

    m_name         = "mimmo.ProjSegmentOnSurface";
    m_topo     = 1;
    m_pointA.fill(0.0);
    m_pointB.fill(1.0);

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.ProjSegmentOnSurface"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Default destructor 
 */
ProjSegmentOnSurface::~ProjSegmentOnSurface(){
    clear();
};

/*!
 * Copy constructor. No resulting projected patch is copied
 */
ProjSegmentOnSurface::ProjSegmentOnSurface(const ProjSegmentOnSurface & other):ProjPrimitivesOnSurfaces(other){
    m_pointA = other.m_pointA;
    m_pointB = other.m_pointB;
};

/*!
 * Assignement operator. No resulting projected patch is copied
 */
ProjSegmentOnSurface & ProjSegmentOnSurface::operator=(ProjSegmentOnSurface other){
    swap(other);
    return *this;
};

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void ProjSegmentOnSurface::swap(ProjSegmentOnSurface & x)   noexcept
{
    std::swap(m_pointA, x.m_pointA);
    std::swap(m_pointB, x.m_pointB);
    ProjPrimitivesOnSurfaces::swap(x);
}

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

    if(norm2(pointB - pointA) < 1.E-18)    return;

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
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ProjSegmentOnSurface::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;
    BaseManipulation::absorbSectionXML(slotXML, name);
    
    if(slotXML.hasOption("Segment")){
        input = slotXML.get("Segment");
        darray3E p1,p2;
        p1.fill(0.0);
        p2.fill(1.0);
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> p1[0]>>p1[1]>>p1[2]>>p2[0]>>p2[1]>>p2[2];
        }
        setSegment(p1,p2);
    };
    if(slotXML.hasOption("nCells")){
        input = slotXML.get("nCells");
        int value = 1000;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setProjElementTargetNCells(value);
    };


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
void ProjSegmentOnSurface::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);
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

    if(m_buildSkdTree){
        output = std::to_string(m_buildSkdTree);
        slotXML.set("SkdTree", output);
    }

    if(m_buildKdTree){
        output = std::to_string(m_buildKdTree);
        slotXML.set("KdTree", output);
    }
};


/*!
 * Core engine for projection.
 */
void
ProjSegmentOnSurface::projection(){

    int counter = 0;
    std::unique_ptr<MimmoObject> dum(new MimmoObject(4));
    //reserving memory
    dum->getPatch()->reserveVertices(m_nC+1);
    dum->getPatch()->reserveCells(m_nC);

    //start filling connectivity of your object.
    bitpit::ElementType eltype = bitpit::ElementType::LINE;

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
        auto dx = (m_pointB - m_pointA);
        dx /= double(m_nC);
        counter = 0;
        for( auto && ele : verts){
            ele  = m_pointA + double(counter)*dx;
            ++counter;
        }
    }

    //...and projecting them onto target surface
    if(!getGeometry()->isSkdTreeSync())    getGeometry()->buildSkdTree();
    counter = 0;
    projs.resize(verts.size());
    for(auto &val : verts){
        projs[counter]= skdTreeUtils::projectPoint(&val, getGeometry()->getSkdTree());
        ++counter;
    }

    //storing the projected points in the MImmoObject:
    long idS = 0;
    for(const auto & vv : projs){
        dum->addVertex(vv, idS);
        ++idS;
    }

    dum->cleanGeometry();
    m_patch = std::move(dum);

};

}
