/*----------------------------------------------------------------------------*\
 *
 *  mimmo
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
#include <SkdTreeUtils.hpp>
#if MIMMO_ENABLE_MPI
    #include "Partition.hpp"
#endif
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
        std::stringstream ss;
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

    //create points on segment first ...
    dvecarr3E points;
#if MIMMO_ENABLE_MPI
    if(m_rank == 0)
#endif
    {
        //create the vertices array, ordered from pointA to pointB
        points.resize(m_nC+1);
        darray3E dx = (m_pointB - m_pointA);
        dx /= double(m_nC);
        int counter = 0;
        for( darray3E & ele : points){
            ele  = m_pointA + double(counter)*dx;
            ++counter;
        }
    }

    //...and projecting them onto target surface
    getGeometry()->buildSkdTree();

    std::size_t npoints = points.size();
    dvecarr3E projs(npoints);
    livector1D ids(npoints);
    ivector1D ranks(npoints,0);
#if MIMMO_ENABLE_MPI
    double radius =  std::numeric_limits<double>::max();
    skdTreeUtils::projectPointGlobal(npoints, points.data(), getGeometry()->getSkdTree(), projs.data(), ids.data(), ranks.data(), radius, false);
#else
    skdTreeUtils::projectPoint(npoints, points.data(), getGeometry()->getSkdTree(), projs.data(), ids.data());
#endif

    //you have now on master 0 ranks where points project, and projected values.
    //istantiate object and fill data of 3D curve
    {
        MimmoSharedPointer<MimmoObject> dum(new MimmoObject(4));
        m_patch = dum;
    }

#if MIMMO_ENABLE_MPI
    //take track of a partition map using ranks of projection.
    std::unordered_map<long, int> partMap;
    if (m_rank ==0) //with rank 0 only in MPI
#endif
    {
        //reserving memory
        m_patch->getPatch()->reserveVertices(m_nC+1);
        m_patch->getPatch()->reserveCells(m_nC);

        //storing the projected points
        long id(0);
        for(darray3E & vv : projs){
            m_patch->addVertex(vv, id);
            ++id;
        }


        //start filling connectivity of your object.
        bitpit::ElementType eltype = bitpit::ElementType::LINE;
        id = 0;
        for(int i=0; i<m_nC; ++i){
            m_patch->addConnectedCell(livector1D({{i, i+1}}), eltype, 0, id, 0);
#if MIMMO_ENABLE_MPI
            partMap[id]= std::min(ranks[i], ranks[i+1]);
#endif
            ++id;
        }
    }
    m_patch->update();

#if MIMMO_ENABLE_MPI
    std::unique_ptr<mimmo::Partition> part(new mimmo::Partition);
    part->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    part->setGeometry(m_patch);
    part->setPartition(partMap);
    part->execute();
    //m_patch will return partitioned and operative.
#endif

};

}
