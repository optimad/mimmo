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

#include "Proj3DCurveOnSurface.hpp"
#include <SkdTreeUtils.hpp>
#include <queue>

namespace mimmo{

/*!
 * Default constructor
 */
Proj3DCurveOnSurface::Proj3DCurveOnSurface(){
    m_name         = "mimmo.Proj3DCurveOnSurface";
    m_topo     = 1;
    m_closed = false;
    m_cobj = NULL;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
Proj3DCurveOnSurface::Proj3DCurveOnSurface(const bitpit::Config::Section & rootXML){

    m_name         = "mimmo.Proj3DCurveOnSurface";
    m_topo     = 1;
    m_closed = false;
    m_cobj = NULL;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.Proj3DCurveOnSurface"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Default destructor
 */
Proj3DCurveOnSurface::~Proj3DCurveOnSurface(){
    clear();
};

/*!Copy constructor. No resulting projected patch is copied
 */
Proj3DCurveOnSurface::Proj3DCurveOnSurface(const Proj3DCurveOnSurface & other):ProjPrimitivesOnSurfaces(other){
    m_cpoints = other.m_cpoints;
    m_closed = other.m_closed;
    m_cobj = other.m_cobj;

};

/*!
 * Assignement operator. No resulting projected patch is copied
 */
Proj3DCurveOnSurface & Proj3DCurveOnSurface::operator=(Proj3DCurveOnSurface other){
    swap(other);
    return *this;
};

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void Proj3DCurveOnSurface::swap(Proj3DCurveOnSurface &x) noexcept
{
    std::swap(m_cpoints, x.m_cpoints);
    std::swap(m_closed , x.m_closed);
    std::swap(m_cobj, x.m_cobj);
    ProjPrimitivesOnSurfaces::swap(x);
}

/*!
 * Clear all elements on your current class
 */
void
Proj3DCurveOnSurface::clear(){
    ProjPrimitivesOnSurfaces::clear();
    m_cpoints.clear();
    m_cobj = NULL;
}

/*!
 * Adding a point defining your curve. Consecutive points will describe it.
 * The method works appending the point to the stored list of points.
 * If this method is selected, whatever is stored with
 * Proj3DCurveOnSurface::setConnectedPoints will be erased
 * \param[in] point  coordinates of the point
 */
void
Proj3DCurveOnSurface::addPoint(darray3E point){
    m_cpoints.push_back(point);
    m_cobj = NULL;
}

/*!
 * Set a list of consecutive points defining your curve.
 * The method works substituting the current list to the stored list of points.
 * If this method is selected, whatever is stored with
 * Proj3DCurveOnSurface::setConnectedPoints will be erased.
 * \param[in] points  3D points list
 */
void
Proj3DCurveOnSurface::setPoints(dvecarr3E points){
    m_cpoints.clear();
    m_cpoints = points;
    m_cobj = NULL;
}

/*!
 * Link as reference 3D curve an external 3D curve passed in a MimmoObject container.
 * If the MimmoObject is not a 3D curve/cloud Point or is empty, nothing will be linked.
 * If a valid MimmoObject is linked, all points already stored will be erased.
 * \param[in] geo pointer to geometry container
 */
void
Proj3DCurveOnSurface::setConnectedPoints(MimmoObject * geo){

    if(geo == NULL)    return;
    if(geo->isEmpty()) return;
    auto type = geo->getType();
    if(type != 3 && type != 4)    return;
    m_cobj = geo;

    if(type == 4){
        m_closed = geo->isClosedLoop();
    }
    m_cpoints.clear();
}


/*!
 * Return if the 3D Curve is defined as a closed loop
 * \return is the 3D Curve a closed loop?
 */
bool
Proj3DCurveOnSurface::isClosedLoop(){
    if(m_cobj != NULL && m_cobj->getType() == 4)    return m_cobj->isClosedLoop();
    return m_closed;
}

/*!
 * Set your 3D curve as open/closed. If a 3D Curve valid geometry is connected through the
 * method Proj3DCurveOnSurface::setConnectedPoints the method ignores the input parameters
 * and retrieve the information directly from 3D connected tessellation.
 * \param[in] flag false for open curve, true for closed one.
 */
void
Proj3DCurveOnSurface::setClosedLoop(bool flag){
    if(m_cobj != NULL && m_cobj->getType() == 4) m_closed = m_cobj->isClosedLoop();
    m_closed = flag;
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void Proj3DCurveOnSurface::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("ClosedLoop")){
        input = slotXML.get("ClosedLoop");
        bool flag = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> flag;
        }
        setClosedLoop(flag);
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
void Proj3DCurveOnSurface::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("ClosedLoop", m_topo);
    std::string output;

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
 * Building ports of the class
 */
void
Proj3DCurveOnSurface::buildPorts(){

    ProjPrimitivesOnSurfaces::buildPorts();
    bool built = m_arePortsBuilt;

    built = (built && createPortIn<dvecarr3E, Proj3DCurveOnSurface>(this, &mimmo::Proj3DCurveOnSurface::setPoints, M_COORDS, true,1));
    built = (built && createPortIn<darray3E, Proj3DCurveOnSurface>(this, &mimmo::Proj3DCurveOnSurface::addPoint, M_POINT,true, 1));
    built = (built && createPortIn<MimmoObject*, Proj3DCurveOnSurface>(this, &mimmo::Proj3DCurveOnSurface::setConnectedPoints, M_GEOM2, true, 1));
    m_arePortsBuilt = built;

}


/*!
 * Core engine for projection.
 */
void
Proj3DCurveOnSurface::projection(){

    std::unique_ptr<MimmoObject> dum(new MimmoObject(4));

    dvecarr3E points;
    std::unordered_map<long, std::array<long,2> > connectivity;
    int furtherCells = fillPreliminaryStructure(points, connectivity);
    refineObject(points, connectivity, furtherCells);


    //...and projecting them onto target surface
    if(!getGeometry()->isSkdTreeSync())    getGeometry()->buildSkdTree();
    int counter=0;
    dvecarr3E projs(points.size());
    for(auto &val : points){
        projs[counter]= skdTreeUtils::projectPoint(&val, getGeometry()->getSkdTree());
        ++counter;
    }

    dum->getPatch()->reserveVertices(points.size());
    dum->getPatch()->reserveCells(connectivity.size());

    //start filling connectivity of your object.
    bitpit::ElementType eltype = bitpit::ElementType::LINE;

    for(auto & con: connectivity){
        livector1D con2(2);
        con2[0] = con.second[0];
        con2[1] = con.second[1];
        dum->addConnectedCell(con2, eltype);
    }

    //storing the projected points in the MimmoObject:
    long idS = 0;
    for(const auto & vv : projs){
        dum->addVertex(vv, idS);
        ++idS;
    }

    dum->cleanGeometry();
    m_patch = std::move(dum);
};


/*!
 * Fill a points/connectivity structure with points data set by the User. Two situation are
 * possible. If a point cloud not connected is available, use it as it is, supposing that
 * each point is displaced consecutively to form the final curve. If a connected curve is provided,
 * copy as it is in the target container.
 * \param[in]    points  reference to structure of 3D points
 * \param[in]    connectivity  reference to connectivity 3D curve
 * \return number of further cells to add required by the class;
 */
int
Proj3DCurveOnSurface::fillPreliminaryStructure(dvecarr3E & points, std::unordered_map<long, std::array<long,2> > &connectivity){
    int nVerts, nCells, res;
    long id;
    liimap mapDataInv;
    if (m_cobj) {
        nVerts = m_cobj->getNVertices();
        mapDataInv = m_cobj->getMapDataInv();
        if(m_cobj->getType() == 4){
            nCells = m_cobj->getNCells();
            res = std::max(0, (m_nC - nCells));
            auto conns= m_cobj->getCompactConnectivity(mapDataInv);
            id=0;
            for(auto & val: conns){
                for(int i=0; i<2; ++i){
                    connectivity[id][i] = val[i];
                }
                ++id;
            }
        }else{
            nCells = nVerts-1+(int)m_closed;
            res = std::max(0, (m_nC - nVerts));
            for(int i=0; i<(nCells -(int)m_closed); ++i){
                connectivity[i][0] = i;
                connectivity[i][1] = i+1;
            }
            if(m_closed){
                connectivity[nCells-1][0] = nCells-1;
                connectivity[nCells-1][1] = 0;
            }
        }
        points = m_cobj->getVerticesCoords();

    }else{
        nVerts = (int) m_cpoints.size();
        nCells = nVerts-1+(int)m_closed;
        res = std::max(0, (m_nC - nVerts));
        points = m_cpoints;
        for(int i=0; i<(nCells -(int)m_closed); ++i){
            connectivity[i][0] = i;
            connectivity[i][1] = i+1;
        }
        if(m_closed){
            connectivity[nCells-1][0] = nCells-1;
            connectivity[nCells-1][1] = 0;
        }
    }

    return (res);
}

/*!
 * Refine 3D curve object adding a number of given cells.
 * Each segment cell will be splitted in two new cells. Priority will be given
 * to longest segments
 * * \param[in]    points  reference to structure of 3D points
 * \param[in]    connectivity  reference to connectivity 3D curve
 * \param[in]   fCells number of cells to add.
 */
void
Proj3DCurveOnSurface::refineObject(dvecarr3E & points, std::unordered_map<long, std::array<long,2> > &connectivity, int fCells){
    if ( fCells == 0) return;

    greatDist mycomp;

    double dist;
    int idV = (int)points.size() -1;
    long idC = 0;
    std::vector<std::pair<double, long>> values;
    values.reserve(int(points.size()) + fCells);
    for (auto &conn : connectivity){
        dist = norm2(points[conn.second[1]] - points[conn.second[0]]);
        values.push_back(std::make_pair(dist, conn.first));
        idC = std::max(idC, conn.first);
    }

    //prep your priority queue
    std::priority_queue<std::pair<double,long>, std::vector<std::pair<double,long>>, greatDist > pQue(mycomp, values);

    int counterCell = 0;
    darray3E point;
    while (counterCell < fCells){

        //extract the longest segment
        auto pp = pQue.top();
        pQue.pop();

        idV++;
        idC++;
        //divide it by two and create a new point
        point = 0.5*(points[connectivity[pp.second][1]] + points[connectivity[pp.second][0]]);
        points.push_back(point);
        //update connectivity
        connectivity[idC]= {{ idV, connectivity[pp.second][1] }};
        connectivity[pp.second][1] = idV;

        //update the priority queue
        pQue.push(std::make_pair(0.5*pp.first, pp.second));
        pQue.push(std::make_pair(0.5*pp.first, idC));

        //adding a new cell to the count;
        ++counterCell;
    }
};


}
