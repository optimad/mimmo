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

#include "MimmoCGUtils.hpp"
#include <CG.hpp>
#include <bitpit_operators.hpp>

namespace mimmo{

namespace mimmoCGUtils{

/*!
 * Custom method to compute distance from a given polygon of vertices V > 3.
 * The polygon is subdivided in triangles using the polygon barycenter, then distances
 * from each sub-triangle are computed. Return the minimun distance between those calculated.
 * \param[in] p target point coordinates
 * \param[in] vertCoords vertices of the polygon
 * \return distance from the polygon
 */
double distancePointPolygon(const darray3E & p,const dvecarr3E & vertCoords){

    //calculate barycenter
    darray3E barycenter = {{0.0,0.0,0.0}};
    for(const auto & vv : vertCoords){
        barycenter += vv;
    }
    barycenter /= double(vertCoords.size());

    //split it in subtriangles - fixed barycenter, run 2-by-2 along vertCoords to get the other two vertices
    double distance(std::numeric_limits<double>::max());
    std::size_t sizeV = vertCoords.size();
    for(std::size_t i=0; i<sizeV-1; ++i){
        distance = std::min(distance, bitpit::CGElem::distancePointTriangle(p, barycenter, vertCoords[i], vertCoords[i+1]));
    }
    return distance;
}


/*!
 * \return true if the target point is inside the segment V0-V1
 * \param[in] p target point coordinates
 * \param[in] V0 segment first vertex
 * \param[in] V1 segment second vertex
 * \param[in] tol geometric tolerance
 */
bool isPointInsideSegment(const darray3E & p, const darray3E & V0,const darray3E & V1, double tol ){
    darray3E Pv      = p - V0;
    darray3E segment = V1 - V0;
    double normSegment = norm2(segment);
    if(normSegment > std::numeric_limits<double>::min()) segment /= normSegment;
    double V = dotProduct(Pv, segment);
    if(norm2(Pv - V*segment)> tol)  return false;
    return ( V >=0.0 && V <= normSegment);
}

/*!
 * \return true if the target point is inside the triangle V0-V1-V2
 * \param[in] p target point coordinates
 * \param[in] V0 triangle first vertex
 * \param[in] V1 triangle second vertex
 * \param[in] V2 triangle third vertex
 * \param[in] tol geometric tolerance
 */
bool isPointInsideTriangle(const darray3E & p, const darray3E & V0,const darray3E & V1,const darray3E & V2, double tol){

    darray3E normal  = crossProduct(V1 - V0,V2 - V0);
    double normNormal = norm2(normal);
    if(normNormal > std::numeric_limits<double>::min()) normal /= normNormal;
    else return false;

    if(std::abs(dotProduct((p-V0), normal)) > tol){ return false;}

    double a1 = bitpit::CGElem::areaTriangle(p,V0,V1)/(0.5*normNormal);
    double a2 = bitpit::CGElem::areaTriangle(p,V1,V2)/(0.5*normNormal);
    double a3 = bitpit::CGElem::areaTriangle(p,V2,V0)/(0.5*normNormal);
    return ( (a1+a2+a3 -1.0) <= tol );
}

/*!
 * \return true if the target point is inside a given polygon of vertices V > 3.
 * The polygon is subdivided in triangles using the polygon barycenter, then check if inside
 * in each sub-triangle is computed.
 * \param[in] p target point coordinates
 * \param[in] vertCoords vertices of the polygon
 * \param[in] tol geometric tolerance
 */
bool isPointInsidePolygon(const darray3E & p,const dvecarr3E & vertCoords, double tol){

    //calculate barycenter
    darray3E barycenter = {{0.0,0.0,0.0}};
    for(const auto & vv : vertCoords){
        barycenter += vv;
    }
    barycenter /= double(vertCoords.size());

    //split it in subtriangles - fixed barycenter, run 2-by-2 along vertCoords to get the other two vertices
    bool check = false;
    std::size_t sizeV = vertCoords.size();
    std::size_t i = 0;
    while(i<sizeV-1 && !check ){
        check = mimmoCGUtils::isPointInsideTriangle(p, barycenter, vertCoords[i], vertCoords[i+1], tol);
        ++i;
    }
    return check;
}


}

}; // end namespace mimmo
