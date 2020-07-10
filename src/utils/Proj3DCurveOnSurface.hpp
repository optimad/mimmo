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
#ifndef __PROJ3DCURVEONSURFACE_HPP__
#define __PROJ3DCURVEONSURFACE_HPP__

#include "ProjPrimitivesOnSurfaces.hpp"

namespace mimmo{

/*!
 * \class Proj3DCurveOnSurface
 * \ingroup utils
 * \brief It is an executable block class capable of projecting a 3D curve defined by multiple
 *  elemental connected segments on a 3D surface mesh defined by a MimmoObject.
 *
 * Proj3DCurveOnSurface project a 3D Curve tessellation on a given 3D surface mesh and return it in a
 * in MimmoObject container. 3DCurve can be set as consecutive list of points, or given
 * as a MimmoObject itself.
 *
 * Ports available in Proj3DCurveOnSurface Class :
 *
 *    =========================================================
 *
     |                   Port Input       ||                               |
     |----------------|----------------------|-------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_COORDS | setPoints                          | (MC_VECARR3, MD_FLOAT)      |
     | M_POINT  | addPoint                           | (MC_ARRAY3, MD_FLOAT)       |
     | M_GEOM2  | setConnectedPoints                 | (MC_SCALAR, MD_MIMMO_)      |


     |             Port Output     ||                                      |
     |----------------|---------------------|--------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


     Inheriting ports from base class ProjPrimitivesOnSurfaces:

     |                   Port Input       ||                               |
     |----------------|----------------------|-------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM   | setGeometry                        | (MC_SCALAR, MD_MIMMO_)      |


     |             Port Output     ||                                      |
     |----------------|---------------------|--------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM    | getProjectedElement                | (MC_SCALAR, MD_MIMMO_)      |


 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B> : name of the class as <tt>mimmo.Proj3DCurveOnSurface</tt>;
 * - <B>Priority</B> : uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B> : boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B> : target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>ClosedLoop</B> : 0/1 to set current 3D input curve as open/closed;
 * - <B>nCells</B> : number of discrete cells of projected 3D curve;
 * - <B>SkdTree</B> : evaluate skdTree true 1/false 0;
 * - <B>KdTree</B> : evaluate kdTree true 1/false 0;
 *
 *
 * Geometry has to be mandatorily passed through port.
 *
 */

class Proj3DCurveOnSurface: public ProjPrimitivesOnSurfaces{

protected:
    dvecarr3E                                       m_cpoints;     /**<points given by the user to define curve*/
    bool                                            m_closed;      /**<flag to identify closed open loops */
    mimmo::MimmoSharedPointer<mimmo::MimmoObject>   m_cobj;        /**<object to deal with connected points */

public:
    Proj3DCurveOnSurface();
    Proj3DCurveOnSurface(const bitpit::Config::Section & rootXML);
    virtual ~Proj3DCurveOnSurface();

    Proj3DCurveOnSurface(const Proj3DCurveOnSurface & other);
    Proj3DCurveOnSurface & operator=(Proj3DCurveOnSurface other);

    void         clear();

    void         addPoint(darray3E point);
    void        setPoints(dvecarr3E points);
    void         setConnectedPoints(MimmoSharedPointer<MimmoObject> geo);

    bool        isClosedLoop();
    void         setClosedLoop(bool flag);

    void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

    void buildPorts();

protected:
    void projection();
    void swap(Proj3DCurveOnSurface & x) noexcept;

private:
    int fillPreliminaryStructure(dvecarr3E & points, std::unordered_map<long, std::array<long,2> > &connectivity);
    void refineObject(dvecarr3E & points, std::unordered_map<long, std::array<long,2> > &connectivity, int fCells);
    /*!
     * Auxiliary structure for sorting std::pair\<double,long\>
     */
    struct greatDist{
        /*!
         * Compare operator
         */
        bool operator()(const std::pair<double,long> & n1, const std::pair<double,long> & n2) const{
            return n1.first < n2.first;
        };
    };
};

REGISTER_PORT(M_COORDS, MC_VECARR3, MD_FLOAT, __PROJ3DCURVEONSURFACE_HPP__)
REGISTER_PORT(M_POINT, MC_ARRAY3, MD_FLOAT, __PROJ3DCURVEONSURFACE_HPP__)
REGISTER_PORT(M_GEOM2, MC_SCALAR, MD_MIMMO_, __PROJ3DCURVEONSURFACE_HPP__)

REGISTER(BaseManipulation, Proj3DCurveOnSurface, "mimmo.Proj3DCurveOnSurface")

};

#endif /* __PROJ3DCURVEONSURFACE_HPP__ */