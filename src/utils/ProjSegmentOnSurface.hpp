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

#ifndef __PROJSEGMENTONSURFACE_HPP__
#define __PROJSEGMENTONSURFACE_HPP__

#include "ProjPrimitivesOnSurfaces.hpp"

namespace mimmo{

/*!
 *  \class ProjSegmentOnSurface
 * \ingroup utils
 *  \brief Executable block class capable of projecting an elemental segment on
           a 3D surface mesh defined by a MimmoObject.
 *
 *  ProjSegmentOnSurface project a given segment on a given 3D surface mesh and
    return a discrete 3D curve mesh in MimmoObject container
 *
 * Ports available in ProjSegmentOnSurface Class :
 *
 *    =========================================================
 *
     |                   Port Input       ||                               |
     |----------------|----------------------|-------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


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
 *  - <B>ClassName</B> : name of the class as <tt>mimmo.ProjSegmentOnSurface</tt>;
 *  - <B>Priority</B>  : uint marking priority in multi-chain execution;
 *  - <B>PlotInExecution</B> : boolean 0/1 print optional results of the class;
 *  - <B>OutputPlot</B> : target directory for optional results writing.
 *
 * Proper of the class
 *  - <B>Segment</B> : pass extremal points of the segmentSegmentPoints (6 coordinates, space separated);
 *  - <B>nCells</B> : number of discrete cells of projected 3D curve;
 *  - <B>SkdTree</B> : evaluate skdTree true 1/false 0;
 *  - <B>KdTree</B> : evaluate kdTree true 1/false 0;

 *
 * Geometry has to be mandatorily passed through port.
 *
 */
class ProjSegmentOnSurface: public ProjPrimitivesOnSurfaces{

private:
    darray3E                 m_pointA;        /**<origin of your segment*/
    darray3E                 m_pointB;        /**<end of your segment*/

public:
    ProjSegmentOnSurface();
    ProjSegmentOnSurface(const bitpit::Config::Section & rootXML);
    virtual ~ProjSegmentOnSurface();

    ProjSegmentOnSurface(const ProjSegmentOnSurface & other);
    ProjSegmentOnSurface & operator=(ProjSegmentOnSurface other);

    void         clear();

    void         setSegment(darray3E pointA, darray3E pointB);
    void        setSegment(darray3E origin, darray3E dir, double length);

    void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void projection();
    void swap(ProjSegmentOnSurface & x) noexcept;
};

REGISTER(BaseManipulation, ProjSegmentOnSurface, "mimmo.ProjSegmentOnSurface")

};

#endif /* __PROJSEGMENTONSURFACE_HPP__ */
