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
#ifndef __PROJPATCHONSURFACE_HPP__
#define __PROJPATCHONSURFACE_HPP__

#include "ProjPrimitivesOnSurfaces.hpp"

namespace mimmo{

/*!
 * \class ProjPatchOnSurface
 * \ingroup utils
 * \brief Executable block class capable of projecting a 3D surface patch
 *  on another 3D surface, both defined as MimmoObject.
 *
 * ProjPatchOnSurface project a 3D surface tessellation on a given 3D surface mesh and return it in a
 * in MimmoObject container.
 *
 * Ports available in ProjPatchOnSurface Class :
 *
 *    =========================================================
 *
     |  Port Input       |||
     |----------------|----------------------|-------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM2  | setPatch   | (MC_SCALAR, MD_MIMMO_)      |

     | Port Output     |||
     |----------------|---------------------|--------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


 Inheriting ports from base class ProjPrimitivesOnSurfaces:
     | Port Input|||
     |----------------|----------------------|-------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM   | setGeometry                        | (MC_SCALAR, MD_MIMMO_)      |

     |Port Output  |||
     |----------------|---------------------|--------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM    | getProjectedElement                | (MC_SCALAR, MD_MIMMO_)      |


 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B> : name of the class as <tt>mimmo.ProjPatchOnSurface</tt>;
 * - <B>Priority</B> : uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B> : boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B> : target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>SkdTree</B> : evaluate skdTree true 1/false 0;
 * - <B>KdTree</B> : evaluate kdTree true 1/false 0;
 *
 *
 * Geometry has to be mandatorily passed through port.
 *
 */

class ProjPatchOnSurface: public ProjPrimitivesOnSurfaces{

protected:
    MimmoObject*             m_cobj;        /**<object to deal with patch */

public:
    ProjPatchOnSurface();
    ProjPatchOnSurface(const bitpit::Config::Section & rootXML);
    virtual ~ProjPatchOnSurface();

    ProjPatchOnSurface(const ProjPatchOnSurface & other);
    ProjPatchOnSurface & operator=(ProjPatchOnSurface other);

    void         clear();

    void         setPatch(MimmoObject * geo);

    void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

    void buildPorts();

protected:
    void projection();
    void swap(ProjPatchOnSurface & x) noexcept;

};

REGISTER_PORT(M_GEOM2, MC_SCALAR, MD_MIMMO_, __PROJPATCHONSURFACE_HPP__)

REGISTER(BaseManipulation, ProjPatchOnSurface, "mimmo.ProjPatchOnSurface")

};

#endif /* __PROJPATCHONSURFACE_HPP__ */
