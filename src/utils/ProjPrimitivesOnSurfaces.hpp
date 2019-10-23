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

#ifndef __PROJPRIMITIVESONSURFACES_HPP__
#define __PROJPRIMITIVESONSURFACES_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \class ProjPrimitivesOnSurfaces
 * \ingroup utils
 * \brief Class for projecting 1D/2D primitives on a target 3D surface mesh
 *
 * Abstract class used as interface for objects of projecting elemental 1D or 2D primitives
 * such as segments, circumferences or circles, triangles etc... on a 3D surface mesh defined
 * by a MimmoObject.
 *
 * Ports available in ProjPrimitivesOnSurfaces Class :
 *
 *    =========================================================
 *
     |                   Port Input       ||                               |
     |----------------|----------------------|-------------------|
     |<B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM   | setGeometry                        | (MC_SCALAR, MD_MIMMO_)      |


     |             Port Output     ||                                      |
     |----------------|---------------------|--------------------|
     |<B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM    | getProjectedElement                | (MC_SCALAR, MD_MIMMO_)      |


 *    =========================================================
 *
 */
class ProjPrimitivesOnSurfaces: public BaseManipulation{

protected:
    int                 m_topo;            /**<Mark topology of your primitive element 1-one dimensional, 2- bi-dimensional*/
    int                 m_nC;           /**<Number of target elements of your 3D curve discrete projection */
    bool                m_buildSkdTree;    /**<If true build SkdTree of of the projected element mesh. */
    bool                m_buildKdTree;    /**<If true build KdTree of the projected element mesh. */
    std::unique_ptr<MimmoObject> m_patch;            /**< resulting projected elements stored as MimmoObject */


public:
    ProjPrimitivesOnSurfaces();
    virtual ~ProjPrimitivesOnSurfaces();

    ProjPrimitivesOnSurfaces(const ProjPrimitivesOnSurfaces & other);
    ProjPrimitivesOnSurfaces & operator=(const ProjPrimitivesOnSurfaces &other);

protected:
    virtual void buildPorts();

public:
    int                             getTopology();
    virtual int                     getProjElementTargetNCells();
    MimmoObject *                   getProjectedElement();

    void        setGeometry(MimmoObject * geo);

    void        setBuildSkdTree(bool build);
    void        setBuildKdTree(bool build);
    void        setProjElementTargetNCells(int nC);

    bool         isEmpty();

    void         clear();
    void         execute();

    void        plotOptionalResults();

protected:
    /*!
     * Projection function - Pure virtual method
     */
    virtual void projection() = 0;
    void swap(ProjPrimitivesOnSurfaces & x) noexcept;
};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_, __PROJPRIMITIVESONSURFACES_HPP__)
};

#endif /* __PROJPRIMITIVESONSURFACES_HPP__ */
