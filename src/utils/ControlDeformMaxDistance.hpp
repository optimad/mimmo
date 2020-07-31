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
#ifndef __CONTROLDEFORMMAXDISTANCE_HPP__
#define __CONTROLDEFORMMAXDISTANCE_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \class ControlDeformMaxDistance
 * \ingroup utils
 * \brief ControlDeformMaxDistance is a class that check a deformation field
    associated to a MimmoObject surface geometry,once a maximum limit distance
    of deformation is fixed, w.r.t. the undeformed state.
 *
 * ControlDeformMaxDistance is derived from BaseManipulation class.
 * It needs a maximum, isotropic limit distance d w.r.t. the undeformed state of
   the target geometry (isolevel at distance d).
 * It returns a double value V, namely the maximum signed distance from
   the constraint iso-level amongst all field points, reporting how much the
   current deformation field violates the constraint itself.
 * if V>0 a violation occurs. if V=0, a contact occurs, otherwise if V<0 no violation occurs.
 * Deformation field must be linked to target geometry and defined on points/nodes.
 * Class absorbs/flushes its parameters from/to xml dictionaries
 *
 * \n
 * Ports available in ControlDeformMaxDistance Class :
 *
 *    =========================================================


     |Port Input  | | |
     |-|-|-|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GDISPLS| setDefField       | (MC_SCALAR, MD_MPVECARR3FLOAT_)       |
     | M_VALUED | setLimitDistance  | (MC_SCALAR, MD_FLOAT)         |
     | M_GEOM   | setGeometry       | (MC_SCALAR, MD_MIMMO_)        |

     |Port Output | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B> |
     | M_SCALARFIELD | getViolationField | (MC_SCALAR, MD_MPVECFLOAT_)             |
     | M_VALUED      | getViolation      | (MC_SCALAR, MD_FLOAT)             |

 *    =========================================================
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.ControlDeformMaxDistance</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>LimitDistance</B>: constraint surface distance from target geometry;
 *
 * Geometry and deformation field have to be mandatorily passed through port.
 *
 */
class ControlDeformMaxDistance: public BaseManipulation{
private:
    double                         m_maxDist;        /**<Limit Distance*/
    dmpvector1D                    m_violationField;    /**<Violation Distance Field */
    dmpvecarr3E                    m_defField;     /**<Deformation field*/

public:
    ControlDeformMaxDistance();
    ControlDeformMaxDistance(const bitpit::Config::Section & rootXML);
    ~ControlDeformMaxDistance();

    ControlDeformMaxDistance(const ControlDeformMaxDistance & other);
    ControlDeformMaxDistance & operator=(ControlDeformMaxDistance other);

    void    buildPorts();

    double                                     getViolation();
    dmpvector1D  *                             getViolationField();

    void    setDefField(dmpvecarr3E *field);
    void    setLimitDistance(double dist);
    void    setGeometry(MimmoSharedPointer<MimmoObject> geo);

    void     execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    void plotOptionalResults();
    void swap(ControlDeformMaxDistance & x) noexcept;
};

REGISTER_PORT(M_GDISPLS, MC_SCALAR, MD_MPVECARR3FLOAT_,__CONTROLDEFORMMAXDISTANCE_HPP__)
REGISTER_PORT(M_VALUED, MC_SCALAR, MD_FLOAT,__CONTROLDEFORMMAXDISTANCE_HPP__)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_,__CONTROLDEFORMMAXDISTANCE_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_SCALAR, MD_MPVECFLOAT_,__CONTROLDEFORMMAXDISTANCE_HPP__)


REGISTER(BaseManipulation, ControlDeformMaxDistance, "mimmo.ControlDeformMaxDistance")
}

#endif /* __CONTROLDEFORMMAXDISTANCE_HPP__ */
