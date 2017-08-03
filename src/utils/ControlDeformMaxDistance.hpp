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
#include "MimmoObject.hpp"

/*!
 * \ingroup Ports
 * \{
 */

//PORTS DEFINITION AS CONSTANTS

#ifndef M_VALUED 
#define M_VALUED "M_VALUED" /**< Port dedicated to communication  of a float/double scalar value */
#endif

#ifndef M_GEOM 
#define M_GEOM "M_GEOM" /**< Port dedicated to communication  of a MimmoObject pointer*/
#endif

#ifndef M_VIOLATION 
#define M_VIOLATION "M_VIOLATION" /**< Port dedicated to communication of a double value(violation) as std::pair< BaseManipulation*, double >. The BaseManipulation pointer retrack the Manipulation Block who trigger the violation.*/
#endif

#ifndef M_GDISPLS 
#define M_GDISPLS "M_GDISPLS" /**< Port dedicated to communication of a a list of displacements referred to geometry vertices*/
#endif

#ifndef M_SCALARFIELD 
#define M_SCALARFIELD "M_SCALARFIELD" /**< Port dedicated to communication of a generic scalar field of floats/doubles [std::vector< double >] */
#endif


/*!
 * \}
 */

/*!
 * \ingroup PortContainers
 * \{
 */

#ifndef MC_SCALAR
#define MC_SCALAR "MC_SCALAR" /**< Single value container identifier*/
#endif

#ifndef MC_VECTOR
#define MC_VECTOR "MC_VECTOR" /**< std::vector< . > container identifier*/
#endif

#ifndef MC_VECARR3 
#define MC_VECARR3 "MC_VECARR3" /**< std::vector< std::array< ., 3 > container identifier*/
#endif

#ifndef MC_PAIR 
#define MC_PAIR "MC_PAIR" /**< std::pair< . , .> container identifier. Needs tuple data*/
#endif

/*!
 * \}
 */

/*!
 * \ingroup PortData
 * \{
 */

#ifndef MD_FLOAT
#define MD_FLOAT "MD_FLOAT" /**< float/double data identifier*/
#endif

#ifndef MD_MIMMO_
#define MD_MIMMO_ "MD_MIMMO_" /**< MimmoObject pointer data identifier*/
#endif

#ifndef MD_PAIRMIMMO_OBJFLOAT_
#define MD_PAIRMIMMO_OBJFLOAT_ "MD_PAIRMIMMO_OBJFLOAT_" /**< tuple (BaseManipulation*, double) data identifier*/
#endif



/*!
 * \}
 */


namespace mimmo{

/*!
 * \class ControlDeformMaxDistance
 * \ingroup utils
 * \brief ControlDeformMaxDistance is a class that check a deformation field associated to a MimmoObject geometry,
 *  once a maximum limit distance of deformation is fixed, w.r.t. the undeformed state.
 *
 * ControlDeformMaxDistance is derived from BaseManipulation class.
 * It needs a maximum, isotropic limit distance d w.r.t. geometry undeformed state, which is used to evaluate the isolevel d
 * of the target geometry. 
 * Returns a double value V, namely the maximum signed distance from constraint iso-level amongst all field points, 
 * reporting how much the current deformation field violate the constraint itself.
 * if V >0 a violation occurs. if V=0, a contact occurs, otherwise if V<0 no violation occurs. 
 * No optional result are plot. 
 * Class absorbs/flushes its parameters from/to xml dictionaries
 *
 * \n
 * Ports available in ControlDeformMaxDistance Class :
 *
 *    =========================================================


     |Port Input  | | |
     |-|-|-|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GDISPLS| setDefField       | (MC_VECARR3, MD_FLOAT)       |
     | M_VALUED | setLimitDistance  | (MC_SCALAR, MD_FLOAT)         |
     | M_GEOM   | setGeometry       | (MC_SCALAR, MD_MIMMO_)        |

     |Port Output | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B> |
     | M_SCALARFIELD | getViolationField | (MC_VECTOR, MD_FLOAT)             |
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
    double                        m_maxDist;        /**<Limit Distance*/
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
    dmpvector1D                                getViolationField();

    void    setDefField(dmpvecarr3E field);
    void    setLimitDistance(double dist);

    void     execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    void plotOptionalResults();
    void swap(ControlDeformMaxDistance & x) noexcept;

};

//Ports
REGISTER_PORT(M_GDISPLS, MC_VECARR3, MD_FLOAT)
REGISTER_PORT(M_VALUED, MC_SCALAR, MD_FLOAT)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_)
REGISTER_PORT(M_SCALARFIELD, MC_VECTOR, MD_FLOAT)
REGISTER_PORT(M_VIOLATION, MC_PAIR, MD_PAIRMIMMO_OBJFLOAT_)


//ManipBlocks
REGISTER(BaseManipulation, ControlDeformMaxDistance, "mimmo.ControlDeformMaxDistance")
}

#endif /* __CONTROLDEFORMMAXDISTANCE_HPP__ */
