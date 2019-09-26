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
#ifndef __REFINEGEOMETRY_HPP__
#define __REFINEGEOMETRY_HPP__

#include "MimmoObject.hpp"
#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \enum RefineMode
 * \ingroup geohandlers
 * \brief Methods available for refining surface geometry.
 */
enum class RefineMode{
    GLOBAL = 0, /**< Refine the whole geometry*/
    	    SELECTION = 1, /**< Refine by giving one ore more external input sub-selection of the target geometry*/
			LINE = 2	/**< Refine by following the projection of one or more external input line on the target geometry*/
};

/*!
 *    \class RefineGeometry
 * \ingroup geohandlers
 *    \brief RefineGeometry is an executable block class capable of
 *         refine a surface geometry
 *
 *    RefineGeometry is the object to refine a surface MimmoObject.
 *    The refining can be performed by using an external geometry, for instance a sub-selection of the surface patch
 *    or a line MimmoObject.
 * 
 * Ports available in RefineGeometry Class :
 *
 *    =========================================================

     |Port Input | | |
     |-|-|-|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM   | setGeometry                     | (MC_SCALAR, MD_MIMMO_)      |
     | M_GEOM2   | addExternalGeometry                     | (MC_SCALAR, MD_MIMMO_)      |


     |Port Output | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>              |
     | M_GEOM    | getGeometry                        | (MC_SCALAR, MD_MIMMO_)         |

 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as "mimmo.RefineGeometry"
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>RefineMode</B>: Refine method.
 
 * Geometries have to be mandatorily passed through port.
 *
 */
class RefineGeometry: public BaseManipulation{

private:
	std::vector<MimmoObject*>   m_extgeo;   /**< Pointers to external geometries*/
	RefineMode					m_mode;		/**< Refine mode. Default 2 - Line Path Refinement*/

//    std::unique_ptr<MimmoObject> m_erasepatch;    /**< erased patch of the original geometry */

public:
    RefineGeometry();
    RefineGeometry(const bitpit::Config::Section & rootXML);
    ~RefineGeometry();

    RefineGeometry(const RefineGeometry & other);
    RefineGeometry & operator=(RefineGeometry other);

    void buildPorts();

    RefineMode	getRefineMode();

    void	addExternalGeometry(MimmoObject * geo);
    void	setRefineMode(RefineMode mode);
    void	setRefineMode(int mode);
    void	clear();
    void	execute();
    
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

    void     plotOptionalResults();
protected:
    void swap(RefineGeometry & x) noexcept;
};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_, __REFINEGEOMETRY_HPP__)
REGISTER_PORT(M_GEOM2, MC_SCALAR, MD_MIMMO_, __REFINEGEOMETRY_HPP__)

REGISTER(BaseManipulation, RefineGeometry, "mimmo.RefineGeometry")

};

#endif /* __REFINEGEOMETRY_HPP__ */
