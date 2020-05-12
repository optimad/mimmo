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
 \ *---------------------------------------------------------------------------*/

#ifndef __RECONSTRUCTSCALAR_HPP__
#define __RECONSTRUCTSCALAR_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \ingroup geohandlers
 * \brief class for setting overlapping criterium for two different fields
 */
enum class OverlapMethod{
    MAX = 1, /**< take max value (with sign) between overlapping fields*/
    MIN = 2, /**< take min value (with sign) between overlapping fields*/
    AVERAGE = 3, /**< take an average value between overlapping fields*/
    SUM = 4 /**< take sum of values between overlapping fields*/
};

/*!
 * \class ReconstructScalar
 * \ingroup geohandlers
 * \brief Reconstruct a scalar field from daughter meshes to mother mesh
 *
 * Class/BaseManipulation Object reconstructing a scalar field on a mimmo::MimmoObject mesh, from several
 * scalar fields defined on sub-patches of the target mesh, where for sub-patches is meant portion
 * of target mesh, preserving their vertex/cell-ids, as in the case of mimmo::Selection Blocks.
 * Field values can be defined on nodes or cells. No interfaces are supported up to now.
 * Reconstructed field on the whole geometry is provided as result as well as
 * the reconstructed fields on the input sub-patches separately.
 *
 *
 * Ports available in ReconstructScalar Class :
 *
 *    =========================================================
 *
     |                   Port Input    ||                                               |
     |----------------|--------------------|----------------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_SCALARFIELD  | addData             | (MC_SCALAR, MD_MPVECFLOAT_)          |
     | M_GEOM         | m_geometry         | (MC_SCALAR, MD_MIMMO_)                 |


     |             Port Output   ||                                          |
     |----------------|--------------------|-----------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_SCALARFIELD  | getResultField      | (MC_SCALAR, MD_MPVECFLOAT_)       |
     | M_VECSFIELDS    | getResultFields   | (MC_VECTOR, MD_MPVECFLOAT_)       |
     | M_GEOM         | getGeometry        | (MC_SCALAR, MD_MIMMO_)      |

 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation
 * - <B>ClassName</B>: name of the class as <tt>mimmo.ReconstructScalar</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing;
 *
 * Proper of the class:
 * - <B>DataLocation</B>: data location of fields 1-POINT, 2-CELL, 3-INTERFACE;
 * - <B>OverlapCriterium</B>: set how to treat fields in the overlapped region 1-MaxVal, 2-MinVal, 3-AverageVal, 4-Summing;
 *
 * Fields and Geometry have to be mandatorily passed through port.
 *
 */
class ReconstructScalar: public BaseManipulation {

private:
    MPVLocation m_loc;                  /**<Data reference Location */
    OverlapMethod m_overlapCriterium;   /**<Overlap Method */
    std::vector<dmpvector1D> m_subpatch;     /**<Vector of input fields on sub-patches. */
    std::vector<dmpvector1D> m_subresults;   /**<Vector of processed/overlapped fields on sub-patches. */
    dmpvector1D m_result;               /**<Output reconstructed field. */

public:
    ReconstructScalar(MPVLocation loc = MPVLocation::POINT);
    ReconstructScalar(const bitpit::Config::Section & rootXML);
    virtual ~ReconstructScalar();
    ReconstructScalar(const ReconstructScalar & other);

    //get-set methods
    OverlapMethod           getOverlapCriteriumENUM();
    int                     getOverlapCriterium();
    int                     getNData();
    dmpvector1D*             getResultField();
    std::vector<dmpvector1D*>     getResultFields();

    void        setOverlapCriteriumENUM( OverlapMethod);
    void        setOverlapCriterium( int);
    void        addData( dmpvector1D * );
    void        removeData(mimmo::MimmoSharedPointer<MimmoObject> );
    void        removeAllData();
    void        buildPorts();
    //cleaners
    void clear();

    //plotting

    void    plotData();
    void    plotSubData(int i);

    //execute
    void        execute();

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
protected:
    virtual void plotOptionalResults();
    void swap(ReconstructScalar &) noexcept;

private:
    void     overlapFields(long int ID, double & locField);
    livector1D idsGeoDataLocation(mimmo::MimmoSharedPointer<MimmoObject>);
};

/*!
 * \class ReconstructVector
 * \ingroup geohandlers
 * \brief Reconstruct a vector field from daughter mesh to mother mesh
 *
 * Class/BaseManipulation Object reconstructing a vector field on a mimmo::MimmoObject mesh, from several
 * vector fields defined on sub-patches of the target mesh, where for sub-patches is meant portion
 * of target mesh, preserving their vertex/cell-ids, as in the case of mimmo::Selection Blocks.
 * Field values can be defined on nodes or cells. No interfaces are supported up to now.
 * Reconstructed field on the whole geometry is provided as result as well as
 * the reconstructed fields on the input sub-patches separately.
 *
 * Ports available in ReconstructVector Class :
 *
 *    =========================================================
 *
     |                   Port Input   ||                                                  |
     |----------------|--------------------|------------------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_VECTORFIELD       | addData     | (MC_SCALAR, MD_MPVECARR3FLOAT_)           |
     | M_GEOM         | m_geometry         | (MC_SCALAR, MD_MIMMO_)                   |

     |             Port Output  ||                                                |
     |----------------|--------------------|----------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_VECTORFIELD       | getResultField     | (MC_SCALAR, MD_MPVECARR3FLOAT_)           |
     | M_VECVFIELDS   | getResultFields    | (MC_VECTOR, MD_MPVECARR3FLOAT_)       |
     | M_GEOM         | getGeometry        | (MC_SCALAR, MD_MIMMO_)           |

 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation
 * - <B>ClassName</B>: name of the class as <tt>mimmo.ReconstructVector</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing;
 *
 * Proper of the class:
 * - <B>DataLocation</B>: data location of fields 1-POINT, 2-CELL, 3-INTERFACE;
 * - <B>OverlapCriterium</B>: set how to treat fields in the overlapped region 1-MaxVal, 2-MinVal, 3-AverageVal, 4-Summing;
 *
 * Fields and Geometry have to be mandatorily passed through port.
 *
 */
class ReconstructVector: public mimmo::BaseManipulation {

private:
    MPVLocation m_loc;                  /**<Data reference Location */
    OverlapMethod m_overlapCriterium;   /**<Overlap Method */
    std::vector<dmpvecarr3E> m_subpatch;     /**<Vector of input fields on sub-patches. */
    std::vector<dmpvecarr3E>  m_subresults;  /**<Vector of processed/overlapped fields on sub-patches. */
    dmpvecarr3E m_result;               /**<Output reconstructed field. */

public:
    ReconstructVector(MPVLocation loc = MPVLocation::POINT);
    ReconstructVector(const bitpit::Config::Section & rootXML);
    virtual ~ReconstructVector();
    ReconstructVector(const ReconstructVector & other);

    //get-set methods
    OverlapMethod               getOverlapCriteriumENUM();
    int                         getOverlapCriterium();
    int                         getNData();
    dmpvecarr3E*                getResultField();
    std::vector<dmpvecarr3E*>    getResultFields();

    void        setOverlapCriteriumENUM( OverlapMethod);
    void        setOverlapCriterium(int );
    void        addData( dmpvecarr3E * );
    void        removeData(mimmo::MimmoSharedPointer<MimmoObject> );
    void        removeAllData();
    void        buildPorts();
    //cleaners
    void clear();


    //plotting
    void    plotData();
    void    plotSubData(int i);

    //execute
    void    execute();

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    virtual void plotOptionalResults();
    void swap(ReconstructVector &) noexcept;
private:
    void    overlapFields(long int ID, darray3E & locField);
    livector1D idsGeoDataLocation(mimmo::MimmoSharedPointer<MimmoObject>);
};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_, __RECONSTRUCTSCALAR_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_SCALAR, MD_MPVECFLOAT_, __RECONSTRUCTSCALAR_HPP__)
REGISTER_PORT(M_VECTORFIELD, MC_SCALAR, MD_MPVECARR3FLOAT_, __RECONSTRUCTSCALAR_HPP__)
REGISTER_PORT(M_VECSCALARFIELDS, MC_VECTOR, MD_MPVECFLOAT_, __RECONSTRUCTSCALAR_HPP__)
REGISTER_PORT(M_VECVECTORFIELDS, MC_VECTOR, MD_MPVECARR3FLOAT_, __RECONSTRUCTSCALAR_HPP__)


REGISTER(BaseManipulation, ReconstructScalar,"mimmo.ReconstructScalar")
REGISTER(BaseManipulation, ReconstructVector,"mimmo.ReconstructVector")

};

#endif /* __RECONSTRUCTFIELDS_HPP__ */
