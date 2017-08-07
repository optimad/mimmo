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
#include <unordered_map>
#include <utility>

namespace mimmo{

/*!
 * @enum OverlapMethod
 * \ingroup geohandlers
 * @brief class for setting overlapping criterium for two different scalar fields:
 *             1) MAX = get max(with sign) between available field values
 *             2) MIN = get min(with sign) between available field values
 *             3) AVERAGE = get simple average between available values
 *             4) SUM = take sum of both values between overlapped fields
 */
enum class OverlapMethod{
    MAX = 1 /**< take max values between overlapped fields*/,
            MIN = 2 /**< take min values between overlapped fields*/,
            AVERAGE = 3 /**< take averaged values between overlapped fields*/,
            SUM = 4 /**< take sum of both values between overlapped fields*/
};

/*!
 * \class ReconstructScalar
 * \ingroup geohandlers
 * \brief Reconstruct a scalar field from daughter meshes to mother mesh
 * 
 * Class/BaseManipulation Object reconstructing a scalar field on a mimmo::MimmoObject mesh, from several
 * scalar fields defined on sub-patches of the target mesh.
 * Field values are defined on nodes.
 * Reconstructed field on the whole geometry is provided as result as well as
 * the reconstructed fields on the input sub-patches separately.
 *
 * 
 * Ports available in ReconstructScalar Class :
 * 
 *    =========================================================
 *
     |                   Port Input    |||                                               |
     |-------|----------------|--------------------|----------------------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | 18    | M_SCALARFIELD  | addData             | (MC_MPVECTOR, MD_FLOAT)          |
    | 99    | M_GEOM         | m_geometry         | (MC_SCALAR, MD_MIMMO_)                 |


     |             Port Output   |||                                          |
     |-------|----------------|--------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | 18    | M_SCALARFIELD  | getResultField      | (MC_MPVECTOR, MD_FLOAT)       |
    | 60    | M_VECSFIELDS    | getResultFields   | (MC_VECTOR, MD_MPVECFLOAT)       |
    | 99    | M_GEOM         | getGeometry        | (MC_SCALAR, MD_MIMMO_)      |

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
 * - <B>OverlapCriterium</B>: set how to treat fields in the overlapped region 1-MaxVal, 2-MinVal, 3-AverageVal, 4-Summing;
 *
 * Fields and Geometry have to be mandatorily passed through port.
 *
 */
class ReconstructScalar: public mimmo::BaseManipulation {

private:

    OverlapMethod m_overlapCriterium;                                      /**<Overlap Method */
    vector<dmpvector1D> m_subpatch;                                         /**<Vector of pointers to fields of sub-patches. */
    vector<dmpvector1D> m_subresults;                                         /**<Vector of overlapped fields of sub-patches. */
    dmpvector1D m_result;                                                   /**<Output reconstructed field. */

public:
    ReconstructScalar();
    ReconstructScalar(const bitpit::Config::Section & rootXML);
    virtual ~ReconstructScalar();
    ReconstructScalar(const ReconstructScalar & other);

    //get-set methods
    OverlapMethod           getOverlapCriteriumENUM();
    int                     getOverlapCriterium();
    int                     getNData();
    dmpvector1D             getResultField();
    vector<dmpvector1D>     getResultFields();

    void        setOverlapCriteriumENUM( OverlapMethod);
    void        setOverlapCriterium( int);
    void        addData( dmpvector1D );
    void        removeData(mimmo::MimmoObject* );
    void        removeAllData();
    void        buildPorts();
    //cleaners
    void clear();

    //plotting

    void    plotData(std::string dir, std::string name, bool flag, int counter);
    void    plotSubData(std::string dir, std::string name, int i, bool flag, int counter);

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
};

/*!
 * \class ReconstructVector
 * \ingroup geohandlers
 * \brief Reconstruct a vector field from daughter mesh to mother mesh
 * 
 * Class/BaseManipulation Object reconstructing a vector field on a mimmo::MimmoObject mesh, from several
 * vector fields defined on sub-patches of the target mesh. Field values are defined on nodes.
 * Reconstructed field on the whole geometry is provided as result as well as
 * the reconstructed fields on the input sub-patches separately.
 * 
 * Ports available in ReconstructVector Class :
 * 
 *    =========================================================
 *
     |                   Port Input   |||                                                  |
     |-------|----------------|--------------------|------------------------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | 19    | M_VECTORFIELD       | addData     | (MC_MPVECARR3, MD_FLOAT)           |
    | 99    | M_GEOM         | m_geometry         | (MC_SCALAR, MD_MIMMO_)                   |

     |             Port Output  |||                                                |
     |-------|----------------|--------------------|----------------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | 19    | M_VECTORFIELD       | getResultField     | (MC_MPVECARR3, MD_FLOAT)           |
    | 61    | M_VECVFIELDS   | getResultFields    | (MC_VECTOR, MD_MPVECARR3)       |
    | 99    | M_GEOM         | getGeometry        | (MC_SCALAR, MD_MIMMO_)           |

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
 * - <B>OverlapCriterium</B>: set how to treat fields in the overlapped region 1-MaxVal, 2-MinVal, 3-AverageVal, 4-Summing;
 *
 * Fields and Geometry have to be mandatorily passed through port.
 *
 */
class ReconstructVector: public mimmo::BaseManipulation {

private:

    OverlapMethod m_overlapCriterium;   /**<Overlap Method */
    vector<dmpvecarr3E> m_subpatch;    /**<Vector of pointers to fields of sub-patches. */
    vector<dmpvecarr3E>  m_subresults;   /**<Vector of overlapped fields of sub-patches. */
    dmpvecarr3E m_result;               /**<Output reconstructed field. */

public:
    ReconstructVector();
    ReconstructVector(const bitpit::Config::Section & rootXML);
    virtual ~ReconstructVector();
    ReconstructVector(const ReconstructVector & other);

    //get-set methods
    OverlapMethod               getOverlapCriteriumENUM();
    int                         getOverlapCriterium();
    int                         getNData();
    dmpvecarr3E                 getResultField();
    std::vector<dmpvecarr3E>    getResultFields();

    void        setOverlapCriteriumENUM( OverlapMethod);
    void        setOverlapCriterium(int );
    void        addData( dmpvecarr3E );
    void        removeData(mimmo::MimmoObject* );
    void        removeAllData();
    void        buildPorts();
    //cleaners
    void clear();


    //plotting
    void    plotData(std::string dir, std::string name, bool flag, int counter);
    void    plotSubData(std::string dir, std::string name, int i, bool flag, int counter);

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
};

REGISTER(BaseManipulation, ReconstructScalar,"mimmo.ReconstructScalar")
REGISTER(BaseManipulation, ReconstructVector,"mimmo.ReconstructVector")

};

#endif /* __RECONSTRUCTSCALAR_HPP__ */
