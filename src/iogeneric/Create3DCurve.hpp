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
#ifndef __CREATE3DCURVE_HPP__
#define __CREATE3DCURVE_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \class Create3DCurve
 * \ingroup iogeneric
 * \brief Create a 3DCurve from a point cloud.
 *
 The class takes point data from a raw list and reverse them in a MimmoObject 3D Curve
 container, using the consecutive ordering of the list to create connectivity of the curve.
 List can be set as pure vector or as a MimmoPiercedVector. Please be aware
 one set method exclude the other.
 Optionally, If any scalar/vector field are assigned to raw points, they will be
 trasformed in data containers associated to curve, with location on POINT.
 Please note if the raw points are passed with a specific container list
 (vector or MimmoPiercedVector), only data passed with the same kind of
 container will be taken into account.
*
 MPI version retains data only on the master rank (rank 0).
 For now, MPI Setting of raw data is supposed to happen in this two configuration:
     1) only master rank (0) retains the data
     2) every rank has the same identical data set.
 So, each time the class will consider only data from master rank.
 No 3DCurve distribution among ranks is perfomed. Mesh and Field Data will be retained only on the
 master rank(0).

* \n
* Ports available in Create3DCurve Class :
*
*    =========================================================

  |                 Port Input   ||                                     |
  |---------------|-------------------|-----------------------|
  | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
  | M_COORDS      | setRawPoints        | (MC_VECARR3, MD_FLOAT)      |
  | M_DISPLS      | setRawVectorField   | (MC_VECARR3, MD_FLOAT)      |
  | M_DATAFIELD   | setRawScalarField   | (MC_VECTOR, MD_FLOAT)       |
  | M_VECTORFIELD | setRawPoints        | (MC_SCALAR, MD_MPVECARR3FLOAT_) |
  | M_VECTORFIELD2| setRawVectorField   | (MC_SCALAR, MD_MPVECARR3FLOAT_) |
  | M_SCALARFIELD | setRawScalarField   | (MC_SCALAR,MD_MPVECFLOAT_)     |


  |              Port Output  ||                                        |
  |---------------|-------------------|-----------------------|
  | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
  | M_GEOM        | getGeometry       | (MC_SCALAR,MD_MIMMO_)     |
  | M_SCALARFIELD | getScalarField | (MC_SCALAR,MD_MPVECFLOAT_)     |
  | M_VECTORFIELD | getVectorField | (MC_SCALAR,MD_MPVECARR3FLOAT_) |

*    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B> : name of the class as <tt>mimmo.Create3DCurve</tt>;
 * - <B>Priority</B> : uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B> : boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B> : target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>ClosedLoop</B> : 0/1 to connect the point cloud list as open/closed 3D curve;
 * - <B>nCells</B> : refine curve tessellation up to a target number of cells (>= number of cloud points);
 *
 * Geometry has to be mandatorily passed through port.
 *
 */

class Create3DCurve: public BaseManipulation{

protected:
    //members
    dmpvector1D     m_scalarfield;  /**< MimmoPiercedVector scalar field */
    dmpvecarr3E     m_vectorfield;  /**< MimmoPiercedVector vector field */

    dmpvecarr3E       m_rawpoints;   /**< input cloud points list, in raw format. */
    dmpvector1D       m_rawscalar;   /**< input scalar attached to cloud points, in raw format. */
    dmpvecarr3E       m_rawvector;   /**< input vector attached to cloud points, in raw format */

    bool              m_closed;      /**<flag to create closed/open loops */
    int               m_nCells;      /**<flag to set a target number of cells in the curve */

public:
    Create3DCurve();
    Create3DCurve(const bitpit::Config::Section & rootXML);
    virtual ~Create3DCurve();

    Create3DCurve(const Create3DCurve & other);
    Create3DCurve & operator=(Create3DCurve other);

    dmpvector1D* getScalarField();
    dmpvecarr3E* getVectorField();

    void setRawPoints(dvecarr3E rawPoints);
    void setRawScalarField(dvector1D rawScalarField);
    void setRawVectorField(dvecarr3E rawVectorField);

    void setRawPoints(dmpvecarr3E *rawPoints);
    void setRawScalarField(dmpvector1D *rawScalarField);
    void setRawVectorField(dmpvecarr3E *rawVectorField);

    bool         isClosedLoop();
    void         setClosedLoop(bool flag);

    void         setNCells(int nCells);

    void clear();

    void execute();

    void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");


protected:
    void swap(Create3DCurve & x) noexcept;
    void buildPorts();
    void plotOptionalResults();

private:
    int fillPreliminaryStructure(dmpvecarr3E & points, std::unordered_map<long, std::array<long,2> > &connectivity);
    void refineObject(dmpvecarr3E & points, dmpvector1D & scalarf, dmpvecarr3E & vectorf, std::unordered_map<long, std::array<long,2> > &connectivity, int fCells);
    //disabling interface method
    void setGeometry(MimmoSharedPointer<MimmoObject> geometry){BaseManipulation::setGeometry(geometry);};

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

REGISTER_PORT(M_COORDS, MC_VECARR3, MD_FLOAT,__CREATE3DCURVE_HPP__)
REGISTER_PORT(M_DISPLS, MC_VECARR3, MD_FLOAT,__CREATE3DCURVE_HPP__)
REGISTER_PORT(M_DATAFIELD, MC_VECTOR, MD_FLOAT,__CREATE3DCURVE_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_SCALAR, MD_MPVECFLOAT_,__CREATE3DCURVE_HPP__)
REGISTER_PORT(M_VECTORFIELD, MC_SCALAR, MD_MPVECARR3FLOAT_,__CREATE3DCURVE_HPP__)
REGISTER_PORT(M_VECTORFIELD2, MC_SCALAR, MD_MPVECARR3FLOAT_,__CREATE3DCURVE_HPP__)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_ ,__CREATE3DCURVE_HPP__)


REGISTER(BaseManipulation, Create3DCurve, "mimmo.Create3DCurve")

};

#endif /* __CREATE3DCURVE_HPP__ */
