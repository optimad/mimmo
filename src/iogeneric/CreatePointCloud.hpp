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

#ifndef __CREATEPOINTCLOUD_HPP__
#define __CREATEPOINTCLOUD_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \class CreatePointCloud
 * \ingroup iogeneric
 * \brief CreatePointCloud manages cloud point data in raw format to create
    a MimmoObject Point Cloud container.

    The class takes point data from a raw list and reverse them in a MimmoObject Point
    Cloud container. List can be set as pure vector or as a MimmoPiercedVector. Please be aware
    one set method exclude the other.
    Optionally, If any scalar/vector field are assigned to raw points, they will be
    trasformed in data containers associated to the point cloud.
    Please note if the raw points are passed with a specific container list
    (vector or MimmoPiercedVector), only data passed with the same kind of
    container will be taken into account.
 *
    MPI version retains data only on the master rank (rank 0).
    For now, MPI Setting of raw data is supposed to happen in this two configuration:
        1) only master rank (0) retains the data
        2) every rank has the same identical data sets.
    So, each time the class will consider only data from master rank.
    No Point Cloud distribution among ranks is perfomed. Mesh and Field Data, even if are parallel,
    are filled only for the partition of master rank(0). The other ranks have empty partition
    for Mesh and Field Data.

 * \n
 * Ports available in CreatePointCloud Class :
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
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.CreatePointCloud</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.

 */
class CreatePointCloud: public BaseManipulation{
protected:
    //members
    dmpvector1D     m_scalarfield;  /**< MimmoPiercedVector scalar field */
    dmpvecarr3E     m_vectorfield;  /**< MimmoPiercedVector vector field */

    dmpvecarr3E       m_rawpoints;   /**< input Cloud points list, in raw format. */
    dmpvector1D       m_rawscalar;   /**< input scalar attached to Cloud points, in raw format. */
    dmpvecarr3E       m_rawvector;   /**< input vector attached to Cloud points, in raw format */

public:
    CreatePointCloud();
    virtual ~CreatePointCloud();
    CreatePointCloud(const bitpit::Config::Section & rootXML);
    CreatePointCloud(const CreatePointCloud & other);
    CreatePointCloud & operator=(CreatePointCloud other);

    dmpvector1D* getScalarField();
    dmpvecarr3E* getVectorField();

    void setRawPoints(dvecarr3E rawPoints);
    void setRawScalarField(dvector1D rawScalarField);
    void setRawVectorField(dvecarr3E rawVectorField);

    void setRawPoints(dmpvecarr3E *rawPoints);
    void setRawScalarField(dmpvector1D *rawScalarField);
    void setRawVectorField(dmpvecarr3E *rawVectorField);

    void clear();

    void execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    void swap(CreatePointCloud & x) noexcept;
    void buildPorts();
    void plotOptionalResults();

private:
    //disabling interface method
    void setGeometry(MimmoSharedPointer<MimmoObject> geometry){BaseManipulation::setGeometry(geometry);};
};

REGISTER_PORT(M_COORDS, MC_VECARR3, MD_FLOAT,__CREATEPOINTCLOUD_HPP__)
REGISTER_PORT(M_DISPLS, MC_VECARR3, MD_FLOAT,__CREATEPOINTCLOUD_HPP__)
REGISTER_PORT(M_DATAFIELD, MC_VECTOR, MD_FLOAT,__CREATEPOINTCLOUD_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_SCALAR, MD_MPVECFLOAT_,__CREATEPOINTCLOUD_HPP__)
REGISTER_PORT(M_VECTORFIELD, MC_SCALAR, MD_MPVECARR3FLOAT_,__CREATEPOINTCLOUD_HPP__)
REGISTER_PORT(M_VECTORFIELD2, MC_SCALAR, MD_MPVECARR3FLOAT_,__CREATEPOINTCLOUD_HPP__)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_ ,__CREATEPOINTCLOUD_HPP__)

REGISTER(BaseManipulation, CreatePointCloud, "mimmo.CreatePointCloud")
}

#endif /* __CREATEPOINTCLOUD_HPP__ */
