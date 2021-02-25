/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

#ifndef __CREATESEEDSONSURFACE_HPP__
#define __CREATESEEDSONSURFACE_HPP__

#include "BaseManipulation.hpp"
#include "OBBox.hpp"



namespace mimmo{

/*!
 * \ingroup utils
 * \brief Enum class for engine choice to set up initial points on a 3D surface.
 */
enum class CSeedSurf{
    RANDOM = 0 /**< Engine type, sows randomly points on surface */,
    LEVELSET = 1 /**< Engine type, sows points around a seed on surface,using geodesic distance between points */,
    CARTESIANGRID=2 /**< Engine type, sows points projecting a 3D Cartesian grid on surface */

};

/*!
 * \class CreateSeedsOnSurface
 * \ingroup utils
 * \brief Distribute a raw list of points on a target 3D surface
 *
 * Class/BaseManipulation Object to create a raw set of 3D points lying on a 3D surface.

 * MPI version will treat this raw list as shared among the ranks (every rank will retain
   exactly the same list)

 * Three type of engines to compute point positions are available: \n
 * - CSeedSurf::RANDOM : sows points randomly on surface, trying to displace them
 * at maximum euclidean distance; \n
 * - CSeedSurf::CARTESIANGRID    evaluate points by projection of a volumetric
 * cartesian grid of surface and decimating the list up the desired value of points,
 * trying to displace them at maximum euclidean distance possible. \n
 * - CSeedSurf::LEVELSET : starting from an initial seed, sows points around it, trying
 * to displace them at maximum geodesic distance possible on the surface; BEWARE this option
   is not available in MPI version with np processors > 1. \n
 *
 * Default engine is CARTESIANGRID.

 * A Sensitivity field, defined on the target 3D surface vertices can be linked,
   to drive the seeding procedure, i.e. to displace points in the most sensible
   location according to the field.
 * \n
 * Ports available in CreateSeedsOnSurface Class :
 *
 *    =========================================================
 *
     | Port Input | | |
     |-|-|-|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_POINT        | setSeed                               | (MC_ARRAY3, MD_FLOAT)       |
     | M_GEOM         | setGeometry                           | (MC_SCALAR, MD_MIMMO_)      |
     | M_FILTER       | setSensitivityField                   | (MC_SCALAR, MD_MPVECFLOAT_)       |

     |Port Output  | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>              |
     | M_COORDS       | getPoints          | (MC_VECARR3, MD_FLOAT)       |
     | M_VALUEI       | getRandomSignature | (MC_SCALAR, MD_INT)          |

 *    =========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.CreateSeedsOnSurface</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: plot optional result during object execution;
 * - <B>OutputPlot</B>: set path to store optional result.
 *
 * Proper of the class:
 * - <B>NPoints</B>: total points to distribute;
 * - <B>Engine</B>: type of distribution engine 0:Random,2:CartesianGrid,1:Levelset(NOT AVAILABLE IN MPI np>1 VERSION);
 * - <B>Seed</B>: initial seed point coordinates (space separated);
 * - <B>MassCenterAsSeed</B>: boolean (0/1), if true use geometry mass center as seed;
 * - <B>RandomFixed</B>: boolean(0/1), if active it fixes distribution pattern when 0:RANDOM engine is selected,
                         otherwise leave it changing as it likes;
 * - <B>RandomSignature</B>: specify unsigned integer seed when RandomFixed is active.
                            This will reproduce the same random distribution each run;
 *
 *
 */
class CreateSeedsOnSurface: public mimmo::BaseManipulation {

private:

    dvecarr3E   m_points;         /**< resulting points of computation */
    int         m_nPoints;        /**< total number of desired points */
    double      m_minDist;        /**< minimum distance of tolerance */
    darray3E    m_seed;           /**< inital seed point */
    CSeedSurf   m_engine;         /**< choose kernel type for points positioning computation */
    bool        m_seedbaricenter; /**< bool activate mass center as starting seed */
    bool        m_randomFixed;    /**< true if User want to reproduce always the same random seed distribution*/
    uint32_t    m_randomSignature;/**< signature for freezing random engine result*/
    dmpvector1D m_sensitivity;    /**< sensitivity map, defined on target geometry to drive placement of seeds*/

    //utility members
    std::unique_ptr<mimmo::OBBox> m_bbox;     /**<pointer to an oriented Bounding box */
    livector1D m_deads;                     /**< inactive ids */
    dmpvector1D m_sensitivity_triangulated; /**! internal member, adjust point-field if re-triangulation occurs*/
    dvector1D m_final_sensitivity;          /**! raw list of interpolated sensitivities on shared points*/

public:
    CreateSeedsOnSurface();
    CreateSeedsOnSurface(const bitpit::Config::Section & rootXML);
    virtual ~CreateSeedsOnSurface();
    CreateSeedsOnSurface(const CreateSeedsOnSurface & other);
    CreateSeedsOnSurface & operator=(CreateSeedsOnSurface other);

    void    buildPorts();

    //get methods
    int          getNPoints();
    dvecarr3E    getPoints();
    CSeedSurf    getEngineENUM();
    int          getEngine();
    darray3E     getSeed();
    bool         isSeedMassCenter();
    double       getMinDistance();
    bool         isRandomFixed();
    long         getRandomSignature();

    //set methods
    void         setNPoints( int);
    void         setEngineENUM( CSeedSurf);
    void         setEngine(int);
    void         setSeed(darray3E);
    void         setMassCenterAsSeed(bool );
    void         setGeometry(MimmoSharedPointer<MimmoObject>);
    void         setRandomFixed(bool fix);
    void         setRandomSignature(uint32_t signature);
    void         setSensitivityMap(dmpvector1D *field);

    void         clear();

    void         solve(bool debug = false);
    void         plotCloud(std::string dir, std::string file, int counter, bool binary);

    void         execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    virtual void plotOptionalResults();
    void    swap(CreateSeedsOnSurface & x) noexcept;
    void    checkField();
    void    normalizeField();
    bool    checkTriangulation();
    MimmoSharedPointer<MimmoObject> createTriangulation();

private:
    void    solveLSet( bool debug = false);
    void    solveGrid(bool debug = false);
    void    solveRandom(bool debug = false);

    void decimatePoints(dvecarr3E & rawpoint, dvector1D & rawsensi);

    void   solveEikonal(double g, double s, bitpit::PatchKernel &tri,std::unordered_map<long,long> & invConn, dmpvector1D & field);
    double updateEikonal(double g, double s, long tVert,long tCell,bitpit::PatchKernel &tri, std::unordered_map<long int, short int> &flag, dmpvector1D & field);

    std::unordered_map<long,long>    getInverseConn(bitpit::PatchKernel &);
    bool                             isDeadFront(const long int label);
    std::set<long>    findVertexVertexOneRing(bitpit::PatchKernel &, const long &, const long & );
};

REGISTER_PORT(M_POINT, MC_ARRAY3, MD_FLOAT,__CREATESEEDSONSURFACE_HPP__)
REGISTER_PORT(M_VALUELI, MC_SCALAR, MD_LONG,__CREATESEEDSONSURFACE_HPP__)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_,__CREATESEEDSONSURFACE_HPP__)
REGISTER_PORT(M_FILTER, MC_SCALAR, MD_MPVECFLOAT_,__CREATESEEDSONSURFACE_HPP__)
REGISTER_PORT(M_COORDS, MC_VECARR3, MD_FLOAT,__CREATESEEDSONSURFACE_HPP__)


REGISTER(BaseManipulation, CreateSeedsOnSurface, "mimmo.CreateSeedsOnSurface")
};

#endif /* __CREATESEEDSONSURFACE_HPP__ */
