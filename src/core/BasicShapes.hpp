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

#ifndef __MIMMOBASICSHAPES_HH
#define __MIMMOBASICSHAPES_HH

#include "MimmoObject.hpp"
#include <mimmo_binary_stream.hpp>
#include <unordered_set>

namespace mimmo{

/*!
 * \ingroup core
 * \brief Identifies the type of elemental shape supported by BasicShape class
 */
enum class ShapeType{
    CUBE    =0, /**< Cubic or generic voxel-shaped objects.*/
    CYLINDER=1, /**< Cylindrical objects.                  */
    SPHERE  =3, /**< Spherical objects.                    */
    WEDGE   =4  /**< Base triangular prism objects         */
};

/*!
 * \ingroup core
 * \brief Specify type of conditions to distribute NURBS nodes in a given coordinate of the shape
 */
enum class CoordType{

    UNCLAMPED =0,  /**<free clamping conditions*/
    CLAMPED   =1,  /**<force nurbs to pass on the extremal point of the interval (clamping)*/
    PERIODIC  =2,  /**<provide periodic condition on the interval extrema*/
    SYMMETRIC =3  /**<provide simmetry condition on the interval extrema*/
};

}
/*!
 * \ingroup binaryStream
 * \{
 */
//Binary stream occurring to Classes memebers for port exchange
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buf, mimmo::CoordType& element);
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buf, const mimmo::CoordType& element);

mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buf, mimmo::ShapeType& element);
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buf, const mimmo::ShapeType& element);
/*!
 *\}
 */

 namespace mimmo{

/*!
 *\class BasicShape
 *\ingroup core
 *\brief Abstract Interface class for Elementary Shape Representation
 *
 * Interface class for Volumetric Core Element, suitable for interaction with
   Data Structure stored in a MimmoObject class.
 * Object orientation in 3D space can be externally manipulated with dedicated
   transformation blocks. Class internally implements transformation to/from
   local sdr to/from world sdr, that can be used in derived objects from it.
 *
 * Class works with three reference systems:
 * + Global Absolute SDR: is the external World reference system
 * + Local Relative SDR: is the local reference system, not affected by Rigid Transformations as RotoTranslations or Scalings
 * + basic SDR: local system remapping to unitary cube, not accounting of the shape type.
 *
 */
class BasicShape {

protected:
    ShapeType	m_shape;        /**< shape identifier, see BasicShape::ShapeType enum */
    darray3E	m_origin;       /**< origin of your shape. May change according to the shape type*/
    darray3E	m_span;         /**< coordinate span of your current shape, in its local reference system*/
    darray3E    m_infLimits;    /**< inferior limits of coordinate of your current shape */
    dmatrix33E	m_sdr;          /**< axis position of the local reference system w.r.t absolute one*/
    dmatrix33E	m_sdr_inverse;  /**< inverse rotating sdr */
    darray3E 	m_scaling;      /**< scaling vector of dimensional coordinates */

    std::array<CoordType,3>     m_typeCoord;    /**< identifiers for coordinate type definition.DEFAULT is clamped*/

    dvecarr3E m_bbox;   /**<points of bounding box of the shape */  //temporary not copied build with get bounding box !!!

public:

    //Costructors, Destructor
    BasicShape();
    virtual ~BasicShape();

    //set-get methods
    void        setOrigin(darray3E);
    void        setSpan(double, double, double);
    void        setSpan(darray3E);
    void        setInfLimits(double val, int dir);
    void        setInfLimits(darray3E val);
    void        setRefSystem(darray3E, darray3E, darray3E);
    void        setRefSystem(int, darray3E);
    void        setRefSystem(dmatrix33E axes);
    void        setCoordinateType(CoordType, int dir);

    darray3E    getOrigin();
    darray3E    getSpan();
    darray3E    getInfLimits();
    dmatrix33E  getRefSystem();
    CoordType   getCoordinateType(int dir);
    ShapeType   getShapeType();
    ShapeType   getShapeType() const;

    darray3E    getScaling();
    darray3E    getLocalSpan();


    //functionalities
    livector1D  includeGeometry(mimmo::MimmoObject * );
    livector1D  excludeGeometry(mimmo::MimmoObject * );

    livector1D  includeGeometry(bitpit::PatchKernel * );
    livector1D  excludeGeometry(bitpit::PatchKernel * );

    livector1D  includeCloudPoints(const dvecarr3E &);
    livector1D  excludeCloudPoints(const dvecarr3E &);

    livector1D  includeCloudPoints(bitpit::PatchKernel * );
    livector1D  excludeCloudPoints(bitpit::PatchKernel * );

    livector1D  includeCloudPoints(mimmo::MimmoObject * );
    livector1D  excludeCloudPoints(mimmo::MimmoObject * );

    bool        isSimplexIncluded(const dvecarr3E &);
    bool        isSimplexIncluded(bitpit::PatchKernel * , const long int &indexT);
    bool        isPointIncluded(const darray3E &);
    bool        isPointIncluded(bitpit::PatchKernel * , const long int &indexV);

    /*!
     * Pure virtual method to get if the current shape an a given Axis Aligned Bounding Box intersects
     * \param[in] bMin min point of AABB
     * \param[in] bMax max point of AABB
     */
    virtual bool    intersectShapeAABBox(const darray3E &bMin, const darray3E &bMax)=0;


    /*!
     * Pure virtual method to convert from Local to World Coordinates
     * \param[in] point 3D point in local shape coordinates
     * \return point in world global coordinate
     */
    virtual	darray3E    toWorldCoord(const darray3E & point)=0;

    /*!
     * Pure virtual method to convert from World to Local Coordinates
     * \param[in] point 3D point in world global coordinates
     * \return point in local world coordinate
     *
     */
    virtual	darray3E    toLocalCoord(const darray3E & point)=0;

    /*!
     * Pure virtual method to get local Coordinate inferior limits of primitive shape
     * \return local origin
     */
    virtual darray3E    getLocalOrigin()=0;

protected:
    uint32_t    intersectShapePlane(int level, const darray3E & target);
    darray3E    checkNearestPointToAABBox(const darray3E &point, const darray3E &bMin, const darray3E &bMax);
    void swap(BasicShape & ) noexcept;
    static dmatrix33E transpose(const dmatrix33E & mat);
    static dmatrix33E inverse(const dmatrix33E & mat);
    static darray3E   matmul(const darray3E & vec, const dmatrix33E & mat);
    static darray3E   matmul(const dmatrix33E & mat,const darray3E & vec);

private:
    /*!
     * Pure Virtual method to transform a 3D point from basic elemental shape coordinate in local ones
     * \param[in] point mapped in basic elemental shape coordinates
     * \return point in local reference system coordinates
     */
    virtual darray3E    basicToLocal(const darray3E &point)=0;

    /*!
     * Pure Virtual method to transform a 3D point from local shape coordinate in basic elemental ones
     * \param[in] point mapped in local shape coordinates
     * \return point in basic elemental shape reference system coordinate
     */
    virtual darray3E    localToBasic(const darray3E &point)=0;

    /*!
     * Pure virtual method to check if your new span values fit your current shape set up
     * and eventually return correct values.
     * \param[in] s0 first span dimension
     * \param[in] s1 second span dimension
     * \param[in] s2 third span dimension
     */
    virtual void    checkSpan(double &s0, double &s1, double &s2)=0;

    /*!
     * Pure virtual method to check if your inf limits origin values fit your current shape set up
     * and eventually return correct values.
     * \param[in] orig inf limits origin array
     * \param[in] dir coordinate index
     * \return true, if valid new value is set.
     */
    virtual bool    checkInfLimits(double & orig, int &dir)=0;

    /*!
     * Pure virtual method to set local span & scaling vectors of your shape
     * \param[in] s0 first span dimension
     * \param[in] s1 second span dimension
     * \param[in] s2 third span dimension
     */
    virtual void    setScaling(const double &s0, const double &s1, const double &s2)=0;

    void        searchKdTreeMatches(bitpit::KdTree<3,bitpit::Vertex,long> & tree, livector1D & result, bool squeeze = true);
    void        searchBvTreeMatches(bitpit::PatchSkdTree & tree, bitpit::PatchKernel * geo, livector1D & result, bool squeeze = true);

    /*!
     * Pure virtual method to get the Axis Aligned Bounding Box of the current shape
     * and return it in m_bbox member.
     */
    virtual void    getTempBBox()=0;
};

/*!
 * \class Cube
 * \ingroup core
 * \brief Elementary Shape Representation of a Cube
 *
 * Volumetric Core Element, shaped as a cube, directly derived from BasicShape class.
 */

class Cube: public BasicShape {

public:

    //Costructors, Destructor, Copy/Assignment operators
    Cube();
    Cube(const darray3E &origin, const darray3E &span);
    ~Cube();

    Cube(const Cube &);

    //reimplementing pure virtuals

    darray3E    toWorldCoord(const darray3E &point);
    darray3E    toLocalCoord(const darray3E &point);
    darray3E    getLocalOrigin();
    bool    intersectShapeAABBox(const darray3E &bMin, const darray3E &bMax);

private:
    darray3E    basicToLocal(const darray3E &point);
    darray3E    localToBasic(const darray3E &point);
    void        checkSpan(double &, double &, double &);
    bool        checkInfLimits(double &, int & dir);
    void        setScaling(const double &, const double &, const double &);
    void        getTempBBox();
};

/*!
 *\class Cylinder
 *\ingroup core
 *\brief Elementary Shape Representation of a Cylinder or portion of it
 *
 * Volumetric Core Element, shaped as a cylinder, directly derived from BasicShape class.
 */

class Cylinder: public BasicShape {

public:
    //Costructors, Destructor, Copy/Assignment operators
    Cylinder();
    Cylinder(const darray3E &origin, const darray3E &span);
    ~Cylinder();

    Cylinder(const Cylinder &);

    //reimplementing pure virtuals
    darray3E	toWorldCoord(const darray3E &point);
    darray3E	toLocalCoord(const darray3E &point);
    darray3E	getLocalOrigin();
    bool		intersectShapeAABBox(const darray3E &bMin, const darray3E &bMax);

private:
    darray3E	basicToLocal(const darray3E &point);
    darray3E	localToBasic(const darray3E &point);
    void 		checkSpan(double &, double &, double &);
    bool 		checkInfLimits( double &, int &);
    void 		setScaling(const double &, const double &, const double &);
    void		getTempBBox();
};


/*!
 *\class Sphere
 *\ingroup core
 *\brief Elementary Shape Representation of a Sphere or portion of it
 *
 * Volumetric Core Element, shaped as a sphere, directly derived from BasicShape class.
 */

class Sphere: public BasicShape {

public:
    //Costructors, Destructor, Copy/Assignment operators
    Sphere();
    Sphere(const darray3E &origin, const darray3E &span);
    ~Sphere();

    Sphere(const Sphere &);

    //reimplementing pure virtuals
    darray3E    toWorldCoord(const darray3E &point);
    darray3E    toLocalCoord(const darray3E &point);
    darray3E    getLocalOrigin();
    bool        intersectShapeAABBox(const darray3E &bMin, const darray3E &bMax);

private:
    darray3E    basicToLocal(const darray3E &point);
    darray3E    localToBasic(const darray3E &point);
    void        checkSpan(double &, double &, double &);
    bool        checkInfLimits(double &, int &);
    void        setScaling(const double &, const double &, const double &);
    void        getTempBBox();
};

/*!
 * \class Wedge
 * \ingroup core
 * \brief Elementary Shape Representation of a Prism with triangular basis
 *
 * Volumetric Core Element, shaped as a triangular prism, directly derived from BasicShape class.
 * Local reference system is an elemental cube with x[0,1], y[0,1], z[-0.5,0.5] obtained by Duffy transformation.
 * Basic reference system is the elemental cube x[0,1], y[0,1], z[0,1];
 */

class Wedge: public BasicShape {

public:

    //Costructors, Destructor, Copy/Assignment operators
    Wedge();
    Wedge(const darray3E &origin, const darray3E &span);
    ~Wedge();

    Wedge(const Wedge &);

    //reimplementing pure virtuals
    darray3E    toWorldCoord(const darray3E &point);
    darray3E    toLocalCoord(const darray3E &point);
    darray3E    getLocalOrigin();
    bool        intersectShapeAABBox(const darray3E &bMin, const darray3E &bMax);

private:
    darray3E    basicToLocal(const darray3E &point);
    darray3E    localToBasic(const darray3E &point);
    void        checkSpan( double &, double &, double &);
    bool        checkInfLimits(double &, int &);
    void        setScaling(const double &, const double &, const double &);
    void        getTempBBox();

};



}

#endif //  __MIMMOBASICSHAPES_HH
