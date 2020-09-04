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

#include "Create3DCurve.hpp"
#include <queue>

namespace mimmo{

/*!
 * Default constructor
 */
Create3DCurve::Create3DCurve(){
    m_name         = "mimmo.Create3DCurve";
    m_closed = false;
    m_nCells = 0;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
Create3DCurve::Create3DCurve(const bitpit::Config::Section & rootXML){

    m_name         = "mimmo.Create3DCurve";
    m_closed = false;
    m_nCells = 0;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.Create3DCurve"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Default destructor
 */
Create3DCurve::~Create3DCurve(){};

/*!
    Copy constructor.
    /param[in] other to be copied
 */
Create3DCurve::Create3DCurve(const Create3DCurve & other):BaseManipulation(other){
    m_closed = other.m_closed;
    m_nCells = other.m_nCells;
    m_rawpoints    = other.m_rawpoints;
    m_rawvector    = other.m_rawvector;
    m_rawscalar    = other.m_rawscalar;
};

/*!
 * Assignement operator.
  /param[in] other to be copied
 */
Create3DCurve & Create3DCurve::operator=(Create3DCurve other){
    swap(other);
    return *this;
};

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void Create3DCurve::swap(Create3DCurve &x) noexcept
{
    std::swap(m_nCells, x.m_nCells);
    std::swap(m_closed , x.m_closed);
    m_rawpoints.swap(x.m_rawpoints);
    m_rawscalar.swap(x.m_rawscalar);
    m_rawvector.swap(x.m_rawvector);
    m_scalarfield.swap(x.m_scalarfield);
    m_vectorfield.swap(x.m_vectorfield);

    BaseManipulation::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
Create3DCurve::buildPorts(){
    bool built = true;

    built = (built && createPortIn<dvecarr3E, Create3DCurve>(this, &mimmo::Create3DCurve::setRawPoints, M_COORDS, true,1));
    built = (built && createPortIn<dvecarr3E, Create3DCurve>(this, &mimmo::Create3DCurve::setRawVectorField, M_DISPLS));
    built = (built && createPortIn<dvector1D, Create3DCurve>(this, &mimmo::Create3DCurve::setRawScalarField, M_DATAFIELD));

    built = (built && createPortIn<dmpvecarr3E*, Create3DCurve>(this, &mimmo::Create3DCurve::setRawPoints, M_VECTORFIELD, true,1));
    built = (built && createPortIn<dmpvecarr3E*, Create3DCurve>(this, &mimmo::Create3DCurve::setRawVectorField, M_VECTORFIELD2));
    built = (built && createPortIn<dmpvector1D*, Create3DCurve>(this, &mimmo::Create3DCurve::setRawScalarField, M_SCALARFIELD));

    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, Create3DCurve>(this, &mimmo::Create3DCurve::getGeometry, M_GEOM));
    built = (built && createPortOut<dmpvector1D*, Create3DCurve>(this, &mimmo::Create3DCurve::getScalarField, M_SCALARFIELD));
    built = (built && createPortOut<dmpvecarr3E*, Create3DCurve>(this, &mimmo::Create3DCurve::getVectorField, M_VECTORFIELD));

    m_arePortsBuilt = built;
}

/*!
 * Return the scalar field stored in the class as pointer to MimmoPiercedVector object.
 * \return scalar field stored in the class
 */
dmpvector1D*
Create3DCurve::getScalarField(){
    return &m_scalarfield;
};

/*!
 * Return the vector field stored in the class as pointer to MimmoPiercedVector object.
 * \return vector field stored in the class
 */
dmpvecarr3E*
Create3DCurve::getVectorField(){
    return &m_vectorfield;
};


/*!
 * It sets the point cloud raw list
 * \param[in] rawPoints raw list of points.
 */
void
Create3DCurve::setRawPoints(dmpvecarr3E * rawPoints){
    if(rawPoints)    m_rawpoints = *rawPoints;
};

/*!
 * It sets the point cloud raw list
 * \param[in] rawPoints raw list of points.
 */
void
Create3DCurve::setRawPoints(dvecarr3E rawPoints){
    m_rawpoints.clear();
    m_rawpoints.reserve(rawPoints.size());
    long count(0);
    for(std::array<double,3> &val : rawPoints){
        m_rawpoints.insert(count, val);
        ++count;
    }
};

/*!
 * It sets the vector field optionally attached to point cloud raw list
 * \param[in] rawVectorField raw vector field.
 */
void
Create3DCurve::setRawVectorField(dmpvecarr3E *rawVectorField){
    if(rawVectorField)  m_rawvector = *rawVectorField;
};

/*!
 * It sets the vector field optionally attached to point cloud raw list
 * \param[in] rawVectorField raw vector field.
 */
void
Create3DCurve::setRawVectorField(dvecarr3E rawVectorField){
    m_rawvector.clear();
    m_rawvector.reserve(rawVectorField.size());
    long count(0);
    for(std::array<double,3> &val : rawVectorField){
        m_rawvector.insert(count, val);
        ++count;
    }
};

/*!
 * It sets the scalar field optionally attached to point cloud raw list
 * \param[in] rawScalarField raw scalar field.
 */
void
Create3DCurve::setRawScalarField(dmpvector1D * rawScalarField){
    if(rawScalarField)  m_rawscalar = *rawScalarField;
};

/*!
 * It sets the scalar field optionally attached to point cloud raw list
 * \param[in] rawScalarField raw scalar field.
 */
void
Create3DCurve::setRawScalarField(dvector1D rawScalarField){
    m_rawscalar.clear();
    m_rawscalar.reserve(rawScalarField.size());
    long count(0);
    for(double &val : rawScalarField){
        m_rawscalar.insert(count, val);
        ++count;
    }
};

/*!
 * Clear all data stored in the class
 */
void
Create3DCurve::clear(){
    m_vectorfield.clear();
    m_scalarfield.clear();
    m_rawscalar.clear();
    m_rawvector.clear();
    m_rawpoints.clear();
    m_nCells = 0;
    m_closed = false;
    BaseManipulation::clear();
}


/*!
 * \return true if the class is set to create a final closed 3D Curve
*/
bool
Create3DCurve::isClosedLoop(){
    return m_closed;
}

/*!
 * Set your final 3D curve as open/closed.
 * \param[in] flag false to get an open curve, true for a closed one.
 */
void
Create3DCurve::setClosedLoop(bool flag){
    m_closed = flag;
}

/*!
 * Set target cells of the 3D curve final tessellation, in order to refine/regularize it.
   This number must be >= of the number of raw point list. If it is,
   the class will split in half the longest segments of tessellation, until the target
   number of cells is achieved.
 * \param[in] ncells number of target cells.
 */
void
Create3DCurve::setNCells(int ncells){
    m_nCells = ncells;
}


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void Create3DCurve::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("ClosedLoop")){
        input = slotXML.get("ClosedLoop");
        bool flag = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> flag;
        }
        setClosedLoop(flag);
    };

    if(slotXML.hasOption("nCells")){
        input = slotXML.get("nCells");
        int value = 0;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setNCells(value);
    };

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void Create3DCurve::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("ClosedLoop", std::to_string(int(m_closed)));
    slotXML.set("nCells", std::to_string(int(m_nCells)));
};

/*!
 * Plot 3DCurve and store it in *.vtu file
 */
void
Create3DCurve::plotOptionalResults(){

    if (getGeometry() == nullptr) return;
    BaseManipulation::write(getGeometry(), m_scalarfield, m_vectorfield);
}


/*!
 * Execute command.
 */
void
Create3DCurve::execute(){

    int np = int(m_rawpoints.size());

// rawpoints are shared among all procs or retained by 0 only.
#if MIMMO_ENABLE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &np, 1 , MPI_INT, MPI_MAX, m_communicator);
#endif

    if(np == 0){
        (*m_log)<< "Warning in "<<m_name<<" : no raw points to work on"<<std::endl;
    }

    m_nCells = std::max(m_nCells, (np - 1 + int(m_closed)));

    //create a brand new 3DCurve container
    {
        MimmoSharedPointer<MimmoObject> dum(new MimmoObject(3));
        m_geometry = dum;
    }

    // prepare connectivity of curve and handling refinement if needed
    // then fill the curve container
    dmpvecarr3E points, vector;
    dmpvector1D scalar;
    std::unordered_map<long, std::array<long,2> > connectivity;
#if MIMMO_ENABLE_MPI
    //leave the job to the master rank
    if(m_rank == 0)
#endif
    {
        points = m_rawpoints;
        scalar.reserve(points.size());
        vector.reserve(points.size());
        for(auto it= points.begin(); it!= points.end(); ++it){

            auto itscalar = scalar.insert(it.getId(), 0.0);
            if(m_rawscalar.exists(it.getId())){
                *itscalar = m_rawscalar[it.getId()];
            }
            auto itvector = vector.insert(it.getId(), {{0., 0.,0.}});
            if(m_rawvector.exists(it.getId())){
                *itvector = m_rawvector[it.getId()];
            }

        }

        int furtherCells = fillPreliminaryStructure(points, connectivity);
        refineObject(points, scalar, vector, connectivity, furtherCells);

        m_geometry->getPatch()->reserveVertices(points.size());
        m_geometry->getPatch()->reserveCells(connectivity.size());

        long id;
        for(auto it= points.begin(); it != points.end(); ++it){
            id = it.getId();
            m_geometry->addVertex(*it, id);
        }

        for(auto & tuple : connectivity){
            m_geometry->addConnectedCell(livector1D(tuple.second.begin(),tuple.second.end()), bitpit::ElementType::LINE, 0, tuple.first, 0);
        }

    }

    //initialize and fill the final fields
    m_scalarfield.initialize(m_geometry,MPVLocation::POINT, 0.);
    m_scalarfield.setName("3DCurveScalar");

    m_vectorfield.initialize(m_geometry,MPVLocation::POINT, {{0.,0.,0.}});
    m_vectorfield.setName("3DCurveVector");

    for(auto it=scalar.begin(); it!=scalar.end(); ++it){
        if(m_scalarfield.exists(it.getId()))    m_scalarfield[it.getId()] = *it;
    }
    for(auto it=vector.begin(); it!=vector.end(); ++it){
        if(m_vectorfield.exists(it.getId()))    m_vectorfield[it.getId()] = *it;
    }

    m_geometry->buildAdjacencies();
#if MIMMO_ENABLE_MPI
    m_geometry->update();
#endif


};


/*!
 * Fill a connectivity structure with points data set provided by the User.
 * \param[in]    points  reference to structure of 3D points
 * \param[out]    connectivity  reference to connectivity 3D curve
 * \return number of further cells to add required by the class;
 */
int
Create3DCurve::fillPreliminaryStructure(dmpvecarr3E & points, std::unordered_map<long, std::array<long,2> > &connectivity){
    int nVerts = (int) points.size();
    int nCells = nVerts -1 + (int) m_closed;
    int res = std::max(0, (m_nCells - nCells));

    auto itP = points.begin();
    int countCell(0);
    while (countCell < (nCells -(int)m_closed)){
        connectivity[countCell][0] = itP.getId();
        ++itP;
        connectivity[countCell][1] = itP.getId();
        ++countCell;
    }

    if(m_closed){
        itP = points.begin();
        auto itE = points.end();
        --itE;
        connectivity[nCells-1][0] = itE.getId();
        connectivity[nCells-1][1] = itP.getId();
    }

    return (res);
}

/*!
 * Refine 3D curve object adding a number of given cells.
 * Each segment cell will be splitted in two new cells. Priority will be given
 * to longest segments.
   For each new point pushed, data fields attached will be interpolated onto new points.
 * \param[in,out]    points  reference to structure of 3D points
 * \param[in,out]    scalarf  reference to scalar field attached
 * \param[in,out]    vectorf  reference to vector field attached
 * \param[in,out]    connectivity  reference to connectivity 3D curve
 * \param[in]       fCells number of cells to add.
 */
void
Create3DCurve::refineObject(dmpvecarr3E & points, dmpvector1D & scalarf, dmpvecarr3E & vectorf, std::unordered_map<long, std::array<long,2> > &connectivity, int fCells){
    if ( fCells == 0) return;

    greatDist mycomp;
    double dist;

    long idV(0);
    {
        livector1D ids = points.getIds();
        if(!ids.empty())    idV = *(std::max_element(ids.begin(), ids.end()));
    };
    long idC(0);
    std::vector<std::pair<double, long>> values;
    values.reserve(int(connectivity.size()) + fCells);
    for (auto &conn : connectivity){
        dist = norm2(points[conn.second[1]] - points[conn.second[0]]);
        values.push_back(std::make_pair(dist, conn.first));
        idC = std::max(idC, conn.first);
    }

    //prep your priority queue
    std::priority_queue<std::pair<double,long>, std::vector<std::pair<double,long>>, greatDist > pQue(mycomp, values);

    int counterCell = 0;
    darray3E point;
    while (counterCell < fCells){

        //extract the longest segment
        auto pp = pQue.top();
        pQue.pop();

        idV++;
        idC++;
        //divide it by two and create a new point
        point = 0.5*(points[connectivity[pp.second][1]] + points[connectivity[pp.second][0]]);
        points.insert(idV, point);
        //update scalar and vector
        scalarf.insert(idV, 0.5*(scalarf[connectivity[pp.second][1]] + scalarf[connectivity[pp.second][0]]));
        vectorf.insert(idV, 0.5*(vectorf[connectivity[pp.second][1]] + vectorf[connectivity[pp.second][0]]));

        //update connectivity
        connectivity[idC]= {{ idV, connectivity[pp.second][1] }};
        connectivity[pp.second][1] = idV;

        //update the priority queue
        pQue.push(std::make_pair(0.5*pp.first, pp.second));
        pQue.push(std::make_pair(0.5*pp.first, idC));

        //adding a new cell to the count;
        ++counterCell;
    }
};


}
