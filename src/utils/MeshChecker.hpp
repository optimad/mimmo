#ifndef __MESH_CHECKER_HPP__
#define __MESH_CHECKER_HPP__

#include <BaseManipulation.hpp>


namespace mimmo{

/*!
 * \class MeshChecker
 * \ingroup utils
 * \brief MeshChecker is the class to evaluate the quality of a mesh
 *
 * The available quality indices computed are :
 * - minimum and maximum cell volume
 * - maximum skweness angle inside the mesh and at boundaries
 * - minimum and maximum face validity
 * - minimum volume change ratio
 *
 * It writes a resume of the quality mesh check directly on the mimmo::Logger.
 * Optional results writes sick elements on file vtu.

 * \n
 *
 * Ports available in MeshChecker Class :
 *
 *    =========================================================

   |                     Port Input   ||                                     |
   |------------------|---------------------|----------------------|
   | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | M_GEOM           | setGeometry         | (MC_SCALAR, MD_MIMMO_)     |


   |               Port Output    ||                                         |
   |------------------|--------------------|-----------------------|
   | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | M_VALUEB          | isGood                 | (MC_SCALAR, MD_BOOL)     |
   | M_VALUEI          | getQualityStatusInt    | (MC_SCALAR, MD_INT)     |

 * =========================================================
 * \n
 *
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation/MimmoFvMesh:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.MeshChecker</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 *
 * Proper of the class:
 * - <B>MinVolTOL</B>: tolerance for minimum volume allowable.
 * - <B>MaxVolTOL</B>: tolerance for maximum volume allowable.
 * - <B>MaxSkewTOL</B>: tolerance for maximum skewness allowable.
 * - <B>MaxBoundarySkewTOL</B>: tolerance for maximum skewness on boundary allowable.
 * - <B>MinFaceValidityTOL</B>: tolerance for maximum skewness on boundary allowable.
 * - <B>MinVolChangeTOL</B>: tolerance for maximum skewness on boundary allowable.
 * - <B>ResumeFile</B>: boolean 1-print Resume file, 0 do nothing.

 * Geometry has to be mandatorily passed by port.
 *
 */
class MeshChecker: public BaseManipulation {

public:

    /*!
     * \enum CMeshOutput
     * \ingroup utils
     * Check error flag given by MehsChecker class.
     */
    enum CMeshOutput{
        /*!Check mesh not run.*/               NOTRUN = -1,
    	/*!All checks passed.*/                GOOD = 0,
    	/*!Face validity check failed.*/       FACEVALIDITY = 1,
    	/*!Volume change ratio check failed.*/ VOLUMECHANGERATIO = 2,
    	/*!Boundary Skewness check failed.*/   BOUNDARYSKEWNESS = 3,
    	/*!Skewness check failed.*/            SKEWNESS = 4,
    	/*!Minimum Volume check failed.*/      MINIMUMVOLUME = 5,
    	/*!Maximum Volume check failed.*/      MAXIMUMVOLUME = 6
    };


	MeshChecker();
	MeshChecker(const bitpit::Config::Section & rootXML);
	virtual ~MeshChecker();

	MeshChecker(const MeshChecker & other);
	MeshChecker & operator=(MeshChecker other);

    void setGeometry(MimmoObject * obj);
    void setMinimumVolumeTolerance(double tol);
	void setMaximumVolumeTolerance(double tol);
	void setMaximumSkewnessTolerance(double tol);
	void setMaximumBoundarySkewnessTolerance(double tol);
	void setMinimumFaceValidityTolerance(double tol);
	void setMinimumVolumeChangeTolerance(double tol);
    void setPrintResumeFile(bool flag);

    bool isGood();
    CMeshOutput getQualityStatus();
    int  getQualityStatusInt();
	void execute();
	void buildPorts();

	void plotOptionalResults();
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");


protected:

    void swap(MeshChecker & x) noexcept;
	void setDefault();
    CMeshOutput  checkMeshQuality();
    bool checkVolume();
	bool checkSkewness();
	bool checkFaceValidity();
	void initializeVolumes();
	void clear();
    void printResumeFile();

protected:

	// check values
	double	m_minVolume;            /**< minVolume variable*/
	double	m_maxVolume;            /**< maxVolume variable*/
	double	m_maxSkewness;          /**< maxSkewness variable*/
	double	m_maxSkewnessBoundary;  /**< maxSkewness on Boundary variable*/
	double	m_minFaceValidity;      /**< minFaceValidity variable*/
	double	m_maxFaceValidity;      /**< maxFaceValidity variable*/
	double	m_minVolumeChange;      /**< minVolumeChange variable*/

	// tolerance values
	double	m_minVolumeTol;             /**< check tolerance for minVolume variable*/
	double	m_maxVolumeTol;             /**< check tolerance for maxVolume variable*/
	double	m_maxSkewnessTol;           /**< check tolerance for maxSkewness variable*/
	double	m_maxSkewnessBoundaryTol;   /**< check tolerance for maxSkewness on Boundary variable*/
	double	m_minFaceValidityTol;       /**< check tolerance for minFaceValidity variable*/
	double	m_minVolumeChangeTol;       /**< check tolerance for minVolumeChange variable*/

	bool		m_isGood;		/**< true is good, false there are some errors*/
    CMeshOutput	m_qualityStatus;/**< Quality check flag of the mesh related to tolerance exceeded (first error encountered with this order): -1, not computed; 0, mesh good with input tolerances;
                                    1, minimum face validity too low; 2, volume change error; 3, skewness on boundary error; 4, skewness error;
                                    5, minimum volume error; 6, max volume error.*/

	//Temporary aux variables
	MimmoPiercedVector<double>	m_volumes; /**<Cell volumes.*/

	std::unique_ptr<MimmoObject>	m_volume; /**<Cells with poor volume.*/
	std::unique_ptr<MimmoObject>	m_skewness; /**<Cells with poor skewness.*/
	std::unique_ptr<MimmoObject>	m_facevalidity; /**<Cells with poor face validity.*/
	std::unique_ptr<MimmoObject>	m_volumechange; /**<Cells with poor volume change ratio.*/

    bool m_printResumeFile; /**< true, print a resume file after checking */

};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_,__MESH_CHECKER_HPP__)
REGISTER_PORT(M_VALUEB, MC_SCALAR, MD_BOOL,__MESH_CHECKER_HPP__)
REGISTER_PORT(M_VALUEI, MC_SCALAR, MD_INT,__MESH_CHECKER_HPP__)


REGISTER(BaseManipulation, MeshChecker, "mimmo.MeshChecker")
}

#endif
