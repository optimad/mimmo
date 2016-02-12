#ifndef MIMMOTYPEDEF_HH
#define MIMMOTYPEDEF_HH

#include <vector>
#include <array>

typedef std::vector<double> 	 dvector1D;
typedef std::vector<float>	 fvector1D;
typedef std::vector<int>	 ivector1D;
typedef std::vector<short int> 	 shivector1D;
typedef std::vector<bool> 	 bvector1D;
typedef std::vector<char> 	 cvector1D;
typedef std::vector<std::string> svector1D;

typedef std::array<double,2>    darray2E;
typedef std::array<double,3>	darray3E;
typedef std::array<double,4>	darray4E;
typedef std::array<float,3>	farray3E;

typedef std::vector<darray2E>	dvecarr2E;
typedef std::vector<darray3E>	dvecarr3E;
typedef std::vector<darray4E>	dvecarr4E;
typedef std::vector<farray3E>	fvecarr3E;

typedef std::vector<dvector1D>	dvector2D;
typedef std::vector<fvector1D>	fvector2D;


typedef std::vector<ivector1D>		ivector2D;
typedef std::vector<shivector1D>	shivector2D;
typedef std::vector<bvector1D>		bvector2D;
typedef std::vector<cvector1D>		cvector2D;
typedef std::vector<svector1D>		svector2D;

typedef std::vector< bvector2D >            bvector3D;
typedef std::vector< bvector3D >            bvector4D;

typedef std::vector< cvector2D >            cvector3D;
typedef std::vector< cvector3D >            cvector4D;

typedef std::vector< ivector2D >            ivector3D;
typedef std::vector< ivector3D >            ivector4D;

typedef std::vector< dvector2D >            dvector3D;
typedef std::vector< dvector3D >            dvector4D;

typedef std::vector< svector2D >            svector3D;
typedef std::vector< svector3D >            svector4D;


#endif //MIMMOTYPEDEF_HH
