/*#---------------------------------------------------------------------------
#
#  mimmo
#
#  Copyright (C) 2015-2016 OPTIMAD engineering Srl
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of mimmo.
#
#  MIMMO is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License v3 (LGPL)
#  as published by the Free Software Foundation.
#
#  MIMMO is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with MIMMO. If not, see <http://www.gnu.org/licenses/>.
#
#---------------------------------------------------------------------------*/

using namespace std;
using namespace bitpit;
using namespace mimmo;

// =================================================================================== //
/*!
	\example MIMMO_example_00001.cpp

	\brief some test that still need to be provided

	<b>To run</b>: ./MIMMO_example_00001 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
*/

int main(int argc, char *argv[]) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);
	
#if BITPIT_ENABLE_MPI==1
	MPI_Init(&argc, &argv);

	{
#endif
		
		
#if BITPIT_ENABLE_MPI==1
	}

	MPI_Finalize();
#endif
	
	return (1);
}
