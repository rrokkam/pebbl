/*  _________________________________________________________________________
 *
 *  Acro: A Common Repository for Optimizers
 *  Copyright (c) 2008 Sandia Corporation.
 *  This software is distributed under the BSD License.
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *  For more information, see the README.txt file in the top Acro directory.
 *  _________________________________________________________________________
 */

#include <pebbl_config.h>
#include <pebbl/example/parKnapsack.h>

using namespace pebbl;
using namespace std;


// If not parallel, make a dummy definition of the parallelBranching class

#ifndef ACRO_HAVE_MPI
typedef void parallelBinaryKnapsack;
#endif


int main(int argc, char* argv[])
{
#ifdef ACRO_HAVE_MPI
	MPI_Init(&argc, &argv);
	CommonIO::begin();
	CommonIO::setIOFlush(1);
	CommonIO::begin_tagging();
	MPI_Comm pebblComm, boundingComm;
	pebbl::setupBoundingCommunicators(&argc, argv, 
	  MPI_COMM_WORLD, &pebblComm, &boundingComm);


	bool pebblProc = false;
	int pebblCommSize;
	if (pebblComm != MPI_COMM_NULL) {
		pebblProc = true;
		MPI_Comm_size(pebblComm, &pebblCommSize);
	}
	if (pebblProc && pebblCommSize == 1) {
		binaryKnapsack instance;
		instance.boundComm = boundingComm;
		bool flag = instance.setup(argc,argv);
		if (flag)
		{
			instance.reset();
			instance.solve();
			instance.endBounders();
		}
	}
	else if (pebblProc) {
		parallelBinaryKnapsack instance(pebblComm);
		instance.boundComm = boundingComm;
		bool flag = instance.setup(argc,argv);
		if (flag)
		{
			instance.reset();
			instance.printConfiguration();
			instance.solve();
			instance.endBounders();
		}
	}
	else {
		binaryKnapsack instance;
		instance.boundComm = boundingComm;
		instance.doBoundWork();
	}
	CommonIO::end();
	uMPI::done();
#endif
	return 0;
}
