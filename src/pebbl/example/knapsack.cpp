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
  int size, rank;
  MPI_Init(&argc, &argv);
  MPI_Comm pebblComm, boundComm;
  pebbl::setupBoundingCommunicators(&argc, argv, MPI_COMM_WORLD, &pebblComm, &boundComm);
  pebbl::binKnapSub::setBoundComm(boundComm); 
  if (pebblComm == MPI_COMM_NULL){
    pebbl::binKnapSub::doBoundWork();
    MPI_Finalize();
  } 
  else {
    return driver<binaryKnapsack,parallelBinaryKnapsack>(argc,argv,pebblComm);
    pebbl::binKnapSub::finish();
    MPI_Finalize();
  }
}


