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
#include <pebbl/pbb/parTeamBranching.h>

using namespace pebbl;
using namespace std;


// If not parallel, make a dummy definition of the parallelBranching class

#ifndef ACRO_HAVE_MPI
typedef void parallelBinaryKnapsack;
#endif

// These aren't defined yet

class teamBinaryKnapsack : public virtual binaryKnapsack, public virtual teamBranching {


  void teamOrganize() { alertBound(); }

  void minionBound() {
    std::cout << "Minion " << teamComm.myRank() << " bounding\n" << std::endl;
  }

  void minionSplit() {
    std::cout << "Minion " << teamComm.myRank() << " separating\n" << std::endl;
  }

  void minionMakeChild() {
    std::cout << "Minion " << teamComm.myRank() << " making child\n" << std::endl;
  }

  public:

  teamBinaryKnapsack(MPI_Comm _comm) :
    binaryKnapsack(),
    teamBranching(_comm) { }
};

class parallelTeamBinaryKnapsack : public virtual parallelBinaryKnapsack, public virtual parallelTeamBranching {

  void teamOrganize() { alertBound(); }


  void minionBound() {
    std::cout << "Minion " << teamComm.myRank() << " bounding\n" << std::endl;
  }

  void minionSplit() {
    std::cout << "Minion " << teamComm.myRank() << " separating\n" << std::endl;
  }

  void minionMakeChild() {
    std::cout << "Minion " << teamComm.myRank() << " making child\n" << std::endl;
  }

  public:
    parallelTeamBinaryKnapsack(MPI_Comm _comm) :
      parallelTeamBranching(_comm) { }


  void reset(bool VBflag=true)
    {
      binaryKnapsack::reset();
      registerFirstSolution(new binKnapSolution(this));
      parallelTeamBranching::reset();
    }

    bool setup(int& argc, char**& argv){
      parallelTeamBranching::setup(argc, argv);
    }

    virtual void printSolution(const char* header,
		      const char* footer,
		      std::ostream& outStream){
      parallelBinaryKnapsack::printSolution(header, footer, outStream);
    }

};

int main(int argc, char* argv[])
{
  return driver<binaryKnapsack, parallelBinaryKnapsack, teamBinaryKnapsack, parallelTeamBinaryKnapsack>(argc, argv);
}
