#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <sstream>
#include <iostream>

using namespace std;
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

/**
 * \file clustering.h
 * \author Jonathan Eckstein
 *
 * Logic for organizing things (typically processors) into clusters.
 */

#ifndef pebbl_clustering_h
#define pebbl_clustering_h

//#include <pebbl_config.h>

//namespace pebbl {


class clusterObj
{
public:
  
  // Initialize things to avoid valgrind complaints when registering
  // parameters for parMIP.
  // Meaningful initialization is by reset(...) with real arguments.

  clusterObj() { /* reset(0,1,1,1,0);*/ };

  void reset(int rank_,
	     int size_,
	     int sizeWanted,
	     int clustersWanted,
	     int forceSeparateSize,
		 int boundingGroupSize);

  int rank;
  int size;

  bool separateFunctions;
  bool iAmLeader;
  bool iAmFollower;

  int clusterNumber;
  int typicalSize;
  int clusterSize;
  int lastClusterSize;
  int leader;
  int numPureLeaders;
  int numFollowers;
  int positionInCluster;
  int followerPosInCluster;
  int numClusters;
  int totalNumFollowers;

  bool typicallySeparated;
  bool lastSeparated;

  int typicalFollowers;
  int myGlobalFollowerRank;

  int follower(int whichFollower)
    {
      return leader + whichFollower + separateFunctions;
    };

  int whichFollower(int n)
    {
      return n - leader - separateFunctions;
    };

  int leaderOfCluster(int whichCluster)
    {
      return whichCluster*typicalSize;
    };

  int whichCluster(int n)
    {
      return n/typicalSize;
    };

  int hasPureLeader(int cNum)
    {
      if (cNum + 1 == numClusters)
	return lastSeparated;
      else
	return typicallySeparated;
    };

  // Translates an index into the total collection of followers
  // into an index into the total collection of things.

  int globalFollower(int n)
    {
      int targetCluster = n/typicalFollowers;
      int targetPos     = n - targetCluster*typicalFollowers;
      return (targetCluster*typicalSize + 
	      hasPureLeader(targetCluster) + targetPos);
    };

  // Translates an index into the total collection of things into
  // an index into the total collection of followers.

  int followerRank(int globalRank)
    {
      int targetCluster = whichCluster(globalRank);
      int relPos        = globalRank - targetCluster*typicalSize;
      return (targetCluster*typicalFollowers + 
	      relPos - hasPureLeader(targetCluster));
    };

  // Returns whether a given global index corresponds to a leader

  bool isLeader(/*int n*/)
  {
    return (/*n*/ rank % typicalSize) == 0;
  }

  // Returns whether a given global index corresponds to a follower

  bool isFollower(/*int n*/)
  {
    int n = rank;
    if (!isLeader(/*n*/))
      return true;
    if (whichCluster(n) == numClusters - 1)
      return !lastSeparated;
    return !typicallySeparated;
  }

};
  
//}

#endif


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
//
// clustering.cpp
//
// Logic for organizing things (typically processors) into clusters.
//
// Jonathan Eckstein
//


//#include <pebbl_config.h>
//#include <pebbl/misc/fundamentals.h>
//#include <pebbl/misc/clustering.h>

//namespace pebbl {


void clusterObj::reset(int rank_,
		       int size_,
		       int sizeWanted,
		       int clustersWanted,
		       int forceSeparateSize,
			   int boundingGroupSize)
{
  rank = rank_;
  size = size_;

  typicalSize = (int) boundingGroupSize * ceil(((double) size)/(boundingGroupSize * std::max(clustersWanted,1)));
  if (typicalSize > sizeWanted)
    typicalSize = sizeWanted;
  if (typicalSize < 1)
    typicalSize = 1;

  numClusters = (int) ceil(((double) size)/typicalSize);

  clusterNumber   = rank/typicalSize;
  leader          = clusterNumber*typicalSize;
  clusterSize     = std::min(typicalSize,size - leader);
  lastClusterSize = size - (numClusters - 1)*typicalSize;

  separateFunctions  = (clusterSize     >= forceSeparateSize);
  typicallySeparated = (typicalSize     >= forceSeparateSize);
  lastSeparated      = (lastClusterSize >= forceSeparateSize);

  typicalFollowers = typicalSize - typicallySeparated;

  numFollowers         = clusterSize - separateFunctions;
  positionInCluster    = rank - leader;
  followerPosInCluster = positionInCluster - separateFunctions;
  iAmLeader            = (positionInCluster == 0);
  iAmFollower          = !iAmLeader || !separateFunctions;

  numPureLeaders       = (numClusters - 1)*typicallySeparated + lastSeparated;
  totalNumFollowers    = size - numPureLeaders;
  myGlobalFollowerRank = clusterNumber*typicalFollowers + followerPosInCluster;
}

//} // namespace pebbl


MPI_Comm boundComm = MPI_COMM_NULL;
MPI_Comm headComm = MPI_COMM_NULL;
int worldSize = -1;
int worldRank = -1;
int headRank = -1;
int boundRank = -1;
int isHub = -1;
int isWorker = -1;

void setupCommunicators(int hubsDontWorkSize,
			int clusterSize,
			int boundingGroupSize,
			int minClusters)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  
  clusterObj worldCluster; // use to find pure hubs, then throw away
  
  bool hubsWork = hubsDontWorkSize > clusterSize;
  int procsPerGroup = !hubsWork + (clusterSize - !hubsWork) * boundingGroupSize;
  int clustersWanted = worldSize / procsPerGroup;
  int forceSeparateSize = 1 + (hubsDontWorkSize - 1) * boundingGroupSize;
  worldCluster.reset(worldRank, worldSize, procsPerGroup, minClusters, 
		  forceSeparateSize, boundingGroupSize);
 
  isHub = worldCluster.isLeader();
  isWorker = worldCluster.isFollower();
  
  MPI_Comm dummy; // communicator containing all bounding processors 
  int dummyRank;
  
  MPI_Comm_split(MPI_COMM_WORLD, isWorker, worldRank, &dummy);
  MPI_Comm_rank(dummy, &dummyRank);
  
  if (isWorker)
    MPI_Comm_split(dummy, dummyRank / boundingGroupSize, worldRank, &boundComm); 
  MPI_Comm_free(&dummy);
  
  int bounder = !isHub && (dummyRank % boundingGroupSize != 0);
  MPI_Comm_split(MPI_COMM_WORLD, bounder, worldRank, &headComm);
  
  if (bounder) {
    MPI_Comm_free(&headComm);
    headComm = MPI_COMM_NULL;
  } 

}

void computeRanks() {
  if (headComm != MPI_COMM_NULL)
    MPI_Comm_rank(headComm, &headRank);
  if (boundComm != MPI_COMM_NULL)
    MPI_Comm_rank(boundComm, &boundRank);
}

void printRanks() {
  stringstream stream;
  stream << worldRank << " " << headRank << " " << boundRank << 
  	  " " << isHub << " " << isWorker << endl;
  for (int i = 0; i < worldSize; i++) {
    if (i == worldRank)
      cout << stream.str();
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

int main (int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  if (argc < 5){
    cout << "USAGE: hubsDontWorkSize clusterSize boundingGroupSize numClusters" << endl;
    return 1;
  }
  
  int hdwSize = strtol(argv[1], NULL, 10);
  int clusterSize = strtol(argv[2], NULL, 10);
  int bgSize = strtol(argv[3], NULL, 10);
  int numClusters = strtol(argv[4], NULL, 10);
  
  setupCommunicators(hdwSize, clusterSize, bgSize, numClusters);
  computeRanks();
  printRanks();

  MPI_Finalize();
  return 0;
}
