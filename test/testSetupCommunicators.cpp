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
	     int forceSeparateSize);

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
		       int forceSeparateSize)
{
  rank = rank_;
  size = size_;

  typicalSize = (int) ceil(((double) size)/std::max(clustersWanted,1));
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


namespace uMPI {
  MPI_Comm boundComm;
  MPI_Comm headComm;
}

void setupCommunicators(MPI_Comm comm_,
		int hubsDontWorkSize,
		int clusterSize,
		int boundingGroupSize)
{
  int headRank = -1;
  int boundRank = -1;
  
  clusterObj worldCluster; // use to find pure hubs, then throw away
  int worldRank, worldSize;
  MPI_Comm_rank(comm_, &worldRank);
  MPI_Comm_size(comm_, &worldSize);
  
  bool hubsWork = hubsDontWorkSize > clusterSize;
  int fullClusterSize = !hubsWork + 
	  (clusterSize - !hubsWork) * boundingGroupSize;
  int clustersWanted = worldSize / fullClusterSize;
  int forceSeparateSize = 1 + (hubsDontWorkSize - 1) * boundingGroupSize;
  worldCluster.reset(worldRank, worldSize, fullClusterSize, clustersWanted,
		  forceSeparateSize);
 
  int inBoundingGroup = worldCluster.isFollower();
  MPI_Comm boundingProcessors; // one for all pure hubs, one for people in bounding groups
  MPI_Comm_split(comm_, inBoundingGroup, worldRank, &boundingProcessors);
  int boundingProcessorsRank;
  MPI_Comm_rank(boundingProcessors, &boundingProcessorsRank);
  if (inBoundingGroup)
  {
    int groupNum = boundingProcessorsRank / boundingGroupSize;
    MPI_Comm_split(boundingProcessors, groupNum, worldRank, &uMPI::boundComm); 
  }
  else
  {
    MPI_Comm_free(&boundingProcessors);
    boundingProcessors = MPI_COMM_NULL;
    boundingProcessorsRank = -1;
  }
  int isWorker = (boundingProcessorsRank % boundingGroupSize == 0); 
  int notMinion = worldCluster.isLeader() || isWorker;
  MPI_Comm_split(comm_, notMinion, worldRank, &uMPI::headComm);
  if (!notMinion)
  {
    MPI_Comm_free(&uMPI::headComm);
    uMPI::headComm = MPI_COMM_NULL;
  } 

  // might not need isLeader if !isFollower => isLeader
  bool pureHub = !worldCluster.isFollower() && worldCluster.isLeader();
  // set variables for printing
  if (uMPI::headComm != MPI_COMM_NULL)
  {
    MPI_Comm_rank(uMPI::headComm, &headRank);
  }
  if (uMPI::boundComm != MPI_COMM_NULL && !pureHub)
  {
    MPI_Comm_rank(uMPI::boundComm, &boundRank);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  // printouts for debugging
  std::stringstream stream;
  stream << worldRank << " " << headRank << " " << boundRank << 
	  " " << worldCluster.isFollower() << " " << worldCluster.isLeader() << std::endl;
  std::cout << stream.str();
}

int main (int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  if (argc < 4)
  {
    std::cout << "USAGE: hubsDontWorkSize clusterSize boundingGroupSize" << std::endl;
    return 1;
  }
  
  int hubsDontWorkSize = strtol(argv[1], NULL, 10);
  int clusterSize = strtol(argv[2], NULL, 10);
  int boundingGroupSize = strtol(argv[3], NULL, 10);
  
  setupCommunicators(MPI_COMM_WORLD, hubsDontWorkSize, clusterSize, boundingGroupSize);
  
  MPI_Finalize();
  return 0;
}
