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


#include <pebbl_config.h>
#include <pebbl/misc/fundamentals.h>
#include <pebbl/misc/clustering.h>

using namespace std;

namespace pebbl {


void clusterObj::reset(int rank_,
		       int size_,
		       int sizeWanted,
		       int clustersWanted,
		       int forceSeparateSize,
		       int boundingGroupSize)
{
  rank = rank_;
  size = size_;

  // workersPerCluster includes both the followers and the bounding minions
  int workersPerCluster = std::max(clustersWanted,1) * boundingGroupSize; 
  typicalSize = boundingGroupSize * (int) ceil(((double) size)/workersPerCluster);
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

} // namespace pebbl

