#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <sstream>
#include <iostream>

MPI_Comm headComm = MPI_COMM_NULL;
MPI_Comm boundComm = MPI_COMM_NULL;
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

  bool hubsWork = hubsDontWorkSize > clusterSize;
  int procsPerCluster = !hubsWork + (clusterSize - !hubsWork) * boundingGroupSize;
  // todo: adjust procsPerCluster based on minClusters parameter

  int numClusters = worldSize / procsPerCluster;
  int lastClusterNonHubs = worldSize - (numClusters - 1) * procsPerCluster - !hubsWork;
  int lastClusterWorkers = 1 + (lastClusterNonHubs - 1) / boundingGroupSize; // ceiling division
  int lastSeparated = lastClusterWorkers >= hubsDontWorkSize; 
  int inLastCluster = (worldRank / procsPerCluster) == numClusters;

  isHub = (worldRank % procsPerCluster) == 0;
  isWorker = !isHub || hubsWork || (inLastCluster && !lastSeparated);

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
  std::stringstream stream;
  stream << worldRank << " " << headRank << " " << boundRank << 
          " " << isHub << " " << isWorker << std::endl;
  for (int i = 0; i < worldSize; i++) {
    if (i == worldRank)
      std::cout << stream.str();
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  if (argc < 5) {
    std::cout << "USAGE: hubsDontWorkSize clusterSize boundingGroupSize numClusters" << std::endl;
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

