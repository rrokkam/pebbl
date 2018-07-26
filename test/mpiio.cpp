#include <mpi.h>
#include <iostream>

// With MPICH, running this with any number of processors gives a return value
// of 6, but still successfully prints to stdout. 
//
// Another question: which process' return value is given to the terminal?

void syncPrint(int val, int rank, int size) {
	for(int i = 0; i < size; i++) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (i == rank) {
			std::cout << "Proc " << i << ": " << val << std::endl;
		}
	}
}

int main(int argc, char *argv[]) {
	std::cout << "We can print stuff" << std::endl;
	MPI_Init(&argc, &argv);
	int error, flag, rank, haveIO, ret, size;
	MPI_Aint *result;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	error = MPI_Comm_get_attr(MPI_COMM_SELF, MPI_IO, &result, &flag);
	if      (error != MPI_SUCCESS)      ret = 1;  // failed to get MPI_IO attribute
	else if (!flag)                     ret = 2;  // some other MPI failure?
	else if (*result == rank)           ret = 3;  // We are an IO processor
	else if (*result == MPI_ANY_SOURCE) ret = 4;  // Anyone can do IO
	else if (*result == MPI_PROC_NULL)  ret = 5;  // No one can do IO
	else                                ret = 6;  // We cannot do IO

	syncPrint(ret, rank, size);
	syncPrint((int)*result, rank, size);

	MPI_Finalize();
	return ret;
}

