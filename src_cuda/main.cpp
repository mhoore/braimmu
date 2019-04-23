#include <mpi.h>
#include "math.h"

#include "pointers.h"
#include "brain.h"

using namespace std;
using namespace brain_NS;

/* ----------------------------------------------------------------------*/
int main(int narg, char **arg) {

  MPI_Init(NULL,NULL);
  int me,nproc;


  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  Brain *brain = new Brain(narg,arg,me,nproc);

  if (!brain->me && narg > 2 || narg < 2) {
    printf("ERROR: enter arg1=input_file_name. \n");
    exit (EXIT_FAILURE);
  }

  delete brain;

  MPI_Finalize();

  if (!me)
    printf("Finished! \n");

  return 0;
}
