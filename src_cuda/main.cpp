#include "scenario_geometry.h"
#include "scenario_connectome.h"
#include "pointers.h"

using namespace std;

/* ----------------------------------------------------------------------*/
int main(int narg, char **arg) {

  MPI_Init(NULL,NULL);
  int me,nproc;

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  if (me && narg > 3 || narg < 3) {
    printf("ERROR: enter arg1=scenario_type, arg2=input_file.\n"
           "Possible scenarios: geometry, connectome \n");
    exit (EXIT_FAILURE);
  }

  string scenario = arg[1];

  if (!scenario.compare("connectome"))
    new ScenarioConnectome(narg,arg,me,nproc);
  if (!scenario.compare("geometry"))
    new ScenarioGeometry(narg,arg,me,nproc);

  MPI_Finalize();

  if (!me)
    printf("Finished! \n");

  return 0;
}
