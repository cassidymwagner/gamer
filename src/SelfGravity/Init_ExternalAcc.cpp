#include "GAMER.h"
#include <stdlib.h>
#include <stdio.h>

#ifdef GRAVITY


#include "CUPOT.h"
double ExtAcc_AuxArray[EXT_ACC_NAUX_MAX];

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Init_ExternalAcc();

// this function pointer may be overwritten by various test problem initializers
void (*Init_ExternalAcc_Ptr)() = Init_ExternalAcc;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExternalAcc
// Description :  Set the array "ExtAcc_AuxArray" used by the external acceration routines
//                "CUPOT_ExternalAcc.cu / CPU_ExternalAcc.cpp"
//
// Note        :  1. Invoked by "Init_GAMER" using the function pointer "Init_ExternalAcc_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Enabled by the runtime option "OPT__GRAVITY_TYPE == 2/3"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_ExternalAcc()
{

// ExtAcc_AuxArray has the size of EXT_ACC_NAUX_MAX defined in CUPOT.h (default = 10)
// --> by default we set
//     ExtAcc_AuxArray[0] = x coordinate of the external acceleration center
//     ExtAcc_AuxArray[1] = y ...
//     ExtAcc_AuxArray[2] = z ..
//     ExtAcc_AuxArray[3] = gravitational_constant*point_source_mass
//     ExtAcc_AuxArray[4] = soften_length (<=0.0 --> disable)
// --> to change the this default behavior, please edit "GPU_Gravity/CUPOT_ExternalAcc.cu"

  size_t nread = 0;
  FILE *bin_data = fopen("UM_IC", "r");

  float *density = (float *) calloc(256*256*256, sizeof(float));
  nread = fread(density, sizeof(float), 256*256*256, bin_data);

  float *field_x = (float *) calloc(256*256*256, sizeof(float));
  nread = fread(field_x, sizeof(float), 256*256*256, bin_data);

  float *field_y = (float *) calloc(256*256*256, sizeof(float));
  nread = fread(field_y, sizeof(float), 256*256*256, bin_data);

  float *field_z = (float *) calloc(256*256*256, sizeof(float));
  nread = fread(field_z, sizeof(float), 256*256*256, bin_data);

  for (int i = 0; i < 256*256*256; i++) {
    field_x[i] /= density[i];
    field_y[i] /= density[i];
    field_z[i] /= density[i];
  }

  fprintf(stderr, "Last row: %0.3g %0.3g %0.3g %0.3g\n",
      density[256*256*256-1],
      field_x[256*256*256-1],
      field_y[256*256*256-1],
      field_z[256*256*256-1]);

  free(density);

  const double M   = 1.0;
  const double GM  = NEWTON_G*M;
  const double Eps = 0.0;

  ExtAcc_AuxArray[0] = field_x;
  ExtAcc_AuxArray[1] = field_y;
  ExtAcc_AuxArray[2] = field_z;
  ExtAcc_AuxArray[3] = GM;
  ExtAcc_AuxArray[4] = Eps;  

  
} // FUNCTION : Init_ExternalAcc



#endif // #ifdef GRAVITY
