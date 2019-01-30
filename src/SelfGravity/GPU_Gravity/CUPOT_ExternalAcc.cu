#ifdef __CUDACC__
#include "Macro.h"
#else
#include "GAMER.h"
#endif
#include "CUPOT.h"

#ifdef GRAVITY

// soften length implementation
//#  define SOFTEN_PLUMMER
//#  define SOFTEN_RUFFERT
#   define DRIV_TURB

#  if   (defined DRIV_TURB)
extern double* ExtAcc_InitialField[3];
#  endif




//-----------------------------------------------------------------------------------------
// Function    :  CUPOT_ExternalAcc / CPU_ExternlAcc
// Description :  Calculate the external acceleration at the given coordinates and time
//
// Note        :  1. This function will be invoked by both CPU and GPU
//                2. "__forceinline__" is required since this device function will be invoked
//                   by more than one kernels (e.g., CUPOT_HydroGravitySolver, CUFLU_ComputeFlux)
//                3. The auxiliary array "UserArray" is set by "Init_ExternalAcc_Ptr", which
//                   points to "Init_ExternalAcc()" by default but may be overwritten by various
//                   test problem initializers
//                4. By default we assume
//                     UserArray[0] = x coordinate of the external acceleration center
//                     UserArray[1] = y ...
//                     UserArray[2] = z ..
//                     UserArray[3] = gravitational_constant*point_source_mass
//                     UserArray[4] = soften_length (<=0.0 --> disable)
//                   --> but one can easily modify this file to change the default behavior
//                5. Two different soften length implementations are supported
//                   --> SOFTEN_PLUMMER & SOFTEN_RUFFERT
//
// Parameter   :  Acc       : Array to store the output external acceleration
//                x/y/z     : Target spatial coordinates
//                Time      : Current physical time
//                UserArray : User-provided auxiliary array (set by "Init_ExternalAcc_Ptr")
//
// Return      :  Acc
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__forceinline__ __device__
void CUPOT_ExternalAcc( real Acc[], const double x, const double y, const double z, const double Time, const double UserArray[] )
#else
void   CPU_ExternalAcc( real Acc[], const double x, const double y, const double z, const double Time, const double UserArray[] )
#endif
{
   const double Cen[3] = { UserArray[0], UserArray[1], UserArray[2] };
   const real GM       = (real)UserArray[3];
   const real eps      = (real)UserArray[4];
   const real dx       = (real)(x - Cen[0]);
   const real dy       = (real)(y - Cen[1]);
   const real dz       = (real)(z - Cen[2]);
   const real r        = SQRT( dx*dx + dy*dy + dz*dz );

// Plummer
#  if   ( defined SOFTEN_PLUMMER )
   const real _r3 = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps), (real)-1.5 );

// Ruffert 1994
#  elif ( defined SOFTEN_RUFFERT )
   const real tmp = EXP( -SQR(r)/SQR(eps) );
   const real _r3 = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps)*tmp, (real)-1.5 )*( (real)1.0 - tmp );

#  else
   const real _r3 = (real)1.0/CUBE(r);
#  endif

   Acc[0] = -GM*_r3*dx;
   Acc[1] = -GM*_r3*dy;
   Acc[2] = -GM*_r3*dz;

#  if   (defined DRIV_TURB)
   int m_temp, ix, iy, iz;
  
   ix = (int) ((x - amr->BoxEdgeL[0])/(amr->BoxEdgeR[0] - amr->BoxEdgeL[0]) * 256);
   iy = (int) ((y - amr->BoxEdgeL[1])/(amr->BoxEdgeR[1] - amr->BoxEdgeL[1]) * 256);
   iz = (int) ((z - amr->BoxEdgeL[2])/(amr->BoxEdgeR[2] - amr->BoxEdgeL[2]) * 256);

   if (ix < 0) ix += 256;
   if (iy < 0) iy += 256;
   if (iy < 0) iz += 256;
   if (ix > 255) ix -= 256;
   if (iy > 255) iy -= 256;
   if (iz > 255) iz -= 256;

   if((ix < 0 || ix > 255) || (iy < 0 || iy > 255) || (iz < 0 || iz > 255))
     Aux_Message(stderr, "At %lf %lf %lf index %d %d %d\n",
        x, y, z, ix, iy, iz);

   m_temp = (iz + 256 * (iy + 256 * ix));

   Acc[0] = ExtAcc_InitialField[0][m_temp];
   Acc[1] = ExtAcc_InitialField[1][m_temp];
   Acc[2] = ExtAcc_InitialField[2][m_temp];
  
   //free( m );
   //if ((ix == iy) && (iy == iz) && (iz == 0)) 
   // Aux_Message(stderr, "At %lf %lf %lf acc %lf %lf %lf\n",
   //     x, y, z, Acc[0], Acc[1], Acc[2]);

#  endif

} // FUNCTION : CUPOT_ExternalAcc / CPU_ExternalAcc



#endif // #ifdef GRAVITY
