#include "GAMER.h"
#include "TestProb.h"

#ifdef GRAVITY

// =====================================================================================================
#include "CUPOT.h"

double *ExtAcc_InitialField[3] = {NULL, NULL, NULL}; 
#endif

// ====================================================================================================

// problem-specific global variables
// =======================================================================================
static double SphCol_Dens_Bg;       // background mass density
static double SphCol_Dens_Delta;    // top-hat mass density --> total density = Dens_Bg*( 1 + Dens_Delta )
static double SphCol_Engy_Bg;       // background energy density
static double SphCol_Radius;        // top-hat radius
static double SphCol_Center[3];     // top-hat center
static double DrivAmp;		    // amplitude for driving term
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

//#  ifndef COMOVING
//   Aux_Error( ERROR_INFO, "COMOVING must be enabled !!\n" );
//#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

   if ( MPI_Rank == 0 )
   {
#     ifdef GRAVITY
      if ( OPT__GRAVITY_TYPE != GRAVITY_SELF )
         Aux_Message( stderr, "WARNING : OPT__GRAVITY_TYPE != GRAVITY_SELF ??\n" );
#     endif

#     ifndef DUAL_ENERGY
         Aux_Message( stderr, "WARNING : it's recommended to enable DUAL_ENERGY for this test\n" );
#     endif

      for (int f=0; f<6; f++)
      if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
         Aux_Message( stderr, "WARNING : non-periodic BC for fluid ??\n" );

#     ifdef GRAVITY
      if ( OPT__BC_POT != BC_POT_PERIODIC )
         Aux_Message( stderr, "WARNING : non-periodic BC for gravity ??\n" );
#     endif

      if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
         Aux_Message( stderr, "WARNING : simulation domain is not cubic ??\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "SphCol_Dens_Bg",    &SphCol_Dens_Bg,        -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "SphCol_Dens_Delta", &SphCol_Dens_Delta,     -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "SphCol_Engy_Bg",    &SphCol_Engy_Bg,        -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "SphCol_Radius",     &SphCol_Radius,         -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "SphCol_Center_X",   &SphCol_Center[0],      -1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "SphCol_Center_Y",   &SphCol_Center[1],      -1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "SphCol_Center_Z",   &SphCol_Center[2],      -1.0,          NoMin_double,     NoMax_double      );

   ReadPara->Add( "DrivAmp",   	       &DrivAmp,      		1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Read( FileName );

   delete ReadPara;

// set the default values
   for (int d=0; d<3; d++)
      if ( SphCol_Center[d] < 0.0 )    SphCol_Center[d] = 0.5*amr->BoxSize[d];


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 5.0e-2;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID           = %d\n",     TESTPROB_ID       );
      Aux_Message( stdout, "  background mass density   = %13.7e\n", SphCol_Dens_Bg    );
      Aux_Message( stdout, "  top-hat over-density      = %13.7e\n", SphCol_Dens_Delta );
      Aux_Message( stdout, "  background energy density = %13.7e\n", SphCol_Engy_Bg    );
      Aux_Message( stdout, "  top-hat radius            = %13.7e\n", SphCol_Radius     );
      Aux_Message( stdout, "  top-hat center x          = %13.7e\n", SphCol_Center[0]  );
      Aux_Message( stdout, "  ...            y          = %13.7e\n", SphCol_Center[1]  );
      Aux_Message( stdout, "  ...            z          = %13.7e\n", SphCol_Center[2]  );
      Aux_Message( stdout, "  driv amplitude            = %13.7e\n", DrivAmp	       );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   const double r = sqrt( SQR(x-SphCol_Center[0]) + SQR(y-SphCol_Center[1]) + SQR(z-SphCol_Center[2]) );

   fluid[DENS] = SphCol_Dens_Bg;
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;
   fluid[ENGY] = SphCol_Engy_Bg;

   if ( r <= SphCol_Radius )  fluid[DENS] *= ( 1.0 + SphCol_Dens_Delta );

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == HYDRO )



#ifdef GRAVITY
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExternalAcc
// Description :  Set the array "ExtAcc_AuxArray" used by the external acceration routines
//                "CUPOT_ExternalAcc.cu / CPU_ExternalAcc.cpp"
//
// Note        :  0. ***TO BE FINISHED***
//                1. Linked to the function pointer "Init_ExternalAcc_Ptr"
//                2. Enabled by the runtime option "OPT__GRAVITY_TYPE == 2/3"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_ExternalAcc()
{

  int dim = NX0[0];

  size_t nread = 0;
  FILE *bin_data = fopen("UM_IC", "r");

  float *density = (float *) calloc(dim*dim*dim, sizeof(float));
  nread = fread(density, sizeof(float), dim*dim*dim, bin_data);

  float *field_x = (float *) calloc(dim*dim*dim, sizeof(float));
  nread = fread(field_x, sizeof(float), dim*dim*dim, bin_data);

  float *field_y = (float *) calloc(dim*dim*dim, sizeof(float));
  nread = fread(field_y, sizeof(float), dim*dim*dim, bin_data);

  float *field_z = (float *) calloc(dim*dim*dim, sizeof(float));
  nread = fread(field_z, sizeof(float), dim*dim*dim, bin_data);

  if (ExtAcc_InitialField[0] == NULL) { 
    ExtAcc_InitialField[0] = (double*)calloc(dim*dim*dim, sizeof(double));
  }
   
  if (ExtAcc_InitialField[1] == NULL) { 
    ExtAcc_InitialField[1] = (double*)calloc(dim*dim*dim, sizeof(double));
  }

  if (ExtAcc_InitialField[2] == NULL) { 
    ExtAcc_InitialField[2] = (double*)calloc(dim*dim*dim, sizeof(double));
  }

  for( int i = 0; i < dim*dim*dim; i++){

    ExtAcc_InitialField[0][i] = DrivAmp * field_x[i] / density[i];
    ExtAcc_InitialField[1][i] = DrivAmp * field_y[i] / density[i];
    ExtAcc_InitialField[2][i] = DrivAmp * field_z[i] / density[i];
  }

  free(density);
  free(field_x);
  free(field_y);
  free(field_z);

  
  //free(ExtAcc_InitialField[0]);
  //free(ExtAcc_InitialField[1]);
  //free(ExtAcc_InitialField[2]);
  //free(ExtAcc_InitialField);
} // FUNCTION : Init_ExternalAcc



#endif // #ifdef GRAVITY



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_SGS
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_SGS()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   assert(0);

// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
   Output_User_Ptr          = NULL;
   Flag_User_Ptr            = NULL;
   Mis_GetTimeStep_User_Ptr = NULL;
   Aux_Record_User_Ptr      = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = NULL;
   End_User_Ptr             = NULL;
#  ifdef GRAVITY
   Init_ExternalAcc_Ptr     = Init_ExternalAcc;
   Init_ExternalPot_Ptr     = NULL;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_SGS
