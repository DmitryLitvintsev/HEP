//==============================================================================
//
//                        Logon file (important to pick up shlib's ! )
//
//==============================================================================
{
#include <iomanip.h>
#include <time.h>

  printf("\n");
  printf("    *=======================================*\n");
  printf("    *                                       *\n");
  printf("    *             Hello Dmitri!             *\n");
  printf("    *                                       *\n");
  printf("    *       Setting up for stripping.       *\n");
  printf("    *                                       *\n");
  printf("    *=======================================*\n");
  printf("\n");

  // Look for include files
  gSystem->SetIncludePath("-I/cdf/atom/home/litvinse/Q5/code/include -I$PROJECT_DIR/include");

  // Load in ROOT physics vectors and event generator libraries
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  
  // Load all shared libraries (USESHLIBS=1)
  gSystem->Load("/cdf/atom/home/litvinse/Q5/code/shlib/$BFARCH/libStntuple_base.so");
  gSystem->Load("/cdf/atom/home/litvinse/Q5/code/shlib/$BFARCH/libStntuple_obj.so");
  gSystem->Load("/cdf/atom/home/litvinse/Q5/code/shlib/$BFARCH/libStntuple_loop.so");
  gSystem->Load("/cdf/atom/home/litvinse/Q5/code/shlib/$BFARCH/libBottomMods_tools.so");
  gSystem->Load("/cdf/atom/home/litvinse/Q5/code/shlib/$BFARCH/libBottomMods_stn.so");  
  gSystem->Load("/cdf/atom/home/litvinse/Q5/code/shlib/$BFARCH/libBottomMods_ana.so");
  gSystem->Load("/cdf/atom/home/litvinse/Q5/code/shlib/$BFARCH/libAna.so");

 END:;
}
//==============================================================================
