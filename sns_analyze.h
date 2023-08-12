#include "sns_header.h"
#include "TOV_lib.h"
#include "physics_star_lib.h"

/* useful macro */
#define MAX_STR_LEN (400)
#define LINE_STR    "-------------------------------------------------------------------------"


void sns_print_physical_system_properties(Physics_T *const phys,
                                          FILE *const file,
                                          const int iteration,
                                          const int pr_screen);

void sns_analyze(Physics_T *const phys,const int iteration);
static void compute_properties(Physics_T *const phys/* sns */);



