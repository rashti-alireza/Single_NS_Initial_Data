/*
// Alireza Rashti
// August 2023
*/

/* analyzing initial data such as mass, momentum, constraints etc.  */

#include "sns_analyze.h"

/* analyzing physics properties, constraints etc */
void sns_analyze(Physics_T *const phys,const int iteration)
{
  if (!phys) return;

  FUNC_TIC
  
  const char *const properties_file_name = P_"properties.txt";
  FILE *file = 0;
  char str[MAX_STR_LEN];
   
  /* compute properties and constraints */ 
  physics(phys,ADM_COMPUTE_CONSTRAINTS);
  
  /* compute various properties */
  compute_properties(phys);
  
  /* open properties file in "my_directory" and save */
  sprintf(str,"%s/%s",Pgets(P_"my_directory"),properties_file_name);
  file = Fopen(str,"w");
  sns_print_physical_system_properties(phys,file,iteration,0);
  Fclose(file);

  /* open properties file in "Diagnostics" and save */
  sprintf(str,"%s/%s",Pgets(P_"Diagnostics"),properties_file_name);
  file = Fopen(str,"a");
  sns_print_physical_system_properties(phys,file,iteration,0);
  Fclose(file);
  
  /* prints */
  print_fields_0D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  print_fields_1D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  print_fields_2D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  print_fields_3D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  
  FUNC_TOC
}

/* print physical system properties such as mass, spin etc in the given
// file, if pr_screen is 1, it also prints in stdout */
void sns_print_physical_system_properties(Physics_T *const phys,
                                          FILE *const file,
                                          const int iteration,
                                          const int pr_screen)
{
  Physics_T *const ns = init_physics(phys,NS);

  if (pr_screen)
  {
    printf(Pretty0"iteration = %d:\n",iteration);
  }
  fprintf(file,"%s\n",LINE_STR);
  fprintf(file,"# iteration = %d\n",iteration);
  fprintf(file,"\n");
  
  star_print_properties(ns,Pgets(P_"NS_properties"),file,pr_screen);
  sys_print_properties(phys,Pgets(P_"SNS_properties"),file,pr_screen);
  
  free_physics(ns);
}

/* compute variety of properties.
// NOTE: order of parameter calculations matter. 
// NOTE: if there is a confusion between target params and current
//       params, "current" suffix added to the latter. */
static void compute_properties(Physics_T *const phys/* sns */)
{
  Physics_T *const ns = init_physics(phys,NS);
  TOV_T *tov          = 0;
  
  const double x_CM = Pgetd(P_"x_CM");
  const double y_CM = Pgetd(P_"y_CM");
  const double z_CM = Pgetd(P_"z_CM");
  double p[3]  = {0.};
  double j[3]  = {0.};
  double s[3]  = {0.};
  double cm[3] = {0.};
  double m     = 0.;
  
  /* NS: */
  observe(ns,"ADM(M)",Pgets("NS_Observe_ADM_M"),&m);
  Psetd("NS_ADM_mass",m);
  
  observe(ns,"Komar(M)",Pgets("NS_Observe_Komar_M"),&m);
  Psetd("NS_Komar_mass",m);
  
  observe(ns,"Baryonic(M)",Pgets("NS_Observe_baryonic_M"),&m);
  Psetd("NS_baryonic_mass_current",m);
  
  tov = TOV_init();
  tov->exit_if_error = 0;
  tov->phys  = ns;
  tov->bar_m = Pgetd("NS_baryonic_mass_current");
  tov = TOV_solution(tov);
  if (tov->status == 0)
  {
    Psetd("NS_TOV_ADM_mass",tov->ADM_m);
    /* Note: compactness = adm_mass/Schwarzschild_radius 
      (not isotropic radius) */
    Psetd("NS_TOV_compactness",tov->ADM_m/tov->r[tov->N-1]);
    Psetd("NS_TOV_radius",tov->rbar[tov->N-1]);
  }
  TOV_free(tov);
  
  Psetd("NS_mass_shedding_indicator",star_NS_mass_shedding_indicator(ns));
  
  observe(ns,"CM",Pgets("NS_Observe_CM"),cm);
  Psetd("NS_x_CM",cm[0]+x_CM);
  Psetd("NS_y_CM",cm[1]+y_CM);
  Psetd("NS_z_CM",cm[2]+z_CM);

  observe(ns,"ADM(P)",Pgets("NS_Observe_ADM_P"),p);
  Psetd("NS_Px_ADM",p[0]);
  Psetd("NS_Py_ADM",p[1]);
  Psetd("NS_Pz_ADM",p[2]);

  observe(ns,"ADM(J)",Pgets("NS_Observe_ADM_J"),j);
  Psetd("NS_Jx_ADM",j[0]);
  Psetd("NS_Jy_ADM",j[1]);
  Psetd("NS_Jz_ADM",j[2]);
  
  observe(ns,"spin",Pgets("NS_Observe_spin"),s);
  Psetd("NS_Spin_x",s[0]);
  Psetd("NS_Spin_y",s[1]);
  Psetd("NS_Spin_z",s[2]);
  
  /* SNS: */
  observe(phys,"ADM(M)",Pgets(P_"Observe_ADM_M"),&m);
  Psetd(P_"adm_mass",m);
  /* for NS chi we need adm_mass. Note: s is calculated a few lines above. */
  Psetd("NS_chi_x",s[0]/Pow2(m));
  Psetd("NS_chi_y",s[1]/Pow2(m));
  Psetd("NS_chi_z",s[2]/Pow2(m));
  
  observe(phys,"Komar(M)",Pgets(P_"Observe_Komar_M"),&m);
  Psetd(P_"Komar_mass",m);
  
  observe(phys,"ADM(P)",Pgets(P_"Observe_ADM_P"),p);
  Psetd(P_"Px_ADM",p[0]);
  Psetd(P_"Py_ADM",p[1]);
  Psetd(P_"Pz_ADM",p[2]);

  observe(phys,"ADM(J)",Pgets(P_"Observe_ADM_J"),j);
  Psetd(P_"Jx_ADM",j[0]);
  Psetd(P_"Jy_ADM",j[1]);
  Psetd(P_"Jz_ADM",j[2]);

  /* virial error */
  double v_e = fabs(1.-Pgetd(P_"adm_mass")/Pgetd(P_"komar_mass"));
  Psetd(P_"virial_error",v_e);
  
  free_physics(ns);
}
