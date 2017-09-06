#include "header.h"

params 		state;
pade		pade_data;
map 		conf;
map 		alt_map;

int main( int argc, char* argv[]) {
  load_parameters(argc, argv);		// read configuration from file	 //
  init_memory();			// initialize memory unit 	 //
  init_timemarching();			// initialize timemarching 	 //
  init_pade();				// initialize rational approx	 //
  allocate_memory();			// allocate all memory		 //

  unsigned int format_flag = 0;
  if (strcmp(state.txt_format,"ascii") == 0 )  		format_flag = 1;
  if (strcmp(state.txt_format,"binary") == 0 ) 		format_flag = 2;
  if (strcmp(state.txt_format,"pade") == 0 )   		format_flag = 3;
  if (strcmp(state.txt_format,"none") == 0 )   		format_flag = 4;
  if (strcmp(state.txt_format,"jon_wilkening") == 0 )   format_flag = 5;
 
  switch (format_flag) {

    case 1:
      printf("Reading ASCII file\n");
      load_ascii();
      set_mapping();
      break;

    case 2:
      printf("Reading binary save\n");
      printf("Placeholder\n");
      exit(1);
      break;

    case 3:
      printf("Reading Pade data\n");
      set_mapping();
      load_pade();
      break;

    case 4:
      printf("Starting new simulation\n");
      set_mapping();
      set_initial_data();
      convertQtoZ(data[0], tmpc[5]);  
      restore_potential(data[0], data[1], tmpc[2]);
      break;

    case 5: 
      printf("Setting data from Jon Wilkening file format.\n");
      set_mapping();
      set_initial_JW();
      convertQtoZ(data[0], tmpc[5]);  
      restore_potential(data[0], data[1], tmpc[2]);
      break;

    default:
      printf("Unknown text format\n");
      exit(1);
  }
  //complex_array_out("zt-original.txt", data[0]);
  //convertZtoQ(data[0], data[0]);
  //complex_array_out("q-original.txt", data[0]);
  //long double c = 0.5L;
  //for (long int j = 0; j < state.number_modes; j++) {
    //data[1][j] = 1.0IL*c*(1.L - data[0][j]*data[0][j]);
  //} 
  //complex_array_out("v-original.txt", data[1]);
  //complex_array_out("pre.Phi.ph.txt", tmpc[2]); 

  //new_map.scaling 	= 0.025L;
  //new_map.image_offset 	= 0.0L;
  //remap(&new_map, 2048); 
  //convertQtoZ(data[0], tmpc[5]);  
  //complex_array_out("zt-recovered-2.txt", tmpc[5]);
  //restore_potential(data[0], data[1], tmpc[3]);  
  print_constants();
  //allocate_pade(3);
  //exit(1);
  //rk6_step(data[0], data[1], 0.05);
  //convertQtoZ(data[0], tmpc[5]);  
  //complex_array_out("zt-after-step.txt", tmpc[5]);
  //restore_potential(data[0], data[1], tmpc[3]);  
  //print_constants();
  //convertQtoZ(data[0], tmpc[5]);  
  //complex_array_out("zt-recovered.txt", tmpc[5]);
  evolve_rk6(data[0], data[1]);

}
