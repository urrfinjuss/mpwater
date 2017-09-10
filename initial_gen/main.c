#include "header.h"

params 		state;
pade		pade_data;
ic		gen_ic;
map 		conf;
map 		alt_map;
mpfr_t 		Pie;

int main( int argc, char* argv[]) {
  load_parameters(argc, argv);		// read configuration from file	 //
  init_memory();			// initialize memory unit 	 //
  allocate_memory();			// allocate all memory		 //
  
  unsigned int format_flag = 0;
  if (strcmp(state.txt_format,"generate_ic") == 0 ) format_flag = 1;
  
  switch (format_flag) {

    case 1:
      printf("Starting new simulation\n");
      set_mapping();
      set_initial_data();
      break;

    default:
      printf("Unknown text format\n");
      exit(1);
  }
  convertQtoZ(data[0], tmpc[5]);
  complex_array_out("z.txt", tmpc[5]);     
 
}
