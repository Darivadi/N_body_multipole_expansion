/***************************************************************
              INCLUDING LIBRARIES AND MODULES
***************************************************************/

/*+++++ Libraries +++++*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/*+++++ Module files +++++*/
#include "vars_and_structs.c"
#include "read_write.c"
#include "functions.c"


int main( int argc, char *argv[] )
{
  char *infile=NULL;
  FILE *outFile=NULL, *outFile1=NULL, *outFile2=NULL; 
  int i, j, k, l, nread, index, Np, idPart, n_res, tag;
  double time_step=0.0, aux_x, aux_y, aux_z, aux_factor;
  double xp, yp, zp, xc, yc, zc, Window_fn, aux_pos[3], aux_acc[3];
  int resolutions[3];
  char buffer[100], buffer1[100];
  
  //For time measurement 
  clock_t start, end;
  double time_used;
  

  if(argc < 2)
    {
      printf("Error: Incomplete number of parameters. Execute as follows:\n");
      printf("%s Parameters_file\n", argv[0]);
      exit(0);
    }//if                                                                                                      
  
  /*+++++ Reading parameters file +++++*/
  infile = argv[1]; // Parameters file 
  printf("Reading parameters file\n");
  read_parameters( infile );
  printf("Reading parameters file: finished!\n");
  printf("-------------------------------------------------------------------------------------\n");


  /*+++++ Allocating memory for particle's structure +++++*/
  printf("Allocating memory for particle structure\n");
  part = (struct Particle *) calloc( (size_t) GV.NPARTS, sizeof(struct Particle) ); 


  /*+++++ Reading particles +++++*/
  printf("Reading particle's data\n");
  read_parts();
  init_conds();

  //exit(0);
  
  /*+++++ Computing forces by direct sum +++++*/
  /*----- Initializing time -----*/
  start = clock();

  forces_calc_direct_sum( GV.NPARTS );

  end = clock();
  time_used = (double) (1.0*end - 1.0*start )/CLOCKS_PER_SEC;
  printf("Computing forces by direct sum for %ld particles took %10.5lf seconds to execute\n", 
	 GV.NPARTS, time_used);
  printf("-------------------------------------------------------------------------------------\n");
  
  snprintf(buffer, sizeof(char)*100, "./../../Data/Execution_time_NParts%ld_NGrid%d.txt", GV.NPARTS, GV.NCELLS);
  outFile = fopen(buffer, "w");  
  fprintf(outFile, "Computing forces by direct sum for %lu particles took %10.5lf seconds to execute\n", 
	  GV.NPARTS, time_used);
  

  /*----- Saving data file with forces by direct sum -----*/  
  snprintf(buffer1, sizeof(char)*100, "./../../Data/Forces_direct_sum_NParts%ld_NGrid%d.bin", GV.NPARTS, GV.NCELLS);
  //snprintf(buffer1, sizeof(char)*100, "./../../Data/Forces_direct_sum_NParts%ld_NGrid%d.txt", GV.NPARTS, GV.NCELLS);
  outFile1 = fopen(buffer1, "w"); 
  
  
  fwrite(&(GV.BoxSize),  sizeof(double),            1, outFile1);  //Box Size                                           
  fwrite(&(GV.NPARTS),   sizeof(long int), 1, outFile1);  //Number of particles     
  fwrite(&(GV.NCELLS),   sizeof(int),               1, outFile1);  //Grid resolution
  
  
  for(i=0; i<GV.NPARTS; i++)
    {
      
      aux_pos[X] = part[i].pos[X];
      aux_pos[Y] = part[i].pos[Y];
      aux_pos[Z] = part[i].pos[Z];
      fwrite(&(aux_pos[X]), sizeof(double), 3, outFile1);
      
      aux_acc[X] = part[i].acc[X];
      aux_acc[Y] = part[i].acc[Y];
      aux_acc[Z] = part[i].acc[Z];
      fwrite(&(aux_acc[X]), sizeof(double), 3, outFile1);
      
      /*
      fprintf(outFile1, "%16.8lf %16.8lf %16.8lf\n", 
	      part[i].acc[X], part[i].acc[Y], part[i].acc[Z]);
      */
    }//for i

  fclose(outFile1);


  /*+++++ Computing forces: 
    Tag = 0: Monopole approximation 
    Tag = 1: Monopole + dipole 
    Tag = 2: Monopole + dipole + quadrupole +++++*/  
  for(tag=0; tag<3; tag++)
    {
      
      if(tag==0)
	printf("Forces using monopole approximation\n");
      else if (tag==1)
	printf("Forces using monopole + dipole\n");
      else if (tag==2)
	printf("Forces using monopole + dipole + quadrupole\n");
      
      
      start = clock();
      
      printf("P3M code: Allocating memory\n");
      gp = (struct grid *) calloc((size_t) GV.NTOTALCELLS, sizeof(struct grid));
      
      GV.CellSize[X] = GV.BoxSize[X] / (1.0*GV.NCELLS);
      GV.CellSize[Y] = GV.BoxSize[Y] / (1.0*GV.NCELLS);
      GV.CellSize[Z] = GV.BoxSize[Z] / (1.0*GV.NCELLS);
      
      printf("CellSize_x = %lf, CellSize_y = %lf, CellSize_z = %lf\n", 
	     GV.CellSize[X], GV.CellSize[Y], GV.CellSize[Z]);
      
      for(i=0; i<GV.NTOTALCELLS; i++)
	{
	  gp[i].Np_cell   = 0;  
	  gp[i].mass      = 0.0; 
	  gp[i].pos_cm[X] = 0.0;
	  gp[i].pos_cm[Y] = 0.0;
	  gp[i].pos_cm[Z] = 0.0;	  
	  gp[i].mass_pos[X] = 0.0;
	  gp[i].mass_pos[Y] = 0.0;
	  gp[i].mass_pos[Z] = 0.0;
	  gp[i].mass_pos2[X] = 0.0; 
	  gp[i].mass_pos2[Y] = 0.0; 
	  gp[i].mass_pos2[Z] = 0.0; 
	  gp[i].mass_posxy = 0.0; 
	  gp[i].mass_posxz = 0.0; 
	  gp[i].mass_posyz = 0.0;
	}//for i
      
      /*----- Locating cells -----*/
      printf("Locating particles in the grid\n");
      for(i=0; i<GV.NPARTS; i++)
	{      
	  locateCell(part[i].pos[X], part[i].pos[Y], part[i].pos[Z], i, gp);           	
	}//for i
      printf("Particles located in the grid\n");
      
      
      /*----- Applying MAS -----*/  
      printf("Starting NGP MAS\n");
      
      for(i=0; i<GV.NCELLS; i++)
	{
	  for(j=0; j<GV.NCELLS; j++)
	    {
	      for(k=0; k<GV.NCELLS; k++)
		{	      
		  /* Index of the cell  */
		  index = INDEX(i,j,k); // C-order                                                                        	      	     
		  
		  /* Number of particles in the cell */
		  Np = gp[index].Np_cell;
		  
		  for(l=0; l<Np; l++)
		    {
		      /* Particle ID */
		      idPart = gp[index].id_part[l];
		      
		      /* Coordinates of the particle  */
		      xp = part[idPart].pos[X];
		      yp = part[idPart].pos[Y];
		      zp = part[idPart].pos[Z];
		  
		      //--- Mass and momentum assignment to neighbour cells (CIC) ---                                       
		      //indexaux = INDEX(i, j, k );
		      xc = GV.CellSize[X]*(0.5 + i);
		      yc = GV.CellSize[Y]*(0.5 + j);
		      zc = GV.CellSize[Z]*(0.5 + k);
		      
		      //--- Mass with CIC assignment scheme ---                                                       
		      Window_fn = W(xc-xp, yc-yp, zc-zp, GV.CellSize[X], GV.CellSize[Y], GV.CellSize[Z]);
		      gp[index].mass += GV.mass * Window_fn;

		      /*
			if(tag == 1 || tag == 2)
			{
		      */
		      
		      gp[index].pos_cm[X] /= gp[index].Np_cell;
		      gp[index].pos_cm[Y] /= gp[index].Np_cell;
		      gp[index].pos_cm[Z] /= gp[index].Np_cell;
		      
		      gp[index].mass_pos[X] += (gp[index].pos_cm[X] - part[idPart].pos[X]);
		      gp[index].mass_pos[Y] += (gp[index].pos_cm[Y] - part[idPart].pos[Y]);
		      gp[index].mass_pos[Z] += (gp[index].pos_cm[Z] - part[idPart].pos[Z]);
		      
		      gp[index].mass_pos2[X] += POW2( (gp[index].pos_cm[X] - part[idPart].pos[X]) ); 
		      gp[index].mass_pos2[Y] += POW2( (gp[index].pos_cm[Y] - part[idPart].pos[Y]) ); 
		      gp[index].mass_pos2[Z] += POW2( (gp[index].pos_cm[Z] - part[idPart].pos[Z]) ); 
		      
		      gp[index].mass_posxy += (gp[index].pos_cm[X] - part[idPart].pos[X]) * (gp[index].pos_cm[Y] - part[idPart].pos[Y]); 
		      gp[index].mass_posxz += (gp[index].pos_cm[X] - part[idPart].pos[X]) * (gp[index].pos_cm[Z] - part[idPart].pos[Z]); 
		      gp[index].mass_posyz += (gp[index].pos_cm[Y] - part[idPart].pos[Y]) * (gp[index].pos_cm[Z] - part[idPart].pos[Z]); 
		      //}//if
		      
		    }//for l 
		  
		}//for k                                                                                                  
	    }//for j                                                                                                    
	}//for i 

      
      /*if(tag == 1 || tag == 2)
	{*/
      for(i=0; i<GV.NTOTALCELLS; i++)
	{
	  gp[i].mass_pos[X]  *= GV.mass;
	  gp[i].mass_pos[Y]  *= GV.mass;
	  gp[i].mass_pos[Z]  *= GV.mass;	  
	  gp[i].mass_pos2[X] *= GV.mass;
	  gp[i].mass_pos2[Y] *= GV.mass;
	  gp[i].mass_pos2[Z] *= GV.mass;
	  gp[i].mass_posxy   *= GV.mass;
	  gp[i].mass_posxz   *= GV.mass;
	  gp[i].mass_posyz   *= GV.mass;
	}//for i
      //}//if
      
      
      printf("NGP MAS finished. Let's begin with the computing of forces\n");
      printf("-------------------------------------------------------------------------------------\n");
      
      /*----- Computing forces -----*/
      p3m_forces(tag);  
      
      end = clock();
      time_used = (double) (1.0*end - 1.0*start )/CLOCKS_PER_SEC;
      
      printf("Computing forces with P3M scheme between %ld particles and a grid with %d ^3 cells and size (%lf,%lf,%lf) units took %10.5lf seconds to execute. Tag = %d\n",
	     GV.NPARTS, GV.NCELLS, GV.BoxSize[X], GV.BoxSize[Y], GV.BoxSize[Z], time_used, tag);
      printf("-------------------------------------------------------------------------------------\n");
      
      //fprintf(outFile, "Computing forces by P3M scheme took %10.5lf seconds to execute\n", time_used);    
      fprintf(outFile, 
	      "Computing forces with P3M scheme between %ld particles and a grid with %d ^3 cells and size (%lf,%lf,%lf) units took %10.5lf seconds to execute. Tag = %d\n"
	      , GV.NPARTS, GV.NCELLS, GV.BoxSize[X], GV.BoxSize[Y], GV.BoxSize[Z], time_used, tag);
      
      
      /*----- Saving data file with forces in the approximations -----*/
      snprintf(buffer1, sizeof(char)*100, "./../../Data/Forces_P3M_NParts%ld_NGrid%d_tag%d.bin", GV.NPARTS, GV.NCELLS, tag);
      //snprintf(buffer1, sizeof(char)*100, "./../../Data/Forces_P3M_NParts%ld_NGrid%d.txt", GV.NPARTS, GV.NCELLS);
      outFile1 = fopen(buffer1, "w"); 
      
      
      fwrite(&(GV.BoxSize),  sizeof(double),            1, outFile1);  //Box Size                                           
      fwrite(&(GV.NPARTS),   sizeof(long int), 1, outFile1);  //Number of particles     
      fwrite(&(GV.NCELLS),   sizeof(int),               1, outFile1);  //Grid resolution
      
      //printf("Size of unsigned long int %zu bytes\n", sizeof(unsigned long int));
      
      for(i=0; i<GV.NPARTS; i++)
	{
	  
	  aux_pos[X] = part[i].pos[X];
	  aux_pos[Y] = part[i].pos[Y];
	  aux_pos[Z] = part[i].pos[Z];
	  fwrite(&(aux_pos[X]), sizeof(double), 3, outFile1);
	  
	  aux_acc[X] = part[i].acc[X];
	  aux_acc[Y] = part[i].acc[Y];
	  aux_acc[Z] = part[i].acc[Z];
	  fwrite(&(aux_acc[X]), sizeof(double), 3, outFile1);
	  
	  /*
	    fprintf(outFile1, "%16.8lf %16.8lf %16.8lf\n", 
	    part[i].acc[X], part[i].acc[Y], part[i].acc[Z]);
	  */
	}//for i
      
      fclose(outFile1);
      printf("File writen for Nparts = %ld, NCells=%d an tag = %d\n", 
	     GV.NPARTS, GV.NCELLS, tag);
      printf("Freeing up memory\n");
      
      free(gp);
      
    }//for tag
      
  //}//for n_res

  fclose(outFile);

  printf("Code finished successfully!\n");
  printf("See ya! Mr. Barbarian\n");
  printf("-------------------------------------------------------------------------------------\n");  
  
  return 0;
}//main
