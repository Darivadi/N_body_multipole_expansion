/***************************************************************************************************
NAME: conf2dump
FUNCTION: Reads the input file with parameters
INPUT: Parameters file
RETURN: 0
***************************************************************************************************/

int conf2dump( char filename[] )
{
  char cmd[1000];
  int nread;
  
  sprintf( cmd, "grep -v \"#\" %s | grep -v \"^$\" | gawk -F\"=\" '{print $2}' > %s.dump", 
	     filename, filename );
  /*
    sprintf( cmd, "grep -v \"#\" %s | grep -v \"^$\" | awk -F\"=\" '{print $2}' > %s.dump", 
    filename, filename );
  */
  nread = system( cmd );
  
  return 0;
}


/***************************************************************************************************
NAME: read_parameters
FUNCTION: Reads the parameters
INPUT: Parameters file
RETURN: 0
***************************************************************************************************/
int read_parameters( char filename[] )
{
  int nread;
  char cmd[1000], filenamedump[1000];
  FILE *file;
  
  /*+++++ Loading the file +++++*/
  file = fopen( filename, "r" );
  if( file==NULL )
    {
      printf( "  * The file '%s' doesn't exist!\n", filename );
      return 1;
    }
  fclose(file);
  
  /*+++++ Converting to plain text +++++*/
  conf2dump( filename );
  sprintf( filenamedump, "%s.dump", filename );
  file = fopen( filenamedump, "r" );


  /*+++++ Reading parameters +++++*/
  nread = fscanf(file, "%lf", &GV.mass); //Mass of particles

  nread = fscanf(file, "%s", GV.FILENAME); //File with data

  nread = fscanf(file, "%ld", &GV.NP_TOT); //Total number of particles

  nread = fscanf(file, "%ld", &GV.NPARTS); //Number of particles in the sample

  nread = fscanf(file, "%d", &GV.NCELLS); //Number of cells in the grid

  //nread = fscanf(file, "%lf", &GV.BoxSize); //Number of cells in the grid

  fclose( file );

  GV.NTOTALCELLS = POW3(GV.NCELLS);
  
  printf( "  * The file '%s' has been loaded!\n", filename );
  //printf("with integration time %lf steps\n", GV.Total_time);
  printf("with %ld particles, a sample of %ld and %d cells per axis with mass of %lf\n", 
	 GV.NP_TOT, GV.NPARTS, GV.NCELLS, GV.mass);
  
  
  sprintf( cmd, "rm -rf %s.dump", filename );
  nread = system( cmd );
  
  return 0;
}


/***************************************************************************************************
NAME: read_parts
FUNCTION: Reads the particles positions
INPUT: Parameters file
RETURN: 0
***************************************************************************************************/
int read_parts()
{
  int count, n, nread;
  FILE *inFile=NULL;
  double max_x = 0.0, max_y = 0.0, max_z = 0.0;
  double min_x = 1e10, min_y = 1e10, min_z = 1e10;
  time_t t;
  int *rand_nums;
  double dummy_pos[3], aux_num;

  printf("Ensayando... %ld, %ld\n", 
	     GV.NPARTS, GV.NP_TOT);
  printf("Reading %s data file\n", GV.FILENAME);
  inFile = fopen(GV.FILENAME, "r");
  printf("It passed!\n");
  fflush(stdout);

  /*+++++ Generating random numbers to choose the GV.NPARTS FROM THE GV.NP_TOT +++++*/  
  printf("Let's comparing the number of particles\n");
  
  if( GV.NPARTS < GV.NP_TOT ) 
    {
      printf("Generating %lu random indices between 0 and %lu\n", 
	     GV.NPARTS, GV.NP_TOT);      
      
      rand_nums = (int *) calloc( (size_t) GV.NPARTS, sizeof(int) );

      printf("Memory allocated\n");
      
      srand((unsigned) time(&t));
      
      
      for(n=0; n<GV.NPARTS; n++)
	{

	  if(n==0)
	    {
	      do
		{
		  rand_nums[n] = (int) rand()%(GV.NP_TOT);
		}
	      while(rand_nums[n] > GV.NP_TOT/GV.NPARTS);
	    }//if
	  else
	    {
	      do
		{
		  rand_nums[n] = (int) rand()%(GV.NP_TOT);
		}
	      //while(rand_nums[n] <= rand_nums[n-1]*GV.NPARTS/GV.NP_TOT);
 	      while( rand_nums[n] <= rand_nums[n-1] || rand_nums[n] > n*GV.NP_TOT/GV.NPARTS);
	    }//else	  	  
	  
	  if(n%2000==0)
	    printf("Generated random number #%d with result %d\n", n, rand_nums[n]);
	  //fflush(stdout);
	}//For n
      
      printf("Random numbers generated\n");

      /*
      printf("Reordering random data in ascending order\n");
      
      for(n=0; n<GV.NPARTS; n++)
	{
	  for(count=n+1; count<GV.NPARTS; count++)
	    {
	      if(rand_nums[n] > rand_nums[count])
		{
		  aux_num = rand_nums[n];
		  rand_nums[n] = rand_nums[count];
		  rand_nums[count] = aux_num;
		}//if
	    }//for count
	}//for n
     
      printf("Random data sorted\n");            
      
      for(n=0; n<GV.NPARTS; n++)
	{
	  printf("Generated random number #%d with result %d\n", n, rand_nums[n]);
	}      
      */
    }//if
  else if(GV.NPARTS == GV.NP_TOT)
    {
      printf("Number of sample particles %ld is the same than the total particles %ld\n", 
	     GV.NPARTS, GV.NP_TOT);
    }//else
  else
    {
      printf("Number of sample particles must be less than total particles\n");
      exit(0);
    }//else
  
   
  /*+++++ Reading data +++++*/ 
  if(GV.NPARTS < GV.NP_TOT)
    {
      count = 0;
      printf("Reading data for nparts < np_tot\n");
      
      for(n=0; n<GV.NP_TOT; n++ )
	{             
	  
	  if( n == rand_nums[count] )
	    {
	      //printf("Entered n=%d, count=%d, rand=%d\n", n, count, rand_nums[count]);
	      
	      nread = fscanf(inFile, "%lf %lf %lf", 
			     &part[count].pos[X], &part[count].pos[Y], &part[count].pos[Z]);

	      /*----- Looking for the maximum -----*/
	        if(part[count].pos[X] > max_x)
		    max_x = part[count].pos[X];
		  
		  if(part[count].pos[Y] > max_y)
		    max_y = part[count].pos[Y];
		  
		  if(part[count].pos[Z] > max_z)
		    max_z = part[count].pos[Z];
		    
		  /*----- Looking for the minimum -----*/
		  if(part[count].pos[X] < min_x)
		    min_x = part[count].pos[X];
		  
		  if(part[count].pos[Y] < min_y)
		    min_y = part[count].pos[Y];
		  
		  if(part[count].pos[Z] < min_z)
		    min_z = part[count].pos[Z];
		  
		  /*
		  printf("Reading count=%d x=%lf y=%lf z=%lf\n", 
			 count, 
			 part[count].pos[X], part[count].pos[Y], part[count].pos[Z]);
		  */		  
		  count++;
	    }//if
	  else
	    {
	      //printf("Not entered n=%d, count=%d, rand=%d\n", n, count, rand_nums[count]);
	      nread = fscanf(inFile, "%lf %lf %lf", 
			     &dummy_pos[X], &dummy_pos[Y], &dummy_pos[Z]);
	    }//else	 	  
	}//for n      
      
        printf("counter must be Nparts: %d\n", count);
      
    }//if nparts z nptot   
  else
    {
      for(n=0; n<GV.NPARTS; n++ )
	{             
	  nread = fscanf(inFile, "%lf %lf %lf", 
			 &part[n].pos[X], &part[n].pos[Y], &part[n].pos[Z]);
	  
	  /*----- Looking for the maximum -----*/
	  if(part[n].pos[X] > max_x)
	    max_x = part[n].pos[X];
	  
	  if(part[n].pos[Y] > max_y)
	    max_y = part[n].pos[Y];
	  
	  if(part[n].pos[Z] > max_z)
	    max_z = part[n].pos[Z];
	  
	  /*----- Looking for the minimum -----*/
	  if(part[n].pos[X] < min_x)
	    min_x = part[n].pos[X];
	  
	  if(part[n].pos[Y] < min_y)
	    min_y = part[n].pos[Y];
	  
	  if(part[n].pos[Z] < min_z)
	  min_z = part[n].pos[Z];
	}//for n
    }//else
  
  

  fclose(inFile);
  
  printf("File read\n");

  
  printf("Maximum values: x=%lf, y=%lf, z=%lf\n", max_x, max_y, max_z);
  printf("Minimum values: x=%lf, y=%lf, z=%lf\n", min_x, min_y, min_z);
   
  GV.BoxSize[X] = ceil( (fabs(min_x) + fabs(max_x))*10.0)/10.0;
  GV.BoxSize[Y] = ceil( (fabs(min_y) + fabs(max_y))*10.0)/10.0;
  GV.BoxSize[Z] = ceil( (fabs(min_z) + fabs(max_z))*10.0)/10.0;
  
  //GV.BoxSize[X] = GV.BoxSize[Y] = GV.BoxSize[Z] = 450;

  printf("Box_x = %lf, Box_y = %lf, Box_z = %lf\n", 
	 GV.BoxSize[X], GV.BoxSize[Y], GV.BoxSize[Z]);

  printf("Number of particles %ld\n", GV.NPARTS);


  for(n=0; n<GV.NPARTS; n++)
    {
      part[n].pos[X] += fabs(min_x) + 1e-8; //+ 0.5 * GV.BoxSize[X];
      part[n].pos[Y] += fabs(min_y) + 1e-8; //+ 0.5 * GV.BoxSize[Y];
      part[n].pos[Z] += fabs(min_z) + 1e-8; //+ 0.5 * GV.BoxSize[Z];

      /*
      printf("n=%d x=%lf, y=%lf, z=%lf\n", n,  part[n].pos[X], part[n].pos[Y], part[n].pos[Z]);
      fflush(stdout);
           
      if( part[n].pos[X] < 0.0 || part[n].pos[Y] < 0.0 || part[n].pos[Z] < 0.0 )
	printf("Be careful with ID=%d\n", n);
      
      if( part[n].pos[X] > GV.BoxSize[X] || part[n].pos[Y] > GV.BoxSize[Y] || part[n].pos[Z] > GV.BoxSize[X] )
	printf("Be careful with ID=%d\n", n);
      */
    }//for n
  
  
  return 0;
}//read_parts
