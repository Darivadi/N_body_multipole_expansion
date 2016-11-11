/****************************************************************************************************
NAME: init_conds
FUNCTION: Sets initial conditions in positions, velocities and accelerations
INPUT: None
RETURN: 0
****************************************************************************************************/
int init_conds( void )
{
  int i;
  double aux_mass, aux_mass2, aux_factor, aux_factor1;
  double MSun, Gmks, AU_mks;

  /*----- Units -----*/
  GV.Grav = 43.0071;//1.0; //Internal Units of (10AU)^3 MSun^-1 (Unit_t)^-2    
  GV.ZERO = 1.0e-300;
}//init_conds


/****************************************************************************************************
NAME: forces_c_calc
FUNCTION: 
INPUT: 
RETURN: 0
****************************************************************************************************/
int forces_calc_direct_sum( long int NPARTS )
{
  int i, j;
  double aux_mass, aux_factor, dist, aux_x, aux_y, aux_z, aux_pos[3];
  
  for(i=0; i<NPARTS; i++)
    {
      part[i].acc[X] = 0.0;
      part[i].acc[Y] = 0.0;
      part[i].acc[Z] = 0.0;

      for(j=0; j<NPARTS; j++)
	{
	  if(i != j)
	    {
	      /*----- Computing distance -----*/
	      aux_pos[X] = part[i].pos[X] - part[j].pos[X];
	      aux_pos[Y] = part[i].pos[Y] - part[j].pos[Y];  
	      aux_pos[Z] = part[i].pos[Z] - part[j].pos[Z];
	      
	      /*
	      if(aux_pos[X] < GV.ZERO || aux_pos[Y] < GV.ZERO || aux_pos[Z] < GV.ZERO)
		printf("Be carefull with parts i=%d, j=%d\n", i, j);
	      */

	      aux_factor = POW2(aux_pos[X]) + POW2(aux_pos[Y]) + POW2(aux_pos[Z]);
	      dist       = pow(aux_factor, 1.5);
	      
	      /*----- Acceleration in X -----*/
	      aux_factor       = aux_pos[X] / dist;
	      part[i].acc[X] += -1.0 * GV.mass * aux_factor;
	      
	      /*----- Acceleration in Y -----*/
	      aux_factor       = aux_pos[Y] / dist;
	      part[i].acc[Y] += -1.0 * GV.mass * aux_factor;
	      
	      /*----- Acceleration in Z -----*/
	      aux_factor       = aux_pos[Z] / dist;
	      part[i].acc[Z] += -1.0 * GV.mass * aux_factor;
	      
	    }//if
	}//for j
      
      part[i].acc[X] *= GV.Grav * GV.mass;
      part[i].acc[Y] *= GV.Grav * GV.mass;
      part[i].acc[Z] *= GV.Grav * GV.mass;
      
      if(i%2000==0)
	printf("i=%d a_x=%lf a_y=%lf, a_z=%lf\n", i, part[i].acc[X], part[i].acc[Y], part[i].acc[Z]);
      
      /*
      if(i%100==0)
	{
	  printf("Ready for i=%d\n", i);
	}
      */
    }//for i

  return 0;  
}//forces_c_calc


/****************************************************************************************************
NAME: locateCell
FUNCTION: locates the particle in the grid point
INPUT: 
RETURN: 0
****************************************************************************************************/
void locateCell(double xp, double yp, double zp, int indexPartArray, struct grid *gp)
{  
  //Numbers in order to find the cell
  //Counters in x,y and z
  int i,j,k,n;
  
  // x-axis
  i = floor( (xp / GV.BoxSize[X]) * GV.NCELLS );
  // Y-axis
  j = floor( (yp / GV.BoxSize[Y]) * GV.NCELLS );
  // Z-axis
  k = floor( (zp / GV.BoxSize[Z]) * GV.NCELLS );
  
  n = INDEX(i,j,k);

  //printf("For part=%d n=%d, i=%d, j=%d, k=%d\n", indexPartArray, n, i, j, k);
  //fflush(stdout);
  
  // Storing particles in cell
  gp[n].Np_cell = gp[n].Np_cell + 1;
  //printf("Particle added\n");
  gp[n].id_part = (long int*) realloc(gp[n].id_part, gp[n].Np_cell*sizeof(long int)); 
  gp[n].id_part[gp[n].Np_cell-1] = indexPartArray;
  
  gp[n].pos_cm[X] += part[indexPartArray].pos[X];
  gp[n].pos_cm[Y] += part[indexPartArray].pos[Y]; 
  gp[n].pos_cm[Z] += part[indexPartArray].pos[Z]; 

  //gp[n].ID = n; //ID of the cell
  part[indexPartArray].GridID = n; //ID of the cell in which the particle is
  

  //printf("Ready for part=%d in cell n=%d, i=%d, j=%d, k=%d\n", indexPartArray, n, i, j, k);
  //fflush(stdout);

}//locateCells


/****************************************************************************************************
NAME: W
FUNCTION: Window function for NGP MAS
INPUT: 
RETURN: 0
****************************************************************************************************/
double W(double x, double y, double z, double Hx, double Hy, double Hz)
{
  double Wx, Wy, Wz;
  //if( fabs(x) <= H*0.5 )
  if( fabs(x) <= Hx*0.5 )
    {
      Wx = 1.0;
    }
  /*else if( fabs(x) == H*0.5  )
    { 
    Wx = 0.5;
    }*/
  else
    {
      Wx = 0.0;
    }


  //if( fabs(y) < H*0.5 )
  if( fabs(y) <= Hy*0.5 )
    {
      Wy = 1.0;
    }  
  /*else if( fabs(y) == H*0.5  )
    { 
      Wy = 0.5;
      }*/
  else
    {
      Wy = 0.0;
    }


  //if( fabs(z) < H*0.5 )
  if( fabs(z) <= Hz*0.5 )
    {
    Wz = 1.0;
    }
  /*else if( fabs(z) == H*0.5  )
    { 
      Wz = 0.5;
      }*/
  else
    {
      Wz = 0.0;
    }

}//W




/****************************************************************************************************
NAME: two_part_force
FUNCTION: computes the force between 
INPUT: 
RETURN: 0
****************************************************************************************************/
int two_part_force(int i, int j)
{
  double aux_mass, aux_factor, dist, aux_x, aux_y, aux_z, aux_pos[3];
  
  /*----- Computing distance -----*/
  aux_pos[X] = part[i].pos[X] - part[j].pos[X];
  aux_pos[Y] = part[i].pos[Y] - part[j].pos[Y];  
  aux_pos[Z] = part[i].pos[Z] - part[j].pos[Z];
  
  aux_factor = POW2(aux_pos[X]) + POW2(aux_pos[Y]) + POW2(aux_pos[Z]);
  dist       = pow(aux_factor, 1.5);
  
  /*----- Acceleration in X -----*/
  aux_factor       = aux_pos[X] / dist;
  part[i].acc[X] += -1.0 * GV.mass * aux_factor;
  //part[i].acc[X] += -1.0 * GV.Grav * GV.mass * GV.mass * aux_factor;
  
  /*----- Acceleration in Y -----*/
  aux_factor       = aux_pos[Y] / dist;
  part[i].acc[Y] += -1.0 * GV.mass * aux_factor;
  
  /*----- Acceleration in Z -----*/
  aux_factor       = aux_pos[Z] / dist;
  part[i].acc[Z] += -1.0 * GV.mass * aux_factor;
  
  return 0;
}//two_part_force



/****************************************************************************************************
NAME: pc_force_monopole
FUNCTION: computes the force between the particle and a cell
INPUT: 
RETURN: 0
****************************************************************************************************/
int pc_force_monopole(int m, int i, int j, int k)
{
  double aux_mass, aux_factor, dist, aux_x, aux_y, aux_z, aux_pos[3];
  double xc, yc, zc;
  int n;
  
  /*----- Computing positions of the cell according to index-----*/
  n = INDEX(i,j,k);

  /*----- Computing distances -----*/  
  aux_pos[X] = gp[n].pos_cm[X] - part[m].pos[X];
  aux_pos[Y] = gp[n].pos_cm[Y] - part[m].pos[Y];
  aux_pos[Z] = gp[n].pos_cm[Z] - part[m].pos[Z];

  /*
  aux_pos[X] = part[m].pos[X] - xc;
  aux_pos[Y] = part[m].pos[Y] - yc;
  aux_pos[Z] = part[m].pos[Z] - zc;
  */
  aux_factor = POW2(aux_pos[X]) + POW2(aux_pos[Y]) + POW2(aux_pos[Z]);
  dist       = pow(aux_factor, 1.5);
  
  /*----- Acceleration in X -----*/
  aux_factor     = aux_pos[X] / dist;  
  //part[i].acc[X] += -1.0 * GV.Grav * GV.mass * gp[n].mass * aux_factor;
  part[m].acc[X] += -1.0 * gp[n].mass * aux_factor;
  
  /*----- Acceleration in Y -----*/
  aux_factor     = aux_pos[Y] / dist;
  part[m].acc[Y] += -1.0 * gp[n].mass * aux_factor;
  
  /*----- Acceleration in Z -----*/
  aux_factor     = aux_pos[Z] / dist;
  part[m].acc[Z] += -1.0 * gp[n].mass * aux_factor;
  
  return 0;
}//pc_force_monopole



/****************************************************************************************************
NAME: pc_force_mono_dip
FUNCTION: computes the force between the particle and a cell
INPUT: 
RETURN: 0
****************************************************************************************************/
int pc_force_mono_dip(int m, int i, int j, int k)
{
  double aux_mass, aux_factor, dist, aux_x, aux_y, aux_z, aux_pos[3];
  double xc, yc, zc;
  double dist5, aux_factor1;
  int n;
  
  /*----- Computing positions of the cell according to index-----*/
  n = INDEX(i,j,k);

  /*----- Computing distances -----*/  
  aux_pos[X] = gp[n].pos_cm[X] - part[m].pos[X];
  aux_pos[Y] = gp[n].pos_cm[Y] - part[m].pos[Y];
  aux_pos[Z] = gp[n].pos_cm[Z] - part[m].pos[Z];

  /*
  aux_pos[X] = part[l].pos[X] - xc;
  aux_pos[Y] = part[l].pos[Y] - yc;
  aux_pos[Z] = part[l].pos[Z] - zc;
  */
  aux_factor  = POW2(aux_pos[X]) + POW2(aux_pos[Y]) + POW2(aux_pos[Z]);
  dist        = pow(aux_factor, 1.5);
  aux_factor1 = pow(aux_factor, 2.5);
  aux_factor1 = 1.0 / aux_factor1;

  /*----- Acceleration in X -----*/   
  /*..... Monopole .....*/
  aux_factor     = aux_pos[X] / dist;
  part[m].acc[X] += -1.0 * gp[n].mass * aux_factor;  
  
  /*..... Dipole .....*/
  aux_factor     = 1.0 / dist;
  aux_factor     += -3.0 * POW2(aux_pos[X]) * aux_factor1;
  part[m].acc[X] += aux_factor * gp[n].mass_pos[X];

  aux_factor     = -3.0 * aux_pos[X] * aux_pos[Y] * aux_factor1;
  part[m].acc[X] += aux_factor * gp[n].mass_pos[Y];

  aux_factor     = -3.0 * aux_pos[X] * aux_pos[Z] * aux_factor1;
  part[m].acc[X] += aux_factor * gp[n].mass_pos[Z];
 
  
  /*----- Acceleration in Y -----*/
  /*..... Monopole .....*/
  aux_factor     = aux_pos[Y] / dist;
  part[m].acc[Y] += -1.0 * gp[n].mass * aux_factor;

  /*..... Dipole .....*/
  aux_factor     = 1.0 / dist;
  aux_factor     += -3.0 * POW2(aux_pos[Y]) * aux_factor1;
  part[m].acc[Y] += aux_factor * gp[n].mass_pos[Y];

  aux_factor     = -3.0 * aux_pos[X] * aux_pos[Y] * aux_factor1;
  part[m].acc[Y] += aux_factor * gp[n].mass_pos[X];

  aux_factor     = -3.0 * aux_pos[Y] * aux_pos[Z] * aux_factor1;
  part[m].acc[Y] += aux_factor * gp[n].mass_pos[Z];
 
  
  /*----- Acceleration in Z -----*/
  /*..... Monopole .....*/
  aux_factor     = aux_pos[Z] / dist;
  part[m].acc[Z] += -1.0 * gp[n].mass * aux_factor;

  /*..... Dipole .....*/
  aux_factor     = 1.0 / dist;
  aux_factor     += -3.0 * POW2(aux_pos[Z]) * aux_factor1;
  part[m].acc[Z] += aux_factor * gp[n].mass_pos[Z];

  aux_factor     = -3.0 * aux_pos[X] * aux_pos[Z] * aux_factor1;
  part[m].acc[Z] += aux_factor * gp[n].mass_pos[X];

  aux_factor     = -3.0 * aux_pos[Y] * aux_pos[Z] * aux_factor1;
  part[m].acc[Z] += aux_factor * gp[n].mass_pos[Z];
 
  return 0;
}//pc_force_mono_dip



/****************************************************************************************************
NAME: pc_force_mono_dip_quad
FUNCTION: computes the force between the particle and a cell
INPUT: 
RETURN: 0
****************************************************************************************************/
int pc_force_mono_dip_quad(int m, int i, int j, int k)
{
  double aux_mass, aux_factor, dist, aux_x, aux_y, aux_z, aux_pos[3];
  double xc, yc, zc;
  double dist5, dist7, aux_factor1, aux_factor2;
  int n;
  
  /*----- Computing positions of the cell according to index-----*/
  n = INDEX(i,j,k);

  /*----- Computing distances -----*/  
  aux_pos[X] = gp[n].pos_cm[X] - part[m].pos[X];
  aux_pos[Y] = gp[n].pos_cm[Y] - part[m].pos[Y];
  aux_pos[Z] = gp[n].pos_cm[Z] - part[m].pos[Z];

  /*
  aux_pos[X] = part[l].pos[X] - xc;
  aux_pos[Y] = part[l].pos[Y] - yc;
  aux_pos[Z] = part[l].pos[Z] - zc;
  */
  aux_factor  = POW2(aux_pos[X]) + POW2(aux_pos[Y]) + POW2(aux_pos[Z]);
  dist        = pow(aux_factor, 1.5);
  dist5       = pow(aux_factor, 2.5);
  dist5       = 1.0 / dist5;
  dist7       = pow(aux_factor, 3.5);
  dist7       = 1.0 / dist7;

  /*----- Acceleration in X -----*/   
  /*..... Monopole .....*/
  aux_factor     = aux_pos[X] / dist;
  part[m].acc[X] += -1.0 * gp[n].mass * aux_factor;  
  
  /*..... Dipole .....*/
  aux_factor     = 1.0 / dist;
  aux_factor     += -3.0 * POW2(aux_pos[X]) * dist5;
  part[m].acc[X] += aux_factor * gp[n].mass_pos[X];

  aux_factor     = -3.0 * aux_pos[X] * aux_pos[Y] * dist5;
  part[m].acc[X] += aux_factor * gp[n].mass_pos[Y];

  aux_factor     = -3.0 * aux_pos[X] * aux_pos[Z] * dist5;
  part[m].acc[X] += aux_factor * gp[n].mass_pos[Z];

  /*..... Quadrupole .....*/
  aux_factor     = 15.0 * POW3( aux_pos[X] ) * dist7;
  aux_factor     += -9.0 * aux_pos[X] * dist5;
  part[m].acc[X] += -0.5 * aux_factor * gp[n].mass_pos2[X]; 

  aux_factor = 15.0 * aux_pos[X] * POW2(aux_pos[Y]) * dist7;
  aux_factor += -3.0 * aux_pos[X] * dist5;
  part[m].acc[X] += -0.5 * aux_factor * gp[n].mass_pos2[Y];

  aux_factor = 15.0 * aux_pos[X] * POW2(aux_pos[Z]) * dist7;
  aux_factor += -3.0 * aux_pos[X] * dist5;
  part[m].acc[X] += -0.5 * aux_factor * gp[n].mass_pos2[Z];

  aux_factor = 15.0 * POW2( aux_pos[X] ) * aux_pos[Y] * dist7;
  aux_factor += -3.0 * aux_pos[X] * dist5;
  part[m].acc[X] += -1.0 * aux_factor * gp[n].mass_posxy;
  
  part[m].acc[X] += -15.0 * aux_pos[X] * aux_pos[Y] * aux_pos[Z] * dist7 * gp[n].mass_posyz;

  aux_factor = 15.0 * POW2( aux_pos[X] ) * aux_pos[Z] * dist7;
  aux_factor += -3.0 * aux_pos[Z] * dist5;
  part[m].acc[X] += -1.0 * aux_factor * gp[n].mass_posxz;
  
  
  /*----- Acceleration in Y -----*/
  /*..... Monopole .....*/
  aux_factor     = aux_pos[Y] / dist;
  part[m].acc[Y] += -1.0 * gp[n].mass * aux_factor;

  /*..... Dipole .....*/
  aux_factor     = 1.0 / dist;
  aux_factor     += -3.0 * POW2(aux_pos[Y]) * dist5;
  part[m].acc[Y] += aux_factor * gp[n].mass_pos[Y];

  aux_factor     = -3.0 * aux_pos[X] * aux_pos[Y] * dist5;
  part[m].acc[Y] += aux_factor * gp[n].mass_pos[X];

  aux_factor     = -3.0 * aux_pos[Y] * aux_pos[Z] * dist5;
  part[m].acc[Y] += aux_factor * gp[n].mass_pos[Z];
 
  /*..... Quadrupole .....*/  
  aux_factor = 15.0 * POW2(aux_pos[X]) * aux_pos[Y] * dist7;
  aux_factor += -3.0 * aux_pos[Y] * dist5;
  part[m].acc[Y] += -0.5 * aux_factor * gp[n].mass_pos2[X];

  aux_factor     = 15.0 * POW3( aux_pos[Y] ) * dist7;
  aux_factor     += -9.0 * aux_pos[Y] * dist5;
  part[m].acc[Y] += -0.5 * aux_factor * gp[n].mass_pos2[Y]; 

  aux_factor = 15.0 * aux_pos[X] * POW2(aux_pos[Z]) * dist7;
  aux_factor += -3.0 * aux_pos[Y] * dist5;
  part[m].acc[Y] += -0.5 * aux_factor * gp[n].mass_pos2[Z];

  aux_factor = 15.0 * aux_pos[X] * POW2(aux_pos[Y]) * dist7;
  aux_factor += -3.0 * aux_pos[X] * dist5;
  part[m].acc[Y] += -1.0 * aux_factor * gp[n].mass_posxy;

  aux_factor = 15.0 * POW2( aux_pos[Y] ) * aux_pos[Z] * dist7;
  aux_factor += -3.0 * aux_pos[Z] * dist5;
  part[m].acc[Y] += -1.0 * aux_factor * gp[n].mass_posyz;

  part[m].acc[Y] += -15.0 * aux_pos[X] * aux_pos[Y] * aux_pos[Z] * dist7 * gp[n].mass_posxz;


  /*----- Acceleration in Z -----*/
  /*..... Monopole .....*/
  aux_factor     = aux_pos[Z] / dist;
  part[m].acc[Z] += -1.0 * gp[n].mass * aux_factor;

  /*..... Dipole .....*/
  aux_factor     = 1.0 / dist;
  aux_factor     += -3.0 * POW2(aux_pos[Z]) * dist5;
  part[m].acc[Z] += aux_factor * gp[n].mass_pos[Z];

  aux_factor     = -3.0 * aux_pos[X] * aux_pos[Z] * dist5;
  part[m].acc[Z] += aux_factor * gp[n].mass_pos[X];

  aux_factor     = -3.0 * aux_pos[Y] * aux_pos[Z] * dist5;
  part[m].acc[Z] += aux_factor * gp[n].mass_pos[Z];


  /*..... Quadrupole .....*/
  aux_factor = 15.0 * POW2(aux_pos[X]) * aux_pos[Z] * dist7;
  aux_factor += -3.0 * aux_pos[Z] * dist5;
  part[m].acc[Z] += -0.5 * aux_factor * gp[n].mass_pos2[X];

  aux_factor = 15.0 * POW2(aux_pos[Y]) * aux_pos[Z] * dist7;
  aux_factor += -3.0 * aux_pos[Z] * dist5;
  part[m].acc[Z] += -0.5 * aux_factor * gp[n].mass_pos2[Y];

  aux_factor     = 15.0 * POW3( aux_pos[Z] ) * dist7;
  aux_factor     += -9.0 * aux_pos[Z] * dist5;
  part[m].acc[Z] += -0.5 * aux_factor * gp[n].mass_pos2[Z]; 
  
  part[m].acc[Z] += -15.0 * aux_pos[X] * aux_pos[Y] * aux_pos[Z] * dist7 * gp[n].mass_posxy;

  aux_factor = 15.0 * aux_pos[Y] * POW2(aux_pos[Z]) * dist7;
  aux_factor += -3.0 * aux_pos[Y] * dist5;
  part[m].acc[Z] += -1.0 * aux_factor * gp[n].mass_posyz;
  
  aux_factor = 15.0 * aux_pos[X] * POW2(aux_pos[Z]) * dist7;
  aux_factor += -3.0 * aux_pos[X] * dist5;
  part[m].acc[Z] += -1.0 * aux_factor * gp[n].mass_posxz;
 
  return 0;
}//pc_force_mono_dip_quad





/****************************************************************************************************
NAME: p3m_forces
FUNCTION: 
INPUT: 
RETURN: 0
****************************************************************************************************/
int p3m_forces( int tag )
{
  int i, j, k, m, n, p;
  int aux_index;


  if(tag == 0) //Monopole
    {
      /*+++++  First iteration over all particles+++++*/
      for(m=0; m<GV.NPARTS; m++)
	{
	  part[m].acc[X] = 0.0;
	  part[m].acc[Y] = 0.0;
	  part[m].acc[Z] = 0.0;
	  
	  /*----- For each particle, I must iterate over all cells -----*/
	  for(i=0; i<GV.NCELLS; i++)
	    {
	      for(j=0; j<GV.NCELLS; j++)
		{
		  for(k=0; k<GV.NCELLS; k++)
		    {
		      
		      /*----- For each cell I must know its index -----*/
		      n = INDEX(i,j,k);
		      
		      /*----- For each cell, I need to know if the particle is inside it or not -----*/
		      if(n != part[m].GridID)
			{
			  /*---- If the particle is outside this cell, I must compute the force between the 
			    particle[m] and all the mass inside the cell[n] -----*/			  			  
			  pc_force_monopole(m, i, j, k);		      								      
			}//if n
		      else
			{		      		      		      
			  /*----- If the particle is inside this cell, I need to iterate over the 
			    particles inside this cell in order to compute a direct-sum-force -----*/
			  for(p=0; p<gp[n].Np_cell; p++)
			    {
			      aux_index = gp[n].id_part[p];
			      
			      if(aux_index != m)
				{
				  /*----- Computing forces between the 2 particles: m and aux_index -----*/
				  two_part_force(m, aux_index);
				}//if aux_index			  
			    }//for p (particles inside the cell)		      		      
			  
			}//else n
		      
		    }//for k (grid-z)
		}//for j (grid-y)
	    }//for i (grid-x)
	  
	  part[m].acc[X] *= GV.Grav * GV.mass;
	  part[m].acc[Y] *= GV.Grav * GV.mass;
	  part[m].acc[Z] *= GV.Grav * GV.mass;

	  if(m%2000==0)
	    printf("i=%d a_x=%lf a_y=%lf, a_z=%lf\n", m, part[m].acc[X], part[m].acc[Y], part[m].acc[Z]);
	  
	  
	  
	}//for m (particles)            
    }//Monopole
  else if(tag == 1) //Dipole
    {			  
      /*+++++  First iteration over all particles+++++*/
      for(m=0; m<GV.NPARTS; m++)
	{
	  part[m].acc[X] = 0.0;
	  part[m].acc[Y] = 0.0;
	  part[m].acc[Z] = 0.0;
	  
	  /*----- For each particle, I must iterate over all cells -----*/
	  for(i=0; i<GV.NCELLS; i++)
	    {
	      for(j=0; j<GV.NCELLS; j++)
		{
		  for(k=0; k<GV.NCELLS; k++)
		    {
		      
		      /*----- For each cell I must know its index -----*/
		      n = INDEX(i,j,k);
		      
		      /*----- For each cell, I need to know if the particle is inside it or not -----*/
		      if(n != part[m].GridID)
			{
			  /*---- If the particle is outside this cell, I must compute the force between the 
			    particle[m] and all the mass inside the cell[n] -----*/		    
			  pc_force_mono_dip(m, i, j, k);		      
			}//if n
		      else
			{		      		      		      
			  /*----- If the particle is inside this cell, I need to iterate over the 
			    particles inside this cell in order to compute a direct-sum-force -----*/
			  for(p=0; p<gp[n].Np_cell; p++)
			    {
			      aux_index = gp[n].id_part[p];
			      
			      if(aux_index != m)
				{
				  /*----- Computing forces between the 2 particles: m and aux_index -----*/
				  two_part_force(m, aux_index);
				}//if aux_index			  
			    }//for p (particles inside the cell)		      		      
			  
			}//else n
		      
		    }//for k (grid-z)
		}//for j (grid-y)
	    }//for i (grid-x)
	  
	  part[m].acc[X] *= GV.Grav * GV.mass;
	  part[m].acc[Y] *= GV.Grav * GV.mass;
	  part[m].acc[Z] *= GV.Grav * GV.mass;
	  

	  if(m%2000==0)
	    printf("i=%d a_x=%lf a_y=%lf, a_z=%lf\n", m, part[m].acc[X], part[m].acc[Y], part[m].acc[Z]);
	  	  
	  
	}//for m (particles)                        
    }//Dipole
  else if(tag == 2) //Quadrupole
    {
      /*+++++  First iteration over all particles+++++*/
      for(m=0; m<GV.NPARTS; m++)
	{
	  part[m].acc[X] = 0.0;
	  part[m].acc[Y] = 0.0;
	  part[m].acc[Z] = 0.0;
	  
	  /*----- For each particle, I must iterate over all cells -----*/
	  for(i=0; i<GV.NCELLS; i++)
	    {
	      for(j=0; j<GV.NCELLS; j++)
		{
		  for(k=0; k<GV.NCELLS; k++)
		    {
		      
		      /*----- For each cell I must know its index -----*/
		      n = INDEX(i,j,k);
		      
		      /*----- For each cell, I need to know if the particle is inside it or not -----*/
		      if(n != part[m].GridID)
			{
			  /*---- If the particle is outside this cell, I must compute the force between the 
			    particle[m] and all the mass inside the cell[n] -----*/
			  pc_force_mono_dip_quad(m, i, j, k);		      
			}//if n
		      else
			{		      		      		      
			  /*----- If the particle is inside this cell, I need to iterate over the 
			    particles inside this cell in order to compute a direct-sum-force -----*/
			  for(p=0; p<gp[n].Np_cell; p++)
			    {
			      aux_index = gp[n].id_part[p];
			      
			      if(aux_index != m)
				{
				  /*----- Computing forces between the 2 particles: m and aux_index -----*/
				  two_part_force(m, aux_index);
				}//if aux_index			  
			    }//for p (particles inside the cell)		      		      
			  
			}//else n
		      
		    }//for k (grid-z)
		}//for j (grid-y)
	    }//for i (grid-x)
	  
	  part[m].acc[X] *= GV.Grav * GV.mass;
	  part[m].acc[Y] *= GV.Grav * GV.mass;
	  part[m].acc[Z] *= GV.Grav * GV.mass;

	  if(m%2000==0)
	    printf("i=%d a_x=%lf a_y=%lf, a_z=%lf\n", m, part[m].acc[X], part[m].acc[Y], part[m].acc[Z]);	  	 	  
	  
	}//for m (particles)
    }//Quadrupole
  
  return 0;
}//p3m_forces

