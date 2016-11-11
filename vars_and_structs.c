/***************************************************************
               DEFINITIONS AND PREPROCESOR DIRECTIVES
***************************************************************/
#define X 0
#define Y 1
#define Z 2

#define INDEX(i,j,k) (k)+GV.NCELLS*((j)+GV.NCELLS*(i))
#define POW2(x) ((x)*(x))
#define POW3(x) ((x)*(x)*(x))
#define DIST2(x,y,z) ((x)*(x)+(y)*(y)+(z)*(z))


/***************************************************************
                        STRUCTURES
***************************************************************/

/*+++++ Global Variables+++++*/
struct GlobalVariables
{   
  char               FILENAME[1000]; //Path of the data file
  double             Grav; //Gravitational constant
  long int  NPARTS; //Number of particles of the sample
  long int  NP_TOT; //Total number of particles
  double             mass; //Mass of each particle
  int                NCELLS; //Number of cells in the grid
  long int           NTOTALCELLS; //Total Number of cells in the grid
  double             unit_t; //Time unit
  double             dt; //Time step
  double             Total_time; /*Total time of integration 
			  (in units of unit_t i.e. Total_time = 2 will be 2*unit_t)*/
  double             ZERO; // Zero for the computer  
  double             Tot_energy; //Total energy from analytical solution
  double             Sim_energy; //Total energy from simulation
  double             Poten;
  double             dist_poten;
  double             BoxSize[3];
  double             CellSize[3]; //size of each cell
}GV;//globalVariables


/*+++++ Particle +++++*/
struct Particle
{
  /*----- General variables for correction -----*/
  //unsigned long int     ID; //ID of the particle
  unsigned long int     GridID; //ID of the grid in which the particle is. If Part[i].GridID == Part[j].GridID => Direct sum. Else: P3M scheme
  double   pos[3]; //Positions 
  //double   vel[3]; //Velocities  
  double   acc[3]; //Accelerations
  double   dipole[3];
  double   quadrupole[3];
  //double    mass; //Mass  
}*part=NULL;


/*+++++ Grid +++++*/
struct grid
{
  int       Np_cell; // Number of particles in the Cell.                                                        
  long int* id_part; // Array with the ID of the particles inside the cell.
  //long int    ID; //ID of the grid cell
  //double  pos[3]; //Position in the center of the grid
  double    mass; //Total mass assigned into the grid point.
  double  pos_cm[3]; //Position of the CM of the cell
  double  mass_pos[3]; //Multiplication sum_i{m_i * x_i}
  double  mass_pos2[3]; //Multiplication sum_i{m_i * x_i^2}
  double  mass_posxy; //Multiplication sum_i{m_i * x_i*y_i}
  double  mass_posxz; //Multiplication sum_i{m_i * x_i*z_i}
  double  mass_posyz; //Multiplication sum_i{m_i * y_i*z_i}
}*gp;
