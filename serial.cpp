#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

//
//  benchmarking program
//

/************************* Data Structures for Barnes-Hut Tree *************************/
// Object Quad: Represent a quadrant in domain
class Quad {
// Constructor
    Quad(double, double, double);
private:
    double x, y; // x, y coordinate of the lower left point of the square
    double length; // length of one side of the square
public:
// Member functions
    void copy(Quad q);
    bool contains(double, double); // function that returns if a point is in quadrant
    Quad NW();// These four methods create and return a new Quad representing a sub-quadrant of the invoking quadrant.

    Quad NE();
    Quad SW();
    Quad SE();
};

Quad::Quad(double x_in, double y_in, double length_in){
    x = x_in;
    y = y_in;
    length = length_in;
};

void copy(Quad &q){
    x = q.x;
    y = q.y;
    length = q.length;
};

bool Quad::contains(double x_q, double y_q){
    return (x_q >= x) && (x_q <= (x + length)) && (y_q >= y) && (y_q <= y + length);
};

Quad Quad::NW(){
    return Quad(x, y+(length/2), length/2);
};

Quad Quad::NE(){
    return Quad(x+(length/2), y+(length/2), length/2);
};

Quad Quad::SW(){
    return Quad(x, y, length/2);
};

Quad Quad::SE(){
    return Quad(x+(length/2), y, length/2);
};

// Object Body: Represent a body in BHTree
class Body{
// Constructor
    Body(particle_t); // Construct body from particle_t
    Body(int, double, double); // Construct clusters by providing the attributes
private:
    int n_part = 0; // Number of particles in body
    double px; // x position of cluster or particle
    double py; // y position of cluster or particle
    particle_t* p_particle; // pointer to particle (for external node)

// Member function
public:
    bool in(Quad); // Function to test if body is in a certain quad domain
    Body add(Body a, Body b); // Return a new Body that represents the center-of-mass of the two bodies a and b.
};

Body::Body(particle_t &p){
    n_part = 1;
    px = pt.x;
    py = pt.y;
    p_particle = &p;
};

Body::(int n_part_in, double px_in, double py_in){
    n_part = n_part_in;
    px = px_in;
    py = py_in;
}

bool Body::in(Quad& q){
    return q.contains(px, py);
}

Body Body::add(Body a, Body, b){
    int n_part_new = a.n_part + b.n_part;
    double px_new = (a.n_part * a.px + b.n_part * b.px ) / n_part_new;
    double py_new = (a.n_part * a.py + b.n_part * b.py ) / n_part_new;
    return Body(n_part_new, px_new, py_new);
}

// Object BHTree: Barnes-Hut tree structure
class BHTree{
// Constructor
    BHTree(Quad); // create a Barnes-Hut tree with no bodies, representing the given quadrant.
private:
    Body body; // body or aggregate body stored in this node
    Quad quad; // square region that the tree represents
    BHTree* NW; // tree representing northwest quadrant
    BHTree* NE; // tree representing northeast quadrant
    BHTree* SW; // tree representing southwest quadrant
    BHTree* SE; // tree representing southeast quadrant
    
public:
    void insert(Body); // add the body to the involking Barnes-Hut tree
    void totalForce(particle_t*); // apply force on particle from all bodies in the invoking Barnes-Hut tree
};

BHTree::BHTree(Quad &q){
    quad.copy(q);
}




int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
	
    for( int step = 0; step < NSTEPS; step++ )
    {
	navg = 0;
        davg = 0.0;
	dmin = 1.0;
        //
        //  compute forces
        //
        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            for (int j = 0; j < n; j++ )
				apply_force( particles[i], particles[j],&dmin,&davg,&navg);
        }
 
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );		

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
