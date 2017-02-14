#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

#define DELTA 0.5
#define nullptr NULL

#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

double SIZE;
//
//  benchmarking program
//

/******************** Data Structures for Barnes-Hut Tree ********************/
/* ---------------------Class and function declarations--------------------- */
// Object Quad: Represent a quadrant in domain
class Quad {
public:
    double x, y; // x, y coordinate of the lower left point of the square
    double length; // length of one side of the square
// Constructor
    Quad(double, double, double);
// Member functions
    bool contains(double, double); // returns whether a point is in quadrant
};

// Object Body: Represent a body in BHTree
class Body{
public:
    int n_part; // Number of particles in body
    double px; // x position of cluster or particle
    double py; // y position of cluster or particle
    particle_t* p_particle; // pointer to particle (for external node)

    // Member function
    // Constructor
    Body(particle_t &); // Construct body from particle_t
    Body(int, double, double); // Construct clusters by providing the attributes
    bool in(Quad &); // Function to test if body is in a certain quad domain
};

// Object BHTree: Barnes-Hut tree structure
class BHTree{
private:
    Body* body; // body or aggregate body stored in this node
    Quad* quad; // square region that the tree represents
    BHTree* NW; // tree representing northwest quadrant
    BHTree* NE; // tree representing northeast quadrant
    BHTree* SW; // tree representing southwest quadrant
    BHTree* SE; // tree representing southeast quadrant
    void fork(); // fork into 4 child nodes
    void insertChild(Body &);

public:
    // Constructor
    BHTree(Quad *); // create a Barnes-Hut tree with no bodies, representing
                    // the given quadrant.
    ~BHTree(){
        delete body;
        delete quad;
        delete NW;
        delete NE;
        delete SW;
        delete SE;
    };
    void insert(Body &); // add the body to the involking Barnes-Hut tree
    void totalForce(particle_t*, double *, double *, int *); // apply force on
                                                             // particle from
                                                             // all bodies in
                                                             // the invoking
                                                             // Barnes-Hut tree
};

// Auxilliary Functions
Body* addBody(const Body &, const Body &);
BHTree* buildTree(int n, particle_t* particles);
/* ---------------------Function Implementations--------------------- */
/* Implementations for Quad */
Quad::Quad(double x_in, double y_in, double length_in){
    x = x_in;
    y = y_in;
    length = length_in;
};

bool Quad::contains(double x_q, double y_q){
    return (x_q >= x) && (x_q <= (x + length)) &&
        (y_q >= y) && (y_q <= y + length);
};

/* Implementations for Body */
Body::Body(particle_t &p){
    n_part = 1;
    px = p.x;
    py = p.y;
    p_particle = &p;
};

Body::Body(int n_part_in, double px_in, double py_in){
    n_part = n_part_in;
    px = px_in;
    py = py_in;
};

bool Body::in(Quad & q){
    return q.contains(px, py);
};

/* Implementation for BHTree */
BHTree::BHTree(Quad* q){
    quad = q;
    body = nullptr;
    NW = nullptr;
    NE = nullptr;
    SW = nullptr;
    SE = nullptr;
};

void BHTree::fork(){
    double hl = (quad->length)/2;
    Quad* q1 = new Quad(quad->x, (quad->y)+hl, hl);
    Quad* q2 = new Quad(quad->x + hl, quad->y + hl, hl);
    Quad* q3 = new Quad(quad->x, quad->y, hl);
    Quad* q4 = new Quad(quad->x + hl, quad->y, hl);
    NW = new BHTree(q1);
    NE = new BHTree(q2);
    SW = new BHTree(q3);
    SE = new BHTree(q4);
};

void BHTree::insertChild(Body &b){    // insert body into a child node
    Quad q1 = *(NW->quad);
    Quad q2 = *(NE->quad);
    Quad q3 = *(SW->quad);
    Quad q4 = *(SE->quad);
    if (b.in(q1))
    NW->insert(b);
    else if (b.in(q2))
    NE->insert(b);
    else if (b.in(q3))
    SW->insert(b);
    else if (b.in(q4))
    SE->insert(b);
    else
    perror("Could not locate a quadrant for body.");
};

void BHTree::insert(Body & b){
    if (body == nullptr){ // empty node
        body = &b;
    }
    else if (body->n_part > 1){ // internal node
        // update center of mass
        body->px = (body->px * body->n_part + b.px * b.n_part)
            / (body->n_part+b.n_part);
        body->py = (body->py * body->n_part + b.py * b.n_part)
            / (body->n_part+b.n_part);
        body->n_part += b.n_part;
        // insert body into a child node
        insertChild(b);
    }
    else if (body->n_part == 1){ // external node
        // create new combined body
        Body* new_body = addBody(*body, b);
        // fork this node
        fork();
        // insert body into a child node
        insertChild(*body);
        insertChild(b);
        body = new_body;
    }
    else {
        perror("Error while inserting body into BHTree.");
    }
};

void BHTree::totalForce(particle_t* ptc, double* dmin, double* davg, int* navg){
    if (body == nullptr){ // empty node
    }
    else if (body->n_part > 1){ // internal node
        // check distance
        double dx = body->px - ptc->x;
        double dy = body->py - ptc->y;
        double r = sqrt(dx * dx + dy * dy);
        if (r - 1.414*(quad->length) < cutoff){ // not entire cluster beyond
                                                // cutoff
            NW->totalForce(ptc, dmin, davg, navg);
            NE->totalForce(ptc, dmin, davg, navg);
            SW->totalForce(ptc, dmin, davg, navg);
            SE->totalForce(ptc, dmin, davg, navg);
        };
    }
    else if (body->n_part == 1){ // external node
        apply_force(*ptc, *(body->p_particle), dmin, davg, navg);
    }
};

/* Auxilliary Function Implementations */
Body* addBody(const Body & a, const Body & b){
    // Return a new Body that represents the center-of-mass
    // of the two bodies a and b.
    int n_part_new = a.n_part + b.n_part;
    double px_new = (a.n_part * a.px + b.n_part * b.px ) / n_part_new;
    double py_new = (a.n_part * a.py + b.n_part * b.py ) / n_part_new;
    Body * bp = new Body(n_part_new, px_new, py_new);
    return bp;
};

BHTree* buildTree(int n, particle_t* particles, double SIZE){
    // build root node
    Quad * quad_root = new Quad(0, 0, SIZE);
    BHTree* Tp = new BHTree(quad_root);
    // insert particles into the root
    for (int i = 0; i < n; i++){
        // pack particle into body
        Body* b = new Body(particles[i]);
        Tp->insert(*b);
    };
    return Tp;
};

/* Main Function */
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
    SIZE = sqrt(density * n);

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
        BHTree* tree_ptr = buildTree(n, particles, SIZE); // build tree
        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            tree_ptr->totalForce(&particles[i], &dmin, &davg, &navg);
        };
//        delete tree_ptr; // chop tree

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
    //  -The minimum distance absmin between 2 particles during the run of the
    //   simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of
    //   cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than
    //   0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting
    //   correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that\
             some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that\
             most particles are not interacting");
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
