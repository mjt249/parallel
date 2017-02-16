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

#define CHILDREN 4
#define NO -1
#define NW 0
#define NE 1
#define SW 2
#define SE 3

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
    BHTree* children[CHILDREN];
    BHTree* parent;
    int child_index;
    void branch(); // branch into 4 child nodes
    void insertChild(Body &);
    void emancipateChild(int);
    void updateCenterOfMass();

public:
    // Constructor
    BHTree(Quad *, BHTree*, int); // create a Barnes-Hut tree with no bodies,
                             // representing the given quadrant.
    ~BHTree(){
        delete body;
        delete quad;
        for (int i = 0; i < CHILDREN; i++) {
            delete children[i];
        }
    };
    void insert(Body &, int); // add the body to the involking Barnes-Hut tree
    void totalForce(particle_t*, double *, double *, int *); // apply force on
                                                             // particle from
                                                             // all bodies in
                                                             // the invoking
                                                             // Barnes-Hut tree
    void moveBodies();
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
    bool ret = q.contains(px, py);
    // if (ret)
    //     printf("(%f, %f) is between (%f, %f) and (%f,%f)\n", px, py, q.x, q.y, q.x + q.length, q.y + q.length);
    return ret;
};

/* Implementation for BHTree */
BHTree::BHTree(Quad* q, BHTree* rent, int child_order=NO){
    quad = q;
    body = nullptr;
    for (int i = 0; i < CHILDREN; i++) {
        children[i] = nullptr;
    }
    if (rent == nullptr)
        parent = this;
    else
        parent = rent;
    child_index = child_order;
};

void BHTree::branch(){
    printf("\tBranching %p\n", this);
    double hl = (quad->length)/2;
    Quad* q1 = new Quad(quad->x, (quad->y)+hl, hl);
    Quad* q2 = new Quad(quad->x + hl, quad->y + hl, hl);
    Quad* q3 = new Quad(quad->x, quad->y, hl);
    Quad* q4 = new Quad(quad->x + hl, quad->y, hl);
    children[NW] = new BHTree(q1, this, NW);
    children[NE] = new BHTree(q2, this, NE);
    children[SW] = new BHTree(q3, this, SW);
    children[SE] = new BHTree(q4, this, SE);
    // printf("Finished branching\n");
};

// returns whether this subtree has changed
void BHTree::moveBodies(){
    if (body == nullptr) {
    } else if (body->n_part > 1) { //internal node
        for (int i = 0; i < CHILDREN; i++) {
            if (children[i] != nullptr) {
                children[i]->moveBodies();
                children[i]->updateCenterOfMass();
                if ((children[i]->body != nullptr) &&
                        !(children[i]->body->in(*(children[i]->quad)))) {
                    // we're not in the right place anymore
                    emancipateChild(child_index);
                    insertChild(*body);
                }
            }
        }
        // update center of mass
    } else if (body->n_part == 1) {
        // external node
        move(*(body->p_particle));
    } else { perror("moveBodies error"); }
}

void BHTree::insertChild(Body &b){    // insert body into a child node
    for (int i = 0; i < CHILDREN; i++) {
        Quad q = *(children[i]->quad);
        if (b.in(q)) {
            if (children[i]->body != nullptr) {
                printf("\tchild node %p (index %d of %p) occupied\n",
                        children[i], i, this);
            }
            children[i]->insert(b, i);
            // printf("\tinserted particle %p into %p (child %d of %p)\n",
            //         b->p_particle, children[i], i, this);
            return;
        }
    }
    // couldn't find a quadrant
    if (parent != this) {
        printf("\t\tre-inserting child %p into %p\n",this, parent);
        parent->emancipateChild(child_index);
        parent->insert(b, NO);
    }
    else
        perror("Could not locate a quadrant for body.");
};

void BHTree::emancipateChild(int index) {
    Body *b = children[index]->body;
    // update center of mass
    body->px = (body->px * body->n_part - b->px * b->n_part)
        / (body->n_part-b->n_part);
    body->py = (body->py * body->n_part - b->py * b->n_part)
        / (body->n_part-b->n_part);
    body->n_part -= b->n_part;
    // detach child
    children[index] = nullptr;

    if (body->n_part == 0) {
        // now empty
        parent->emancipateChild(child_index);
        delete this;
    }
}

void BHTree::insert(Body & b, int ci){
    if (body == nullptr){ // empty node
        printf("\tInserted particle %p into %p (child %d of %p)\n",
                b.p_particle, this, child_index, parent);
        // printf("body of %p is null, b is %p\n", this, &b);
        child_index = ci;
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
        // branch this node
        // printf("Branch me once, shame on you\n");
        branch();
        // insert body into a child node
        insertChild(*body);
        // printf("Inserted *body (%p)\n", body);
        insertChild(b);
        // printf("Inserted b (%p)\n", &b);
        body = new_body;
    }
    else {
        perror("Error while inserting body into BHTree.");
    }
};

void BHTree::totalForce(particle_t* ptc, double* dmin, double* davg, int* navg){
    if (body == nullptr){ // empty node
        // printf("body of %p is null\n", this);
    }
    else if (body->n_part > 1){ // internal node
        // printf("body of %p has npart %d > 1\n", this, body->n_part);
        // check distance
        double dx = body->px - ptc->x;
        double dy = body->py - ptc->y;
        double r = sqrt(dx * dx + dy * dy);
        if (r - 1.414*(quad->length) < cutoff){ // not entire cluster beyond
                                                // cutoff
            children[NW]->totalForce(ptc, dmin, davg, navg);
            children[NE]->totalForce(ptc, dmin, davg, navg);
            children[SW]->totalForce(ptc, dmin, davg, navg);
            children[SE]->totalForce(ptc, dmin, davg, navg);
        };
    }
    else if (body->n_part == 1){ // external node
        // printf("body of %p has npart %d == 1\n", *this, body->n_part);
        apply_force(*ptc, *(body->p_particle), dmin, davg, navg);
        // printf("applied force\n");
        // parent->insert(*body, NO);
        // printf("insert complete\n");
    }
};

void BHTree::updateCenterOfMass() {
    if (body != nullptr && body->n_part > 1) { //internal node
        body->px = 0;
        body->py = 0;
        for (int i = 0; i < CHILDREN; i++) {
            // children[i]->updateCenterOfMass();
            if (children[i]->body != nullptr) {
                body->px += children[i]->body->px * children[i]->body->n_part;
                body->py += children[i]->body->py * children[i]->body->n_part;
            }
        }
        body->px /= body->n_part;
        body->py /= body->n_part;
    }
}

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
    BHTree* Tp = new BHTree(quad_root, nullptr);
    printf("built root\n");
    // insert particles into the root
    for (int i = 0; i < n; i++){
        printf("inserting particle %d\n", i);
        // pack particle into body
        Body* b = new Body(particles[i]);
        // printf("body created\n");
        Tp->insert(*b, NO);
        // printf("inserted particle\n");
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

    printf("building tree\n");
    BHTree* tree_ptr = buildTree(n, particles, SIZE); // build tree
    printf("///////////////////////////////DONE BUILDING TREE\n");
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
            tree_ptr->totalForce(&particles[i], &dmin, &davg, &navg);
            // printf("called total force\n");
        };
        // delete tree_ptr; // chop tree

        //
        //  move particles
        //
        tree_ptr->moveBodies();
        // for( int i = 0; i < n; i++ ) {
        //     move( particles[i] );
        //     tree_ptr->moveBody();
        // }

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
