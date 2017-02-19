#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <assert.h>

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
    bool contains(double, double) const; // returns whether a point is in quadrant
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
    bool in(const Quad &) const; // Function to test if body is in a certain quad domain
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
    void collapse();
    void insertChild(Body *);
    Body* emancipateChild(int);

public:
    // Constructor
    BHTree(Quad *, BHTree*, int); // create a Barnes-Hut tree with no bodies,
                             // representing the given quadrant.
    ~BHTree(){
        // delete body;
        delete quad;
        for (int i = 0; i < CHILDREN; i++) {
            delete children[i];
        }
    };
    void insert(Body *, int); // add the body to the involking Barnes-Hut tree
    void totalForce(particle_t*, double *, double *, int *); // apply force on
                                                             // particle from
                                                             // all bodies in
                                                             // the invoking
                                                             // Barnes-Hut tree
    void updateBodies();
    void updateCenterOfMass();
    int bodyCount(bool, int);
    bool consistent();
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

bool Quad::contains(double x_q, double y_q) const {
    return (x_q >= x) && (x_q < (x + length)) &&
        (y_q >= y) && (y_q < y + length);
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
    p_particle = nullptr;
};

bool Body::in(const Quad & q) const{
    bool ret = q.contains(px, py);
    // if (ret)
    //     printf("(%f, %f) is between (%f, %f) and (%f,%f)\n", px, py, q.x, q.y, q.x + q.length, q.y + q.length);
    // else
    //     printf("(%f, %f) is not between (%f, %f) and (%f,%f)\n", px, py, q.x, q.y, q.x + q.length, q.y + q.length);
    return ret;
};

/* Implementation for BHTree */
BHTree::BHTree(Quad* q, BHTree* rent, int child_order=NO){
    quad = q;
    body = new Body(0,0,0);
    for (int i = 0; i < CHILDREN; i++) {
        children[i] = nullptr;
    }
    if (rent == nullptr)
        parent = this;
    else
        parent = rent;
    child_index = child_order;
    printf("\t\tCreated node %p as child %d of %p\n",this, child_index, parent);
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

void BHTree::collapse() {
    bool didSomething = false;
    printf("Collapsing %p. Before:\n", this);
    bodyCount(true,0);
    int childcount = 0;
    for (int i = 0; i < CHILDREN; i++) {
        if (children[i] != nullptr && children[i]->body->n_part != 0) {
            childcount++;
        }
    }
    if (childcount > 1) {
        printf("error collapsing %p\n", this);
        parent->bodyCount(true, 0);
    }
    assert (childcount <= 1);
    body->n_part = 0;
    for (int i = 0; i < CHILDREN; i++) {
        if (children[i] != nullptr) {
            // children[i]->collapse();
            if (children[i]->body->n_part != 0) {
                if (children[i]->body->p_particle == nullptr) {
                    children[i]->collapse();
                }
                body = children[i]->body;
                // children[i]->body = nullptr;
                printf("Collapsed node %p into %p\n", children[i], this);
                // else {
                // }
            }
            delete children[i];
            children[i] = nullptr;
            didSomething = true;
        }
    }
    if (didSomething) {
        printf("After:\n");
        bodyCount(true, 0);
    }
}

Body* BHTree::emancipateChild(int index) {
    printf("Emancipating %p (index %d, particle %p) from %p\n",
            children[index],index, children[index]->body->p_particle, this);
    assert (!(children[index]->body->in(*(children[index]->quad))));
    Body* b = children[index]->body;
    // update center of mass
    body->px = (body->px * body->n_part - b->px * b->n_part)
        / (body->n_part-b->n_part);
    body->py = (body->py * body->n_part - b->py * b->n_part)
        / (body->n_part-b->n_part);
    body->n_part -= b->n_part;
    printf("Changed %p's npart from %d to %d\n",
            this, body->n_part + b->n_part, body->n_part);
    // detach child
    children[index]->body = new Body(0,0,0);

    if (body->n_part <= 1) {
        // now empty
        // printf("Deleting empty node %p\n", this);
        // parent->collapse();
        printf("Collapsing %p from emancipateChild\n", this);
        collapse();
    }

    assert (b != nullptr);
    return b;
}

// returns whether this subtree has changed
void BHTree::updateBodies(){
    if (body->p_particle == nullptr) { //internal node
        for (int i = 0; i < CHILDREN; i++) {
            if (children[i] != nullptr) {
                // printf("updating %p\n", children[i]);
                children[i]->updateBodies();
                if ((children[i]->body->p_particle != nullptr) &&
                        (children[i]->body->n_part != 0) &&
                        !(children[i]->body->in(*(children[i]->quad)))) {
                    // we're not in the right place anymore
                    // printf("Calling emancipate from %p.updateBodies (child %d, %p). Particle is %p\n",
                    //         this, i, children[i], children[i]->body->p_particle);
                    Body* b = emancipateChild(i);
                    body->n_part -= b->n_part;
                    printf("Changed %p's n_part from %d to %d\n", this,
                            body->n_part + b->n_part, body->n_part);
                    insert(b, NO);
                    printf("moved particle %p, now n_part is %d\n",
                            b->p_particle, body->n_part);
                    parent->bodyCount(true, 0);
                    // if (body->n_part == 1) {
                    //     collapse();
                    // }
                }
            }
        }
        // update center of mass
    } else if (body->n_part > 0) {
        // external node
        double oldx = body->p_particle->x;
        double oldy = body->p_particle->y;
        move(*(body->p_particle));
        body->px = body->p_particle->x;
        body->py = body->p_particle->y;
        // printf("Moved %p's %p from (%f,%f) to (%f,%f)\n", this, body->p_particle,
        //         oldx, oldy, body->p_particle->x, body->p_particle->y);
        // printf("Quad is (%f, %f) to (%f,%f)\n", quad->x, quad->y, quad->x + quad->length, quad->y + quad->length);
    } else {
        // printf("Attempt to update %p which has %d particles\n",
        //         this, body->n_part);
        // assert(false);
    }
    // if (body->n_part == 1 && body->p_particle == nullptr) {
    //     printf("Collapsing %p from updateBodies\n", this);
    //     collapse();
    // }
}

void BHTree::insertChild(Body * b){    // insert body into a child node
    assert (b->p_particle != nullptr);
    for (int i = 0; i < CHILDREN; i++) {
        Quad q = *(children[i]->quad);
        if (b->in(q)) {
            if (children[i]->body->n_part != 0) {
                // printf("\tchild node %p (index %d of %p) occupied\n",
                //         children[i], i, this);
            }
            children[i]->insert(b, i);
            // printf("\tinserted particle %p into %p (child %d of %p)\n",
            //         b->p_particle, children[i], i, this);
            return;
        }
    }
    // couldn't find a quadrant
    if (parent != this) {
        // printf("\t\tre-inserting particle %p into %p\n",b->p_particle, parent);
        // parent->emancipateChild(child_index);
        body->n_part -= b->n_part;
        parent->insert(b, NO);
    }
    else
        perror("Could not locate a quadrant for body.");
};

void BHTree::insert(Body * b, int ci){
    // printf("Inserting particle %p into %p\n", b->p_particle, this);
    assert(b->p_particle != nullptr);
    if (body->n_part == 0 || body == b){ // empty node
        printf("\tInserted particle %p into %p (child %d of %p). Now n_part is %d\n",
                b->p_particle, this, child_index, parent, b->n_part);
        // printf("body of %p is null, b is %p\n", this, &b);
        child_index = ci;
        delete body;
        body = b;
    }
    // else if (body->n_part > 1){ // internal node
    else if (body->p_particle == nullptr){ // internal node
        // update center of mass
        // body->px = (body->px * body->n_part + b->px * b->n_part)
        //     / (body->n_part+b->n_part);
        // body->py = (body->py * body->n_part + b->py * b->n_part)
        //     / (body->n_part+b->n_part);
        body->n_part += b->n_part;
        // insert body into a child node
        insertChild(b);
    }
    else { // external node
        assert (body->p_particle != b->p_particle);
        // create new combined body
        Body* new_body = addBody(*body, *b);
        // branch this node
        branch();
        // insert body into a child node
        insertChild(body);
        body = new_body;
        // printf("Inserted *body (%p)\n", body);
        insertChild(b);
        // printf("Inserted b (%p)\n", &b);
    }
};

void BHTree::totalForce(particle_t* ptc, double* dmin, double* davg, int* navg){
    if (body == nullptr){ // empty node
        // printf("body of %p is null\n", this);
    }
    else if (body->p_particle == nullptr && body->n_part != 0){ // internal node
        // printf("body of %p has npart %d > 1\n", this, body->n_part);
        // check distance
        double dx = body->px - ptc->x;
        double dy = body->py - ptc->y;
        double r = sqrt(dx * dx + dy * dy);
        if (r - 1.414*(quad->length) < cutoff){ // not entire cluster beyond
                                                // cutoff
            for (int i = 0; i < CHILDREN; i++) {
                // printf("Child %d of %p is %p\n", i, this, children[i]);
                if (children[i] == nullptr) {
                    printf("Error updating forces: child %d of %p is null.\n", i, this);
                    parent->bodyCount(true, 0);
                    assert(false);
                }
                children[i]->totalForce(ptc, dmin, davg, navg);
            }
        };
    }
    else if (body->p_particle != nullptr){ // external node
        // printf("body of %p has npart %d == 1\n", *this, body->n_part);
        apply_force(*ptc, *(body->p_particle), dmin, davg, navg);
        // printf("Applied force from %p to %p\n", ptc, this);
        // printf("applied force\n");
        // parent->insert(*body, NO);
        // printf("insert complete\n");
    }
};

void BHTree::updateCenterOfMass() {
    // printf("Updating centers of mass: %p. n_part = %d\n", this, body->n_part);
    if (body->n_part > 1) { //internal node
        body->px = 0;
        body->py = 0;
        body->n_part = 0;
        for (int i = 0; i < CHILDREN; i++) {
            children[i]->updateCenterOfMass();
            body->px += children[i]->body->px * children[i]->body->n_part;
            body->py += children[i]->body->py * children[i]->body->n_part;
            body->n_part += children[i]->body->n_part;
            // printf("%p now has n_part %d\n", this, body->n_part);
        }
        if (body->n_part == 1) {
            collapse();
        } else {
            body->px /= body->n_part;
            body->py /= body->n_part;
        }
    } else if (body->p_particle != nullptr && body->n_part == 1) {
        body->px = body->p_particle->x;
        body->py = body->p_particle->y;
    } else if (body->p_particle != nullptr) {
        printf("Node %p is broken: Has particle %p but n_part %d\n",
                this, body->p_particle, body->n_part);
        assert(false);
    }
}

int BHTree::bodyCount(bool verbose, int level) {
    if (verbose) {
        for (int i = 0; i < level; i++) {
            printf("\t");
        }
        printf("%p = %d ", this, body->n_part);
        if (body->p_particle != nullptr)
            printf("(%p)", body->p_particle);
    }
    if (verbose) printf("\n");
    for (int i = 0; i < CHILDREN; i++) {
        if (children[i] != nullptr) {
            children[i]->bodyCount(verbose, level + 1);
        }
    }
    return body->n_part;
}

bool BHTree::consistent() {
    if (body->p_particle != nullptr) { //external node
        if (body->n_part != 1){
            printf("%p is a particle leaf with n_part == %d\n", this, body->n_part);
            return false;
        }
        for (int i = 0; i < CHILDREN; i++) {
            if (children[i] != nullptr){
                printf("%p is a particle leaf with a child %p\n", this, children[i]);
                return false;
            }
        }
    } else { //internal node - expect non-null leaves
        int nulls = 0;
        for (int i = 0; i < CHILDREN; i++) {
            if (children[i] != nullptr) {
                nulls++;
                if (!children[i]->consistent()) return false;
            }
        }
        if (nulls != 0 && nulls != 4) {
            printf("%p has mixed null and non-null chlidren\n", this);
            return false;
        }
    }
    return true;
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
        Tp->insert(b, NO);
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
    printf("///////////////////////////////DONE BUILDING TREE. Particles: %d\n",
            tree_ptr->bodyCount(true, 0));
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
        assert(tree_ptr->bodyCount(false, 0) == n);
        //
        //  move particles
        //
        tree_ptr->updateCenterOfMass();
        tree_ptr->updateBodies();
        tree_ptr->updateCenterOfMass();
        // if (tree_ptr->bodyCount(false) != n)
        if (!(tree_ptr->consistent()) || !(tree_ptr->bodyCount(false, 0) == n)) {
            printf("BODY COUNT: %d\n",tree_ptr->bodyCount(true, 0));
            assert(false);
        } else {
            printf("consistent step %d\n", step);
        }
        assert(tree_ptr->bodyCount(false, 0) == n);
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
