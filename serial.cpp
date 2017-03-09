#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <iostream>

#define K 1.0 // where n/K is the amount of space to allocate per bin

//
//  benchmarking program
//
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
//        for( int i = 0; i < n; i++ )
//        {
//            particles[i].ax = particles[i].ay = 0;
//            for (int j = 0; j < n; j++ )
//				apply_force( particles[i], particles[j],&dmin,&davg,&navg);
//
//
// ________________________________________________________________________________



        double density = .0005;
        double size = sqrt( density * n );

        // Create Bins as a cell_num by cell_num matrix of Vectors
        int cell_num = 16;
        particle_t** bins[cell_num][cell_num];
        int bindices[cell_num][cell_num];
        double cell_size = size / cell_num;

        //initialize bins
        for (int i = 0; i < cell_num; i++)
        {
            for (int j = 0; j < cell_num; j++) {
                bins[i][j] = (particle_t**) malloc( ceil(n/K) * sizeof(particle_t));
                bindices[i][j] = 0;
            }
        }

        // Fill bins matrix with particles
        for (int i = 0; i < n; i++)
        {
            double current_x = particles[i].x;
            int bin_x = floor(current_x / cell_size);

            double current_y = particles[i].y;
            int bin_y = floor(current_y / cell_size);
            int bindex = bindices[bin_x][bin_y];
            bins[bin_x][bin_y][bindex] = &particles[i];
            bindices[bin_x][bin_y] += 1;
        }

        // Compute forces


        for (int current_row = 0; current_row < cell_num; current_row++)
        {
            for (int current_col = 0; current_col < cell_num; current_col++)
            {
                //Add current bin and neighboring bins to close_bins
                particle_t** current_bin = bins[current_row][current_col];
                for (int p = 0; p < bindices[current_row][current_col]; p++ ) {
                    current_bin[p]->ax = current_bin[p]->ay = 0;
                    //Iterate through close_bins
                    for (int i = current_row - 1; i <= current_row + 1; i++) {
                        for (int j = current_col - 1; j <= current_col + 1; j++) {
                            if (i >= 0 && i < cell_num && j >= 0 && j < cell_num) {
                                // this bin is a neighbor.
                                // iterate through particles in this bin.
                                for (int np = 0; np < bindices[i][j]; np++){
                                particle_t** current_neighbor = bins[i][j];
                                apply_force( *current_bin[p],
                                        *current_neighbor[np], &dmin, &davg,
                                        &navg);
                                }
                            }
                        }
                    }
                }
            }
        }






 //____________________________________________________________________________________________________________
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ )  {
            //std::cout << particles[i].ax << std::endl;
            move( particles[i] );
        }

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
