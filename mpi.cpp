#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include <vector>
#include <array>
#include <math.h>
#define m 4  // number of bins in each direction
#define M m*m  // total number of bins to use

//
//  auxilliary function for calculating a list of neighboring bins
int neighbors(const int & this_bin, int* neighbor_bins){
    // input: this_bin, id of this bin (from 0 to M)
    // output:
    // n_neighbors: number of neighbors
    // neighbor_bins: id of all neighbor bins

    int this_col = this_bin % m;
    int this_row = (this_bin - this_col) / m;
    int count = 0;
    for (int i = -1; i < 2; i++){
        for (int j = -1; j < 2; j++){
            int nr = this_row + i;
            int nc = this_col + j;
            if ((nr >= 0) && (nr < m) && (nc >= 0) && (nc < m)){
                neighbor_bins[count] = nr * m + nc;
                count++;
            }
        }
    }
    return count;
}
//
//  benchmarking program
//
int main( int argc, char **argv )
{
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg;
    
    //
    //  process command line parameters
    //
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
    
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;
    
    
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    
    particle_t *particle_binned = (particle_t*) malloc( n * sizeof(particle_t) );
    int *bin_population = (int*) malloc( M * sizeof(int) );
    int *bin_offsets = (int*) malloc( M * sizeof(int) );
    //
    //  set up the data partitioning across processors
    //
    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    for( int i = 0; i < n_proc+1; i++ )
        partition_offsets[i] = min( i * particle_per_proc, n );
    
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    for( int i = 0; i < n_proc; i++ )
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
    
    //
    //  allocate storage for local partition
    //
    int nlocal = partition_sizes[rank];
    particle_t *local = (particle_t*) malloc( nlocal * sizeof(particle_t) );
    int* local_bininfo = (int*) malloc(nlocal * sizeof(int));
    
    std::array<std::vector<particle_t>, M> bins;
    int *bin_count = (int*) malloc( M * sizeof(int) );
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    //  step 1 : root creates N particles
    //  step 2 : root sends N/P particles to each of the P processors
    set_size( n );
    double density = .0005;
    double size = sqrt( density * n );
    double bin_size = size / m;
    if( rank == 0 )
        init_particles( n, particles );
    MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
    
    // auxilliary parameters
    particle_t* sendbuf;
    int* sendcounts = (int*) malloc( n_proc * sizeof(int) );
    int* sdispls = (int*) malloc( n_proc * sizeof(int) );
    int* recvcounts = (int*) malloc( n_proc * sizeof(int) );
    int* rdispls = (int*) malloc( n_proc * sizeof(int) );
    int* rc = (int*) malloc( n_proc * sizeof(int) );
    int* dp = (int*) malloc( n_proc * sizeof(int) );
    for (int pid = 0; pid < n_proc; pid++){
        rc[pid] = 1;
        dp[pid] = pid;
    }
    int this_bin;
    int* neighbor_bins = (int*) malloc( 8 * sizeof(int) );
    int n_neighbors, particle_id;
    int bin_x, bin_y, bin_id;
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;

        //  step 3 & 4: parallel binning
        //  loop through all particles in local partition and sort into M bins
        for (int pid = 0; pid < nlocal; pid++){
            bin_x = ceil(local[pid].x / bin_size);
            bin_y = ceil(local[pid].y / bin_size);
            bin_id = (bin_y - 1) * m + bin_x - 1;
            bins[bin_id].push_back(local[pid]);
            local_bininfo[pid] = bin_id;
        }
        //  save the size of each of the local bins
        //  allreduce bin_count across processors into bin_population
        //  calculate bin_offsets based on rolling sum of bin_population
        for (int bin_id = 0; bin_id < M; bin_id++){
            bin_count[bin_id] = bins[bin_id].size();
            MPI_Allreduce(&bin_count[bin_id], &bin_population[bin_id], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (bin_id == 0)
                bin_offsets[bin_id] = 0;
            else
                bin_offsets[bin_id] = bin_offsets[bin_id-1]+bin_population[bin_id-1];
        }
        MPI_Barrier(MPI_COMM_WORLD);
        //  alltoallv particles in each bin to particle_binned
        for (int bin_id = 0; bin_id < M; bin_id++){
            sendbuf = &bins[bin_id][0];
            for (int pid = 0; pid < n_proc; pid++){
                sendcounts[pid] = bin_count[bin_id];
                sdispls[pid] = 0;
            }
            MPI_Allgatherv(&bin_count[bin_id], 1, MPI_INT, recvcounts, rc, dp, MPI_INT, MPI_COMM_WORLD);  // load recvcounts
            rdispls[0] = 0;
            for (int pid = 1; pid < n_proc; pid++){  // load rdispls
                rdispls[pid] = rdispls[pid-1]+recvcounts[pid-1];
            }
            int bo = bin_offsets[bin_id];
            MPI_Alltoallv(sendbuf, sendcounts, sdispls, PARTICLE, &particle_binned[bo], recvcounts, rdispls, PARTICLE, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if( find_option( argc, argv, "-no" ) == -1 )
            if( fsave && (step%SAVEFREQ) == 0 )
                save( fsave, n, particles );
        
        // step 5, 6 & 7: par_update
        // par_update(local, particle_binned, int* bin_population, bin_offsets, &dmin, &davg, &navg);
        // loop through particles in local to update force
        for( int i = 0; i < nlocal; i++ )
        {
            local[i].ax = local[i].ay = 0;
            this_bin = local_bininfo[i];
            n_neighbors = neighbors(this_bin, neighbor_bins);
            for (int j = 0; j < n_neighbors; j++ ){ // loop through neighbor bins
                int bin_num = neighbor_bins[j];
                for (int k = 0; k < bin_population[bin_num]; k++){
                    particle_id = bin_offsets[bin_num] + k;
                    apply_force( local[i], particle_binned[particle_id], &dmin, &davg, &navg );
                }
            }
        }
    
        
        if( find_option( argc, argv, "-no" ) == -1 )
        {
            
            MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
            MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
            
            
            if (rank == 0){
                //
                // Computing statistical data
                //
                if (rnavg) {
                    absavg +=  rdavg/rnavg;
                    nabsavg++;
                }
                if (rdmin < absmin) absmin = rdmin;
            }
        }
        //
        //  move particles
        //
        for( int i = 0; i < nlocal; i++ )
            move( local[i] );
        for (int bin_id = 0; bin_id < M; bin_id++){
            bins[bin_id].clear();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    simulation_time = read_timer( ) - simulation_time;
    
    if (rank == 0) {
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
            fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }
    
    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );
    free( partition_offsets );
    free( partition_sizes );
    free( local );
    free( particles );
    free( bin_population );
    free( bin_offsets );
    free( local_bininfo );
    free( bin_count );
    free( neighbor_bins );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
