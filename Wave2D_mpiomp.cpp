// Wave.cpp : Defines the entry point for the console application.
//

#include "mpi.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Timer.h"
#include <iostream>
#include <fstream>
#include<string>
#include <omp.h>

using namespace std;

const int left_to_right = 20;
const int right_to_left = 10;
const int worker_to_master = 30;


const int default_size = 100;  // the default system size
const int defaultCellWidth = 8;
const double c = 1.0;      // wave speed
const double dt = 0.1;     // time quantum
const double dd = 2.0;     // change in system

double z[3][600][600] = { 0 };

void calculate(int start, int end, int size, int max_time, int interval, int rank, int processes);

int main(int argc, char* argv[])
{
    int rank;
    int processes;
    double wtime;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);

    if (argc < 5) {
        cerr << "usage: Wave2D size max_time interval threads" << endl;
        return -1;
    }

    int size = atoi(argv[1]);
    int max_time = atoi(argv[2]);
    int interval = atoi(argv[3]);
    int nThreads = atoi(argv[4]);

    // change # of threads
    omp_set_num_threads(nThreads);

    Timer time;
    time.start();
    int start;
    int end;
    int stripe = size / processes;

    start = rank  * stripe;
    end = (rank + 1) * stripe;
    end = size > end ? end : size;
    stripe = end - start;
    cerr << "rank[" << rank << "]'s range = " << start << " ~ " << (end - 1) << endl;
    calculate(start, end, size, max_time, interval, rank, processes);


    if (rank == 0)
    {
        cerr << "Elapsed time = " << time.lap() << endl;
    }

    MPI_Finalize();                
    return 0;
}



void calculate(int start, int end, int size, int max_time, int interval, int rank, int processes)
{
    /*string fileName = "D:\\Feiping\\parallel\\ConsoleApplication1\\x64\\Debug\\";
    fileName += std::to_string(rank);
    fileName += ".txt";
    ofstream file;
    file.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);
    cerr << "open file:" << fileName << endl;*/
	
    int stripe = end - start;
    int weight = size / default_size;


#pragma omp parallel for firstprivate(size) shared(z)
    //initialize z
	for (int p = 0; p < 3; p++)
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                z[p][i][j] = 0.0; // no wave

    int dtime = 2;
    if (interval != 0 && dtime % interval == 0 && rank == 0)
    {
        cout << dtime << endl;
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                cout << z[0][j][i] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    dtime++;
    max_time--;
	
	//time =1
#pragma omp parallel for firstprivate(size, weight) shared(z)
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i > 40 * weight && i < 60 * weight  &&
                j > 40 * weight && j < 60 * weight) {
                z[0][i][j] = 20.0;
            }
            else {
                z[0][i][j] = 0.0;
            }
        }
    }

    if (interval != 0 && dtime % interval == 0 && rank == 0)
    {
        cout << dtime << endl;
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                cout << z[0][i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    dtime++;
    max_time--;
    // time = 1
    // calculate z[1][][] 
    // cells not on edge
#pragma omp parallel for firstprivate(size) shared(z)
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (i == 0 || j == 0 || i == (size - 1) || j == (size - 1))
            {
                z[1][i][j] = 0.0;
            }
            else
            {
                z[1][i][j] = z[0][i][j] + (c * c) / 2 * ((dt / dd)*(dt / dd))*
				(z[0][i + 1][j] + z[0][i - 1][j] + z[0][i][j + 1] + z[0][i][j - 1] - 4.0 * z[0][i][j]);
            }
        }
    }

    if (interval != 0 && dtime % interval == 0 && rank == 0)
    {
        cout << dtime << endl;
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                cout << z[1][i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    dtime++;
    max_time--;
	
	
    // simulate wave diffusion from time = 2
    for (int t = 2; t < max_time; t++)
    {
#pragma omp parallel for firstprivate(start, end, size, t) shared(z)
        for (int i = start; i < end; i++)
        {
            for (int j = 0; j < size; ++j)
            {
                if (i == 0 || j == 0 || i == (size - 1) || j == (size - 1))
                {
                    z[t % 3][i][j] = 0.0;
                }
                else
                {
                    z[(t % 3)][i][j] = 2.0 * z[(t - 1) % 3][i][j] - z[(t - 2) % 3][i][j] + (c*c) * (dt / dd)* (dt / dd)*
					(z[(t - 1) % 3][i + 1][j] + z[(t - 1) % 3][i - 1][j] + z[(t - 1) % 3][i][j - 1] + z[(t - 1) % 3][i][j + 1] - 
					4.0 * z[(t - 1) % 3][i][j]);
                }
            }
        }

		// rank 0 1 2 send data to right neighbors
        if (rank < processes - 1)
        {
            MPI_Send(&z[t % 3][end - 1], size, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            //file << "rank: " << rank << "send to rank" << rank + 1 << " t: " << t << " start: " << end - 1 << " size: " << size << endl;
            //for (int l = 0; l < size; ++l)
            //{
            //    file <<  z[t % 3][end - 1][l] << " ";
            //}
            //file << endl;
        }
		// rank 1 2 3 send data to left neighbors
        if (rank > 0)
        {
            //send to left neighbor
            MPI_Send(&z[t % 3][start], size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            //file << "rank: " << rank << "send to rank" << rank -1 << " t: " << t << " start: " << start << " size: " << size << endl;
            //for (int l = 0; l < size; ++l)
            //{
            //    file << z[t % 3][start][l] << " ";
            //}
            //file << endl;
        }
		// rank 0 1 2 receive data from right neighbors
        if (rank < processes - 1)
        {
            MPI_Recv(&z[t % 3][end], size, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //file << "rank: " << rank << "rev from rank" << rank + 1 << " t: " << t << " start: " << end << " size: " << size << endl;
            //for (int l = 0; l < size; ++l)
            //{
            //    file << z[t % 3][end][l] << " ";
            //}
            //file << endl;*/
        }

		// rank 1 2 3 receive data from left neighbors
        if (rank > 0)
        {
            MPI_Recv(&z[t % 3][start - 1], size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //file << "rank: " << rank << "rev from rank" << rank - 1 << " t: " << t << " start: " << start -1 << " size: " << size << endl;
            //for (int l = 0; l < size; ++l)
            //{
            //    file << z[t % 3][start-1][l] << " ";
            //}
            //file << endl;
        }
        
		
        if (interval != 0)
        {
			if ((dtime % interval == 0) || (t == max_time - 1))
			{
            if (rank > 0)
            {
                for (int k = start; k < end; ++k)
                {
					// slaves send data to master
                    MPI_Send(&z[t % 3][k], size, MPI_DOUBLE, 0, worker_to_master, MPI_COMM_WORLD);
                    //file << "[collection]rank: " << rank << "send to rank" << 0 << " t: " << t << " start: " << k << " size: " << size << endl;
                    //for (int m = start; m < end; ++m)
                    //{
                    //    for (int l = 0; l < size; ++l)
                    //    {
                    //        file << z[t % 3][m][l] << " ";
                    //    }
                    //    file << endl;
                    //}
                    //file << endl;
                }
            }
            else
            {
				//master receives data from slaves.
                for (int i = 1; i < processes; i++)
                {
                    for (int k = 0; k < stripe; ++k)
                    {
                        MPI_Recv(&z[t % 3][i*stripe + k], size, MPI_DOUBLE, i, worker_to_master, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        //file << "[collection]rank: " << rank << "recv from rank" << i << " t: " << t << " start: " << i*stripe + k << " size: " << size<< endl;
                        //for (int m = i*stripe; m < i*stripe + size; ++m)
                        //{
                        //    for (int l = 0; l < size; ++l)
                        //    {
                        //        file << z[t % 3][m][l] << " ";
                        //    }
                        //    file << endl;
                        //}
                    }
                }
				
				//master sends data to slaves
                cout << dtime << endl;
                for (int i = 0; i < size; i++) {
                    for (int j = 0; j < size; j++) {
                        cout << z[t % 3][i][j] << " ";
                    }
                    cout << endl;
                }

                cout << endl;
            }
			}
        }

        dtime++;
    } // end of simulation

   // file.close();
}
