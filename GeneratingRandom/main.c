#include <omp.h>
#include "stdlib.h"
#include "stdio.h"

double scheduleSequential(int N) {
    int *a = malloc(sizeof(int) * N);
    double start_time = omp_get_wtime();

        unsigned int seed = 42 + omp_get_thread_num();

        for (int i = 0; i < N; ++i) {
            a[i] = rand_r(&seed);
        }

    double time = omp_get_wtime() - start_time;
    free(a);
    return time;
}

double scheduleStatic(int N, int chunk) {
    int *a = malloc(sizeof(int) * N);
    double start_time = omp_get_wtime();
#pragma omp parallel
    {
        unsigned int seed = 42 + omp_get_thread_num();
#pragma omp for \
    schedule(static, chunk)
        for (int i = 0; i < N; ++i) {
            a[i] = rand_r(&seed);
        }
    }
    double time = omp_get_wtime() - start_time;
    free(a);
    return time;
}

double scheduleDynamic(int N, int chunk) {
    int *a = malloc(sizeof(int) * N);
    double start_time = omp_get_wtime();
#pragma omp parallel
    {
        unsigned int seed = 42 + omp_get_thread_num();
#pragma omp for \
    schedule(dynamic, chunk)
        for (int i = 0; i < N; ++i) {
            a[i] = rand_r(&seed);
        }
    }
    double time = omp_get_wtime() - start_time;
    free(a);
    return time;
}

double scheduleGuided(int N, int chunk) {
    int *a = malloc(sizeof(int) * N);
    double start_time = omp_get_wtime();
#pragma omp parallel
    {
        unsigned int seed = 42 + omp_get_thread_num();
#pragma omp for \
    schedule(guided, chunk)
        for (int i = 0; i < N; ++i) {
            a[i] = rand_r(&seed);
        }
    }
    double time = omp_get_wtime() - start_time;
    free(a);
    return time;
}

int main() {
    const long int N = 1000 * 1000 * 10;
    const int CH = 100;

    FILE *fp;
    if ((fp = fopen("times.csv", "a")) == NULL) {
        printf("Can't open data.csv in append mode!\n");
        exit(1);
    }
    fprintf(fp, "chunks,static,dynamic,guided\n");
    for (int chunk = 1; chunk <= CH; ++chunk) {
        fprintf(fp, "%d,%lf,%lf,%lf\n",
                chunk,
                scheduleStatic(N, chunk),
                scheduleDynamic(N, chunk),
                scheduleGuided(N, chunk)
                );
    }
    fclose(fp);

    if ((fp = fopen("sequential.csv", "a")) == NULL) {
        printf("Can't open data.csv in append mode!\n");
        exit(1);
    }
    fprintf(fp, "sequential\n");
        fprintf(fp, "%lf\n",
                scheduleSequential(N));
    fclose(fp);

    return 0;
}
