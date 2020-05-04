#include <omp.h>
#include "stdlib.h"
#include "stdio.h"

void checkIfSorted(const int N, const int *A) {
    int sorted = 1;
    for (int i = 0; i < N - 1; i++) {
        if (A[i] > A[i + 1])
            sorted = 0;
    }
    if (!sorted)
        printf("The data is not sorted!!!\n");
}

void printArray(const int N, const int *A) {
    for (int j = 0; j < N; ++j) {
        printf("%i, ", A[j]);
    }
    printf("\n");
}

void fillWithRandoms(const int N, int *A, const int RANGE) {
    const int CH = 100;

#pragma omp parallel
    {
        unsigned int seed = 42 + omp_get_thread_num();
#pragma omp for \
    schedule(guided, CH)
        for (int i = 0; i < N; ++i) {

            A[i] = rand_r(&seed) % RANGE;
        }
    }
}

int my_compare(const void *a, const void *b) {
    int _a = *(int *) a;
    int _b = *(int *) b;
    if (_a < _b) return -1;
    else if (_a == _b) return 0;
    else return 1;
}

int main(int argc, char **argv) {
//    rozmiar tablicy, liczba kubełków, liczba wątków
//    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    int N, k, t;
    N = strtoul(argv[1], NULL, 10);
    k = strtoul(argv[2], NULL, 10);
    t = strtoul(argv[3], NULL, 10);

    omp_set_num_threads(t);

    int n = sizeof(int)*N;
    int *A = malloc(n);
    if (A == NULL) {
        fprintf(stderr, "Fatal: failed to allocate %zu bytes.\n", n);
        exit(1);
    }
    printf("allocate %zu Mbytes.\n", n/1024/1024);

    const int RANGE = 10000;
    const int BUCKETS_PER_THREAD = 1;

    double start_time = omp_get_wtime();
    fillWithRandoms(N, A, RANGE);
    double AfillingWithRandomsTime = omp_get_wtime() - start_time;
    printf("a) filling time: %lf\n", AfillingWithRandomsTime);

    const int nbuckets = 10;

    int global_elementsInBucket[nbuckets];
    int global_startingPosition[nbuckets];

    double parallel_start_time = omp_get_wtime();

#pragma omp parallel
    {
        int num_threads = omp_get_num_threads();
        int my_id = omp_get_thread_num();
        int thread_length = RANGE / num_threads;
        int start = thread_length * my_id;
        int end = thread_length * (my_id + 1) - 1;
        int b_index = 0;
        int *B = malloc(sizeof(int) * N);

        //spread data to threads

        for (int i = 0; i < N; ++i) {
            if (A[i] >= start && A[i] <= end) {
                B[b_index++] = A[i];
            }
        }

        global_elementsInBucket[my_id] = b_index;

        //sort values in bucket
        qsort(B, b_index, sizeof(int), my_compare);



#pragma omp barrier
        double mergingStart = omp_get_wtime();

        //calculating start position for threads
#pragma omp master
        {
            global_startingPosition[0] = 0;
            for (int k = 1; k < num_threads; ++k) {
                global_startingPosition[k] = global_startingPosition[k - 1] + global_elementsInBucket[k - 1];
            }
        };

#pragma omp barrier
//merging data from threads
        for (int l = 0; l < b_index; ++l) {
            A[global_startingPosition[my_id] + l] = B[l];
        }

    };

    double parallelTime = omp_get_wtime() - parallel_start_time;
    printf("parallel time: %lf\n", parallelTime);

    checkIfSorted(N, A);

    free(A);

    return 0;
}
