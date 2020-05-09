#include <omp.h>
#include "stdlib.h"
#include "stdio.h"

// structure to hold exectution times
struct times {
    double spread, sorting, combine;
};

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

struct times bucket_sort(int N, int k, int *A, const int SORTED_NUMBERS_RANGE) {

    int global_elementsInBucket[k];
    int global_startingPosition[k];

    struct times times;

#pragma omp parallel
    {
        int num_threads = omp_get_num_threads();
        int my_id = omp_get_thread_num();
        int thread_length = SORTED_NUMBERS_RANGE / num_threads;

        int sorting_start_number = thread_length * my_id;
        int sorting_end_number = thread_length * (my_id + 1) - 1;
        int b_index = 0;
        int *B = malloc(sizeof(int) * N);

        double spreadStart = omp_get_wtime();
        //spread data to threads

        for (int i = 0; i < N; ++i) {
            if (A[i] >= sorting_start_number && A[i] <= sorting_end_number) {
                B[b_index++] = A[i];
            }
        }

        global_elementsInBucket[my_id] = b_index;

        times.spread = omp_get_wtime() - spreadStart;

        double sortStart = omp_get_wtime();
        //sort values in bucket
        qsort(B, b_index, sizeof(int), my_compare);

        times.sorting = omp_get_wtime() - sortStart;


        double combineStart = omp_get_wtime();

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

        times.combine = omp_get_wtime() - combineStart;

    };

    return times;
}

//program realizuje algorytm sortowania kubełkowego 1
int main(int argc, char **argv) {
//    omp_set_dynamic(0);     // Explicitly disable dynamic teams

//    rozmiar tablicy, liczba kubełków, liczba wątków
    int N, k, t;
    N = strtoul(argv[1], NULL, 10);
    k = strtoul(argv[2], NULL, 10);
    t = strtoul(argv[3], NULL, 10);

    omp_set_num_threads(t);

    double global_start_time = omp_get_wtime();

    int n = sizeof(int) * N;
    int *A = malloc(n);
    if (A == NULL) {
        fprintf(stderr, "Fatal: failed to allocate %zu bytes.\n", n);
        exit(1);
    }

    const int SORTED_NUMBERS_RANGE = 10000;

    double start_time = omp_get_wtime();

    fillWithRandoms(N, A, SORTED_NUMBERS_RANGE);

    double AfillingWithRandomsTime = omp_get_wtime() - start_time;


    struct times times = bucket_sort(N, k, A, SORTED_NUMBERS_RANGE);

    double global_time = omp_get_wtime() - global_start_time;
    printf("allocate %zu Mbytes.\n", n / 1024 / 1024);
    printf("a) filling time: %lf\n", AfillingWithRandomsTime);
    printf("b) spreading time: %lf\n", times.spread);
    printf("c) sorting time: %lf\n", times.sorting);
    printf("d) combining time: %lf\n", times.combine);
    printf("e) global algorithm time: %lf\n", global_time);


    checkIfSorted(N, A);

    free(A);

    return 0;
}
