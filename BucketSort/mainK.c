#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

// structure to hold each bucket,
// index is offset in the main array
// for first element in particular bucket
struct bucket {
    float *data;
    int count;
    int index;
};

// structure to hold exectution times
struct times {
    float placing, sorting, append;
};

// the most simply bubble sort implementation
// for sorting elements in buckets
void bubble_sort(float *data, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (data[i] <= data[j]) {
                float temp = data[i];
                data[i] = data[j];
                data[j] = temp;
            }
        }
    }
}

int compare(const void *first, const void *second) {

    float a = *((float *) first), b = *((float *) second);
    if (a < b) {
        return -1;
    } else if (a > b) {
        return 1;
    } else {
        return 0;
    }
}

// method to check if elements are sorted
int check_sort(float *data, long n) {
    long i;

    for (i = 1; i < n; i++) {
        if (data[i] < data[i - 1]) {
            return 0;
        }
    }
    return 1;
}

// main function for parallel bucket sort,
// returns structure with execution times for every steps
struct times *bucket_sort(float data[], long n, long bucket_size) {

    long i, j;

    // step is value for calculating the range of bucket,
    // every i-bucket will conatain values within [i*step, (i+1)*step)
    // we assume that date will be in range 0 - 1.0
    float start, stop, step = 1.0 / bucket_size;

    struct bucket *buckets = malloc(sizeof(struct bucket) * bucket_size);
    struct times *alg_time = malloc(sizeof(struct times));

    // initialization of all buckets with default values
#pragma omp parallel for default(shared) private(i)
    for (i = 0; i < bucket_size; i++) {
        buckets[i].data = malloc(sizeof(float) * n);
        buckets[i].count = 0;
    }

    // we are going through all data elements
    // and checking for which bucket it suits,
    // if particular number suits for bucket
    // then we synchronise writing this number
    // to global array of buckets
    start = omp_get_wtime();
#pragma omp parallel for default(shared) private(i, j)
    for (i = 0; i < n; i++) {
        for (j = 0; j < bucket_size; j++) {
            if (j * step <= data[i] && data[i] < (j + 1) * step) {
#pragma omp critical
                buckets[j].data[buckets[j].count++] = data[i];
                break;
            }
        }
    }
    stop = omp_get_wtime();
    alg_time->placing = stop - start;

    // loop to calculate proper offset values for each bucket
    buckets[0].index = 0;
    for (j = 1; j < bucket_size; j++) {
        buckets[j].index = buckets[j - 1].index + buckets[j - 1].count;
    }


    // loop for run bubble_sort of bucket on each thread
    // and then place sorted numbers in main data array using
    // offset calculated above
    start = omp_get_wtime();
#pragma omp parallel for default(shared) private(j)
    for (j = 0; j < bucket_size; j++) {
        // printf("%d,", buckets[j].count);
        qsort(buckets[j].data, buckets[j].count, sizeof(float), &compare);
        // bubble_sort(buckets[j].data, buckets[j].count);
    }
    stop = omp_get_wtime();
    alg_time->sorting = stop - start;

    // loop for placing sorted numbers in main data array using
    // offset calculated earlier, this for could be combined into one
    // with loop above but for measuring purposes it is split
    start = omp_get_wtime();
#pragma omp parallel for default(shared) private(i, j)
    for (j = 0; j < bucket_size; j++) {
        for (i = 0; i < buckets[j].count; i++) {
            data[buckets[j].index + i] = buckets[j].data[i];
        }
    }
    stop = omp_get_wtime();
    alg_time->append = stop - start;

    return alg_time;
}


// generate n random date within range [0,1.0]
float *generate_data(long n) {
    long i;
    float *data = (float *) malloc(n * sizeof(float));

#pragma omp parallel
    {
        unsigned int seed = omp_get_thread_num();
#pragma omp parallel for default(shared) private(i)
        for (i = 0; i < n; i++) {
            data[i] = ((float) rand_r(&seed) / RAND_MAX);
        }
    }
    return data;
}


int main(int argc, char *argv[]) {
    // srand((unsigned int)time(NULL));

    long size, bucket_size, threads;
    double start, end;

    if (argc != 4) {
        printf("Error, required flag size\n");
        return -1;
    }

    size = strtol(argv[1], NULL, 10);
    bucket_size = strtol(argv[2], NULL, 10);
    threads = strtol(argv[3], NULL, 10);

    omp_set_num_threads(threads);

    start = omp_get_wtime();
    float *data = generate_data(size);
    end = omp_get_wtime();
    // printf("Genrating %ld elements took %f seconds\n", size, (float)end - start);

    start = omp_get_wtime();
    struct times *alg_time = bucket_sort(data, size, bucket_size);
    end = omp_get_wtime();


    // printf("Placing numbers into buckets took %f seconds\n", alg_time->placing);
    // printf("Sorting numbers with bubble sort in each bucket took %f seconds\n", alg_time->sorting);
    // printf("Placing numbers back in array took %f seconds\n", alg_time->placing);

    printf("%lf\n", (float) end - start);
    // printf("Sorting validation: %d\n", check_sort(data, size));


    return 0;
}
