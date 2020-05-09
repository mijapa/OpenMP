#include <stdio.h>
#include <omp.h>
#include <stdlib.h>

struct node
{
    float val;
    struct node *next;
};

int getMaxValueFromtabrray(int tab[],int n){
  int maximum = tab[0];
  int i = 0;
  for (i = 1; i < n; i++)
  {
    if (tab[i] > maximum)
    {
       maximum  = tab[i];
    }
  }
  return maximum;
}


void randomizeValues(int tab[], int n)
{
  int i = 0;
    #pragma omp parallel
  {
    unsigned int myseed = omp_get_thread_num();
    #pragma omp for
    for (i = 0; i < n; i++) {
      tab[i] = rand_r(&myseed);
    }
  }
}

void sort_bucket(struct node *list)
{
    int i,j,a;

    struct node *temp1;
    struct node *temp2;

    for(temp1=list;temp1!=NULL;temp1=temp1->next)
      {
        for(temp2=temp1->next;temp2!=NULL;temp2=temp2->next)
          {
            if(temp2->val < temp1->val)
              {
                a = temp1->val;
                temp1->val = temp2->val;
                temp2->val = a;
              }
           }
       }
}


void bucketSort(int *tab, unsigned int size, int count, int threads) {

        double bucketsStart,bucketsEnd,sortStart, sortEnd, resStart,resEnd;


    struct node **buckets = (struct node **) calloc(sizeof(struct node*), count* threads);
    int max = getMaxValueFromtabrray(tab, size);
    int bucket = 0;
    double bucketSize = (double) max / (double) count;

    bucketsStart = omp_get_wtime();
    #pragma omp parallel num_threads(threads)
    {
            #pragma omp for private(bucket)
        for (int i = 0; i < size; ++i) {
            int thread = omp_get_thread_num();

            bucket = (double)tab[i] / bucketSize;
            if (bucket == count) bucket--;

            struct node *tmp = (struct node *) malloc(sizeof(struct node));
            tmp->val = tab[i];
            tmp->next = buckets[bucket + thread * count];
            buckets[bucket + thread * count] = tmp;
        }

        #pragma omp for
        for(int i = 0; i< count; i++) {
            for(int thread = 1; thread < threads; thread++){
                struct node *head = buckets[i];
                        while(head != NULL && head->next != NULL){
                                head = head->next;
                }
                if(head != NULL){
                                head->next = buckets[i + thread*count];
                }
                else{
                                head = buckets[i + thread*count];
                            buckets[i]=head;
                        }
            }
        }

        bucketsEnd = omp_get_wtime();
        sortStart = omp_get_wtime();

        #pragma omp for
        for(int i = 0; i< count; i++){
                sort_bucket(buckets[i]);
        }
        sortEnd = omp_get_wtime();
    }

    int p = 0;
    resStart = omp_get_wtime();

    #pragma omp for
        for (int i = 0; i < count; i++) {
        for (struct node *a = buckets[i]; a != NULL; a = a->next) {
            tab[p] = a->val;
            p++;
        }
    }
    resEnd = omp_get_wtime();

    printf("tworzenie i laczenie kubelkow: %f\n", bucketsEnd - bucketsStart);
    printf("sortowanie: %f\n", sortEnd - sortStart);
    printf("tablica wynikowa: %f\n", resEnd - resStart);
}


int check_order(int *tab, int n) {
  for(int i = 0; i < n - 1; ++i) {
    if (tab[i] > tab[i+1]) return 0;
  }
  return 1;
}

int main(int argc, char **argv)
{
   int size = atoi(argv[1]);
   int buckets = atoi(argv[2]);
   int threads = atoi(argv[3]);

   omp_set_num_threads(threads);

   int tab[size];
   double totalStart = omp_get_wtime();

   double randomizeStart = omp_get_wtime();
   randomizeValues(tab, size);
   double randomizeEnd = omp_get_wtime();

   bucketSort(tab, size, buckets,threads);
   double totalEnd = omp_get_wtime();
   printf("tworzenie tablicy: %f\n", randomizeEnd - randomizeStart);
   printf("czas trwania programu: %f\n", totalEnd - totalStart);
//        for(int i = 0; i< size; i++){printf("%d ",tab[i]);}

  printf("Czy tablica jest posortowana? %s ",check_order(tab,size)? "Tak": "Nie :(");
   return 0;
}