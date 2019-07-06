#include "heapsort.h"
#include <stdio.h>
#include <stdlib.h>
 
#define IS_LESS(v1, v2)  (v1 < v2)
 
void siftDown( int *a, int start, int n);
 
#define SWAP(r,s)  do{int t=r; r=s; s=t; } while(0)
 
void heap_sort( int *a, int n)
{
    int start, end;
 
    /* heapify */
    for (start = (n-2)/2; start >=0; start--) {
        siftDown( a, start, n);
    }
 
    for (end=n-1; end > 0; end--) {
        SWAP(a[end],a[0]);
        siftDown(a, 0, end);
    }
}
 
void siftDown( int *a, int start, int end)
{
    int root = start;
 
    while ( root*2+1 < end ) {
        int child = 2*root + 1;
        if ((child + 1 < end) && IS_LESS(a[child],a[child+1])) {
            child += 1;
        }
        if (IS_LESS(a[root], a[child])) {
            SWAP( a[child], a[root] );
            root = child;
        }
        else
            return;
    }
}
