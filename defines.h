/* from 2darray.c */
/* Define a unique (i.e. random) value that can be used to verify a pointer
   points to an LSRD_2D_ARRAY. This is used to verify the operation succeeds to
   get an LSRD_2D_ARRAY pointer from a row pointer. */
#define SIGNATURE 0x326589ab

/* Given an address returned by the allocate routine, get a pointer to the
   entire structure. */
#define GET_ARRAY_STRUCTURE_FROM_PTR(ptr) \
    ((LSRD_2D_ARRAY *)((char *)(ptr)-offsetof(LSRD_2D_ARRAY, memory_block)))

#define VZA_LOW_BOUND 0
#define VZA_HIGH_BOUND 60
#define MAX_VZA_GROUPS 4
#define N_TIMES 3
#define MIN_NUM_C 4 /* Minimum number of coefficients           */
#define MID_NUM_C 4 /* Medium number of coefficients           */
#define MAX_NUM_C 4 /* Max number of coefficients           */
#define NUM_YEARS 365.25
#define MIN_YEARS 1
#define T_MAX_CG_VSA 99999999
#define LASSO_MIN 6
#define PI 3.1415926535897935
#define NSIGN 45
#define NUM_FC 40
#define NO_CHANGE 0
#define CHANGE_DETECTED 1
#define OUTLIER 2