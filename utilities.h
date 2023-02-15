/* this file manages system functions.
Selected functions from misc.c and utilities.c of pycold package */
#include <stdbool.h>

#include <stdio.h>

#define WARNING_MESSAGE(message, module)          \
    write_message((message), (module), "WARNING", \
                  __FILE__, __LINE__, stdout);

#define ERROR_MESSAGE(message, module)          \
    write_message((message), (module), "ERROR", \
                  __FILE__, __LINE__, stdout);

#define RETURN_ERROR(message, module, status)       \
    {                                               \
        write_message((message), (module), "ERROR", \
                      __FILE__, __LINE__, stdout);  \
        return (status);                            \
    }

#define RETURN_WARNING(message, module, status)       \
    {                                                 \
        write_message((message), (module), "WARNING", \
                      __FILE__, __LINE__, stdout);    \
        return (status);                              \
    }

void write_message(
    const char *message, /* I: message to write to the log */
    const char *module,  /* I: module the message is from */
    const char *type,    /* I: type of the error */
    char *file,          /* I: file the message was generated in */
    int line,            /* I: line number in the file where the message was
                               generated */
    FILE *fd             /* I: where to write the log message */
);

int adjust_median_variogram(
    int *clrx,                 /* I: dates                                          */
    float **array,             /* I: input time series array                        */
    int dim1_len,              /* I: dimension 1 length in input array              */
    int dim2_start,            /* I: dimension 2 start index                        */
    int dim2_end,              /* I: dimension 2 end index                          */
    float *date_vario,         /* I/O: outputted median variogran for dates           */
    float *max_neighdate_diff, /* I/O: maximum difference for two neighbor times       */
    float *output_array,       /* I/O: output array                                   */
    int option                 /* I: option for median variogram: 1 - normal; 2 - adjust (PYCCD version) */
);

int auto_robust_fit2(
    int *clrx,     /* I: the inputted date */
    float *clry,   /* I: the inputted band */
    int start,     /* I: the start index of the input for regression. Note it starts from 0 */
    int end,       /* I: the end index of the input for regression. Note it starts from 0 */
    int df,        /* I: the degree of freedom, should be 4, 6, or 8 */
    double *coefs, /* I/O: the harmonic coefficients */
    float *rmse,   /* I/O: rmse from robust regression */
    float *v_dif   /* I/O: the residue vector at per obs basis */
);

void quick_sort_float(float arr[], int left, int right);

int auto_ts_predict_singleband(
    int *clrx,
    double *coefs,
    int df,
    int start,
    int end,
    float *pred_y);

void matlab_2d_float_median(
    float **array,       /* I: input array */
    int dim1_index,      /* I: 1st dimension index */
    int dim2_len,        /* I: number of input elements in 2nd dim */
    float *output_median /* O: output norm value */
);

float compute_adjust_rmse(
    float *clry,  /* I: input time series array                        */
    int clear_num /* I: the length of clry                             */
);

float MeanAngl_float_1d(
    float *v_diff, // input: a two-dimensional vector of different (i_count)
    int i_count    // input: the number of consecutive observations
);

void update_cft(
    int i_span,
    int n_times,
    int min_num_c,
    int mid_num_c,
    int max_num_c,
    int num_c,
    int *update_number_c);

float median_1d_float(
    float *array, /* I: input array                         */
    int start,    /* I: number of start elements in 2nd dim */
    int end       /* I: number of end elements in 2nd dim */
);

float max_1d_float(float *array, int len, int *max_id);

float max_1d_int(int *array, int len);