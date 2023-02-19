#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <libgen.h>
#include <stdarg.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include "utilities.h"
#include "2d_array.h"
#include "defines.h"

/*****************************************************************************
  NAME:  write_message

  PURPOSE:  Writes a formatted log message to the specified file handle.

  RETURN VALUE:  None

  NOTES:
      - Log Message Format:
            yyyy-mm-dd HH:mm:ss pid:module [filename]:line message
*****************************************************************************/

void write_message(
    const char *message, /* I: message to write to the log */
    const char *module,  /* I: module the message is from */
    const char *type,    /* I: type of the error */
    char *file,          /* I: file the message was generated in */
    int line,            /* I: line number in the file where the message was
                               generated */
    FILE *fd             /* I: where to write the log message */
)
{
    time_t current_time;
    struct tm *time_info;
    int year;
    pid_t pid;

    time(&current_time);
    time_info = localtime(&current_time);
    year = time_info->tm_year + 1900;

    pid = getpid();

    fprintf(fd, "%04d:%02d:%02d %02d:%02d:%02d %d:%s [%s]:%d [%s]:%s\n",
            year,
            time_info->tm_mon,
            time_info->tm_mday,
            time_info->tm_hour,
            time_info->tm_min,
            time_info->tm_sec,
            pid, module, basename(file), line, type, message);
}

/******************************************************************************
MODULE:  quick_sort_float

PURPOSE:  sort the scene_list & sdate based on yeardoy string

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
02/02/2020   Su Ye           original values

NOTES:
******************************************************************************/
void quick_sort_float(float arr[], int left, int right)
{
    int index = partition_float(arr, left, right);

    if (left < index - 1)
    {
        quick_sort_float(arr, left, index - 1);
    }
    if (index < right)
    {
        quick_sort_float(arr, index, right);
    }
}

/******************************************************************************
MODULE:  partition_float

PURPOSE:  partition the sorted list.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
i               partitioned value

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
02/02/2023   Su Ye        Original Development

NOTES:
******************************************************************************/
int partition_float(float arr[], int left, int right)
{
    int i = left, j = right;
    float tmp;
    float pivot = arr[(left + right) / 2];

    while (i <= j)
    {
        while (arr[i] < pivot)
        {
            i++;
        }
        while (arr[j] > pivot)
        {
            j--;
        }
        if (i <= j)
        {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            i++;
            j--;
        }
    }

    return i;
}

/******************************************************************************
MODULE:  dofit_ls

PURPOSE: Declare data type and allocate memory and do multiple linear least-square
         fit used for auto_robust_fit

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
7/20/2019    Su Ye         Original Development

NOTES:
******************************************************************************/
void dofit_ls(const gsl_matrix *X, const gsl_vector *y,
              gsl_vector *c, gsl_matrix *cov)
{
    double chisq;
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(X->size1, X->size2);
    gsl_multifit_linear(X, y, c, cov, &chisq, work);
    gsl_multifit_linear_free(work);
    // work = NULL; // SY 03242019
}

/******************************************************************************
MODULE:  auto_ts_predict

PURPOSE:  Using lasso regression fitting coefficients to predict new values

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
02042023    Su Ye            Original development
NOTES:
******************************************************************************/
int auto_ts_predict_single_band(
    int *clrx,
    double *coefs,
    int df,
    int start,
    int end,
    float *pred_y)
{
    char FUNC_NAME[] = "auto_ts_predict";
    int i;
    int nums = end - start + 1;
    double w, w2, w3;
    w = TWO_PI / AVE_DAYS_IN_A_YEAR;
    w2 = 2.0 * w;
    w3 = 3.0 * w;

    for (i = 0; i < nums; i++)
    {
        if (df == 2)
        {
            pred_y[i] = coefs[0] + coefs[1] *
                                       (double)clrx[i + start];
        }
        else if (df == 4)
        {
            pred_y[i] = coefs[0] + coefs[1] * (double)clrx[i + start] + coefs[2] * cos((double)clrx[i + start] * w) + coefs[3] * sin((double)clrx[i + start] * w);
        }
        else if (df == 6)
        {
            pred_y[i] = coefs[0] + coefs[1] * (double)clrx[i + start] + coefs[2] * cos((double)clrx[i + start] * w) + coefs[3] * sin((double)clrx[i + start] * w) + coefs[4] * cos((double)clrx[i + start] * w2) + coefs[5] * sin((double)clrx[i + start] * w2);
        }
        else if (df == 8)
        {
            pred_y[i] = coefs[0] + coefs[1] * (double)clrx[i + start] + coefs[2] * cos((double)clrx[i + start] * w) + coefs[3] * sin((double)clrx[i + start] * w) + coefs[4] * cos((double)clrx[i + start] * w2) + coefs[5] * sin((double)clrx[i + start] * w2) + coefs[6] * cos((double)clrx[i + start] * w3) + coefs[7] * sin((double)clrx[i + start] * w3);
        }
        else
        {
            RETURN_ERROR("Unsupported df number", FUNC_NAME, ERROR);
        }
    }

    return (SUCCESS);
}

/******************************************************************************
MODULE:  norm_1darray

PURPOSE:  calculate the norm of 1d array

RETURN VALUE:
Type = void
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/9/2015   Song Guo         Original Development

NOTES:
******************************************************************************/

void norm_1darray(
    float *array,      /* I: input array                                   */
    int dim2_len,      /* I: number of input elements in 2nd dim           */
    float *output_norm /* O: output norm value                             */
)
{
    int i;
    float sum = 0.0;

    for (i = 0; i < dim2_len; i++)
    {
        sum += array[i] * array[i];
    }
    *output_norm = sqrtf(sum);
}

/******************************************************************************
MODULE:  auto_robust_fit_singleband

PURPOSE:  Robust fit for one band

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
02/16/2023    Su Ye        Original Development

NOTES:
******************************************************************************/
void auto_robust_fit_singleband(
    double **clrx,
    float *clry,
    int nums,
    int start,
    double *coefs,
    int df)
{
    int i, j;
    const int p = df; /* coefficient number*/
    gsl_matrix *x, *cov;
    gsl_vector *y, *c;

    /******************************************************************/
    /*                                                                */
    /* Defines the inputs/outputs for robust fitting                  */
    /*                                                                */
    /******************************************************************/

    x = gsl_matrix_alloc(nums, p);
    y = gsl_vector_alloc(nums);

    c = gsl_vector_alloc(p);
    cov = gsl_matrix_alloc(p, p);

    /******************************************************************/
    /*                                                                */
    /* construct design matrix x for linear fit                       */
    /*                                                                */
    /******************************************************************/

    for (i = 0; i < nums; ++i)
    {
        for (j = 0; j < p; j++)
        {
            if (j == 0)
            {
                gsl_matrix_set(x, i, j, 1.0);
            }
            else
            {
                gsl_matrix_set(x, i, j, clrx[i][j - 1]);
            }
        }
        gsl_vector_set(y, i, (double)clry[i + start]);
    }

    /******************************************************************/
    /*                                                                */
    /* perform robust fit                                             */
    /*                                                                */
    /******************************************************************/

    dofit(gsl_multifit_robust_bisquare, x, y, c, cov);

    for (j = 0; j < (int)c->size; j++)
    {
        coefs[j] = gsl_vector_get(c, j);
    }

    /******************************************************************/
    /*                                                                */
    /* Free the memories                                              */
    /*                                                                */
    /******************************************************************/

    gsl_matrix_free(x);
    gsl_vector_free(y);
    gsl_vector_free(c);
    gsl_matrix_free(cov);
};

/******************************************************************************
MODULE:  auto_robust_fit2

PURPOSE:  Harmonic Robust fit for one band. Equal to auto_robust_fit2 in MATLAB

RETURN VALUE: bool (success or failure)

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/3/2023    Su Ye        Original Development

NOTES:
******************************************************************************/
int auto_robust_fit2(
    int *clrx,     /* I: the inputted date */
    float *clry,   /* I: the inputted band */
    int start,     /* I: the start index of the input for regression. Note it starts from 0 */
    int end,       /* I: the end index of the input for regression. Note it starts from 0 */
    int df,        /* I: the degree of freedom, should be 4, 6, or 8 */
    double *coefs, /* I/O: the harmonic coefficients */
    float *rmse,   /* I/O: rmse from robust regression */
    float *v_dif   /* I/O: the residue vector at per obs basis */
)
{
    int i, j;
    char FUNC_NAME[] = "auto_robust_fit2";
    int nums = end - start + 1;
    double w = TWO_PI / 365.25;
    float *yhat;
    float v_dif_norm;
    int status;
    gsl_matrix *x, *cov;
    gsl_vector *y, *c;

    /******************************************************************/
    /*                                                                */
    /* Defines the inputs/outputs for robust fitting                  */
    /*                                                                */
    /******************************************************************/

    x = gsl_matrix_alloc(nums, df);
    y = gsl_vector_alloc(nums);

    c = gsl_vector_alloc(df);
    cov = gsl_matrix_alloc(df, df);

    yhat = (float *)malloc(nums * sizeof(float));
    if (yhat == NULL)
    {
        RETURN_ERROR("Allocating yhat memory", FUNC_NAME, ERROR);
    }

    /******************************************************************/
    /*                                                                */
    /* Defines the inputs/outputs for robust fitting                  */
    /*                                                                */
    /******************************************************************/
    if (df != 2 && df != 4 && df != 6 && df != 8)
    {
        RETURN_ERROR("Only supported df = 4, 6 and 8", FUNC_NAME, ERROR);
    }

    /******************************************************************/
    /*                                                                */
    /* construct design matrix x for linear fit                       */
    /*                                                                */
    /******************************************************************/
    if (df == 2)
    {
        for (i = 0; i < nums; i++)
        {
            gsl_matrix_set(x, i, 0, 1.0);
            gsl_matrix_set(x, i, 1, (double)clrx[i + start]);
            gsl_vector_set(y, i, (double)clry[i + start]);
        }
    }

    if (df == 4)
    {
        for (i = 0; i < nums; i++)
        {
            gsl_matrix_set(x, i, 0, 1.0);
            gsl_matrix_set(x, i, 1, (double)clrx[i + start]);
            gsl_matrix_set(x, i, 2, (double)cos(w * (double)clrx[i + start]));
            gsl_matrix_set(x, i, 3, (double)sin(w * (double)clrx[i + start]));
            gsl_vector_set(y, i, (double)clry[i + start]);
        }
    }
    else if (df == 6)
    {
        for (i = 0; i < nums; i++)
        {
            gsl_matrix_set(x, i, 0, 1.0);
            gsl_matrix_set(x, i, 1, (double)clrx[i + start]);
            gsl_matrix_set(x, i, 2, (double)cos(w * (double)clrx[i + start]));
            gsl_matrix_set(x, i, 3, (double)sin(w * (double)clrx[i + start]));
            gsl_matrix_set(x, i, 4, (double)cos(2 * w * (double)clrx[i + start]));
            gsl_matrix_set(x, i, 5, (double)sin(2 * w * (double)clrx[i + start]));
            gsl_vector_set(y, i, (double)clry[i + start]);
        }
    }
    else if (df == 8)
    {
        for (i = 0; i < nums; ++i)
        {
            gsl_matrix_set(x, i, 0, 1.0);
            gsl_matrix_set(x, i, 1, (double)clrx[i + start]);
            gsl_matrix_set(x, i, 2, (double)cos(w * (double)clrx[i + start]));
            gsl_matrix_set(x, i, 3, (double)sin(w * (double)clrx[i + start]));
            gsl_matrix_set(x, i, 4, (double)cos(2 * w * (double)clrx[i + start]));
            gsl_matrix_set(x, i, 5, (double)sin(2 * w * (double)clrx[i + start]));
            gsl_matrix_set(x, i, 6, (double)cos(3 * w * (double)clrx[i + start]));
            gsl_matrix_set(x, i, 7, (double)sin(3 * w * (double)clrx[i + start]));
            gsl_vector_set(y, i, (double)clry[i + start]);
        }
    }

    /******************************************************************/
    /*                                                                */
    /* perform robust fit                                             */
    /*                                                                */
    /******************************************************************/
    dofit(gsl_multifit_robust_bisquare, x, y, c, cov);
    // dofit_ls(x, y, c, cov);

    for (j = 0; j < (int)c->size; j++)
    {
        coefs[j] = gsl_vector_get(c, j);
    }

    /******************************************************************/
    /*                                                                */
    /* predict robust fitting results                                  */
    /*                                                                */
    /******************************************************************/
    // auto_ts_predict_single_band(clrx, coefs, df, start, end, yhat);
    for (i = 0; i < nums; ++i)
    {
        // double xi = gsl_vector_get(x, i);
        // double yi = gsl_vector_get(y, i);
        gsl_vector_view v = gsl_matrix_row(x, i);
        double y_ols, y_rob, y_err;
        // double tt;
        // tt = (float)yhat[i];

        gsl_multifit_robust_est(&v.vector, c, cov, &y_rob, &y_err);
        v_dif[i] = (float)clry[i + start] - y_rob;
        // printf("%g %g \n", tt, y_rob);
        // printf("%g %g \n", y_err, (float)clry[i + start] - yhat[i]);
        // v_dif[i] = (float)y_err;
        // gsl_multifit_robust_est(&v.vector, c_ols, cov, &y_ols, &y_err);
    }
    // for (i = 0; i < nums; i++)
    // {
    //     v_dif[i] = (float)clry[i + start] - yhat[i];
    // }
    norm_1darray(v_dif, nums, &v_dif_norm);
    *rmse = (float)(v_dif_norm / sqrtf((float)(nums - df)));

    /******************************************************************/
    /*                                                                */
    /* Free the memories                                              */
    /*                                                                */
    /******************************************************************/
    free(yhat);

    /******************************************************************/
    /*                                                                */
    /* Free the memories                                              */
    /*                                                                */
    /******************************************************************/

    gsl_matrix_free(x);
    gsl_vector_free(y);
    gsl_vector_free(c);
    gsl_matrix_free(cov);

    return SUCCESS;
}

/******************************************************************************
MODULE:  dofit

PURPOSE: Declare data type and allocate memory and do multiple linear robust
         fit used for auto_robust_fit

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
3/5/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
void dofit(const gsl_multifit_robust_type *T,
           const gsl_matrix *X, const gsl_vector *y,
           gsl_vector *c, gsl_matrix *cov)
{
    gsl_multifit_robust_workspace *work = gsl_multifit_robust_alloc(T, X->size1, X->size2);
    gsl_multifit_robust(X, y, c, cov, work);
    gsl_multifit_robust_free(work);
    // work = NULL; // SY 03242019
}

/******************************************************************************
MODULE:  auto_ts_predict_singleband

PURPOSE:  Using regression fitting coefficients to predict new values for a single band

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
02062023     Su Ye            original
NOTES:
******************************************************************************/
int auto_ts_predict_singleband(
    int *clrx,
    double *coefs,
    int df,
    int start,
    int end,
    float *pred_y)
{
    char FUNC_NAME[] = "auto_ts_predict_singleband";
    int i;
    int nums = end - start + 1;
    double w, w2, w3;
    w = TWO_PI / AVE_DAYS_IN_A_YEAR;
    w2 = 2.0 * w;
    w3 = 3.0 * w;

    for (i = 0; i < nums; i++)
    {
        if (df == 2)
        {
            pred_y[i] = (float)(coefs[0] + coefs[1] *
                                               (double)clrx[i + start]);
        }
        else if (df == 4)
        {
            pred_y[i] = (float)(coefs[0] + coefs[1] * (double)clrx[i + start] + coefs[2] * cos((double)clrx[i + start] * w) + coefs[3] * sin((double)clrx[i + start] * w));
        }
        else if (df == 5)
        {
            pred_y[i] = (float)(coefs[0] + coefs[1] * cos((double)clrx[i + start] * w) + coefs[2] * sin((double)clrx[i + start] * w) + coefs[3] * cos((double)clrx[i + start] * w2) + coefs[4] * sin((double)clrx[i + start] * w2));
        }
        else if (df == 6)
        {
            pred_y[i] = (float)(coefs[0] + coefs[1] * (double)clrx[i + start] + coefs[2] * cos((double)clrx[i + start] * w) + coefs[3] * sin((double)clrx[i + start] * w) + coefs[4] * cos((double)clrx[i + start] * w2) + coefs[5] * sin((double)clrx[i + start] * w2));
        }
        else if (df == 8)
        {
            pred_y[i] = (float)(coefs[0] + coefs[1] * (double)clrx[i + start] + coefs[2] * cos((double)clrx[i + start] * w) + coefs[3] * sin((double)clrx[i + start] * w) + coefs[4] * cos((double)clrx[i + start] * w2) + coefs[5] * sin((double)clrx[i + start] * w2) + coefs[6] * cos((double)clrx[i + start] * w3) + coefs[7] * sin((double)clrx[i + start] * w3));
        }
        else
        {
            RETURN_ERROR("Unsupported df number", FUNC_NAME, ERROR);
        }
    }

    return (SUCCESS);
}

/******************************************************************************
MODULE:  matlab_2d_float_median

PURPOSE:  simulate matlab median function for 1 dimesion in 2d array float point
          number case only

RETURN VALUE:
Type = void
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
6/26/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
void matlab_2d_float_median(
    float **array,       /* I: input array */
    int dim1_index,      /* I: 1st dimension index */
    int dim2_len,        /* I: number of input elements in 2nd dim */
    float *output_median /* O: output norm value */
)
{
    int m = dim2_len / 2;

    if (dim2_len % 2 == 0)
    {
        *output_median = (array[dim1_index][m - 1] + array[dim1_index][m]) / 2.0;
    }
    else
    {
        *output_median = array[dim1_index][m];
    }
}

/******************************************************************************
MODULE:  mean included angle

PURPOSE:  caculated mean included angle of v_diff vector

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/11/2023    Su Ye        Original Development

NOTES:
******************************************************************************/

float MeanAngl_float_1d(
    float *v_diff, // input: a two-dimensional vector of difference normalized by RMSE
    int i_count    // input: the number of consecutive observations
)
{
    float y;
    int i, j;
    float product;
    float norm1;
    float norm2;
    float angl_sum = 0;
    float tmp;
    if (i_count > 1)
    {
        for (i = 0; i < i_count - 1; i++)
        {
            product = 0;
            norm1 = 0;
            norm2 = 0;

            product += v_diff[i] * v_diff[i + 1];
            norm1 += v_diff[i] * v_diff[i];
            norm2 += v_diff[i + 1] * v_diff[i + 1];

            tmp = acosf(product / (sqrtf(norm1) * sqrtf(norm2)));
            angl_sum += (tmp * 180.0) / PI;
        }
        y = angl_sum / (i_count - 1);
    }
    else
    {
        y = 0;
    }
    return y;
}

/******************************************************************************
MODULE:  update_cft

PURPOSE:  determine the number of coefficient use in the time series model

RETURN VALUE:
Type = void
Value           Description
-----           -----------
HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
void update_cft(
    int i_span,
    int n_times,
    int min_num_c,
    int mid_num_c,
    int max_num_c,
    int num_c,
    int *update_number_c)
{
    /* start with 4 coefficients model */
    if (i_span < mid_num_c * n_times)
    {
        *update_number_c = min(min_num_c, num_c);
    }
    /* start with 6 coefficients model */
    else if (i_span < max_num_c * n_times)
    {
        *update_number_c = min(mid_num_c, num_c);
    }
    /* start with 8 coefficients model */
    else
    {
        *update_number_c = min(max_num_c, num_c);
    }
}

/******************************************************************************
MODULE:  median_1d_float

PURPOSE: calculate median value of 1-d float array

RETURN VALUE:
Type = void
Value           Description
-----           -----------

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
02/14/203  Su Ye           Original Development

NOTES:
******************************************************************************/

float median_1d_float(
    float *array, /* I: input array */
    int start,    /* I: number of start elements in 2nd dim */
    int end       /* I: number of end elements in 2nd dim */
)
{
    int m = (end + start + 1) / 2;
    float output_median;

    if ((end - start + 1) % 2 == 0)
    {
        output_median = (array[m - 1] + array[m]) / 2.0;
    }
    else
    {
        output_median = array[m];
    }
    return output_median;
}

/******************************************************************************
MODULE:  max_1d_float

PURPOSE:  calculated max value of 1d float array

RETURN VALUE: max value
Type = float
Value           Description
-----           -----------

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/02/2023   Su Ye         Original Development

NOTES:
******************************************************************************/
float max_1d_float(float *array, int len, int *max_id)
{
    int i;
    float max = -9999999;
    for (i = 0; i < len; i++)
    {
        if (array[i] > max)
        {
            max = array[i];
        }
    }
    *max_id = i;
    return max;
}

/******************************************************************************
MODULE:  max_1d_int

PURPOSE:  calculated max value of 1d int array

RETURN VALUE: max value
Type = int
Value           Description
-----           -----------

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/02/2023   Su Ye         Original Development

NOTES:
******************************************************************************/
float max_1d_int(int *array, int len)
{
    int i;
    int max = -999999999;
    for (i = 0; i < len; i++)
    {
        if (array[i] > max)
        {
            max = array[i];
        }
    }
    return max;
}
