#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defines.h"
#include "2d_array.h"
#include "output.h"
#include "utilities.h"

/******************************************************************************
MODULE: check_datarange
PURPOSE: preprocessing night time dataset by multiple qa. Translated by the MATLAB function datarange_brdf
RETURN VALUE:
Type = int* an list of 0 or 1 representing 'not good' or 'good' obs
HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
02/03/2023   Su Ye           Original Development
******************************************************************************/
void check_datarange(
    long *buf_rad,         /* I: Lunar-BRDF-corrected DNB radiances */
    long *buf_tvza,        /* I:  viewing zenith angles from the raw Black Marble product. */
    long *buf_tmq,         /* I: mandatory quality flag from the lunar-brdf-corrected black marble product */
    long *buf_tm2,         /* I: cloud and snow mask summarized from the cloud mask of the raw Black Marble product. */
    long *buf_QFDNB,       /* I: quality flag from the Lunar-BRDF-corrected Black Marble product.  */
    long *buf_tm_buf,      /* I: cloud and snow buffer mask.  */
    int buf_len,           /* I: the number of valid scenes.  */
    double vza_low_bound,  /* I: the upper limit of vza  */
    double vza_high_bound, /* I: the upper limit of vza  */
    int *id_good           /* O: a id list for indicating 'good' observations  */
)
{
    int i;
    for (i = 0; i < buf_len; i++)
    {
        if ((buf_rad[i] < 65535) && (buf_rad[i] > 0) && (buf_QFDNB[i] == 0) && (buf_tm2[i] == 1) && (buf_tm_buf[i] == 1) && (buf_tmq[i] < 2) && (buf_tvza[i] < vza_high_bound * VSA_SCALE) && (buf_tvza[i] >= vza_low_bound * VSA_SCALE))
        {
            id_good[i] = 1;
        }
        else
        {
            id_good[i] = 0;
        }
    }
}

/******************************************************************************
MODULE: get_last_ibreak
PURPOSE: get the date closest to the t_break within an input array
RETURN VALUE: element id for the break
Type = Int
HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
02/12/2023   Su Ye           Original Development
******************************************************************************/
int get_last_ibreak(
    int num_fc,          /* I: the current number of the curve */
    int total_obs,       /* I: the total clear observation number */
    int *clrx,           /* I: the clear observation dates */
    output_t_vsa *rec_cg /* I: the change records */
)
{
    int i_break = 0;
    int j;
    if (num_fc > 0)
    {
        for (j = 0; j < total_obs; j++)
        {
            if (clrx[j] > rec_cg[num_fc - 1].t_break)
            {
                i_break = j - 1;
                break;
            }
        }
    }
    return i_break;
}

/******************************************************************************
MODULE:
PURPOSE: sorted out clrx, clry, clry_mid. Translated by the MATLAB function find_ids
RETURN VALUE:
Type = bool, success or not
HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
02/03/2023   Su Ye         Original Development
******************************************************************************/
int extract_clearobs(
    long *buf_rad,     /* I: Lunar-BRDF-corrected DNB radiances */
    long *buf_tvza,    /* I:  viewing zenith angles from the raw Black Marble product. */
    long *sdate,       /* I:  date as matlab serial date form (counting from Jan 0, 0000). Note ordinal date in python is from (Jan 1th, 0001) */
    int *id_good,      /* I: a id list for indicating 'good' observations  */
    int buf_len,       /* I: the number of valid scenes.  */
    int *vsa_bin_edge, /* I: the bin edge of defined vsa groups, such as [0, 20, 40, 60]  */
    int vsa_model_num, /* I: bin number of vsa angle, such as 3  */
    int *clear_num,    /* I/O: the number of clear obseration.  */
    int *clrx,         /* I/O: date list for clear obserations.  */
    float *clry,       /* I/O: band value list for clear obserations.  */
    int *clry_mid      /* I/O: view zenith angle group id for clear obs. If not in given range, give 0 */
)
{
    int i;
    int group_id;
    char FUNC_NAME[] = "extract_clearobs"; /* for error messages                   */

    *clear_num = 0;
    for (i = 0; i < buf_len; i++)
    {
        if (id_good[i] == 1)
        {
            clrx[*clear_num] = sdate[i];
            clry[*clear_num] = (float)buf_rad[i];

            /* initialize clry_mid as 0 meaning that not in any vsa_groups */
            clry_mid[*clear_num] = 0;
            for (group_id = 1; group_id < vsa_model_num + 1; group_id++)
            {
                if ((buf_tvza[i] < vsa_bin_edge[group_id] * VSA_SCALE) && (buf_tvza[i] >= vsa_bin_edge[group_id - 1] * VSA_SCALE))
                {
                    clry_mid[*clear_num] = (int)group_id;
                    break;
                }
            }
            if (clry_mid[*clear_num] == 0)
                RETURN_ERROR("clrvza include angle outi_starte the bounds", FUNC_NAME, ERROR);

            *clear_num = *clear_num + 1;
        }
    }
    return SUCCESS;
}

/******************************************************************************
MODULE:  compute_adjust_rmse

PURPOSE:  single-band calculated adjusted RMSE for a given time series

RETURN VALUE: adjusted rmse
Type = float
Value           Description
-----           -----------

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/02/2023   Su Ye         Original Development

NOTES:
******************************************************************************/

float compute_adjust_rmse(
    float *clry,  /* I: input time series array                        */
    int clear_num /* I: the length of clry                             */
)
{
    int j;                                    /* loop indecies                                    */
    float *var;                               /* pointer for allocation variable memory           */
    char FUNC_NAME[] = "compute_adjust_rmse"; /* for error messages                   */
    int m;
    float rmse;

    if (clear_num == 0)
        return 0;

    if (clear_num == 1)
        return clry[0];

    var = (float *)malloc((clear_num - 1) * sizeof(int));
    if (var == NULL)
    {
        RETURN_ERROR("Allocating var memory", FUNC_NAME, ERROR);
    }

    for (j = 0; j < clear_num - 1; j++)
    {
        var[j] = fabsf(clry[j + 1] - clry[j]);
    }

    quick_sort_float(var, 0, clear_num - 2);
    m = (clear_num - 1) / 2;
    if ((clear_num - 1) % 2 == 0)
        rmse = (var[m - 1] + var[m]) / 2.0;
    else
        rmse = var[m];

    free(var);

    return rmse;
}

/******************************************************************************
MODULE:  identify change

PURPOSE:  converted from MATLAB function'confirmChange'

RETURN VALUE: identification result
Type = int
Value               Description
-----               -----------
NO_CHANGE           No change detected
OUTLIER             Outlier identified
CHANGE_DETECTED   Change being identified
HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/14/2023   Su Ye         Original Development

NOTES:
******************************************************************************/
int identifyChange(
    float *v_dif_mag, /* I: 1-band difference between actual and predicted value */
    float rmse,       /* I: minimum rmse  */
    int num_toler,    /* I: outlier tolerance number */
    double t_cg,      /* I: change threshold */
    float nsign,      /* I: threshold for mean included angle */
    int start,        /* I: start id of v_dif_mag */
    int end,          /* I: end id of v_dif_mag  */
    float Tmax_cg)    /* I: outlier threshold */
{
    int i;
    int n_anomaly = 0;
    float mean_angle;
    float *tmp_v_diff_1d;
    int num = end - start + 1;
    tmp_v_diff_1d = (float *)calloc(num, sizeof(float));

    for (i = start; i < end + 1; i++)
    {
        tmp_v_diff_1d[i - start] = fabsf(v_dif_mag[i] / rmse);
    }

    /* if the first obs is smaller than t_cg, it is always no change */
    if (tmp_v_diff_1d[0] < t_cg)
    {
        return NO_CHANGE;
    }

    for (i = 0; i < num; i++)
    {
        if (tmp_v_diff_1d[i] > t_cg)
        {
            n_anomaly = n_anomaly + 1;
        }
    }
    mean_angle = MeanAngl_float_1d(tmp_v_diff_1d, num);

    if ((n_anomaly >= (num - num_toler)) & (mean_angle < nsign))
    {
        free(tmp_v_diff_1d);
        return CHANGE_DETECTED;
    }
    else
    {
        if (tmp_v_diff_1d[0] > Tmax_cg)
        {
            free(tmp_v_diff_1d);
            return OUTLIER;
        }
        else
        {
            free(tmp_v_diff_1d);
            return NO_CHANGE;
        }
    }
}

/******************************************************************************
MODULE:  vsa_cold
PURPOSE:  main function for vsa_cold
RETURN VALUE:
Type = int (SUCCESS OR FAILURE)
HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
02/02/2023   Su Ye         Original Development
The datatype long is used for an easier python wrapper in future.
The use of long won't occupy too much memory as this function is only pixel-based
******************************************************************************/
int vsa_cold_detect(
    long *sdate,       /* I:  date as matlab serial date form (counting from Jan 0, 0000). Note ordinal date in python is from (Jan 1th, 0001) */
    long *buf_rad,     /* I: Lunar-BRDF-corrected DNB radiances */
    long *buf_tmq,     /* I: mandatory quality flag from the lunar-brdf-corrected black marble product */
    long *buf_tvza,    /* I:  viewing zenith angles from the raw Black Marble product. */
    long *buf_tvaa,    /* I: viewing azimuth angles from the raw Black Marble product.  */
    long *buf_tm2,     /* I: cloud and snow mask summarized from the cloud mask of the raw Black Marble product. */
    long *buf_QFDNB,   /* I: quality flag from the Lunar-BRDF-corrected Black Marble product. */
    long *buf_tm_buf,  /* I: cloud and snow buffer mask.  */
    int *vsa_bin_edge, /* I: the bin edge of defined vsa groups, such as [0, 20, 40, 60]  */
    int vsa_model_num, /* I: model number of vsa angle. E.g., '3' models are [0, 20], [20, 40] and [40, 60]  */
    int buf_len,       /* I: number of valid scenes  */
    int pos,           /* I: the position id of pixel */
    int conse,         /* I: consecutive obs number for break identification */
    double t_cg,       /* I: change probaility threshold for break */
    int num_toler,     /* I: outlier number for tolerance */
    int *num_fc,       /* O: number of fitting curves                       */
    output_t_vsa *rec_cg /* O: outputted structure for VSA_COLD results    */)
{
    int status; /* function return result */
    int i, j, m, k, b;
    int id_last = 0; /* anomaly id at the last, used to indicate probability */
    char FUNC_NAME[] = "vsa_cold_detect";
    int result = 0;
    int clear_num = 0; /* clear observations */

    int *clrx;   /* clear observation dates */
    float *clry; /* clear observation radiance */
    int *clrx_tmp;
    float *clry_tmp;
    int *clry_mid;                         /* vsa-based group id, e.g., 1, 2, 3. 0 means outside  */
    int *id_good;                          /* index for clear observations */
    float sub_adjust_rmse[MAX_VZA_GROUPS]; /* sub-group adjust RMSE */
    int sub_istart[MAX_VZA_GROUPS];        /* sub-group i_start */
    int sub_last_ibreak[MAX_VZA_GROUPS];   /* sub-group i_break */
    int sub_bl_train[MAX_VZA_GROUPS];      /* sub-group bl_train. 0 means initialization stage; 1 means continuous monitoring stage */
    int sub_i[MAX_VZA_GROUPS];             /* sub-group i_break */
    int sub_end[MAX_VZA_GROUPS];           /* sub-group i_end */
    float sub_rmse[MAX_VZA_GROUPS];        /* sub-group rmse, i.e., max(adjust_rmse, normal_rmse) */
    int sub_day_span[MAX_VZA_GROUPS];      /* sub-group day span */
    float sub_prob[MAX_VZA_GROUPS];        /* sub-group probability (at the tail) */
    int sub_prob_tend[MAX_VZA_GROUPS];     /* sub-group t_end associated with each probability (at the tail) */
    int sub_prob_tbreak[MAX_VZA_GROUPS];   /* sub-group t_break associated with each probability (at the tail) */
    // int sub_i[MAX_VZA_GROUPS];
    int model_num = vsa_model_num + 1; /* all_model + n vsa_model */
    int model_2process[] = {0, 0};     // the model id list to be processed in one loop, the first element is always 0 (all model)
    // double tmp_rmse;

    // Assigned 2-d array to save the time for re-extracting clry_n and clrx_n when each obs
    // is introduced in MATLAB codes
    int **sub_clrx;
    float **sub_clry;
    double **sub_fit_cft;
    double *tmp_fit_cft;
    float *tmp_v_dif;
    int mid, i_local;
    float r_dif_norm;        /* norm of rei_startue difference values          */
    float r_start;           /* residue at the start of observation(s)    */
    float r_end;             /* residue at the end of observastion(s)     */
    float r_slope;           /* residue for anormalized slope values   */
    int globel_t_start = NA; /* the unique i_start for all_model and 3 vsa model  */
    int i_break;
    int i_start;
    int i_ini, ini_conse;
    float *min_rmse;
    float **v_dif_mag; /* change magnitudes. Note not normalized by RMSE yet */
    float *tmp_v_diff_1d, *tmp_v_dif_mag_1d;
    float *tmp_sub_tstart; /* tmp float istart for sorting  */

    int i_conse;
    float ts_pred_temp;
    int num_c = 4;
    float *tmp_dif_mag;
    int identify_result;
    int max_id;
    float tmp;
    int i_span;
    int tmp_num_c;
    float time_span;
    int curve_indicator; /* indicate if any curve is found at the end of time series */

    if (buf_len == 0)
        return SUCCESS;

    if (model_num > MAX_VZA_GROUPS)
        RETURN_ERROR("vsa_model_num should be no more than the maximum number - 1", FUNC_NAME, ERROR);

    sub_clrx = (int **)allocate_2d_array(model_num, buf_len,
                                         sizeof(int));
    if (sub_clrx == NULL)
    {
        RETURN_ERROR("Allocating sub_clrx memory", FUNC_NAME, FAILURE);
    }

    sub_clry = (float **)allocate_2d_array(model_num, buf_len,
                                           sizeof(int));
    if (sub_clry == NULL)
    {
        RETURN_ERROR("Allocating sub_clry memory", FUNC_NAME, FAILURE);
    }

    sub_fit_cft = (double **)allocate_2d_array(model_num, MAX_NUM_C,
                                               sizeof(double));
    if (sub_fit_cft == NULL)
    {
        RETURN_ERROR("Allocating sub_fit_cft memory", FUNC_NAME, FAILURE);
    }

    id_good = (int *)calloc(buf_len, sizeof(int));
    clrx = (int *)calloc(buf_len, sizeof(int));
    clry = (float *)calloc(buf_len, sizeof(int));
    clrx_tmp = (int *)calloc(buf_len, sizeof(int));
    clry_tmp = (float *)calloc(buf_len, sizeof(float));
    clry_mid = (int *)calloc(buf_len, sizeof(int));
    tmp_v_dif = (float *)calloc(buf_len, sizeof(float));
    tmp_sub_tstart = (float *)calloc(model_num, sizeof(float));
    min_rmse = (float *)calloc(model_num, sizeof(float));

    v_dif_mag = (float **)allocate_2d_array(model_num, conse,
                                            sizeof(float));
    if (v_dif_mag == NULL)
    {
        RETURN_ERROR("Allocating v_dif_mag memory",
                     FUNC_NAME, FAILURE);
    }

    for (j = 0; j < buf_len - 1; j++)
        if (sdate[j] > sdate[j + 1])
            RETURN_ERROR("The inputted data do not follow an ascending order!", FUNC_NAME, ERROR);

    /********************************************************************/
    /*                                                                  */
    /*               select valid observation using qa band             */
    /*                                                                  */
    /********************************************************************/
    check_datarange(buf_rad, buf_tvza, buf_tmq, buf_tm2, buf_QFDNB, buf_tm_buf, buf_len,
                    vsa_bin_edge[0], vsa_bin_edge[vsa_model_num], id_good);

    extract_clearobs(buf_rad, buf_tvza, sdate, id_good, buf_len, vsa_bin_edge, vsa_model_num,
                     &clear_num, clrx, clry, clry_mid);

    /********************************************************************/
    /*                                                                  */
    /*               initialize sub_group model parameters              */
    /*                                                                  */
    /********************************************************************/
    for (j = 0; j < MAX_VZA_GROUPS; j++)
    {
        sub_adjust_rmse[j] = 0;
        sub_istart[j] = 0;
        sub_last_ibreak[j] = 0;
        sub_bl_train[j] = 0;
        sub_i[j] = 0;
        sub_end[j] = 0;
        sub_rmse[j] = 0;
        sub_day_span[j] = 0;
        sub_prob[j] = -9999999; // we assigned the highest prob as the final one, so give a large negative value
        sub_prob_tend[j] = 0;
        sub_prob_tbreak[j] = 0;
    }

    /********************************************************************/
    /*                                                                  */
    /*               sorted our observations into sub group             */
    /*                                                                  */
    /********************************************************************/
    for (j = 0; j < clear_num; j++)
    {
        // sub_model plus 1
        if (clry_mid[j] > 0)
        {
            sub_clry[clry_mid[j]][sub_end[clry_mid[j]]] = clry[j];
            sub_clrx[clry_mid[j]][sub_end[clry_mid[j]]] = clrx[j];
            sub_end[clry_mid[j]] = sub_end[clry_mid[j]] + 1;
        }

        // all_model always plus 1
        sub_clry[0][sub_end[0]] = clry[j];
        sub_clrx[0][sub_end[0]] = clrx[j];
        sub_end[0] = sub_end[0] + 1;
    }

    /********************************************************************/
    /*                                                                  */
    /*               compute minimum RMSE                               */
    /*                                                                  */
    /********************************************************************/
    for (m = 0; m < model_num; m++)
    {
        sub_adjust_rmse[m] = compute_adjust_rmse(sub_clry[m], sub_end[m]);
    }

    i = 0; // we start from 0 as we need distribute obs into sub_clrx and sub_clry
    while (i < clear_num - conse)
    {
        model_2process[1] = clry_mid[i]; // decide the sub model to be processed
        /************************************************************/
        /*                                                          */
        /* looping model to be processed: all_model + 1 vsa_model   */
        /*                                                          */
        /************************************************************/

        for (mid = 0; mid < 2; mid++) // 2 = all_model(0) + 1 VSA_model
        {

            /********************************************************************************/
            /*                                                                              */
            /* mid, i_local and i_start are used to shorten line length for code elegance   */
            /*                                                                              */
            /********************************************************************************/
            mid = model_2process[mid]; // current model id, started from 0
            i_local = sub_i[mid];      // current observation id, started from 0
            i_start = sub_istart[mid]; // current i_start id, started from 0
            i_span = sub_i[mid] - sub_istart[mid] + 1;
            time_span = (sub_clrx[mid][i_local] - sub_clrx[mid][i_start]) / NUM_YEARS;

            /**********************************************************/
            /*                                                        */
            /* basic requrirements 1: 1) enough observations;          */
            /*                      2) enough time                    */
            /*                                                        */
            /**********************************************************/
            if ((i_span < N_TIMES * MIN_NUM_C) || (time_span < (float)MIN_YEARS))
            {
                sub_i[mid] = sub_i[mid] + 1;
                continue;
            }

            /**********************************************************/
            /*                                                        */
            /* basic requrirements 2: 1) enough obs left at the end;  */
            /*                                                        */
            /**********************************************************/
            if (sub_end[mid] - 1 - sub_i[mid] + 1 < conse)
            {
                sub_i[mid] = sub_i[mid] + 1;
                continue;
            }

            if (sub_bl_train[mid] == 0)
            {
                // no tmask procedure
                // no i_dense after confirming with Zhe as nighttime data is dense enough. SY 02/03/2023

                status = auto_robust_fit2(sub_clrx[mid], sub_clry[mid], i_start, i_local, MAX_NUM_C,
                                          sub_fit_cft[mid], &sub_rmse[mid], tmp_v_dif);

                if (status != SUCCESS)
                {
                    RETURN_ERROR("Calling auto_robust_fit2 during model initilization\n",
                                 FUNC_NAME, FAILURE);
                }

                /**********************************************/
                /*                                            */
                /*                  Stablity test.             */
                /*                                            */
                /**********************************************/
                r_dif_norm = 0.0;
                min_rmse[mid] = fmax(sub_rmse[mid], sub_adjust_rmse[mid]);
                r_start = tmp_v_dif[0] / min_rmse[mid];
                r_end = tmp_v_dif[i_local - i_start] / min_rmse[mid];
                r_slope = sub_fit_cft[mid][1] *
                          (sub_clrx[mid][i_local] - sub_clrx[mid][i_start] + 1) / min_rmse[mid];

                r_dif_norm += (fabs(r_slope) + fmax(fabs(r_start), fabs(r_end)));

                /**************************************************/
                /*                                                */
                /* if not stable, i and i_start move forward      */
                /*                                                */
                /**************************************************/
                if (r_dif_norm > t_cg)
                {
                    sub_istart[mid] = sub_istart[mid] + 1;
                    sub_i[mid] = sub_i[mid] + 1;
                    continue;
                }

                /**************************************************/
                /*                                                */
                /* Model is ready.                                */
                /*                                                */
                /**************************************************/
                sub_bl_train[mid] = 1;

                /**************************************************/
                /*                                                */
                /* Count difference of i for each iteration.      */
                /*                                                */
                /**************************************************/
                sub_day_span[mid] = 0;

                /**************************************************/
                /*                                                */
                /* Find the previous break point.                 */
                /*                                                */
                /**************************************************/
                i_break = get_last_ibreak(*num_fc, sub_i[mid], sub_clrx[mid], rec_cg);

                if (i_start > i_break)
                {
                    /**********************************************/
                    /*                                            */
                    /* Model fit at the beginning of the time     */
                    /* series.                                    */
                    /*                                            */
                    /**********************************************/

                    for (i_ini = i_start - 1; i_ini >= i_break; i_ini--)
                    {
                        if (i_ini - i_break < conse)
                        {
                            ini_conse = i_ini - i_break + 1;
                        }
                        else
                        {
                            ini_conse = conse;
                        }

                        if (ini_conse == 0)
                        {
                            RETURN_ERROR("No data point for model fit at "
                                         "the begining",
                                         FUNC_NAME, FAILURE);
                        }

                        /******************************************/
                        /*                                        */
                        /* Allocate memory for model_v_dif,       */
                        /* vec_magg for the non-stdin             */
                        /* branch here.                           */
                        /*                                        */
                        /******************************************/
                        tmp_v_diff_1d = (float *)malloc(ini_conse * sizeof(float));
                        if (tmp_v_diff_1d == NULL)
                        {
                            RETURN_ERROR("Allocating tmp_v_diff_1d memory",
                                         FUNC_NAME, FAILURE);
                        }

                        tmp_v_dif_mag_1d = (float *)malloc(ini_conse * sizeof(float));
                        if (tmp_v_dif_mag_1d == NULL)
                        {
                            RETURN_ERROR("Allocating tmp_v_dif_mag_1d memory",
                                         FUNC_NAME, FAILURE);
                        }

                        /******************************************/
                        /*                                        */
                        /* Detect change.                         */
                        /* value of difference for adj_conse      */
                        /* observations                          */
                        /* Record the magnitude of change.        */
                        /*                                        */
                        /******************************************/
                        for (i_conse = 1; i_conse < ini_conse + 1; i_conse++) // SY 09192018
                        {
                            auto_ts_predict_singleband(clrx, sub_fit_cft[mid], MAX_NUM_C, i_ini - i_conse + 1,
                                                       i_ini - i_conse + 1, &ts_pred_temp);
                            tmp_v_dif_mag_1d[i_conse - 1] = (float)(clry[i_ini - i_conse + 1] - ts_pred_temp);
                        }

                        /********************************************************************/
                        /*                                                                  */
                        /* identify change based on t_c, mean_included angle and num_toler  */
                        /*                                                                  */
                        /********************************************************************/
                        min_rmse[mid] = fmax(sub_rmse[mid], sub_adjust_rmse[mid]);
                        identify_result = identifyChange(tmp_v_dif_mag_1d, min_rmse[mid], num_toler,
                                                         (double)t_cg, NSIGN, 0, ini_conse, T_MAX_CG_VSA);

                        if (identify_result == CHANGE_DETECTED) /* change detected */
                        {
                            free(tmp_v_diff_1d);
                            free(tmp_v_dif_mag_1d);
                            break;
                        }
                        else if (identify_result == OUTLIER) /* false change */
                        {
                            for (k = i_ini; k < sub_i[mid]; k++)
                            {
                                sub_clrx[mid][k] = sub_clrx[mid][k + 1];
                                sub_clry[mid][k] = sub_clry[mid][k + 1];
                            }
                            sub_i[mid] = sub_i[mid] - 1;
                            sub_end[mid] = sub_end[mid] - 1;
                        }

                        free(tmp_v_diff_1d);
                        free(tmp_v_dif_mag_1d);

                        /**************************************/
                        /*                                    */
                        /* Update i_start if i_ini is not a   */
                        /* confirmed break.                   */
                        /*                                    */
                        /**************************************/
                        sub_istart[m] = i_ini;

                    } // end for (i_ini = i_start-1; i_ini >= i_break; i_ini--)
                }     // end for if (i_start > i_break)
            }         // if (sub_bl_train[mid] == 0)

            /******************************************************/
            /*                                                    */
            /* Continuous monitoring started!!!                   */
            /*                                                    */
            /******************************************************/
            if (sub_bl_train[mid] == 1)
            {
                /**************************************************/
                /*                                                */
                /* Determine the time series model.               */
                /*                                                */
                /**************************************************/

                update_cft(i_span, N_TIMES, MIN_NUM_C, MID_NUM_C, MAX_NUM_C,
                           num_c, &tmp_num_c);

                /************************************************************/
                /*                                                          */
                /* initial model fit when there are not many observations.  */
                /* if (sub_day_span == 0 || ids_old_len < (N_TIMES * MIN_NUM_C)) */
                /*                                                          */
                /************************************************************/

                if (sub_day_span[mid] == 0 || i_span <= (N_TIMES * MIN_NUM_C))
                {
                    /**********************************************/
                    /*                                            */
                    /* update sub_day_span at each iteration.      */
                    /*                                            */
                    /**********************************************/

                    sub_day_span[mid] = clrx[i_local] - clrx[i_start] + 1;

                    status = auto_robust_fit2(sub_clrx[mid], sub_clry[mid], i_start, i_local, MAX_NUM_C,
                                              sub_fit_cft[mid], &sub_rmse[mid], tmp_v_dif);
                    if (status != SUCCESS)
                    {
                        RETURN_ERROR("Calling auto_ts_fit_float during continuous monitoring\n",
                                     FUNC_NAME, FAILURE);
                    }

                    /**********************************************/
                    /*                                            */
                    /* Updating information for the first         */
                    /* iteration.  Record time of curve start and */
                    /* time of curve end.                         */
                    /*                                            */
                    /**********************************************/

                    rec_cg[*num_fc].t_start = clrx[i_start];
                    rec_cg[*num_fc].t_end = clrx[i_local];

                    /**********************************************/
                    /*                                            */
                    /* No break at the moment.                    */
                    /*                                            */
                    /**********************************************/

                    rec_cg[*num_fc].t_break = 0;

                    /**********************************************/
                    /*                                            */
                    /* Record change probability, number of       */
                    /* observations, fit category.                */
                    /*                                            */
                    /**********************************************/
                    rec_cg[*num_fc].change_prob = 0;
                    rec_cg[*num_fc].num_obs[mid] = i_local - i_start + 1;
                    rec_cg[*num_fc].category = 0 + tmp_num_c;
                    rec_cg[*num_fc].pos = pos;

                    /******************************************/
                    /*                                        */
                    /* Record rmse of the pixel.              */
                    /*                                        */
                    /******************************************/
                    rec_cg[*num_fc].rmse[mid] = sub_rmse[mid];

                    /******************************************/
                    /*                                        */
                    /* Record change magnitude.               */
                    /*                                        */
                    /******************************************/

                    for (k = 0; k < MAX_NUM_C; k++)
                    {
                        /**************************************/
                        /*                                    */
                        /* Record fitted coefficients.        */
                        /*                                    */
                        /**************************************/
                        rec_cg[*num_fc].coefs[mid][k] = sub_fit_cft[mid][k];
                    }

                    /**********************************************/
                    /*                                            */
                    /* Detect change, value of difference for     */
                    /* adj_conse observations.                     */
                    /*                                            */
                    /**********************************************/
                    // for (k = 0; k < conse; k++)
                    // {
                    //     for (b = 0; b < NUM_LASSO_BANDS; b++)
                    //         v_diff[b][k] = 0;
                    //     for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                    //         v_dif_mag[b][k] = 0;
                    // }

                    for (i_conse = 1; i_conse < conse + 1; i_conse++) // SY 09192018
                    {
                        /**************************************/
                        /*                                    */
                        /* Absolute differences.              */
                        /*                                    */
                        /**************************************/
                        status = auto_ts_predict_singleband(clrx, sub_fit_cft[mid], tmp_num_c,
                                                            i_local + i_conse, i_local + i_conse, &ts_pred_temp);
                        v_dif_mag[mid][i_conse - 1] = (float)sub_clry[mid][i_local + i_conse] - ts_pred_temp; // SY 09192018

                        /******************************/
                        /*                            */
                        /* Minimum rmse,              */
                        /* z-scores.                  */
                        /*                            */
                        /******************************/
                    }
                    min_rmse[mid] = (float)fmax(sub_adjust_rmse[mid], sub_rmse[mid]);
                } // end for if (sub_day_span == 0 || i_span <= (N_TIMES * MIN_NUM_C))
                else
                {
                    /*********************************************/
                    /*                                           */
                    /* Update coefficent at each iteration year. */
                    /*                                           */
                    /*********************************************/
                    sub_day_span[mid] = clrx[i_local] - clrx[i_start] + 1;

                    status = auto_robust_fit2(sub_clrx[mid], sub_clry[mid], i_start, i_local,
                                              MAX_NUM_C, sub_fit_cft[mid], &sub_rmse[mid], tmp_v_dif);
                    if (status != SUCCESS)
                    {
                        RETURN_ERROR("Calling auto_robust_fit2 during continuous monitoring\n",
                                     FUNC_NAME, FAILURE);
                    }

                    /***************************************************/
                    /*                                                 */
                    /* Record fitted coefficients and RSME.            */
                    /*                                                 */
                    /***************************************************/
                    for (k = 0; k < MAX_NUM_C; k++)
                    {

                        rec_cg[*num_fc].coefs[mid][k] = sub_fit_cft[mid][k];
                    }
                    rec_cg[*num_fc].rmse[mid] = sub_rmse[mid];

                    /******************************************/
                    /*                                        */
                    /* Record number of observations, fit     */
                    /* category.                              */
                    /*                                        */
                    /******************************************/
                    rec_cg[*num_fc].num_obs[mid] = i_local - i_start + 1;
                    rec_cg[*num_fc].category = 0 + tmp_num_c;

                    /**********************************************/
                    /*                                            */
                    /* Record time of curve end.                  */
                    /*                                            */
                    /**********************************************/
                    rec_cg[*num_fc].t_end = clrx[i_local];

                    /**********************************************/
                    /*                                            */
                    /* Move the ith col to i-1th col.             */
                    /*                                            */
                    /**********************************************/

                    for (k = 0; k < conse - 1; k++)
                    {
                        v_dif_mag[mid][k] = v_dif_mag[mid][k + 1];
                    }
                    v_dif_mag[mid][conse - 1] = 0.0;

                    /************************************************/
                    /*                                              */
                    /*      predicting obs values                   */
                    /*                                              */
                    /************************************************/
                    status = auto_ts_predict_singleband(clrx, sub_fit_cft[mid], MIN_NUM_C,
                                                        i_local + conse, i_local + conse, &ts_pred_temp);
                    v_dif_mag[mid][conse - 1] = (float)sub_clry[mid][i_local + conse] - ts_pred_temp;

                    /**************************************/
                    /*                                    */
                    /*     Minimum rmse.                  */
                    /*                                    */
                    /**************************************/
                    min_rmse[mid] = fmax((double)sub_adjust_rmse[mid], sub_rmse[mid]);

                } // else update frequency

                /********************************************************************/
                /*                                                                  */
                /* identify change based on t_c, mean_included angle and num_toler  */
                /*                                                                  */
                /********************************************************************/
                identify_result = identifyChange(v_dif_mag[mid], min_rmse[mid], num_toler,
                                                 (double)t_cg, NSIGN, 0, conse - 1, T_MAX_CG_VSA);

                /**********************************************/
                /*                                            */
                /* change identified.                         */
                /*                                            */
                /**********************************************/
                if (identify_result == CHANGE_DETECTED)
                {
                    /***************************************************************/
                    /*                                                             */
                    /* check the beginning of the curve. use the latest istart     */
                    /*                                                             */
                    /***************************************************************/
                    for (m = 0; m < model_num; m++)
                    {
                        tmp_sub_tstart[m] = sub_clrx[m][sub_istart[m]];
                    }

                    // get the last ibreak for all_model
                    for (m = 0; m < model_num; m++)
                        sub_last_ibreak[m] = get_last_ibreak(*num_fc, sub_i[m], sub_clrx[m], rec_cg);

                    /**********************************************************************/
                    /*                                                                    */
                    /* look back to decide temporal segment at the beginning of the curve */
                    /*                                                                    */
                    /**********************************************************************/
                    if ((*num_fc == 0) || (sub_clrx[0][sub_istart[0]] - sub_clrx[0][sub_last_ibreak[0]] + 1 > NUM_YEARS)) // we used all_model to compute length
                    {
                        /**********************************************************************/
                        /*                                                                    */
                        /* as long as one model has enough obs, one need to save curve        */
                        /*                                                                    */
                        /**********************************************************************/
                        curve_indicator = FALSE;
                        for (m = 0; m < model_num; m++)
                        {
                            if (sub_istart[m] - sub_last_ibreak[m] >= LASSO_MIN) // actually sub_istart[m] - 1 - sub_last_ibreak[m] + 1 for here
                            {
                                curve_indicator = TRUE;
                                break;
                            }
                        }

                        if (curve_indicator == TRUE)
                        {
                            for (m = 0; m < model_num; m++)
                            {
                                /* not adequate for this sub model */
                                if (sub_istart[m] - sub_last_ibreak[m] < LASSO_MIN)
                                {
                                    break;
                                }

                                tmp_fit_cft = (double *)calloc(MAX_NUM_C, sizeof(double));
                                tmp_dif_mag = (float *)calloc(MAX_NUM_C, sizeof(float));
                                status = auto_robust_fit2(sub_clrx[m], sub_clry[m], sub_last_ibreak[m],
                                                          sub_istart[m] - 1, MIN_NUM_C, tmp_fit_cft, &sub_rmse[m],
                                                          tmp_v_dif);
                                if (status != SUCCESS)
                                {
                                    RETURN_ERROR("Calling auto_robust_fit2 with enough observations\n",
                                                 FUNC_NAME, FAILURE);
                                }

                                for (i_conse = 1; i_conse < conse + 1; i_conse++)
                                {
                                    auto_ts_predict_singleband(sub_clrx[m], sub_fit_cft[m], MIN_NUM_C,
                                                               sub_istart[m] - i_conse, sub_istart[m] - i_conse,
                                                               &ts_pred_temp);
                                    tmp_dif_mag[i_conse - 1] = (float)(sub_clry[m][sub_istart[m] - i_conse] - ts_pred_temp);
                                }

                                quick_sort_float(tmp_dif_mag, 0, conse - 1);
                                rec_cg[*num_fc].magnitude[mid] = -median_1d_float(tmp_dif_mag, 0, conse - 1);
                                // rec_cg[*num_fc].magnitude[m] = -tmp_dif_mag[(int)(conse / 2)];

                                for (k = 0; k < MAX_NUM_C; k++)
                                {
                                    if (k < MIN_NUM_C)
                                        rec_cg[*num_fc].coefs[m][k] = sub_fit_cft[m][k];
                                    else
                                        rec_cg[*num_fc].coefs[m][k] = 0;
                                }

                                /******************************************/
                                /*                                        */
                                /* Record rmse of the pixel.              */
                                /*                                        */
                                /******************************************/
                                rec_cg[*num_fc].rmse[m] = sub_rmse[m];
                                rec_cg[*num_fc].num_obs[m] = sub_istart[m] - 1 - sub_last_ibreak[m] + 1;
                                free(tmp_fit_cft);
                                free(tmp_dif_mag);
                            }

                            rec_cg[*num_fc].t_end = sub_clrx[0][sub_istart[0] - 1];
                            rec_cg[*num_fc].t_break = sub_clrx[0][sub_istart[0]];
                            rec_cg[*num_fc].category = 10 + MIN_NUM_C;
                            rec_cg[*num_fc].change_prob = 100;
                            rec_cg[*num_fc].t_start = sub_clrx[mid][sub_last_ibreak[0]];
                            rec_cg[*num_fc].pos = pos;
                            *num_fc = *num_fc + 1;
                        }
                    }

                    /***********************************************************/
                    /*                                                         */
                    /*              look forward to save change records        */
                    /*                                                         */
                    /***********************************************************/
                    rec_cg[*num_fc].t_start = max_1d_float(tmp_sub_tstart, model_num, &max_id);
                    rec_cg[*num_fc].t_break = sub_clrx[mid][i_local];
                    rec_cg[*num_fc].category = MIN_NUM_C;
                    rec_cg[*num_fc].t_end = sub_clrx[mid][i_local - 1];
                    rec_cg[*num_fc].change_prob = 100;
                    rec_cg[*num_fc].pos = pos;
                    *num_fc = *num_fc + 1;

                    for (b = 0; b < model_num; b++)
                    {
                        quick_sort_float(v_dif_mag[b], 0, conse - 1);
                        rec_cg[*num_fc].magnitude[mid] = median_1d_float(v_dif_mag[b], 0, conse - 1);
                    }

                    /**********************************************/
                    /*                                            */
                    /* update all sub variables                   */
                    /*                                            */
                    /**********************************************/
                    globel_t_start = NA;
                    for (b = 0; b < MAX_VZA_GROUPS; b++)
                    {
                        sub_bl_train[b] = 0;
                        sub_istart[b] = sub_i[b] + 1;
                    }

                    /* break the looping for sub models */
                    break;
                }
                else if (identify_result == OUTLIER) /*false change*/
                {
                    /**********************************************/
                    /*                                            */
                    /* Remove noise.                              */
                    /*                                            */
                    /**********************************************/
                    for (k = i_ini; k < sub_i[mid]; k++)
                    {
                        sub_clrx[mid][k] = sub_clrx[mid][k + 1];
                        sub_clry[mid][k] = sub_clry[mid][k + 1];
                    }
                    sub_i[mid] = sub_i[mid] - 1;
                    sub_end[mid] = sub_end[mid] - 1;
                }
                sub_i[mid] = sub_i[mid] + 1;
            } // (sub_bl_train[mid] == 1)
        }     // for (mid = 0; mid < 2; mid++)
        i++;
    }

    /**************************************************************/
    /*                                                            */
    /* Two ways for processing the end of the time series.        */
    /*                                                            */
    /**************************************************************/
    curve_indicator = FALSE;
    for (mid = 0; mid < model_num; mid++)
    {
        id_last = conse;
        if (sub_bl_train[mid] == 1)
        {
            /**********************************************************/
            /*                                                        */
            /* If no break, find at the end of the time series,       */
            /* define probability of change based on adj_conse.           */
            /*                                                        */
            /**********************************************************/
            for (i_conse = conse - 1; i_conse >= 1; i_conse--)
            {
                identify_result = identifyChange(v_dif_mag[mid], min_rmse[mid], num_toler,
                                                 (double)t_cg, NSIGN, conse - i_conse, conse - 1,
                                                 T_MAX_CG_VSA);
                if (identify_result == NO_CHANGE)
                {
                    /**************************************************/
                    /*                                                */
                    /* The last stable ID.                            */
                    /*                                                */
                    /**************************************************/

                    id_last = i_conse + 1;
                    break;
                }
            }

            /**********************************************************/
            /*                                                        */
            /* Update change probability, end time of the curve.      */
            /*                                                        */
            /**********************************************************/
            if (conse > id_last)
            {
                /******************************************************/
                /*                                                    */
                /* Update magnitude of change.                        */
                /*                                                    */
                /******************************************************/

                quick_sort_float(v_dif_mag[mid], id_last, conse - 1);
                rec_cg[*num_fc].magnitude[mid] = median_1d_float(v_dif_mag[mid], id_last, conse - 1);
                sub_prob[mid] = (conse - id_last) * 1.0 / conse;
                sub_prob_tend[mid] = sub_clrx[mid][sub_end[mid] - 1 - conse + id_last];
                sub_prob_tbreak[mid] = sub_clrx[mid][sub_end[mid] - conse + id_last];
            }
            else
            {
                sub_prob[mid] = 0;
                sub_prob_tend[mid] = sub_clrx[mid][sub_end[mid] - 1];
                sub_prob_tbreak[mid] = 0;
            }

            if (curve_indicator == FALSE)
            {
                curve_indicator = TRUE;
            }
        }
        else if (sub_bl_train[mid] == 0)
        {
            if ((sub_end[mid] - sub_istart[mid]) >= LASSO_MIN) // 09/28/2018 SY delete equal sign //11/15/2018 put back equal sign
            {
                status = auto_robust_fit2(sub_clrx[mid], sub_clry[mid], sub_last_ibreak[mid],
                                          sub_istart[mid], MIN_NUM_C, tmp_fit_cft, &sub_rmse[mid],
                                          tmp_v_dif);
                if (status != SUCCESS)
                {
                    RETURN_ERROR("Calling auto_robust_fit2 with enough observations\n",
                                 FUNC_NAME, FAILURE);
                }

                /******************************************************/
                /*                                                    */
                /* Record fitted coefficients, change magnitudes      */
                /*                                                    */
                /******************************************************/

                for (k = 0; k < MAX_NUM_C; k++)
                {
                    rec_cg[*num_fc].coefs[mid][k] = sub_fit_cft[mid][k];
                }
                rec_cg[*num_fc].rmse[mid] = *tmp_v_dif;
                rec_cg[*num_fc].magnitude[mid] = 0.0;
            }

            /******************************************************/
            /*                                                    */
            /* guarantee only num_fc being plus by 1.             */
            /*                                                    */
            /******************************************************/

            if (curve_indicator == FALSE)
            {
                curve_indicator = TRUE;
                *num_fc = *num_fc + 1;
            }
            sub_prob[mid] = 0;
            sub_prob_tend[mid] = sub_clrx[mid][sub_end[mid] - 1];
            sub_prob_tbreak[mid] = 0;
        }

        for (m = 0; m < model_num; m++)
        {
            tmp_sub_tstart[m] = sub_clrx[m][sub_istart[m]];
        }

        rec_cg[*num_fc].t_start = max_1d_float(tmp_sub_tstart, model_num, &max_id);
        rec_cg[*num_fc].category = MIN_NUM_C;
        tmp = max_1d_float(sub_prob, model_num, &max_id);
        rec_cg[*num_fc].change_prob = (short int)(tmp * 100);
        rec_cg[*num_fc].t_end = sub_prob_tend[max_id];
        rec_cg[*num_fc].t_break = sub_prob_tbreak[max_id];
        rec_cg[*num_fc].pos = pos;
    }

    status = free_2d_array((void **)sub_clrx);
    if (status != SUCCESS)
    {
        RETURN_ERROR("Freeing memory: rec_v_dif\n",
                     FUNC_NAME, FAILURE);
    }
    status = free_2d_array((void **)sub_clry);
    if (status != SUCCESS)
    {
        RETURN_ERROR("Freeing memory: sub_clry\n",
                     FUNC_NAME, FAILURE);
    }
    status = free_2d_array((void **)sub_fit_cft);
    if (status != SUCCESS)
    {
        RETURN_ERROR("Freeing memory sub_fit_cft\n",
                     FUNC_NAME, FAILURE);
    }
    status = free_2d_array((void **)v_dif_mag);
    if (status != SUCCESS)
    {
        RETURN_ERROR("Freeing memory v_dif_mag\n",
                     FUNC_NAME, FAILURE);
    }

    free(id_good);
    free(clrx);
    free(clry);
    free(clrx_tmp);
    free(clry_tmp);
    free(clry_mid);
    free(tmp_v_dif);
    free(tmp_sub_tstart);
    free(min_rmse);

    // for debug
    // printf("free stage 6 \n");

    if (result == SUCCESS)
    {
        return (SUCCESS);
    }
    else
    {
        return (FAILURE);
    }
}