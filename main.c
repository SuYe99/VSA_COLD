#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include "defines.h"
#include "utilities.h"
#include "2d_array.h"
#include "output.h"
#include "vsa_cold.h"

int parse_string_byspace(char *str, int *out_list)
{
    char *token;
    int ref_break_count = 0;
    token = strtok(str, "|");
    out_list[ref_break_count] = atoi(token);
    ref_break_count++;
    token = strtok(NULL, "|");
    while (token != NULL)
    {
        out_list[ref_break_count] = atoi(token);
        ref_break_count++;
        token = strtok(NULL, "|");
    }
    return ref_break_count;
}

int main(int argc, char *argv[])
{
    int row_count = 0;
    char FUNC_NAME[] = "main"; /* For printing error messages           */
    FILE *interpret_file, *sample_file;
    int valid_scene_count = 0;
    int status;
    output_t_vsa *rec_cg;
    long *sdate;
    char *csv_row;
    long *buf_rad;    /* Lunar-BRDF-corrected DNB radiances */
    long *buf_tvza;   /* viewing zenith angles from the raw Black Marble product. */
    long *buf_tvaa;   /* viewing azimuth angles from the raw Black Marble product.  */
    long *buf_tmq;    /* mandatory quality flag from the lunar-brdf-corrected black marble product */
    long *buf_tm2;    /* cloud and snow mask summarized from the cloud mask of the raw Black Marble product. */
    long *buf_QFDNB;  /* quality flag from the Lunar-BRDF-corrected Black Marble product. */
    long *buf_tm_buf; /* cloud and snow buffer mask.  */
    int num_fc = 0;   /* the number of resultant curves */
    int result;
    int i, j, k;
    int id;
    char *date_start;
    char *date_end;
    char *chg_type;
    char *chg_conf;
    char filename[MAX_STR_LEN];
    char sample_path[MAX_STR_LEN];
    int denom_id;
    float tmp;
    float start_value, end_value, rad_bfr, rad_aft, chg_mag;
    int tmp_mapped_break[50];
    int tmp_ref_break_start[50];
    int tmp_ref_break_end[50];
    int tmp_map_count;
    int tmp_match_count;
    int sample_count = 0;
    int MAX_REF_BREAK_COUNT = 15;
    int date_start_list[MAX_REF_BREAK_COUNT];
    int date_end_list[MAX_REF_BREAK_COUNT];
    int chg_type_list[MAX_REF_BREAK_COUNT];
    int chg_conf_list[MAX_REF_BREAK_COUNT];
    int ref_break_count;
    int final_ref_break_count;
    int total_ref = 0;
    int total_map = 0;
    int total_match = 0;
    bool low_conf;

    /******************************************************/
    /*                                                    */
    /*     specify VSA_COLD parameters below              */
    /*                                                    */
    /******************************************************/
    double t_cg = 1.1503; // normal distribution value at 0.75. Matlab function norminv(1-0.25/2)
    int conse = 14;
    int num_toler = 1;
    int vsa_model_num = 3;
    int vsa_bin_edge[4] = {0, 20, 40, 60};
    int pos = 1; // assign 1 to a testing site

    /******************************************************/
    /*                                                    */
    /*     specify validation parameters below            */
    /*                                                    */
    /******************************************************/
    int ntl_lim = 1;
    int low_bound = 735235;   /* 1/1/2013 */
    int upper_bound = 738155; /* 12/31/2020 */
    float tlr_t = 365.25 * 0.5;

    /* python code to generate 'MATLAB' ordinal date. Note that it needs to add 366 offset days*/
    // import datetime as dt
    // import pandas as pd
    // pd.Timestamp.toordinal(dt.date(2020, 1, 1)) + 366

    sdate = (long *)malloc(MAX_SCENE_LIST * sizeof(int));
    buf_rad = (long *)malloc(MAX_SCENE_LIST * sizeof(long));
    buf_tmq = (long *)malloc(MAX_SCENE_LIST * sizeof(long));
    buf_tvza = (long *)malloc(MAX_SCENE_LIST * sizeof(long));
    buf_tvaa = (long *)malloc(MAX_SCENE_LIST * sizeof(long));
    buf_tm2 = (long *)malloc(MAX_SCENE_LIST * sizeof(long));
    buf_QFDNB = (long *)malloc(MAX_SCENE_LIST * sizeof(long));
    buf_tm_buf = (long *)malloc(MAX_SCENE_LIST * sizeof(long));
    rec_cg = (output_t_vsa *)malloc(NUM_FC * sizeof(output_t_vsa));
    csv_row = (char *)malloc(MAX_STR_LEN * sizeof(char));

    /* pixel mode */
    if (argv[1][0] == 'p')
    {
        sample_file = fopen(argv[2], "r");
        /******************************************************/
        /*                                                    */
        /*  initialize rec_cg by NA                           */
        /*                                                    */
        /******************************************************/
        for (i = 0; i < NUM_FC; i++)
        {
            rec_cg[i].t_start = NA;
            rec_cg[i].t_end = NA;
            rec_cg[i].t_break = NA;
            rec_cg[i].pos = NA;
            rec_cg[i].category = NA;
            for (j = 0; j < MAX_VZA_GROUPS; j++)
            {
                rec_cg[i].num_obs[j] = NA;
                for (k = 0; k < MIN_NUM_C; k++)
                {
                    rec_cg[i].coefs[j][k] = NA;
                }
                rec_cg[i].rmse[j] = NA;
                rec_cg[i].magnitude[j] = NA;
            }
            rec_cg[i].change_prob = NA;
        }
        num_fc = 0;
        valid_scene_count = 0;

        /******************************************************/
        /*                                                    */
        /*  read csv into 1d arrays                           */
        /*                                                    */
        /******************************************************/
        while (fgets(csv_row, 255, sample_file) != NULL)
        {
            if (row_count != HEADLINE) // we skip first line because it is a header
            {
                sdate[valid_scene_count] = (long)atoi(strtok(csv_row, ","));   // date list
                buf_rad[valid_scene_count] = (long)atoi(strtok(NULL, ","));    // Lunar-BRDF-corrected DNB radiances.
                buf_tmq[valid_scene_count] = (long)atoi(strtok(NULL, ","));    // mandatory quality flag from the Lunar-BRDF-corrected Black Marble product
                buf_tvza[valid_scene_count] = (long)atoi(strtok(NULL, ","));   // viewing zenith angles from the raw Black Marble product
                buf_tvaa[valid_scene_count] = (long)atoi(strtok(NULL, ","));   // viewing azimuth angles from the raw Black Marble product
                buf_tm2[valid_scene_count] = (long)atoi(strtok(NULL, ","));    // cloud and snow mask summarized from the cloud mask of the raw Black Marble product
                buf_QFDNB[valid_scene_count] = (long)atoi(strtok(NULL, ","));  // quality flag from the Lunar-BRDF-corrected Black Marble product
                buf_tm_buf[valid_scene_count] = (long)atoi(strtok(NULL, ",")); // cloud and snow buffer mask
                valid_scene_count++;
            }
            row_count++;
        }
        result = vsa_cold_detect(sdate, buf_rad, buf_tmq, buf_tvza, buf_tvaa, buf_tm2, buf_QFDNB,
                                 buf_tm_buf, vsa_bin_edge, vsa_model_num, valid_scene_count, pos, conse,
                                 t_cg, num_toler, &num_fc, rec_cg);

        fclose(sample_file);
        for (k = 0; k < num_fc; k++)
        {
            printf("break date for %dth break is %d\n", (k + 1), rec_cg[k].t_break);
        }
    }
    else if (argv[1][0] == 'a') /* accuracy assessment mode */
    {
        interpret_file = fopen(argv[2], "r");

        sample_count = 0;
        int row_count = 0;
        while (fgets(csv_row, 255, interpret_file) != NULL)
        {
            if (row_count == 0)
            {
                row_count++;
                continue;
            }
            /*********************************************/
            /*           read interpretation data        */
            /*********************************************/
            id = (int)atoi(strtok(csv_row, ","));
            date_start = strtok(NULL, ",");
            date_end = strtok(NULL, ",");
            chg_type = strtok(NULL, ",");
            chg_conf = strtok(NULL, ",");
            low_conf = FALSE;
            final_ref_break_count = 0;
            if (date_start[0] != 'N') // nA
            {                         /* stable */
                for (i = 0; i < MAX_REF_BREAK_COUNT; i++)
                {
                    date_start_list[i] = 0;
                    date_end_list[i] = 0;
                    chg_type_list[i] = 0;
                    chg_conf_list[i] = 0;
                }

                ref_break_count = parse_string_byspace(date_start, date_start_list);
                ref_break_count = parse_string_byspace(date_end, date_end_list);
                ref_break_count = parse_string_byspace(chg_type, chg_type_list);
                ref_break_count = parse_string_byspace(chg_conf, chg_conf_list);

                for (i = 0; i < ref_break_count; i++)
                {
                    if (chg_conf_list[i] < 2)
                    {
                        final_ref_break_count = 0;
                        low_conf = TRUE;
                        break;
                    }
                    if ((date_start_list[i] < low_bound) || (date_end_list[i] > upper_bound || (chg_type_list[i] == 7) || (chg_type_list[i] == 8)))
                    {
                        continue;
                    }
                    tmp_ref_break_start[final_ref_break_count] = date_start_list[i];
                    tmp_ref_break_end[final_ref_break_count] = date_end_list[i];
                    final_ref_break_count++;
                }
            }
            if (low_conf == TRUE)
            {
                continue; // we didn't process low confidence
            }
            total_ref = total_ref + final_ref_break_count;

            /*********************************************/
            /*           run pixel-based vsa-cold        */
            /*********************************************/
            for (i = 0; i < NUM_FC; i++)
            {
                rec_cg[i].t_start = NA;
                rec_cg[i].t_end = NA;
                rec_cg[i].t_break = NA;
                rec_cg[i].pos = NA;
                rec_cg[i].category = NA;
                for (j = 0; j < MAX_VZA_GROUPS; j++)
                {
                    rec_cg[i].num_obs[j] = NA;
                    for (k = 0; k < MIN_NUM_C; k++)
                    {
                        rec_cg[i].coefs[j][k] = NA;
                    }
                    rec_cg[i].rmse[j] = NA;
                    rec_cg[i].magnitude[j] = NA;
                }
                rec_cg[i].change_prob = NA;
            }
            num_fc = 0;
            valid_scene_count = 0;

            sprintf(filename, "sample_TS_%03d.csv", id);
            sprintf(sample_path, "%s/%s", argv[3], filename);
            sample_file = fopen(sample_path, "r");
            row_count = 0;
            while (fgets(csv_row, 255, sample_file) != NULL)
            {
                if (row_count != HEADLINE) // we skip first line because it is a header
                {
                    sdate[valid_scene_count] = (long)atoi(strtok(csv_row, ","));   // date list
                    buf_rad[valid_scene_count] = (long)atoi(strtok(NULL, ","));    // Lunar-BRDF-corrected DNB radiances.
                    buf_tmq[valid_scene_count] = (long)atoi(strtok(NULL, ","));    // mandatory quality flag from the Lunar-BRDF-corrected Black Marble product
                    buf_tvza[valid_scene_count] = (long)atoi(strtok(NULL, ","));   // viewing zenith angles from the raw Black Marble product
                    buf_tvaa[valid_scene_count] = (long)atoi(strtok(NULL, ","));   // viewing azimuth angles from the raw Black Marble product
                    buf_tm2[valid_scene_count] = (long)atoi(strtok(NULL, ","));    // cloud and snow mask summarized from the cloud mask of the raw Black Marble product
                    buf_QFDNB[valid_scene_count] = (long)atoi(strtok(NULL, ","));  // quality flag from the Lunar-BRDF-corrected Black Marble product
                    buf_tm_buf[valid_scene_count] = (long)atoi(strtok(NULL, ",")); // cloud and snow buffer mask
                    valid_scene_count++;
                }
                row_count++;
            }
            result = vsa_cold_detect(sdate, buf_rad, buf_tmq, buf_tvza, buf_tvaa, buf_tm2, buf_QFDNB,
                                     buf_tm_buf, vsa_bin_edge, vsa_model_num, valid_scene_count, pos,
                                     conse, t_cg, num_toler, &num_fc, rec_cg);

            fclose(sample_file);

            tmp_map_count = 0;
            tmp_match_count = 0;
            for (i = 0; i < num_fc - 1; i++)
            {
                rad_bfr = 999999;
                rad_aft = 999999;
                tmp = -999999;
                denom_id = 0;
                if ((rec_cg[i].change_prob < 100) || (rec_cg[i].t_break == 0))
                {
                    continue;
                }

                if ((rec_cg[i].t_break > upper_bound) || (rec_cg[i].t_break < low_bound))
                {
                    continue;
                }

                /* get the 'break-denominant' model id*/
                for (j = 0; j < (vsa_model_num + 1); j++)
                {
                    if (rec_cg[i].magnitude[j] > tmp)
                    {
                        denom_id = j;
                        tmp = rec_cg[i].magnitude[j];
                    }
                }

                start_value = 0.1 * (rec_cg[i].coefs[denom_id][0] + rec_cg[i].coefs[denom_id][1] * rec_cg[i].t_start);
                end_value = 0.1 * (rec_cg[i].coefs[denom_id][0] + rec_cg[i].coefs[denom_id][1] * rec_cg[i].t_end);
                rad_bfr = fmax(start_value, end_value);

                start_value = 0.1 * (rec_cg[i + 1].coefs[denom_id][0] + rec_cg[i + 1].coefs[denom_id][1] * rec_cg[i + 1].t_start);
                end_value = 0.1 * (rec_cg[i + 1].coefs[denom_id][0] + rec_cg[i + 1].coefs[denom_id][1] * rec_cg[i + 1].t_end);
                rad_aft = fmax(start_value, end_value);

                chg_mag = 0.1 * rec_cg[i].magnitude[denom_id];

                if ((rad_bfr < ntl_lim) && (rad_aft < ntl_lim) && (chg_mag < ntl_lim))
                    continue;

                tmp_mapped_break[tmp_map_count] = rec_cg[i].t_break;
                tmp_map_count++;
            }
            total_map = total_map + tmp_map_count;

            tmp_match_count = 0;
            for (i = 0; i < tmp_map_count; i++)
            {
                for (j = 0; j < final_ref_break_count; j++)
                {
                    if (((tmp_ref_break_start[j] - tlr_t) <= tmp_mapped_break[i]) && (tmp_mapped_break[i] <= (tmp_ref_break_end[j] + tlr_t)))
                    {
                        tmp_match_count = tmp_match_count + 1;
                        break;
                    }
                }
            }
            total_match = total_match + tmp_match_count;
            printf("total_match  is %d\n", total_match);
            printf("total_map is %d\n", total_map);

            sample_count++;
            row_count++;
            printf("processing finished for sampleid = %d\n", id);
        }
        fclose(interpret_file);

        float omission = (total_ref - total_match) / (total_ref * 1.0);
        float commission = (total_map - total_match) / (total_map * 1.0);
        float f1 = 2 * (1 - omission) * (1 - commission) / (2 - commission - omission);

        printf("Match count is  %d\n", total_match);
        printf("Reference count is  %d\n", total_ref);
        printf("Mapped count is  %d\n", total_map);
        printf("Omission is  %f\n", omission);
        printf("Commission count is  %f\n", commission);
        printf("f1 is  %f\n", f1);
    } // else if (argv[1][0] == 'a')

    free(csv_row);
    free(sdate);
    free(buf_rad);
    free(buf_tvza);
    free(buf_tvaa);
    free(buf_tmq);
    free(buf_tm2);
    free(buf_QFDNB);
    free(buf_tm_buf);
    free(rec_cg);
}