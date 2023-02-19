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

int main(int argc, char *argv[])
{
    int row_count = 0;
    char FUNC_NAME[] = "main"; /* For printing error messages           */
    FILE *sample_file;
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

    /******************************************************/
    /*                                                    */
    /*     specify VSA_COLD parameters here               */
    /*                                                    */
    /******************************************************/
    double t_cg = 1.1503; // normal distribution value at 0.75. Matlab function norminv(1-0.25/2)
    int conse = 14;
    int num_toler = 1;
    int vsa_model_num = 3;
    int vsa_bin_edge[4] = {0, 20, 40, 60};
    int pos = 1; // assign 1 to a testing site

    sample_file = fopen(argv[2], "r");

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