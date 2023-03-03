#ifndef OUTPUT_H
#define OUTPUT_H
#include <stdbool.h>
#include "defines.h"

/* Structure for the 'output' data type */
typedef struct
{
    int t_start;        /* time when series model gets started */
    int t_end;          /* time when series model gets ended */
    int t_break;        /* time when the first break (change) is observed */
    int pos;            /* the location of each time series model */
    short int category; /* the quality of the model estimation (what model
                     is used, what process is used) */
    /*
    The current category in output structure:
    first digit:
    0: normal model (no change)
    1: change at the beginning of time series model
    second digit:
    4: model has 3 coefs + 1 const */
    int num_obs[MAX_VZA_GROUPS]; /* the number of "good" observations used for model
                            estimation; the first element is all_model */
    short int change_prob;       /* the probability of a pixel that have undergone
                                    change (between 0 and 100) */
    double coefs[MAX_VZA_GROUPS][MIN_NUM_C];
    /*  coefficients for each time series model for each
        spectral band*/
    float rmse[MAX_VZA_GROUPS];
    /*  RMSE for each time series model for each
        spectral band*/
    float magnitude[MAX_VZA_GROUPS]; /* the magnitude of change (difference between model
                               prediction and observation for each spectral band)*/
} output_t_vsa;

#endif
