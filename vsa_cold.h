#include "output.h"

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
    long *sdate,         /* I:  date as matlab serial date form (counting from Jan 0, 0000). Note ordinal date in python is from (Jan 1th, 0001) */
    long *buf_rad,       /* I: Lunar-BRDF-corrected DNB radiances */
    long *buf_tmq,       /* I: mandatory quality flag from the lunar-brdf-corrected black marble product */
    long *buf_tvza,      /* I:  viewing zenith angles from the raw Black Marble product. */
    long *buf_tvaa,      /* I: viewing azimuth angles from the raw Black Marble product.  */
    long *buf_tm2,       /* I: cloud and snow mask summarized from the cloud mask of the raw Black Marble product. */
    long *buf_QFDNB,     /* I: quality flag from the Lunar-BRDF-corrected Black Marble product. */
    long *buf_tm_buf,    /* I: cloud and snow buffer mask.  */
    int *vsa_bin_edge,   /* I: the bin edge of defined vsa groups, such as [0, 20, 40, 60]  */
    int vsa_model_num,   /* I: bin number of vsa angle, such as 3  */
    int buf_len,         /* I: number of valid scenes  */
    int pos,             /* I: the position id of pixel */
    int conse,           /* I: consecutive obs number for break identification */
    double t_cg,         /* I: change probaility threshold for break */
    int num_toler,       /* I: outlier number for tolerance */
    int *num_fc,         /* O: number of fitting curves                       */
    output_t_vsa *rec_cg /* O: outputted structure for VSA_COLD results    */
);