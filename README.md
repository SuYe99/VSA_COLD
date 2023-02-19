# VSA_COLD
Continuous monitoring of nighttime light changes based on daily NASA's Black Marble product suite

### Author: Su Ye (remotesensingsuy@gmail.com)

## 1. Preparation
TBD


## 2. Installation
TBD


## 3. Usage
```
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
    int vsa_model_num, /* I: bin number of vsa angle, such as 3  */
    int buf_len,       /* I: number of valid scenes  */
    int pos,           /* I: the position id of pixel */
    int conse,         /* I: consecutive obs number for break identification */
    double t_cg,       /* I: change probaility threshold for break */
    int num_toler,     /* I: outlier number for tolerance */
    int *num_fc,       /* O: number of fitting curves                       */
    output_t_vsa *rec_cg /* O: outputted structure for VSA_COLD results    */);
```
Please check [the example file](main.c) for the details.

## 4. Citation
If you use this repo, please read/cite the publication [VSA_COLD](https://www.sciencedirect.com/science/article/pii/S0034425722003753):

Li, T., Zhu, Z., Wang, Z., Román, M. O., Kalb, V. L., & Zhao, Y. (2022). Continuous monitoring of nighttime light changes based on daily NASA's Black Marble product suite. Remote Sensing of Environment, 282, 113269.