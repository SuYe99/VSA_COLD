# VSA_COLD
Continuous monitoring of nighttime light changes based on daily NASA's Black Marble product suite

### Author: Su Ye (remotesensingsuy@gmail.com)

## 1. Preparation
The GSL libraries are required.

For Ubuntu/Debian systems, GSL can be installed via:
```
sudo apt-get update
sudo apt-get install libgsl-dev -y
```

## 2. Installation
There is a variety of compilation ways. For cmake, you could run
```
/usr/local/bin/cmake --no-warn-unused-cli -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_BUILD_TYPE:STRING=Debug -DCMAKE_C_COMPILER:FILEPATH=/usr/local/bin/gcc-7 -DCMAKE_CXX_COMPILER:FILEPATH=/usr/local/bin/g++-7 -S/Users/coloury/Dropbox/Documents/VSA_COLD -B/Users/coloury/Dropbox/Documents/VSA_COLD/build -G "Unix Makefiles"
```

## 3. Usage

The main function 'vsa_cold_detect' is in [vsa_cold.c](vsa_cold.c). I used the pixel-based [example csv](/data/TrialData_StablePixel.csv) Tian provided. Please inquire Tian for the data preprocessing codes for the inputs.

```
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
    int date_start,      /* I: the low bound of dates to be processed, inclusive. Set 0 if you don't have low bound to set */
    int date_end,        /* I: the low bound of dates to be processed, inclusive. Set 999999, if you don't have upper bound to set */
    int *num_fc,         /* O: number of fitting curves                       */
    output_t_vsa *rec_cg /* O: outputted structure for VSA_COLD results    */
);
```
The output 'rec_cg' is a C structure that stores temporal information.  The rec_cg structure is defined as (see [output.h](output.h)):

```
/* Structure for the 'output' data type */
typedef struct
{
    int t_start;        /* time when series model gets started */
    int t_end;          /* time when series model gets ended */
    int t_break;        /* time when the first break (change) is observed */
    int pos;            /* the location of each time series model */
    short int category; /* the quality of the model estimation (what model
                     is used, what process is used) */
    int num_obs[MAX_VZA_GROUPS]; /* the number of "good" observations used for model
                            estimation; the first element is all_model */
    short int change_prob;       /* the probability of a pixel that have undergone
                                    change (between 0 and 100) */
    double coefs[MAX_VZA_GROUPS][MIN_NUM_C]; /*  coefficients for each time series model for each spectral band*/
    float rmse[MAX_VZA_GROUPS]; /*  RMSE for each time series model for each spectral band*/
    float magnitude[MAX_VZA_GROUPS]; /* the magnitude of change (difference between model
                               prediction and observation for each spectral band)*/
} output_t_vsa;
```

Note that I have made some modifications on the original matlab outputs:

[1] The C version now outputs change magntidue for each band, instead of one change magnitude for only the 'feature band' (i.e., the band detects change) in original matlab version. The intention to do this is to output change magnitudes for multiple VSA angles which should be informative for indicating change types for future;

[2] the probability now is between 0 and 100. I removed the digit that the Matlab version used to indicate the 'feature band', as the feature band could be indicated by the maximum 'change magnitudes'. 

[3] The order of model parameters now (such as magnitudes, rmse) follows the order 'all_model, [0,20], [20,40], [40, 60]', not '[0,20], [20,40], [40, 60], all_model' in the matlab

## 4. Efficiency (what if it is too slow?)

The original VSA_COLD is 50-100 slower than the COLD algorithm. There two parameters in [defines.h](defines.h) you may consider to adjust for an efficiency improvement:

**UPDATE_FREQ**:  the model update interval. The default is 1. You may consider to change it to 3, which could reach 3 times efficiency boost without hurting accuracy (see the below testing results)

**ROBUST_FIT_TIMES**: the iteration times of fitting for the robust fitting algorithm. The default is 50. I haven't tested it yet but it seems that it could be decreased at no cost of accuracy.

## 5. Tests

I gave two testing mode, **P**ixel-based and  **A**ccuracy-assessment mode, in [main.c](main.c): 

If you are interested in a pixel-based csv test (P mode), please run:

```
./vsa_cold p /Users/coloury/Dropbox/UCONN/sample_TS_csv/sample_TS_026.csv
```

If you are interested in the accuracy assessment (A mode), please run:

```
./vsa_cold a /Users/coloury/Dropbox/Documents/VSA_COLD/data/sample_interpretation_clean.csv /Users/coloury/Dropbox/UCONN/sample_TS_csv
```

Where 'sample_interpretation_clean.csv' is the interpretation table Tian provided, and 'sample_TS_csv' is the folder that saves 594 original pixel-based input tables (ask me if you need 'sample_TS_csv').

The current accuracies under different UPDATE_FREQs:

**UPDATE_FREQ = 1 (C language)**: 35.1% omission, 36.3% commission, 0.642 F1

**UPDATE_FREQ = 3 (C language)**: 33.7% omission, 37.2% commission, 0.645 F1

**UPDATE_FREQ = 5 (C language)**: 33.7% omission, 40.1% commission, 0.630 F1

**UPDATE_FREQ = 1     (Matlab)**: 33.3% omission, 31.4% commission, 0.676 F1 (Reported by Tian)

(testing date: 03-05-2023, by Su Ye)

## 6. Citation
If you use this repo, please read/cite the publication [VSA_COLD](https://www.sciencedirect.com/science/article/pii/S0034425722003753):

Li, T., Zhu, Z., Wang, Z., Román, M. O., Kalb, V. L., & Zhao, Y. (2022). Continuous monitoring of nighttime light changes based on daily NASA's Black Marble product suite. Remote Sensing of Environment, 282, 113269.