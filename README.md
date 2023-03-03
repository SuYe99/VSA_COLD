# VSA_COLD
Continuous monitoring of nighttime light changes based on daily NASA's Black Marble product suite

### Author: Su Ye (remotesensingsuy@gmail.com)

## 1. Preparation
The GSL libraries are required.

For Ubuntu/Debian systems, GSL can be installed via:
```
sudo apt-get update
sudo apt-get install build-essential  -y
sudo apt-get install libgsl-dev -y
```

## 2. Installation
TBD

## 3. Usage

The main function 'vsa_cold_detect' is in [vsa_cold.c](vsa_cold.c). I used the pixel-based [example csv](/data/TrialData_StablePixel.csv) Tian provided. Please inquire Tian for the data preprocessing codes for the inputs.

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

[1] C version outputs  change magntidue for each band, instead of one change magnitude for single 'feature band' in original matlab (i.e., the band detects change). The intention to do that is to outputted change magnitudes for multiple VSA angles which should be informative for indicating change types for future;

[2] the probability now is between 0 and 100. I removed the digit that the Matlab version used to indicate the 'feature band', as the feature band could be indicated by 'change magnitudes'

[3] The order of model parameters (such as magnitudes, rmse) follows the order 'all_model, [0,20], [20,40], [40, 60]'

## 4. Tests

I gave two testing mode, **P**ixel-based and  **A**ccuracy-assessment mode, in [main.c](main.c): 

If you are interested in a pixel-based csv test (P mode), please run:

```
./vsa_cold p /Users/coloury/Dropbox/UCONN/sample_TS_csv/sample_TS_026.csv
```

If you are interested in the accuracy assessment (A mode), please run:

```
./vsa_cold a /Users/coloury/Dropbox/Documents/VSA_COLD/data/sample_interpretation_clean.csv /Users/coloury/Dropbox/UCONN/sample_TS_csv
```

Where 'sample_interpretation_clean.csv' is the interpretation table Tian provided, and 'sample_TS_csv' is the folder that saves original pixel-based inputs as tables (ask me if you need 'sample_TS_csv').

The current accuracy is 37.7% omission, 41.0% comission, 0.606 F1 (testing date: 03-03-2023, by Su Ye)

## 5. Citation
If you use this repo, please read/cite the publication [VSA_COLD](https://www.sciencedirect.com/science/article/pii/S0034425722003753):

Li, T., Zhu, Z., Wang, Z., Román, M. O., Kalb, V. L., & Zhao, Y. (2022). Continuous monitoring of nighttime light changes based on daily NASA's Black Marble product suite. Remote Sensing of Environment, 282, 113269.