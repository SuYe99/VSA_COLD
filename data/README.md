# Definitions of variables in the data

## 1. sample_interpretation_clean.csv
**date_start**: Julian start date(s) for the reference change event(s). (In cases where multiple change events occurred from 2013 to 2021, multiple dates will be provided.)

**date_start**: Julian end date(s) for the reference change event(s).

**chg_type**: nighttime light change type(s) for the reference change event(s) (0=no change; 1=abrupt recovered change, 2=cross-year abrupt recovered change, 3=abrupt unrecovered change, 4=cross-year abrupt unrecovered change, 5=transit change, 6=cross-year transit change).

**chg_conf**: confidence level(s) for the manual interpretation result, samples with low confidence change event will be excluded in the accuracy assessment process (1=low; 2=medium; 3=high).


## 2. sample_TS_csv/TrialData_StablePixel.csv
**sdate**: Julian observation date.

**line_t**: Lunar-BRDF-corrected VIIRS DNB radiance.

**line_tmq**: Mandatory quality flag.

**line_tvza**: viewing zenith angle.

**line_tvaa**: viewing azimuth angle.

**line_tm2**: cloud and snow mask (0=clear, 1=cloud, 3=snow).

**line_QFDNB**: DNB quality flag.

**line_tm_buf**: cloud and snow buffer mask with a 5x5 buffer window.

