gcan
====

## Usage

gcan (Growth Curve ANalyzer) helps you to quickly analyze growth curve data. 
It expects your growth curves to be saved in csv (with "*.csv" file ending) or text ("*.txt") format. With the
following layout:

-----------------------------------------------
| Time | Sample 1 | Sample 2 | ... | Sample n |
-----------------------------------------------

For sample repetitions simply give the columns the same name.

The script is run with
```bash
./gcan.R data_file
```

The script will generate some PDF plots with the curves + error bars and the groeth rates and area under the curve as 
boxplots. The results are also generated and put into a results.txt file which has the individual data for each curve
as well as additional information regarding the linear models generated to approximate the growth rates such as
F-test p-value, number of points used for regrssion and so on.

The "linear parts" PDF is a verification plot to see whether the cutoffs for the linear parts of the log-OD curves are
well chosen. If the curves do not behave linear within the grey rectangle you should ajust the OD.min and OD.max parameters
in the script to your data. The given values usually work well for yeast growth curves with initial ODs around 0.2.

## Installation

Just copt the script to any location and make it executable with `chmod +x gcan.R`.

## Requirements

The following R packages have to be installed

* ggplot2
* reshape2
