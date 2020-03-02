Automated Household Water End-Use Disaggregation Through A Rule-Based Automated Methodology
Filippo Mazzoni, Stefano Alvisi, Marco Franchini, Marco Ferraris, and Zoran Kapelan

The folder attached in the repository includes the content briefly described in the following:
1. Water use data are included in a Microsoft Excel® datasheet for each household (namely H1data.xlsx, H2data,xslx, H3data.xlsx, H4data.xlsx). Collected data are grouped in tables, where: the first column includes date-time references (format dd/MM/yyyy HH.mm); the second column includes the aggregate water use time series (L/min); columns 3-8 include the water use time series (L/min) for each end-use category (i.e. dishwasher DW, kitchen sink KS, washing machine WM, shower S, bathroom taps BT, toilet T).
2. The code of the automated methodology for water end-use disaggregation is included in a MATLAB® m-file (namely AutomatedDisaggregationTool.m). The code is provided with comments pointing out the idea behind the method.
3. The functions required by the main code are included in several MATLAB® m-file (namely F_Parameters.m, F_ToiletParameterValuesIdentification, F_DW, F_KS, F_WM, F_S, F_TT). Functions are provided with comments pointing out the idea behind them and are recalled in turn during the main code run.

