# Water end-use disaggregation methodology
The developed methodology for automated water use disaggregation detects, disaggregates and classifies individual water uses one end-use category at a time. Specifically, the disaggregation of the aggregate water use time series collected at the domestic water inlet point is done by means of a set of functions. Such functions are applied in an order, starting with functions aimed at detecting end-uses that, given their very nature, are generally more regular in terms of water use, and ending with the most irregular ones. At first, electronic appliances use is investigated through function dishwasher (_F_DW_) and function washing machine (_F_WM_). Then, shower uses are searched through function shower (_F_S_). At last, toilet and tap uses (which are generally less-recognisable at one-minute resolution or often overlapping in time with other uses) are detected and classified by means of the function tap and toilet (_F_TT_). In conclusion, the water use time series for each end-use category is available. Operatively, the automated methodology for water use disaggregation was developed by using MATLAB® programming software and consists of a main code (_AutomatedDisaggregationTool_), where the Microsoft Excel® datasheet including the collected aggregate water use is loaded and functions for water use disaggregation are applied in turn.

**Content of the folder**
- Water use data are included in a Microsoft Excel® datasheet for each household (namely _H1data.xlsx_, _H2data,xslx_, _H3data.xlsx_, _H4data.xlsx_). Collected data are grouped in tables, where: the first column includes date-time references (format dd/MM/yyyy HH.mm); the second column includes the aggregate water use time series (L/min); columns 3-8 include the water use time series (L/min) for each end-use category (i.e. dishwasher DW, kitchen sink KS, washing machine WM, shower S, bathroom taps BT, toilet T).
- The main code of the automated methodology for water end-use disaggregation is included in a MATLAB® m-file (namely _AutomatedDisaggregationTool.m_). The code is provided with comments pointing out the idea behind the method.
- The function including household specific and general parameters (_F_Parameters.m_), the function including the rules for the automatic identification of toilet parameter values if the specific ones are not available (_F_ToiletParameterValuesIdentification.m_) and the functions including the rules for dishwasher (_F_DW.m_), washing machine (_F_WM.m_), shower (_F_S.m_) and toilet-tap (_F_TT.m_) water use disaggregation are included in several MATLAB® m-file. Functions are provided with comments pointing out the idea behind them and are recalled in turn during the main code run.

**Requirements**
The disaggregation code and its related functions have been developed on MATLAB® R2019a.

**How to run the tool**
The main MATLAB® file to run the automated disaggregation tool is _AutomatedDisaggregationTool.m_. The following input settings are required to the user, after launching the run:
- Household whose water use has to be automatically disaggregated into end-uses (1 for household H1, 2 for household H2, 3 for household H3, 4 for household H4) .
- Period of analysis (1 for calibration period, 2 for validation period).
- Water use parameter values (1 for household H1 specific values, 2 for household H2 specific values, 3 for household H3 specific values, 0 for general values).

**Results provided**
- Volumes automatically disaggregated into end-uses and their corresponding baseline volumes (Command Window).
- Disaggregation performance metrics, i.e. Water Contribution Accuracy (_WCA_) and Normalized Root-Mean-Square Error (_NRMSE_) introduced by Cominola et al. (2018) (Command Window).
- End-use time series obtained through the automated disaggregation, i.e. vector _DW_ for dishwasher, _WM_ for washing machine, _S_ for shower, _KS_ for kitchen sink, _BT_ for bathroom taps, _T_ for toilet. The length of each vector equals the length of the selected period of analysis, in minutes (Workspace).

**Authors**
- **Filippo Mazzoni**, **Stefano Alvisi**, **Marco Franchini** - Department of Engineering | University of Ferrara | Ferrara, Italy
- **Marco Ferraris** - Department for Sustainability | Italian National Agency for New Technologies, Energy and Sustainable Economic Development (ENEA) | Bologna, Italy
- **Zoran Kapelan** - Department of Water Management | Delft University of Technology | Delft, Netherlands

**References**
- Cominola, A., Giuliani, M., Castelletti, A., Rosenberg, D. E., and Abdallah, A. M. (2018). “_Implications of data sampling resolution on water use simulation, end-use disaggregation, and demand management._” Environmental Modelling & Software, 102(Apr), 199–212.
