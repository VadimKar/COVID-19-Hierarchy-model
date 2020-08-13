# COVID-19-Hierarchy-model
Simulation and model-fitting code for manuscript by Karatayev, Anand, and Bauch.

covidHier0.4Github contains functions to fit and simulate the model

covidHierData.rds contains the following objects:

•	Distd50: the proportion of people voluntarily distancing by day 21 of the outbreak, data from https://news.gallup.com/opinion/gallup/298310/americans-step-social-distancing-even- 555 further.aspx

•	testRat: estimate of the ratio of true infections to positive tests made by A Lachmann, Correcting under-reported covid-19 case numbers. medRxiv (2020).

•	Msave: Travel matrix, where entries are numbers of commuters between a home and workplace county. Data from Statistics canada, 2016 census, catalogue no. 98-400-x2016391.

•	storeSatIncTrm.65: Pre-calculated values of baseline transmission probability (across values of c and xi) which leads to 65% of population being infected after 1 year without mitigation

•	caseCtBin: the daily number of reported cases, with smaller counties grouped by population density; data from https://www.publichealthontario.ca/en/data-and- 511 analysis/infectious-disease/covid-19-data-surveillance/covid-19-data-tool

•	testT3.1: function giving the (smoothed) testing intensity on any given day t from the day the 50th case positive case was reported (i.e., t=0 if fewer than 50 cumulative cases reported). These values increase from 0 to 1following the (smoothed) relative increase in number of tests taken each day based on data from https://www.ontario.ca/page/2019-novel- 546 coronavirussection-0

•	CountyPopDensities: population densities in each county, data from: https://www12.statcan.gc.ca/health-sante/82-228/details/page.cfm?Lang=E&Tab=1&Geo1=HR&Code1=3570&Geo2=PR&Code2=35&Data=Rate&SearchText=York%20Regional%20Health%20Unit&SearchType=Contains&SearchPR=01&B1=All&Custom=&B2=All&B3=All
