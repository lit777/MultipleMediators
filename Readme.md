# To link all different data sources together

'Master.RData' is the main dataset for the analysis, which links EGUs data to monitoring sites data.

# Explain the linkage strategy

We obtained data on 400 coal-fired power plants from the EPA Air Markets Program Data and the Energy Information Administration, including plant characteristics, emissions control technologies installed (if any), and emissions of SO2, NOx and CO2.  We also obtain annual average ambient PM2.5 concentrations from over 3000 pollution monitoring stations in the EPA Air Quality System (AQS).  Basic meteorologic conditions such as temperature and barometric pressure are also obtained from AQS.  We use the locations of the power plants and the ambient monitoring stations (defined using latitude and longitude) to link each power plant to all ambient monitoring stations within a 150km radius.  Monitors located within 150km of more than one plant are linked to the closest plant, and remaining plants with no linked monitors are discarded.  The 150km range was chosen both to acknowledge that atmospheric processes carry power plant emissions across distances at least this great, but also to minimize the number of monitoring stations considered within range of more than one power plant. Thus, our investigation is confined to causal effects of SO2 scrubbers on PM2.5 within a 150km range.    

This geographical linkage yields an analysis data set consisting of 249 coal-fired power plants, each linked to at least one ambient PM2.5 monitoring station.  We conduct our analysis using annual data from the year 2005 (5 years after initiation of Phase II of the ARP).  A power plant can consist of multiple EGUs, each of which may or may not be equipped with SO2 scrubber. We regard any power plant to be "treated" with an SO2 scrubber in 2005 if at least 10% of the total heat input for that power plant can be attributed to EGUs that have SO2 scrubbers installed as of January 2005. It can be shown that the vast majority of facilities had nearly all or nearly none of their heat input attributed to EGUs with SO2 scrubbers, indicating robustness to this 10% cutoff.

For each power plant, we also obtain the annual average ambient PM2.5 concentration measured across all monitors within a 150km radius during 2005.  We similarly obtain average temperature and barometric pressure in the surrounding 150km area during 2005. 

# MCMC running

    - Run 'Model_BNP.R' script for obtaining posterior samples.
    - 'MCMC.RData' includes all posterior samples.
    
# Post-processing

   - To estimate the causal mediation effects, run 'Mediation_Summary.R' script.
   - To estimate the principal causal effects, run 'PCE_Summary.R' script.

# Generating outputs
   - To get the outputs, run 'Outputs.R' script.
   - To get surface plots, run 'Finite_Sample.R', 'SurfaceFrame1.R' scripts and then run 'SurfaceFrame2.R' script
      Finally, run 'SurfacePlot.R' script.
   - Running 'MapPlot.R' produces two US maps with treated sites and untreated sites.
