## R code of environmental data (temperature, ndvi, snow) manipulation as described in vanGils et al. XXX

### Temperature Data

Daily mean temperatures were obtained for 1983-2013 across the entire breeding range of the Calidris c. canutus subspecies (northern Taimyr Peninsula as defined by Lappo et al. 2002). As a first step, daily temperature data (May-August) were downloaded from the NOAA National Climatic Data Centre (http://www.ncdc.noaa.gov/oa/climate/climatedata.html#daily), which were then plotted as temperature surface maps3, using ‘gravity’ as the interpolation algorithm, while taking a search radius of 500 km and a maximum of ten weather stations. Subsequently, surface maps were overlaid on the knot’s breeding range, and surface values were averaged across this range, yielding a mean temperature for each day. Next, for each year separately a quadratic model was fitted to these daily mean temperatures (using the lm function in R), and the date at which the increasing part of this fit reached 0 °C was defined as D.T0.

Lappo, E. G., Tomkovich, P. S. & Syroechkovskiy, E. Atlas of Breeding Waders in the Russian Arctic.  (Institute of Geography, Russian Academy of Sciences, 2012).

### Normalised Differenced Vegetation Index Data

Based on MODIS satellite images, weekly NDVI data on a scale of 16 × 16 km grid cells were downloaded (ftp://ftp.orbit.nesdis.noaa.gov/pub/corp/scsb/wguo/GVIx/GVIx_VH_16km/NVI/) for the period 1983-2013. Next, cells located in the subspecies’ entire breeding range were selected for further analysis. Then, for each year separately, a smoother was fitted through the data (using the loess function in R, span set to 0.3), where after values due to the albedo effect of snow cover were removed: in case the smoother initially decreased to a minimum before increasing the annual maximum and increased in late autumn (after a minimum) the values before the minimum in spring and after the minimum (phenology thresholds) in autumn were removed. The smoother was then used to determine D.NDVI, i.e. the date at which the fitted NDVI crosses a threshold value of 0 (before reaching the yearly maximum), and to determine A.NDVI, i.e. the area underneath the smoother from the start (D.NDVI) to the end (D2.NDVI) of the season (with the latter defined as the date at which the fitted NDVI crosses a threshold value of 0.1, after having reached the yearly maximum). A threshold of 0.1 was used for defining the end of the season since vegetation index was in general still above 0 at the start of the snow covered season. Note that data were missing for 2001 and the second half of 1994.

### Snow Cover Data

Based on remote sensed NOAA/NCDC climate data record of Northern Hemisphere snow cover extent (Brodzik and Armstrong 2013), weekly snow and ice cover data for the period 1983-2013 on a scale of 24 × 24 km were downloaded (ftp://sidads.colorado.edu/pub/DATASETS/nsidc0046_weekly_snow_seaice/). Next, grid cells falling in the subspecies’ entire breeding range were extracted. Next, large erratic changes in snow cover during summer were removed as these reflect rather unpredictable incidences and do not reflect the main phenology of the seasonal snowfall-thaw cycle. Then, data were modelled using a maximum likelihood fit (mle2 function from R package bbmle) of the asymmetric Gaussian model function (Jönsson and Eklundh 2002) with a binomial error distribution. Using these year-specific fits, date of snowmelt D.SM was then determined as the date that the fitted snow cover falls 0.66 (1/3 of the area ice free), whereas date of snowfall D.SF was determined as the date that the fitted snow cover exceeds 0.33 (1/3 of the area covered with snow).

Brodzik, M. and R. Armstrong. 2013. Northern Hemisphere EASE-Grid 2.0 Weekly Snow Cover and Sea Ice Extent. Version 4. [indicate subset used]. Boulder, Colorado USA: NASA DAAC at the National Snow and Ice Data Center.

Jönsson, P. & Eklundh, L. Seasonality extraction by function fitting to time-series of satellite sensor data. IEEE Transactions on Geoscience and Remote Sensing 40, 1824-1832 (2002).


### Breeding Area of Calidris c. canutus subspecies

Lappo, E. G., Tomkovich, P. S. & Syroechkovskiy, E. Atlas of Breeding Waders in the Russian Arctic.  (Institute of Geography, Russian Academy of Sciences, 2012).
