# Functions to load and set up data set into indicators for further analysis
# 
# Author: Miqdad Asaria
###############################################################################

memory.limit(22000)

library(dplyr)
library(tidyr)

extract_indicator_data_from_database = function(geo, cached){
	if(geo=="CCG"){
		cache_file_name = "data/indicators_db_data_ccg.RData"
	} else if(geo=="LAT") {
		cache_file_name = "data/indicators_db_data_lat.RData"
	} else if(geo=="LAD") {
		cache_file_name = "data/indicators_db_data_lad.RData"
	} else if(geo=="UTLA") {
		cache_file_name = "data/indicators_db_data_utla.RData"
	}
	
	if(cached){
		load(cache_file_name)
	} else {
		source("../db_connection.R")
		con = get_db_connection()
		if(geo=="CCG"){
			geo_sql = "ccg13nm"
		} else if(geo=="LAT") {
			geo_sql = "nhsat13nm AS"
		} else if(geo=="LAD") {
			geo_sql = "lad11nm AS"
		} else if(geo=="UTLA") {
			geo_sql = "utlanm AS"
		}
		sql = paste("SELECT year, ind.lsoa11cd, " ,geo_sql," GEOGRAPHY, indicator, indicator_value, normalised_rank, \'Q\'||quintile AS quintile, \'D\'||decile AS decile 
						FROM 
						HEPI_INDICATORS_LSOA_2011 ind
						INNER JOIN 
						IMD_2010_LSOA_2011 imd
						ON ind.lsoa11cd=imd.lsoa11cd
						INNER JOIN 
						NHS_GEOGRAPHY geo
						ON ind.lsoa11cd=geo.lsoa11cd
						WHERE YEAR BETWEEN 2004 AND 2011",sep="")
		lsoa_data = dbGetQuery(con,sql)
		
		
		dbDisconnect(con)
		
		indicators_list = c("population",
				"need_adjusted_population",
				"amenable_mortality",
				"amenable_mortality_adjusted",
				"gp_fte",
				"cips_count_ind_inpatient_wait",
				"cips_count_ind",
				"unplanned_hospitalisations",
				"unplanned_hospitalisations_adjusted",
				"PHIS",
				"PHIS_REPORTED")
		
		lsoa_data_wide = lsoa_data %>% filter(INDICATOR %in% indicators_list) %>% spread(INDICATOR, INDICATOR_VALUE)
		
		data = list()
		data[["lsoa_data_wide"]] = lsoa_data_wide
		save(data,file=cache_file_name)
	}
	return(data)
}

aggregate_indicator_data = function(lsoa_data_wide){
	indicators = list()
	
	pat_gp_fte = list()
	pat_gp_fte[["indicator_data"]] = lsoa_data_wide %>% 
			select(YEAR,LSOA11CD,GEOGRAPHY,NORMALISED_RANK,QUINTILE,DECILE, numerator=population, denominator=gp_fte) %>% 
			mutate(exposure=numerator, lsoa_indicator=numerator/denominator) %>%
			filter(denominator>0)
	pat_gp_fte[["indicator_label"]] = "Patients per full time equivalent GP"
	pat_gp_fte[["title"]] = "Primary Care Supply (unadjusted)"
	pat_gp_fte[["footnote"]] = "Primary Care Supply: Patients per full time equivalent GP, excluding registrars and retainers"
	pat_gp_fte[["slope_sign"]] = -1
	indicators[["pat_gp_fte"]] = pat_gp_fte
	
	pat_gp_fte_adj = list()
	pat_gp_fte_adj[["indicator_data"]] = lsoa_data_wide %>% 
			select(YEAR,LSOA11CD,GEOGRAPHY,NORMALISED_RANK,QUINTILE,DECILE, numerator=need_adjusted_population, denominator=gp_fte) %>% 
			mutate(exposure=numerator, lsoa_indicator=numerator/denominator) %>%
			filter(denominator>0)
	pat_gp_fte_adj[["indicator_label"]] = "Patients per full time equivalent GP"
	pat_gp_fte_adj[["title"]] = "Primary Care Supply"
	pat_gp_fte_adj[["footnote"]] = "Primary Care Supply: Patients per full time equivalent GP, excluding registrars and retainers, adjusted for age, sex and health deprivation"
	pat_gp_fte_adj[["slope_sign"]] = -1
	indicators[["pat_gp_fte_adj"]] = pat_gp_fte_adj
	
	PHIS = list()
	PHIS[["indicator_data"]] = lsoa_data_wide %>% 
			select(YEAR,LSOA11CD,GEOGRAPHY,NORMALISED_RANK,QUINTILE,DECILE, numerator=PHIS, denominator=population) %>% 
			mutate(numerator=numerator*denominator, exposure=denominator*100, lsoa_indicator=numerator/denominator) 
	PHIS[["indicator_label"]] = "Public health impact score"
	PHIS[["title"]] = "Primary Care Quality (including exceptions)*"
	PHIS[["footnote"]] = "Primary Care Quality: clinical performance in the quality and outcomes framework including exceptions (weighted by public health impact)"
	PHIS[["slope_sign"]] = 1
	indicators[["PHIS"]] = PHIS
	
	PHIS_REPORTED = list()
	PHIS_REPORTED[["indicator_data"]] = lsoa_data_wide %>% 
			select(YEAR,LSOA11CD,GEOGRAPHY,NORMALISED_RANK,QUINTILE,DECILE, numerator=PHIS_REPORTED, denominator=population) %>% 
			mutate(numerator=numerator*denominator, exposure=denominator*100, lsoa_indicator=numerator/denominator) 
	PHIS_REPORTED[["indicator_label"]] = "Public health impact score"
	PHIS_REPORTED[["title"]] = "Primary Care Quality*"
	PHIS_REPORTED[["footnote"]] = "Primary Care Quality: clinical performance in the quality and outcomes framework as reported (weighted by public health impact)"
	PHIS_REPORTED[["slope_sign"]] = 1
	indicators[["PHIS_REPORTED"]] = PHIS_REPORTED
	
	unplanned_hospitalisations = list()
	unplanned_hospitalisations[["indicator_data"]] = lsoa_data_wide %>% 
			select(YEAR,LSOA11CD,GEOGRAPHY,NORMALISED_RANK,QUINTILE,DECILE, numerator=unplanned_hospitalisations, denominator=population) %>% 
			mutate(exposure=1000, denominator=denominator/exposure, lsoa_indicator=numerator/denominator)
	unplanned_hospitalisations[["indicator_label"]] = "Preventable hospital admissions (per 1,000)"
	unplanned_hospitalisations[["title"]] = "Preventable Hospitalisation (unadjusted)"
	unplanned_hospitalisations[["footnote"]] = "Preventable Hospitalisation: hospitalisations per 1,000 population for conditions amenable to healthcare"
	unplanned_hospitalisations[["slope_sign"]] = -1
	indicators[["unplanned_hospitalisations"]] = unplanned_hospitalisations
	
	unplanned_hospitalisations_adj = list()
	unplanned_hospitalisations_adj[["indicator_data"]] = lsoa_data_wide %>% 
			select(YEAR,LSOA11CD,GEOGRAPHY,NORMALISED_RANK,QUINTILE,DECILE, numerator=unplanned_hospitalisations_adjusted, denominator=population) %>% 
			mutate(exposure=1000, denominator=denominator/exposure, lsoa_indicator=numerator/denominator)
	unplanned_hospitalisations_adj[["indicator_label"]] = "Preventable hospital admissions (per 1,000)"
	unplanned_hospitalisations_adj[["title"]] = "Preventable Hospitalisation"
	unplanned_hospitalisations_adj[["footnote"]] = "Preventable Hospitalisation: hospitalisations per 1,000 population for conditions amenable to healthcare adjusted for age and sex"
	unplanned_hospitalisations_adj[["slope_sign"]] = -1
	indicators[["unplanned_hospitalisations_adj"]] = unplanned_hospitalisations_adj
	
	amenable_mortality = list()
	amenable_mortality[["indicator_data"]] = lsoa_data_wide %>% 
			select(YEAR,LSOA11CD,GEOGRAPHY,NORMALISED_RANK,QUINTILE,DECILE, numerator=amenable_mortality, denominator=population) %>% 
			mutate(exposure=1000, denominator=denominator/exposure, lsoa_indicator=numerator/denominator)
	amenable_mortality[["indicator_label"]] = "Amenable mortality (per 1,000)"
	amenable_mortality[["title"]] = "Amenable Mortality (unadjusted)"
	amenable_mortality[["footnote"]] = "Amenable Mortality: deaths per 1,000 population from causes amenable to health care"
	amenable_mortality[["slope_sign"]] = -1
	indicators[["amenable_mortality"]] = amenable_mortality
	
	amenable_mortality_adj = list()
	amenable_mortality_adj[["indicator_data"]] = lsoa_data_wide %>% 
			select(YEAR,LSOA11CD,GEOGRAPHY,NORMALISED_RANK,QUINTILE,DECILE, numerator=amenable_mortality_adjusted, denominator=population) %>% 
			mutate(exposure=1000, denominator=denominator/exposure, lsoa_indicator=numerator/denominator)
	amenable_mortality_adj[["indicator_label"]] = "Amenable mortality (per 1,000)"
	amenable_mortality_adj[["title"]] = "Amenable Mortality"
	amenable_mortality_adj[["footnote"]] = "Amenable Mortality: deaths per 1,000 population from causes amenable to health care adjusted for age and sex"
	amenable_mortality_adj[["slope_sign"]] = -1
	indicators[["amenable_mortality_adj"]] = amenable_mortality_adj
	
	return(indicators)
}

trim_indicator = function(indicator){
	num_sds = 6
	if(indicator[["indicator_label"]]=="Public health impact score"){
		max=100
		num_sds = 3
	} else if(indicator[["indicator_label"]]=="Full time equivalent GPs (per 100,000)"){
		max=10000
	} else{ 
		max=-1
	}
	untrimmed = data_frame()
	num_sds = 6
	for(year in unique(indicator[["indicator_data"]]$YEAR)){
		year_indicator = indicator[["indicator_data"]][indicator[["indicator_data"]]$YEAR==year,]
		year_indicator = year_indicator %>% filter(!is.infinite(lsoa_indicator))
		year_indicator = year_indicator %>% filter(lsoa_indicator >= 0)
		indicator_mean = mean(year_indicator$lsoa_indicator, na.rm=TRUE)
		indicator_sd = sd(year_indicator$lsoa_indicator, na.rm=TRUE)
		lower_limit = indicator_mean - num_sds*indicator_sd
		upper_limit = indicator_mean + num_sds*indicator_sd
		untrimmed = bind_rows(untrimmed,year_indicator %>% filter(lsoa_indicator>=lower_limit & lsoa_indicator<=upper_limit))
	}
	indicator[["indicator_data"]] = untrimmed
	return(indicator)
}

load_indicators = function(cached, trim){
	data = extract_indicator_data_from_database(geo="CCG",cached=TRUE)
	indicators = aggregate_indicator_data(data[["lsoa_data_wide"]])
	if(trim){
	  indicators = lapply(indicators,trim_indicator)
	}
	
	return(indicators)
}