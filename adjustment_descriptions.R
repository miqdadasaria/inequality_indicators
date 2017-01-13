# TODO: Add comment
# 
# Author: ma725
###############################################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

source("../db_connection.R")
con = get_db_connection()

sql = "SELECT indicator, year, imd_quintile, 
SUM(population) population, 
SUM(observed) observed, 
SUM(expected) expected,
SUM(adjusted_count) adjusted
FROM
HEPI_INDICATORS_ADJUSTED_LSOA
GROUP BY 
indicator, year, imd_quintile
ORDER BY 
indicator, year, imd_quintile"

sql_cips = "SELECT year, imd_group AS imd_quintile, 
SUM(indicator_value) cips 
FROM
HEPI_INDICATORS_LSOA
WHERE indicator='cips_count_ind'
GROUP BY 
indicator, year, imd_group
ORDER BY 
indicator, year, imd_group"

indicators = dbGetQuery(con,sql)
cips = dbGetQuery(con,sql_cips)
dbDisconnect(con)

indicators2 = indicators %>%
			mutate(IMD_QUINTILE=as.factor(paste("Q",IMD_QUINTILE,sep="")))  
cips2 = cips %>%
		mutate(IMD_QUINTILE=as.factor(IMD_QUINTILE))
	
data = inner_join(indicators2, cips2, by=c("YEAR","IMD_QUINTILE"))

imd_labels = c("Q1 (most deprived)","Q2","Q3","Q4","Q5 (least deprived)")
year_labels = data.frame(YEAR=2001:2011,YEAR_LABEL=c("01/02","02/03","03/04","04/05","05/06","06/07","07/08","08/09","09/10","10/11","11/12"))

plot_indicator = function(data, indicator, name, denominator){
	if (denominator=="population") {
		graph_data = data %>% 
				filter(INDICATOR==indicator) %>%
				mutate(OBSERVED=OBSERVED*1000/POPULATION,
						EXPECTED=EXPECTED*1000/POPULATION,
						ADJUSTED=ADJUSTED*1000/POPULATION) %>%
				select(YEAR,IMD_QUINTILE,OBSERVED,EXPECTED,ADJUSTED) %>%
				gather(key,value,-YEAR,-IMD_QUINTILE)
	} else {
		graph_data = data %>% 
				filter(INDICATOR==indicator) %>%
				mutate(OBSERVED=OBSERVED*1000/CIPS,
						EXPECTED=EXPECTED*1000/CIPS,
						ADJUSTED=ADJUSTED*1000/CIPS) %>%
				select(YEAR,IMD_QUINTILE,OBSERVED,EXPECTED,ADJUSTED) %>%
				gather(key,value,-YEAR,-IMD_QUINTILE)
	}

	plot = ggplot(graph_data) + 
			aes(x=YEAR, y=value, group=IMD_QUINTILE, colour=IMD_QUINTILE) + 
			geom_line(aes(linetype=IMD_QUINTILE, size=IMD_QUINTILE)) + 
			geom_point(aes(shape=IMD_QUINTILE, colour=IMD_QUINTILE)) +
			xlab("Year") +
			ylab(paste(name," per 1,000",sep="")) +
			ggtitle(name) + 
			scale_y_continuous(labels = comma) +
			scale_x_continuous(
					breaks=min(graph_data$YEAR):max(graph_data$YEAR),	
					labels=year_labels$YEAR_LABEL) +
			scale_colour_manual(name="IMD Group", values=c("black","lightblue","lightblue","lightblue","darkgrey"), labels=imd_labels) +
			scale_shape_manual(name="IMD Group", values=c(19,21,24,0,15), labels=imd_labels) +
			scale_linetype_manual(name="IMD Group", values=c(1,2,2,2,1), labels=imd_labels) +
			scale_size_manual(name="IMD Group", values=c(1,0.5,0.5,0.5,1), labels=imd_labels) +
			facet_wrap(~key,nrow=2) +
			theme_bw() +
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(lineheight=.8, face="bold", size=rel(1.7)))
	ggsave(paste("output/appendix/",gsub("\\s","_",name),".png",sep=""), plot=plot, dpi=300, width=15, height=10)
	
}

plot_indicator(data, indicator="all_mortality_adjusted", name="All cause mortality", denominator="population")

plot_indicator(data, indicator="amenable_mortality_adjusted", name="Amenable mortality", denominator="population")

plot_indicator(data, indicator="hosp_discharge_failure_no_death_adjusted", name="Repeat hospitalisation", denominator="cips")

plot_indicator(data, indicator="unplanned_hospitalisations_adjusted", name="Preventable hospitalisation", denominator="population")
