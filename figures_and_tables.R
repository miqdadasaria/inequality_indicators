# Figures and tables for main overview paper for HEPI indicators
# 
# Author: Miqdad Asaria
# Date: July 2015 
###############################################################################

memory.limit(22000)

library(MASS)
library(ggplot2)
library(scales)
library(grid)
library(gridExtra)
library(dplyr)
library(tidyr)
library(xtable)

source("load_data.R")

calculate_national_data = function(indicator){
	CI95 = qnorm(0.975)
	
	national_quintiles = indicator[["indicator_data"]] %>% 
			group_by(YEAR,QUINTILE) %>%
			summarise(national_indicator=sum(numerator,na.rm=TRUE)/sum(denominator,na.rm=TRUE)) %>% 
			spread(QUINTILE,national_indicator)

	national_ineq = indicator[["indicator_data"]] %>%
			select(YEAR, NORMALISED_RANK, lsoa_indicator, numerator, denominator) %>%
			group_by(YEAR) %>%
			mutate(OVERALL_MEAN=sum(numerator)/sum(denominator),
					OVERALL_MEAN_SE=sqrt((sd(numerator)/mean(numerator))^2 + (sd(denominator)/mean(denominator))^2 - 2*cov(numerator,denominator)/(mean(numerator)*mean(denominator)))/sqrt(n()),
					denominator = sum(denominator)
			) %>%
			group_by(YEAR, OVERALL_MEAN, OVERALL_MEAN_SE, denominator) %>%
			do(model = lm(lsoa_indicator~NORMALISED_RANK, data=.,na.action=na.omit)) %>%
			mutate(SII=coef(model)[2]*indicator[["slope_sign"]],
					SII_LCI=(coef(model)[2]-CI95*sqrt(vcov(model)[2,2])*indicator[["slope_sign"]])*indicator[["slope_sign"]], 
					SII_UCI=(coef(model)[2]+CI95*sqrt(vcov(model)[2,2])*indicator[["slope_sign"]])*indicator[["slope_sign"]],
					RII=SII/OVERALL_MEAN, 
					RII_LCI=SII_LCI/OVERALL_MEAN, 
					RII_UCI=SII_UCI/OVERALL_MEAN) %>%
			select(YEAR,OVERALL_MEAN,OVERALL_MEAN_SE,SII,SII_LCI,SII_UCI,RII,RII_LCI,RII_UCI,denominator)
	
	national_data = inner_join(national_quintiles,national_ineq,by=c("YEAR"))
	
	if(grepl("Primary Care Quality",indicator[["title"]])){
		national_data$denominator=1
	}
	
	national_data = national_data %>%
			mutate(REAL_GAP=SII*0.5*denominator, 
					REAL_GAP_LCI=SII_LCI*0.5*denominator, 
					REAL_GAP_UCI=SII_UCI*0.5*denominator,
					INDICATOR=indicator[["title"]]) %>% 
			arrange(YEAR)
	return(national_data %>% ungroup())
}

calculate_real_gap_ci = function(model,denominator_2004,denominator_2011,slope_sign,ci){
	set.seed(123)
	psa_iterations = 10000
	cholesky = chol(vcov(model))
	randoms = matrix(rnorm(length(model$coef)*psa_iterations),nrow=psa_iterations)
	betas = matrix(rep(model$coef,psa_iterations),nrow=psa_iterations,byrow=TRUE) + randoms%*%cholesky
	real_gap = ((betas[,3]+betas[,4])*denominator_2011 - betas[,3]*denominator_2004)*slope_sign*0.5
	if(ci=="LCI"){
		return(quantile(real_gap,0.025))
	}else{
		return(quantile(real_gap,0.975))
	}
}

calculate_change_data = function(indicator){
	SII = indicator[["indicator_data"]] %>%
			select(YEAR, NORMALISED_RANK, lsoa_indicator, numerator, denominator) %>%
			filter(YEAR %in% c(2004,2011)) %>%
			do(model = lm(lsoa_indicator~factor(YEAR)*NORMALISED_RANK, data=.,na.action=na.omit)) %>%
			mutate(SII_CHANGE=coef(model[[1]])[4]*indicator[["slope_sign"]],
					SII_LCI_CHANGE=(coef(model[[1]])[4]-qnorm(0.975)*sqrt(vcov(model[[1]])[4,4])*indicator[["slope_sign"]])*indicator[["slope_sign"]], 
					SII_UCI_CHANGE=(coef(model[[1]])[4]+qnorm(0.975)*sqrt(vcov(model[[1]])[4,4])*indicator[["slope_sign"]])*indicator[["slope_sign"]]) %>%
			select(SII_CHANGE,SII_LCI_CHANGE,SII_UCI_CHANGE)
	
	RII = indicator[["indicator_data"]] %>%
			select(YEAR, NORMALISED_RANK, lsoa_indicator, numerator, denominator) %>%
			filter(YEAR %in% c(2004,2011)) %>%
			group_by(YEAR) %>%
			mutate(OVERALL_MEAN=sum(numerator)/sum(denominator),
					OVERALL_MEAN_SE=sqrt((sd(numerator)/mean(numerator))^2 + (sd(denominator)/mean(denominator))^2 - 2*cov(numerator,denominator)/(mean(numerator)*mean(denominator)))/sqrt(n())
			) %>%
			ungroup() %>%
			do(model = lm((lsoa_indicator/OVERALL_MEAN)~factor(YEAR)*NORMALISED_RANK, data=.,na.action=na.omit)) %>%
			mutate(RII_CHANGE=coef(model[[1]])[4]*indicator[["slope_sign"]],
					RII_LCI_CHANGE=(coef(model[[1]])[4]-qnorm(0.975)*sqrt(vcov(model[[1]])[4,4])*indicator[["slope_sign"]])*indicator[["slope_sign"]], 
					RII_UCI_CHANGE=(coef(model[[1]])[4]+qnorm(0.975)*sqrt(vcov(model[[1]])[4,4])*indicator[["slope_sign"]])*indicator[["slope_sign"]]) %>%
			select(RII_CHANGE,RII_LCI_CHANGE,RII_UCI_CHANGE)
	
	means = indicator[["indicator_data"]] %>%
			select(YEAR, NORMALISED_RANK, lsoa_indicator, numerator, denominator) %>%
			filter(YEAR %in% c(2004,2011)) %>%
			group_by(YEAR) %>%
			summarise(OVERALL_MEAN=sum(numerator)/sum(denominator),
					OVERALL_MEAN_SE=sqrt((sd(numerator)/mean(numerator))^2 + (sd(denominator)/mean(denominator))^2 - 2*cov(numerator,denominator)/(mean(numerator)*mean(denominator)))/sqrt(n()),
					denominator=sum(denominator)
				)
	
	mean_2011 = subset(means,YEAR==2011)
	mean_2004 = subset(means,YEAR==2004)
	mean_change_se = sqrt(mean_2011$OVERALL_MEAN_SE^2 + mean_2004$OVERALL_MEAN_SE^2) 

	if(indicator[["title"]]=="Primary Care Quality"){
		mean_2011$denominator=1
		mean_2004$denominator=1
	}
	
	REAL_GAP = indicator[["indicator_data"]] %>%
			select(YEAR, NORMALISED_RANK, lsoa_indicator, numerator, denominator) %>%
			filter(YEAR %in% c(2004,2011)) %>%
			do(model = lm(lsoa_indicator~factor(YEAR)*NORMALISED_RANK, data=.,na.action=na.omit)) %>%
			mutate(REAL_GAP_CHANGE=((coef(model[[1]])[3]+coef(model[[1]])[4])*mean_2011$denominator - coef(model[[1]])[3]*mean_2004$denominator)*indicator[["slope_sign"]]*0.5,
					REAL_GAP_LCI_CHANGE=calculate_real_gap_ci(model[[1]],mean_2004$denominator,mean_2011$denominator,indicator[["slope_sign"]],"LCI"),
					REAL_GAP_UCI_CHANGE=calculate_real_gap_ci(model[[1]],mean_2004$denominator,mean_2011$denominator,indicator[["slope_sign"]],"UCI")) %>%
			select(REAL_GAP_CHANGE,REAL_GAP_LCI_CHANGE,REAL_GAP_UCI_CHANGE)
	
	
	change_data = bind_cols(SII,RII,REAL_GAP) %>%
					mutate(	OVERALL_MEAN_CHANGE = mean_2011$OVERALL_MEAN-mean_2004$OVERALL_MEAN,
							OVERALL_MEAN_UCI_CHANGE = OVERALL_MEAN_CHANGE + qnorm(0.975)*mean_change_se,
							OVERALL_MEAN_LCI_CHANGE = OVERALL_MEAN_CHANGE - qnorm(0.975)*mean_change_se,
							INDICATOR=indicator[["title"]])
	return(change_data)
}

get_legend = function(plot, position){
	g = ggplotGrob(plot + theme(legend.position=position))$grobs
	legend = g[[which(sapply(g, function(x) x$name) == "guide-box")]]
	return(legend)
}

composite_panel_plot = function(national_data, indicator, zero_line=FALSE, reverse=FALSE){
	imd_labels = c("Q1 (most deprived)","Q2","Q3","Q4","Q5 (least deprived)")
	year_labels = data.frame(YEAR=2000:2014,YEAR_LABEL=c("00/01","01/02","02/03","03/04","04/05","05/06","06/07","07/08","08/09","09/10","10/11","11/12","12/13","13/14","14/15"))
	graph_data = national_data %>%
			filter(INDICATOR==indicator[["title"]]) %>%
			inner_join(year_labels,by="YEAR") %>% 
			gather(IMD_QUINTILE,value,Q1:Q5) %>%
			select(YEAR=YEAR_LABEL,IMD_QUINTILE,value)
	main_plot = ggplot(graph_data) + 
			aes(x=YEAR, y=value, group=IMD_QUINTILE, colour=IMD_QUINTILE) + 
			geom_line(aes(linetype=IMD_QUINTILE, size=IMD_QUINTILE)) + 
			geom_point(aes(shape=IMD_QUINTILE, colour=IMD_QUINTILE)) +
			xlab("Year") +
			ylab(indicator[["indicator_label"]]) +
			scale_colour_manual(name="IMD Group", values=c("black","lightblue","lightblue","lightblue","darkgrey"), labels=imd_labels) +
			scale_shape_manual(name="IMD Group", values=c(19,21,24,0,15), labels=imd_labels) +
			scale_linetype_manual(name="IMD Group", values=c(1,2,2,2,1), labels=imd_labels) +
			scale_size_manual(name="IMD Group", values=c(1,0.5,0.5,0.5,1), labels=imd_labels) +
			theme_bw() +
			theme(panel.grid.major = element_blank(), 
					panel.grid.minor = element_blank(), 
					legend.position="none")
	
	if(reverse){
		main_plot = main_plot + scale_y_reverse(labels = comma)
	} else {
		main_plot = main_plot + scale_y_continuous(labels = comma)
	}
	
	index_data = national_data %>% 
			filter(INDICATOR==indicator[["title"]]) %>%
			inner_join(year_labels,by="YEAR") %>% 
			select(YEAR,YEAR_LABEL,RII,RII_LCI,RII_UCI,SII,SII_UCI,SII_LCI)
	
	SII_plot = ggplot(index_data) + 
			aes(x=YEAR, y=SII) + 
			geom_line() + 
			geom_point() +
			geom_errorbar(aes(ymin=SII_LCI, ymax=SII_UCI), width=0.1 ) +
			xlab("Year") +
			ylab("") +
			ggtitle("Slope Index of Inequality") +
			scale_x_continuous(
					breaks=min(index_data$YEAR):max(index_data$YEAR),	
					labels=index_data$YEAR_LABEL)+
			scale_y_continuous(labels = comma) +
			theme_bw() +
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
	RII_plot = ggplot(index_data) + 
			aes(x=YEAR, y=RII) + 
			geom_line() + 
			geom_point() +
			geom_errorbar(aes(ymin=RII_LCI, ymax=RII_UCI), width=0.1 ) +
			xlab("Year") +
			ylab("") +
			ggtitle("Relative Index of Inequality") +
			scale_x_continuous(
					breaks=min(index_data$YEAR):max(index_data$YEAR),	
					labels=index_data$YEAR_LABEL) +
			scale_y_continuous(labels = percent) +
			theme_bw() +
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
	if(zero_line){
		#SII_plot = SII_plot + geom_hline(yintercept=0, colour="darkgrey", linetype=2)
		#RII_plot = RII_plot + geom_hline(yintercept=0, colour="darkgrey", linetype=2)
		main_plot = main_plot + geom_hline(yintercept=0, colour="white", linetype=2)
	}
	
	main_panel = arrangeGrob(
			main_plot,
			RII_plot, 
			SII_plot,
			widths = c(0.4,0.3,0.3),
			nrow=1)
	
	panel_plot = 
			arrangeGrob(
					main=textGrob(indicator[["title"]],gp=gpar(fontsize=20,fontface="bold")),
					main_panel,
					sub = textGrob(indicator[["footnote"]], x = 0, hjust = -0.1, vjust=0.1,	gp = gpar(fontface = "italic", fontsize = 12)),
					heights=c(0.1,0.85,0.05))
	return(panel_plot)
}

panel_plot = function(national_data, indicator, sub, zero_line){
	imd_labels = c("Q1 (most deprived)","Q2","Q3","Q4","Q5 (least deprived)")
	year_labels = data.frame(YEAR=2000:2014,YEAR_LABEL=c("00/01","01/02","02/03","03/04","04/05","05/06","06/07","07/08","08/09","09/10","10/11","11/12","12/13","13/14","14/15"))
	graph_data = national_data %>%
			filter(INDICATOR==indicator[["title"]]) %>%
			inner_join(year_labels,by="YEAR") %>% 
			gather(IMD_QUINTILE,value,Q1:Q5) %>%
			select(YEAR=YEAR_LABEL,IMD_QUINTILE,value)
	main_plot = ggplot(graph_data) + 
			aes(x=YEAR, y=value, group=IMD_QUINTILE, colour=IMD_QUINTILE) + 
			geom_line(aes(linetype=IMD_QUINTILE, size=IMD_QUINTILE)) + 
			geom_point(aes(shape=IMD_QUINTILE, colour=IMD_QUINTILE)) +
			xlab("Year") +
			ylab(indicator[["indicator_label"]]) +
			scale_y_continuous(labels = comma) +
			scale_colour_manual(name="IMD Group", values=c("black","lightblue","lightblue","lightblue","darkgrey"), labels=imd_labels) +
			scale_shape_manual(name="IMD Group", values=c(19,21,24,0,15), labels=imd_labels) +
			scale_linetype_manual(name="IMD Group", values=c(1,2,2,2,1), labels=imd_labels) +
			scale_size_manual(name="IMD Group", values=c(1,0.5,0.5,0.5,1), labels=imd_labels) +
			theme_bw() +
			theme(panel.grid.major = element_blank(), 
					panel.grid.minor = element_blank(), 
					legend.position="none")
	if(zero_line){
		main_plot = main_plot + geom_hline(yintercept=0, colour="white", linetype=2)
	}
	
	panel_plot = 
			arrangeGrob(
					main=textGrob(indicator[["title"]],gp=gpar(fontsize=20,fontface="bold")),
					main_plot,
					sub = textGrob(indicator[["footnote"]], x = 0, hjust = -0.1, vjust=0.1, gp = gpar(fontface = "italic", fontsize = 6)),
					heights=c(0.1,0.85,0.05))
	if(sub){
		return(panel_plot)
	} else {
		return(main_plot)
	}
}

index_plot = function(national_data, indicator){
	imd_labels = c("Q1 (most deprived)","Q2","Q3","Q4","Q5 (least deprived)")
	year_labels = data.frame(YEAR=2000:2014,YEAR_LABEL=c("00/01","01/02","02/03","03/04","04/05","05/06","06/07","07/08","08/09","09/10","10/11","11/12","12/13","13/14","14/15"))
	
	index_data = national_data %>% 
			filter(INDICATOR==indicator[["title"]]) %>%
			inner_join(year_labels,by="YEAR") %>% 
			select(YEAR,YEAR_LABEL,RII,RII_LCI,RII_UCI,SII,SII_UCI,SII_LCI)
	
	SII_plot = ggplot(index_data) + 
			aes(x=YEAR, y=SII) + 
			geom_line() + 
			geom_point() +
			geom_errorbar(aes(ymin=SII_LCI, ymax=SII_UCI), width=0.1 ) +
			xlab("Year") +
			ylab("") +
			ggtitle("Slope Index of Inequality") +
			scale_x_continuous(
					breaks=min(index_data$YEAR):max(index_data$YEAR),	
					labels=index_data$YEAR_LABEL)+
			scale_y_continuous(labels = comma) +
			theme_bw() +
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
	RII_plot = ggplot(index_data) + 
			aes(x=YEAR, y=RII) + 
			geom_line() + 
			geom_point() +
			geom_errorbar(aes(ymin=RII_LCI, ymax=RII_UCI), width=0.1 ) +
			xlab("Year") +
			ylab("") +
			ggtitle("Relative Index of Inequality") +
			scale_x_continuous(
					breaks=min(index_data$YEAR):max(index_data$YEAR),	
					labels=index_data$YEAR_LABEL) +
			scale_y_continuous(labels = percent) +
			theme_bw() +
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

	
	index_plot = arrangeGrob(
					main=textGrob(indicator[["title"]],gp=gpar(fontsize=20,fontface="bold")),
					arrangeGrob(RII_plot, SII_plot,ncol=2),
					heights = c(0.1,0.9))
	
	return(index_plot)
}


plot_scatter = function(year, indicator){
	indicator_data = indicator[["indicator_data"]] %>% 
			filter(YEAR==year) %>% 
			select(NORMALISED_RANK, lsoa_indicator)
	average = mean(indicator_data$lsoa_indicator)
	model = lm(lsoa_indicator ~ NORMALISED_RANK, data=indicator_data)
	max_point = predict(model,data_frame(NORMALISED_RANK=1))
	min_point = predict(model,data_frame(NORMALISED_RANK=0))
	
	trend = predict(model,data_frame(NORMALISED_RANK=seq(0,1,0.1)))
	trend_se = predict(model,data_frame(NORMALISED_RANK=seq(0,1,0.1)),se.fit=TRUE)$se.fit
	CI_95 = qnorm(0.975)
	trend_UCI = trend + CI_95*trend_se
	trend_LCI = trend - CI_95*trend_se
	level = rep(max_point, length(trend))
	area = data.frame(DEP=as.factor(seq(0,1,0.1)),level,trend,trend_LCI,trend_UCI,average)
	averages = indicator_data %>%
			mutate(DEPRIVATION=round(NORMALISED_RANK,1)) %>%
			group_by(DEPRIVATION) %>%
			summarise(lsoa_indicator=mean(lsoa_indicator)) %>%
			mutate(DEP=as.factor(DEPRIVATION))
	graph_data = inner_join(area, averages, by="DEP")
	scatter = ggplot(data=graph_data, aes(x=(1-DEPRIVATION),y=lsoa_indicator)) +
			geom_point(size=2, colour="black") +
			xlab("Small area deprivation rank") + 
			ylab(indicator[["indicator_label"]]) +
			ggtitle(paste(year,"/",substr(as.character(year+1),3,4),sep="")) +
			scale_x_continuous(breaks=seq(0,1,0.2), labels=c("least deprived","","","","","most deprived")) +
			geom_line(aes(y=trend), colour="black") +	
			geom_line(aes(y=average),linetype=2, size=1.5, colour="red") +
			geom_ribbon(aes(ymax=trend, ymin=level), alpha=0.1) +
			theme_bw() +
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

	result = list()
	result[["plot"]] = scatter
	result[["limits"]] = c(min(trend,graph_data$lsoa_indicator),max(trend,graph_data$lsoa_indicator))
	return(result)
}

scatter_panel = function(start_year, end_year, indicator, reverse){
	before = plot_scatter(start_year, indicator)
	after = plot_scatter(end_year, indicator)
	ymin = min(before[["limits"]],after[["limits"]])
	ymax = max(before[["limits"]],after[["limits"]])
	before_plot = before[["plot"]] 
	after_plot = after[["plot"]]
	
	if(reverse){
		before_plot = before_plot + scale_y_reverse(labels=comma, limits=c(ymax,ymin))
		after_plot = after_plot + scale_y_reverse(labels=comma, limits=c(ymax,ymin))
	} else {
		before_plot = before_plot + scale_y_continuous(labels=comma, limits=c(ymin,ymax))
		after_plot = after_plot + scale_y_continuous(labels=comma, limits=c(ymin,ymax))
	}
	scatter_panel = arrangeGrob(
					main=textGrob(indicator[["title"]],gp=gpar(fontsize=20,fontface="bold")),
					arrangeGrob(before_plot, after_plot, ncol=2),
					sub = textGrob(indicator[["footnote"]], x = 0, hjust = -0.1, vjust=0.1,gp = gpar(fontface = "italic", fontsize = 10)),
					heights = c(0.1, 0.85, 0.05)
					)
	
	return(scatter_panel)
	
}
#####################################################################################################
# Run code to load data and produce figures and tables used in the JECH paper
#####################################################################################################

indicators = load_indicators(cached=TRUE, trim=TRUE)
national_data = do.call("rbind",lapply(indicators,calculate_national_data))

spacer = rectGrob(gp=gpar(col="white"))
scatter_composite = arrangeGrob(
		scatter_panel(2004, 2011, indicators[["pat_gp_fte_adj"]], FALSE),
		spacer,
		scatter_panel(2004, 2011, indicators[["PHIS_REPORTED"]], TRUE),
		spacer,
		scatter_panel(2004, 2011, indicators[["unplanned_hospitalisations_adj"]], FALSE),
		spacer,
		scatter_panel(2004, 2011, indicators[["amenable_mortality_adj"]], FALSE),
		nrow=7,
		heights=c(0.24,0.02,0.24,0.02,0.24,0.02,0.24))
ggsave(filename="output/figure_1.png",plot=scatter_composite,width=28,height=45,units="cm",dpi=300)

legend = get_legend(panel_plot(national_data, indicators[["pat_gp_fte_adj"]], sub=FALSE, zero_line=FALSE),"right")
panel = arrangeGrob(
		arrangeGrob(
				panel_plot(national_data, indicators[["pat_gp_fte_adj"]], sub=TRUE, zero_line=FALSE),
				spacer,
				panel_plot(national_data, indicators[["PHIS_REPORTED"]], sub=TRUE, zero_line=FALSE),
				spacer,
				panel_plot(national_data, indicators[["unplanned_hospitalisations_adj"]], sub=TRUE, zero_line=TRUE),
				spacer,
				panel_plot(national_data, indicators[["amenable_mortality_adj"]], sub=TRUE, zero_line=TRUE),
				nrow=7,
				heights=c(0.23,0.027,0.23,0.027,0.23,0.027,0.23)),
		legend,
		ncol = 2,
		widths = unit.c(unit(1, "npc") - sum(legend$width), sum(legend$width)))
ggsave(filename="output/figure_2a.png",plot=panel,width=20,height=40,units="cm",dpi=300)

index = arrangeGrob(
		index_plot(national_data, indicators[["pat_gp_fte_adj"]]),
		spacer,
		index_plot(national_data, indicators[["PHIS_REPORTED"]]),
		spacer,
		index_plot(national_data, indicators[["unplanned_hospitalisations_adj"]]),
		spacer,
		index_plot(national_data, indicators[["amenable_mortality_adj"]]),
		nrow=7,
		heights=c(0.23,0.027,0.23,0.027,0.23,0.027,0.23))
ggsave(filename="output/figure_2b.png",plot=index,width=20,height=40,units="cm",dpi=300)

overall_panel = arrangeGrob(
		composite_panel_plot(national_data, indicators[["pat_gp_fte_adj"]], zero_line=FALSE, reverse=FALSE),
		spacer,
		composite_panel_plot(national_data, indicators[["PHIS_REPORTED"]], zero_line=FALSE, reverse=TRUE),
		spacer,
		composite_panel_plot(national_data, indicators[["unplanned_hospitalisations_adj"]], zero_line=TRUE, reverse=FALSE),
		spacer,
		composite_panel_plot(national_data, indicators[["amenable_mortality_adj"]], zero_line=TRUE, reverse=FALSE),
		nrow=7,
		heights=c(0.23,0.027,0.23,0.027,0.23,0.027,0.23))
composite_panel = arrangeGrob(
		overall_panel,
		legend,
		ncol = 2,
		widths = unit.c(unit(1, "npc") - sum(legend$width), sum(legend$width)))
ggsave(filename="output/figure_2_composite.png",plot=composite_panel,width=40,height=40,units="cm",dpi=300)

phis_sensitivity = arrangeGrob(
		spacer,
		composite_panel_plot(national_data, indicators[["PHIS_REPORTED"]], zero_line=FALSE, reverse=TRUE),
		spacer,
		composite_panel_plot(national_data, indicators[["PHIS"]], zero_line=FALSE, reverse=TRUE),
		nrow=4,
		heights=c(0.05,0.45,0.05,0.45))
phis_panel = arrangeGrob(
		phis_sensitivity,
		legend,
		ncol = 2,
		widths = unit.c(unit(1, "npc") - sum(legend$width), sum(legend$width)))
ggsave(filename="output/figure_a1_phis_sensitivity.png",plot=phis_panel,width=40,height=20,units="cm",dpi=300)

results_tables = national_data %>% 
			filter(YEAR %in% c(2004,2011)) %>%
			select(-(Q1:Q5),-denominator) %>%
			mutate(SII_SE = abs(SII-SII_LCI)/qnorm(0.975), 
					RII_SE = abs(RII-RII_LCI)/qnorm(0.975),
					REAL_GAP_SE = abs(REAL_GAP_UCI-REAL_GAP_LCI)/qnorm(0.975),
					OVERALL_MEAN_UCI = OVERALL_MEAN+qnorm(0.975)*OVERALL_MEAN_SE,
					OVERALL_MEAN_LCI = OVERALL_MEAN-qnorm(0.975)*OVERALL_MEAN_SE) %>%
					arrange(INDICATOR)

indicator_names = results_tables %>% select(INDICATOR) %>% unique() %>% arrange(INDICATOR)			
results_2004 = results_tables %>%
				filter(YEAR==2004)
results_2011 = results_tables %>%
		filter(YEAR==2011)
change_data = do.call("rbind",lapply(indicators,calculate_change_data))
colnames(results_2004) = paste(colnames(results_2004),"2004",sep="_")
colnames(results_2011) = paste(colnames(results_2011),"2011",sep="_")

results = inner_join(results_2004,inner_join(results_2011,change_data,by=c("INDICATOR_2011"="INDICATOR")),by=c("INDICATOR_2004"="INDICATOR_2011")) %>%
		select(INDICATOR=INDICATOR_2004,
				OVERALL_MEAN_2004,OVERALL_MEAN_LCI_2004,OVERALL_MEAN_UCI_2004,
				OVERALL_MEAN_2011,OVERALL_MEAN_LCI_2011,OVERALL_MEAN_UCI_2011,
				OVERALL_MEAN_CHANGE,OVERALL_MEAN_LCI_CHANGE,OVERALL_MEAN_UCI_CHANGE,
				RII_2004,RII_LCI_2004,RII_UCI_2004,
				RII_2011,RII_LCI_2011,RII_UCI_2011,
				RII_CHANGE,RII_LCI_CHANGE,RII_UCI_CHANGE,
				SII_2004,SII_LCI_2004,SII_UCI_2004,
				SII_2011,SII_LCI_2011,SII_UCI_2011,
				SII_CHANGE,SII_LCI_CHANGE,SII_UCI_CHANGE,
				REAL_GAP_2004,REAL_GAP_LCI_2004,REAL_GAP_UCI_2004,
				REAL_GAP_2011,REAL_GAP_LCI_2011,REAL_GAP_UCI_2011,
				REAL_GAP_CHANGE,REAL_GAP_LCI_CHANGE,REAL_GAP_UCI_CHANGE)

results = results %>%
		gather(key, value, -INDICATOR) %>%
		extract(key, c("measure","year"),"(^.*)_(.*)$") %>%
		spread(year, value)

write.csv(results,"output/results_table.csv")

