# Exploring the impact of using different models to calculate inequality
# 
# Author: Miqdad Asaria
###############################################################################
library(MASS)
library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(gridExtra)

source("load_data.R")

get_legend = function(plot, position){
  g = ggplotGrob(plot + theme(legend.position=position))$grobs
  legend = g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  return(legend)
}

indicator_plots = function(indicator){

  scatter_data = data_frame()
  inequality_data = data_frame()
  for(year in min(indicator[["indicator_data"]]$YEAR):max(indicator[["indicator_data"]]$YEAR)){
    indicator_data = indicator[["indicator_data"]] %>% 
      filter(YEAR==year)
    
    overall_mean = sum(indicator_data$numerator)/sum(indicator_data$denominator)
    model_lm = lm(lsoa_indicator~NORMALISED_RANK, data=indicator_data)
    model_negbin = glm.nb(lsoa_indicator~NORMALISED_RANK, data=indicator_data)
    
    trend_lm = predict(model_lm,data_frame(NORMALISED_RANK=seq(0,1,0.1)))
    trend_negbin = predict(model_negbin, data_frame(NORMALISED_RANK=seq(0,1,0.1)), type="response")
    deprivation = 1-seq(0,1,0.1)

    averages = indicator_data %>%
      mutate(deprivation=1-round(NORMALISED_RANK,1)) %>%
      group_by(deprivation) %>%
      summarise(observed=mean(lsoa_indicator))
    
    scatter_data = scatter_data %>% bind_rows(data_frame(year, deprivation, model="OLS", predicted=trend_lm),
                                              data_frame(year, deprivation, model="NegBin", predicted=trend_negbin),
                                              data_frame(year, deprivation,model="Observed", predicted=rev(averages$observed)))
    
    sii_lm = coef(model_lm)[2]*indicator$slope_sign
    sii_lm_confint = confint(model_lm)
    sii_lm_uci = sii_lm_confint[2,1]*indicator$slope_sign
    sii_lm_lci = sii_lm_confint[2,2]*indicator$slope_sign
    
    sii_negbin = (exp(sum(coef(model_negbin))) - exp(coef(model_negbin)[1]))*indicator$slope_sign
    sii_negbin_confint = confint(model_negbin)
    sii_negbin_uci = (exp(sum(sii_negbin_confint[,1])) - exp(sii_negbin_confint[1,1]))*indicator$slope_sign
    sii_negbin_lci = (exp(sum(sii_negbin_confint[,2])) - exp(sii_negbin_confint[1,2]))*indicator$slope_sign
    
    inequality_data = inequality_data %>% 
      bind_rows(data_frame(year,model="OLS",SII=sii_lm,SII_LCI=sii_lm_lci,SII_UCI=sii_lm_uci,RII=sii_lm/overall_mean,RII_LCI=sii_lm_lci/overall_mean,RII_UCI=sii_lm_uci/overall_mean),
                data_frame(year,model="NegBin",SII=sii_negbin,SII_LCI=sii_negbin_lci,SII_UCI=sii_negbin_uci,RII=sii_negbin/overall_mean,RII_LCI=sii_negbin_lci/overall_mean,RII_UCI=sii_negbin_uci/overall_mean))
    
  }
  scatter_data$model = factor(scatter_data$model,c("Observed","OLS","NegBin"))
  inequality_data$model = factor(inequality_data$model,c("OLS","NegBin"))
  
  scatter = ggplot(data=scatter_data, aes(x=deprivation, y=predicted, group=model, colour=model, linetype=model, shape=model)) +
    geom_line(size=1) +
    geom_point(size=2) +
    ggtitle(indicator[["indicator_label"]]) +
    xlab("Small area deprivation rank") + 
    ylab(indicator[["indicator_label"]]) +
    scale_x_continuous(breaks=seq(0,1,0.2), labels=c("least deprived","","","","","most deprived")) +
    scale_linetype_manual(name="Model",labels=c("Observed","OLS","NegBin"),values=c("blank","solid","longdash")) +
    scale_colour_manual(name="Model",labels=c("Observed","OLS","NegBin"),values=c("darkgrey","black","darkred")) +
    scale_shape_discrete(name="Model",labels=c("Observed","OLS","NegBin")) +
    facet_wrap(~year,nrow=3) +
    theme_bw() +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.spacing=unit(1,"lines"))

  sii_plot = ggplot(inequality_data) + 
    aes(x=year, y=SII, group=model, linetype=model, colour=model, shape=model) + 
    geom_line() + 
    geom_point() +
    geom_errorbar(aes(ymin=SII_LCI, ymax=SII_UCI), width=0.1 ) +
    xlab("Year") +
    ylab("Slope Index of Inequality") +
    scale_y_continuous(labels = comma) +
    scale_linetype_manual(name="Model",labels=c("OLS","NegBin"),values=c("solid","longdash")) +
    scale_colour_manual(name="Model",labels=c("OLS","NegBin"),values=c("black","darkred")) +
    scale_shape_discrete(name="Model",labels=c("OLS","NegBin")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none")

  rii_plot = ggplot(inequality_data) + 
    aes(x=year, y=RII, group=model, linetype=model, colour=model, shape=model) + 
    geom_line() + 
    geom_point() +
    geom_errorbar(aes(ymin=RII_LCI, ymax=RII_UCI), width=0.1 ) +
    xlab("Year") +
    ylab("Relative Index of Inequality") +
    scale_y_continuous(labels = percent) +
    scale_linetype_manual(name="Model",labels=c("OLS","NegBin"),values=c("solid","longdash")) +
    scale_colour_manual(name="Model",labels=c("OLS","NegBin"),values=c("black","darkred")) +
    scale_shape_discrete(name="Model",labels=c("OLS","NegBin")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  rii_sii_legend = get_legend(rii_plot,"right")
  rii_plot = rii_plot + theme(legend.position="none")
  
  rii_sii_plot = arrangeGrob(
    main=textGrob(gsub("\\*","",indicator[["title"]]),gp=gpar(fontsize=20,fontface="bold")),
    arrangeGrob(rii_plot, sii_plot, rii_sii_legend, ncol=3, widths=c(0.45,0.45,0.1)),
    heights = c(0.1,0.9))
  
  ggsave(filename=paste0("output/appendix_scatter_panel_",gsub("\\*","",gsub(" ","_",tolower(indicator$title))),".png"),plot=scatter,width=20,height=30,units="cm",dpi=300)

  ggsave(filename=paste0("output/appendix_rii_sii_panel_",gsub("\\*","",gsub(" ","_",tolower(indicator$title))),".png"),plot=rii_sii_plot,width=24,height=10,units="cm",dpi=300)
  
  return(TRUE)
}

indicators = load_indicators(cached=TRUE, trim=TRUE)

indicators = indicators[c("pat_gp_fte_adj","PHIS_REPORTED","unplanned_hospitalisations_adj","amenable_mortality_adj")]

lapply(indicators, indicator_plots)