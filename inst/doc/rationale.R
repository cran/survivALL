## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", warning = FALSE, message = FALSE, fig.align = 'center', results = 'asis', fig.show = 'hold', fig.width = 7, fig.height = 5)

## ------------------------------------------------------------------------
library(survivALL)
library(survival)
library(survcomp)
library(Biobase)
library(magrittr)
library(ggplot2)
library(GGally)
library(ggthemes)
library(cowplot); theme_set(theme_grey())

## ------------------------------------------------------------------------
data(nki_subset)

## ---- eval = FALSE-------------------------------------------------------
#  pData(nki_subset)[1:3, ]

## ---- echo = FALSE-------------------------------------------------------
library(pander)
library(magrittr)
pData(nki_subset)[1:3, ] %>% pandoc.table(caption = "Survival Information")

## ------------------------------------------------------------------------
#ERBB2 expression vector
erbb2_xpr <- exprs(nki_subset)["NM_004448",] 
#convert to binary classifier
erbb2_med <- ifelse(erbb2_xpr >= median(erbb2_xpr), "high", "low") 

## ------------------------------------------------------------------------
srv_obj <- survival::Surv(nki_subset$t.dmfs, nki_subset$e.dmfs)
median_fit <- survival::survfit(srv_obj ~ erbb2_med)

p_med <- GGally::ggsurv(median_fit, surv.col = c("#525252", "#bdbdbd")) +
    ylim(0.4, 1) +
    labs(title = "Figure S1", 
         subtitle = "Median approach", 
         x = "Years") +
    theme_pander() + 
    theme(axis.line = element_line(size = 0.1), 
          legend.position = 'bottom')

p_med <- ggdraw(add_sub(p_med, "p = 0.9", 
                    vpadding=grid::unit(0, "lines"), 
                    y = 17, 
                    x = 0.55, 
                    hjust = 0, 
                    size = 10))
p_med

## ------------------------------------------------------------------------
broom::tidy(survival::coxph(srv_obj ~ erbb2_med)) %>% pandoc.table()

## ------------------------------------------------------------------------
erbb2_hypothesis <- ifelse(erbb2_xpr >= quantile(erbb2_xpr, probs = 0.75), "high", "low") 
hypothesis_fit <- survival::survfit(srv_obj ~ erbb2_hypothesis)

p_hyp <- GGally::ggsurv(hypothesis_fit, surv.col = c("#525252", "#bdbdbd")) +
    ylim(0.4, 1) +
    labs(title = "Figure S2", 
         subtitle = "Hypothesis-driven approach", 
         x = "Years") +
    theme_pander() + 
    theme(axis.line = element_line(size = 0.1), 
          legend.position = 'none')
    
p_hyp <- ggdraw(add_sub(p_hyp, "p = 0.043", 
                        vpadding=grid::unit(0, "lines"), 
                        y = 17, 
                        x = 0.55, 
                        hjust = 0, 
                        size = 10))
p_hyp

## ------------------------------------------------------------------------
broom::tidy(survival::coxph(srv_obj ~ erbb2_hypothesis)) %>% pandoc.table()

## ------------------------------------------------------------------------
plotALL(measure = erbb2_xpr,
        srv = pData(nki_subset), 
        time = "t.dmfs", 
        event = "e.dmfs") +
    labs(title = "Figure S3", subtitle = "survivALL plot output") + 
    theme(plot.title = element_text(hjust = 0))

## ------------------------------------------------------------------------
srvall <- survivALL(measure = erbb2_xpr,
               srv = pData(nki_subset), 
               time = "t.dmfs", 
               event = "e.dmfs", 
               measure_name = "ERBB2") 

srvall[which.min(srvall$p),] %>% pandoc.table()

## ------------------------------------------------------------------------
erbb2_data <- ifelse(srvall$clsf == 0, "low", "high") 
srv_obj <- survival::Surv(as.numeric(srvall$event_time), srvall$event)
data_fit <- survival::survfit(srv_obj ~ erbb2_data)

p_data <- GGally::ggsurv(data_fit, surv.col = c("#525252", "#bdbdbd")) +
    ylim(0.4, 1) +
    labs(title = "Figure S4", subtitle = "Data-drive approach", x = "Years") +
    theme_pander() + 
    theme(axis.line = element_line(size = 0.1), 
          legend.position = 'none')

p_data <- ggdraw(add_sub(p_data, 
                         "p = 0.001", 
                         vpadding=grid::unit(0, "lines"), 
                         y = 17, 
                         x = 0.55, 
                         hjust = 0, 
                         size = 10))
p_data

## ------------------------------------------------------------------------
broom::tidy(survival::coxph(srv_obj ~ erbb2_data)) %>% pandoc.table()

