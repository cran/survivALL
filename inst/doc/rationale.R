## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", warning = FALSE, message = FALSE, fig.align = 'center', results = 'asis', fig.show = 'hold', fig.width = 5, fig.height = 3.5)

## ------------------------------------------------------------------------
library(survivALL)
library(Biobase)
library(ggplot2)
library(magrittr)
library(viridis)
library(survival)
library(survcomp)

## ------------------------------------------------------------------------
data(nki_subset)

## ---- eval = FALSE-------------------------------------------------------
#  pData(nki_subset)[1:3, ]

## ---- echo = FALSE-------------------------------------------------------
library(pander)
library(magrittr)
pData(nki_subset)[1:3, ] %>% pandoc.table(caption = "Survival Information")

## ------------------------------------------------------------------------
erbb2_xpr <- exprs(nki_subset)["NM_004448",] #ERBB2 expression vector
erbb2_med <- ifelse(erbb2_xpr >= median(erbb2_xpr), "high", "low") #convert to binary classifier

## ------------------------------------------------------------------------
srv_obj <- survival::Surv(nki_subset$t.dmfs, nki_subset$e.dmfs)

broom::tidy(survival::coxph(srv_obj ~ erbb2_med)) %>% pandoc.table()

median_fit <- survival::survfit(srv_obj ~ erbb2_med)
GGally::ggsurv(median_fit) + ggtitle("ERBB2 median")

## ------------------------------------------------------------------------
erbb2_hypothesis <- ifelse(erbb2_xpr >= quantile(erbb2_xpr, probs = 0.8), "high", "low") 

hypothesis_fit <- survival::survfit(srv_obj ~ erbb2_hypothesis)
GGally::ggsurv(hypothesis_fit) + ggtitle("ERBB2 hypothesis-driven")

broom::tidy(survival::coxph(srv_obj ~ erbb2_hypothesis)) %>% pandoc.table()

## ------------------------------------------------------------------------
forest_dfr <- rbind(
                    data.frame(survcomp::hazard.ratio(erbb2_med, nki_subset$t.dmfs, nki_subset$e.dmfs)[1:6]),
                    data.frame(survcomp::hazard.ratio(erbb2_hypothesis, nki_subset$t.dmfs, nki_subset$e.dmfs)[1:6])
                    )
forest_dfr$stratification <- c("median", "hypothesis\ndriven")

ggplot(forest_dfr, aes(x = stratification)) + 
    geom_hline(yintercept = 1, linetype = 3) + 
    geom_linerange(aes(ymin = lower, ymax = upper)) + 
    geom_point(aes(y = hazard.ratio, colour = p.value)) + 
    scale_color_viridis(direction = -1, breaks = c(0.05, 0.5, 1), limits = c(0, 1)) +
    coord_flip()

## ------------------------------------------------------------------------
plotALL(measure = erbb2_xpr,
        srv = pData(nki_subset), 
        time = "t.dmfs", 
        event = "e.dmfs", 
        title = "ERBB2 Expression") 

## ------------------------------------------------------------------------
dfr <- survivALL(measure = erbb2_xpr,
               srv = pData(nki_subset), 
               time = "t.dmfs", 
               event = "e.dmfs", 
               measure_name = "ERBB2") 

dfr[which.min(dfr$p),] %>% pandoc.table()

## ------------------------------------------------------------------------
erbb2_data <- ifelse(erbb2_xpr >= dfr[which.min(dfr$p),]$measure, "high", "low") #convert to binary classifier
data_fit <- survival::survfit(srv_obj ~ erbb2_data)
GGally::ggsurv(data_fit) + ggtitle("ERBB2 data-driven")

