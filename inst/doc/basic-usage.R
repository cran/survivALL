## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", warning = FALSE, message = FALSE, fig.align = 'center', results = 'asis', fig.show = 'hold', fig.width = 7, fig.height = 5)

## ------------------------------------------------------------------------
library(survivALL)
library(Biobase)
library(ggplot2)

data(nki_subset)

## ------------------------------------------------------------------------
xpr_vec <- exprs(nki_subset)["NM_020974", ] #expression vector for SCUBE2 (anti-correlated with proliferation)

plotALL(
        measure = xpr_vec, #expression data
        srv = pData(nki_subset), #survival information
        time = "t.dmfs", #time-to-outcome
        event = "e.dmfs", #outcome type
        bs_dfr = c(), #thresholding data would go here
        measure_name = "SCUBE2", #our gene's name
        title = "SCUBE2 prognostic capacity in a mixed\npopulation of invasive breast cancer samples", #plot title
        )

## ------------------------------------------------------------------------
a_random_x_axis_value <- 123

plotALL(measure = xpr_vec, 
        srv = pData(nki_subset), 
        time = "t.dmfs", 
        event = "e.dmfs", 
        bs_dfr = c(),
        measure_name = "SCUBE2", 
        title = "SCUBE2 prognostic capacity in a mixed\npopulation of invasive breast cancer samples") + 
    geom_vline(xintercept = a_random_x_axis_value, linetype = 5)

## ---- fig.width = 8------------------------------------------------------
geneset <- data.frame(refseq_id = c("NM_020974", "NM_002051", "NM_004448"), hgnc_id = c("SCUBE2", "GATA3", "ERBB2"), stringsAsFactors = FALSE)

xpr_lst <- lapply(geneset$refseq_id, function(id){
                exprs(nki_subset)[id,]
        })
names(xpr_lst) <- geneset$hgnc_id

plot_lst <- lapply(geneset$hgnc_id, function(id){
                       plotALL(
                               measure = xpr_lst[[id]], #expression data
                               srv = pData(nki_subset), #survival information
                               time = "t.dmfs", #time-to-outcome
                               event = "e.dmfs", #outcome type
                               bs_dfr = c(), #thresholding data 
                               measure_name = id, #our gene's name
                               title = id #plot title
                               ) + 
                           ylim(-2.5, 2.5)  
        })

## ---- eval = FALSE-------------------------------------------------------
#  cowplot::plot_grid(plotlist = plot_lst, nrow = 1)

## ---- echo = FALSE, fig.height = 5, fig.width = 11-----------------------
p <- cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
print(p)

## ------------------------------------------------------------------------
survivall_out <- survivALL(
                           measure = xpr_vec, #expression data
                           srv = pData(nki_subset), #survival information
                           time = "t.dmfs", #time-to-outcome
                           event = "e.dmfs", #outcome type
                           bs_dfr = c(), #thresholding data
                           measure_name = "SCUBE2" #our gene's name
                           )

## ---- eval = FALSE-------------------------------------------------------
#  head(survivall_out)

## ---- echo = FALSE-------------------------------------------------------
library(pander)
library(magrittr)
head(survivall_out) %>% pandoc.table()

## ------------------------------------------------------------------------
survivall_lst <- lapply(geneset$hgnc_id, function(id){
                            survivALL(
                                      measure = xpr_lst[[id]], #expression data
                                      srv = pData(nki_subset), #survival information
                                      time = "t.dmfs", #time-to-outcome
                                      event = "e.dmfs", #outcome type
                                      bs_dfr = c(), #thresholding data
                                      measure_name = id #our gene's name
                                      )
                           })

survivall_dfr <- do.call(rbind, survivall_lst)

ggplot(survivall_dfr, aes(x = index, y = HR, colour = name)) + 
    geom_hline(yintercept = 0, linetype = 3) + 
    geom_point() + 
    ylim(-2.5, 2.5) + 
    ggthemes::theme_pander()

