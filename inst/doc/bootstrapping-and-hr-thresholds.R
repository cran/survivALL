## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", warning = FALSE, message = FALSE, fig.align = 'center', results = 'asis', fig.show = 'hold', fig.width = 7, fig.height = 5)

## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", warning = FALSE, message = FALSE, fig.align = 'center', fig.show = 'hold', fig.width = 5, fig.height = 3.5)

## ------------------------------------------------------------------------
library(survivALL)
library(Biobase)
library(knitr)

## ------------------------------------------------------------------------
data(nki_subset)

#bootstrapping data should be in the format of 1 repeat per column
bs_mtx <- matrix(nrow = ncol(nki_subset), ncol = 20)

system.time(
            for(i in 1:ncol(bs_mtx)){
                bs_mtx[, i] <- allHR(measure = sample(1:ncol(nki_subset), 
                                                      replace = TRUE),
                                     srv = pData(nki_subset),
                                     time = "t.dmfs",
                                     event = "e.dmfs")
            }
)

kable(bs_mtx[1:20, 1:5])

