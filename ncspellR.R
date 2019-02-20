#!/usr/bin/env Rscript

program_name = 'ncspellR.R'
description = 'Script to calculate the yearly spell count for an index from a monthly netCDF file.
The output will contain the number of events for each year in the input. A spell event starts when two consecutive timesteps go under a given threshold, and stops when another higher threshold is surpassed. The year of any given event is marked as when the event started, not where it ended.'
further_description = 'Input file MUST be monthly. Does everything in memory, so make sure your dataset fits in memory! \nVery few checks are performed, so also make sure you know what you are doing.\n Input files must follow the CF Conventions >= 1.5 (http://cfconventions.org/).'
author = 'Adriano Fantini'
version = '0.1'
contact = 'afantini@ictp.it'
gh_url = 'https://github.com/adrfantini/ncspellR'
required_pkgs = c(
    'pacman',
    'optparse',
    'glue',
    'stars',
    'magrittr',
    'futile.logger',
    'lubridate',
    'RNetCDF',
    'ncdf.tools',
    'compiler'
)

# Second approach, with years
library(compiler)
spi_spell = cmpfun(function(v, years, eventstart = -1, eventend = 0) {
    stopifnot( length(v) == length(years) )
    # Define what the labels mean, for clarity
    A = 1; B = 2; C = 3
    # Categorize input values: lower than eventstart (A), higher than eventend (C) and in-between (B)
    values = as.numeric(cut(v, c(-Inf, eventstart, eventend, Inf), labels = c(A, B, C)))
    # Split values and years into separate groups, divided by a C
    idx <- 1 + cumsum( values == C )
    notC <- !( values == C )
    v_groups = split( values[notC], idx[notC] )
    y_groups = split( years[notC], idx[notC] )
    # Loop over groups, checking if there is any AA*C inside (might also be done with a regex?) while keeping track of the year the event begins
    df = as.data.frame(table(sapply(1:length(v_groups), function(i) {
        a = v_groups[[i]]
        b = y_groups[[i]]
        r = rle(a)
        r$values = (r$values == A) & (r$lengths > 1)
        b[ which(inverse.rle(r))[1] ]
    })))
    df$Var1 = as.numeric(levels(df$Var1))
    # Fill in missing years with 0s
    dd = rbind(df, data.frame(Var1 = base::setdiff(unique(years), df$Var1), Freq=0))
    dd[with(dd, order(Var1)), ]
})

v = runif(48, -3, 1)
years = rep(2000:2100, length.out = length(v), each=12)
data.frame(v, years)
spi_spell(v, years)
