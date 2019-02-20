#!/usr/bin/env Rscript

program_name = 'ncspellR.R'
description = 'Script to calculate the yearly spell count for an index from a monthly netCDF file.
The output will contain the number of events for each year in the input.
A spell event starts when two consecutive timesteps go under a given threshold, and stops when another higher threshold is surpassed.
The year of any given event is marked as when the event started, not where it ended.'
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
    'raster'
)

#============= INITIALIZATION =============

fatalerror = function() stop('fatal error', call.=FALSE)
if (!'pacman' %in% rownames(installed.packages())) {
    cat('Please first install package pacman from the R command line with:
        install.packages(pacman)
    pacman will take care of the installation of other necessary packages
    ')
    fatalerror()
}
suppressPackageStartupMessages(library(pacman))
options("pac_update" = FALSE) # Do not try to update relevant packages on package install with p_load
missing = required_pkgs[!p_isinstalled(required_pkgs)]
if (length(missing) > 0L) {
    missing = paste0(missing, collapse = ', ')
} else {
    missing = 'none, well done!'
}
suppressPackageStartupMessages(p_load(optparse))
suppressPackageStartupMessages(p_load(glue))
suppressPackageStartupMessages(p_load(futile.logger))
flog.fatal = function(...) {futile.logger::flog.fatal(...); fatalerror()}

option_list = list( make_option(c("-n", "--nthreads"),
                                type="integer",
                                default=NA,
                                help="Number of threads to be used. Set to NA for automatic detection (all available threads will be used) [default: %default]"),
                    make_option(c("-v", "--varin"),
                                type="character",
                                default='SPI',
                                help="Variable name to use from input file [default: %default]"),
                    make_option(c("-o", "--varout"),
                                type="character",
                                default='number_of_observations',
                                help="Variable name to use in output file [default: %default]"),
                    make_option("--eventstart",
                                type="double",
                                default=-1,
                                help="Threshold for starting an event. [default: %default]"),
                    make_option("--eventend",
                                type="double",
                                default=0,
                                help="Threshold for ending an event. [default: %default]"),
                    make_option(c('-l', "--logfile"),
                                type="character",
                                default=NULL,
                                help="Optional file to write logs to. [default: %default]"),
                    make_option("--debug",
                                action="store_true",
                                help="Print additional debug output. This flag is also useful if you want to check that the options were correctly understood [default: %default]")
                    )
parser = OptionParser(
    usage = "%prog [options] INPUT OUTPUT",
    option_list=option_list,
    epilogue=glue("

    #================= DESCRIPTION =================#
    {program_name} version {version} from {author} ({contact})

    {description}
    {further_description}

    REQUIRED PACKAGES: {paste(required_pkgs, collapse=', ')}
    MISSING: {missing}

    Note that {program_name} will *attempt* to install any missing packages by itself.

    #================= GET IN TOUCH ================#
    For feature requests and bug reports, please get in touch on github:
    {gh_url}/issues

    ")
)

#============= INPUT ARGUMENTS =============

arguments = parse_args(parser, positional_arguments = 2) # , args = c('input/example_input_std_short.nc', 'output/example_output_std_short.nc'))
opt = arguments$options

debug = isTRUE(opt$debug)
if (debug) invisible(flog.threshold(DEBUG))

logfile = opt$logfile
if (!is.null(logfile)) {
    if (!file.exists(logfile)) ftry(file.create(logfile))
    if (file.access(logfile, mode=2) == -1) {
        flog.fatal('cannot write to log file, check permissions and path')
    }
    flog.appender(appender.tee(logfile))
}

flog.debug("Parsing input arguments")

fn_in = arguments$args[1]
flog.debug("Input file: %s", fn_in)
if (!file.exists(fn_in)) flog.fatal('input file does not exist')
if (file.access(fn_in, mode=4) == -1) flog.fatal('input file cannot be read, check permissions and path')

fn_out = arguments$args[2]
flog.debug("Output file: %s", fn_out)
if (file.exists(fn_out)) flog.fatal('input file already exists')
if (file.access(dirname(fn_in), mode=2) == -1) flog.fatal('output file cannot be created, check permissions and path')

var_in = opt$varin
flog.debug("Input variable name: %s", var_in)
var_out = opt$varout
flog.debug("Output variable name: %s", var_out)

nthreads = opt$nthreads
if (is.na(nthreads)) nthreads = parallel::detectCores()
if (!is.integer(nthreads)) flog.fatal('Number of threads must be integer, got "%s".', nthreads)
if (nthreads < 1) flog.fatal('Number of threads must be greater than 1, got %d.', nthreads)
flog.debug("Number of threads: %d", nthreads)

estart = opt$eventstart
eend   = opt$eventend
if (eend <= estart) flog.fatal('The event start threshold must be lower than the event end threshold, got %g and %g respectively.', estart, eend)
flog.debug("Event start threshold: %g", estart)
flog.debug("Event end threshold: %g", eend)

#============= PACKAGE LOAD =============

flog.info("Loading packages")
suppressPackageStartupMessages(p_load(required_pkgs, character.only = TRUE))
flog.debug("Loaded packages")

#============= READ INPUT DATA =============

flog.info("Reading metadata")
nc_in = suppressWarnings(fn_in %>% read_stars(sub = var_in, proxy = TRUE))

if (!st_dimensions(nc_in)$time$refsys %in% c('PCICt', 'POSIXct', 'POSIXlt')) {
    flog.fatal('Input calendar not understood. Are you sure your input file follows the CF Conventions?')
}
times = suppressWarnings(st_get_dimension_values(nc_in, 'time') %>% as.POSIXct %>% round('day'))
years = year(times)
r_years = years %>% rle
lr_years = r_years$lengths %>% length
check12 = all( r_years$lengths[c(-1, -lr_years)] == 12 )
if (!check12) flog.fatal('Some internal years were detected to not have exactly 12 timesteps. Missing any timesptes?')
first_year = years[1]
last_year = years[length(years)]
if (any( unique(years) != (first_year:last_year) )) flog.fatal('Years appear to not be consecutive. Missing any years?')
flog.debug("First year in the file: %d", first_year)
flog.debug("Last  year in the file: %d", last_year)

flog.info("Reading data")
nc_in = suppressWarnings(fn_in %>% read_stars(sub = var_in, proxy = FALSE))

#============= COMPUTE =============

spi_spell = function(v, yrs, eventstart = -1, eventend = 0) {
    stopifnot( length(v) == length(yrs) )
    # Define what the labels mean, for clarity
    A = 1; B = 2; C = 3
    # Categorize input values: lower than eventstart (A), higher than eventend (C) and in-between (B)
    values = as.numeric(cut(v, c(-Inf, eventstart, eventend, Inf), labels = c(A, B, C)))
    # Split values and years into separate groups, divided by a C
    idx <- 1 + cumsum( values == C )
    notC <- !( values == C )
    v_groups = split( values[notC], idx[notC] )
    y_groups = split( yrs[notC], idx[notC] )
    if (length(v_groups) == 0L) return(rep(0, length(unique(yrs))))
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
    dd = rbind(df, data.frame(Var1 = base::setdiff(unique(yrs), df$Var1), Freq=0))
    dd[with(dd, order(Var1)), ]$Freq
}

# TEST constants
# v = runif(60, -3, 1)
# years = rep(2000:2100, length.out = length(v), each=12)
# data.frame(v, years)
# spi_spell(v, years)

if (nthreads > 1 ) {
    p_load(parallel)
    cluster = makeCluster(nthreads)
} else {
    cluster = NULL
}

flog.info("Starting computation using %d threads", nthreads)
spell_res = st_apply(nc_in,
    which((st_dimensions(nc_in) %>% names) != 'time'),
    spi_spell, CLUSTER = cluster, PROGRESS = TRUE,
    yrs = years, eventstart = estart, eventend = eend
)
if (nthreads > 1 ) stopCluster(cluster)
names(spell_res) = var_out
flog.info("Ended computation")

# For some reason it is necessary to permutate the output. Further checks are necessary to verify this is the correct permutation # TODO
shift_vec = function(v) {
    len = length(v)
    return(c(v[2:len], v[1]))
}
spell_res = aperm(spell_res, shift_vec(1:length(dim(spell_res))))

#============= WRITE OUTPUT DATA =============

# Better output writing # TODO

flog.info('Writing to file %s', fn_out)
spell_res %>%
    as('Raster') %>%
    setZ(ymd(paste0(unique(years), '06', '15'))) %>%
    writeRaster(fn_out, varname = var_out, varunit = '1', zname = 'time')

flog.info('All done')
