#!/usr/bin/env Rscript

program_name = 'ncspellR.R'
description = 'Script to calculate the yearly spell count for an index from a monthly netCDF file.
The output will contain the number of events for each year in the input.
A spell event starts when two consecutive timesteps go under a given threshold, and stops when another higher threshold is surpassed.
The year of any given event is marked as when the event started, not when it ended.'
further_description = 'Input file MUST be monthly. Does everything in memory, so make sure your dataset fits in memory! \nVery few checks are performed, so also make sure you know what you are doing.\n Input files must follow the CF Conventions >= 1.5 (http://cfconventions.org/).\n This program is parallel by default, but is not capable of crossing node boundaries (cannot currently run on multiple nodes).'
author = 'Adriano Fantini'
version = '0.2.1'
contact = 'afantini@ictp.it'
gh_url = 'https://github.com/adrfantini/ncspellR'
required_pkgs = c(
    'pacman',
    'optparse',
    'glue',
    'magrittr',
    'futile.logger',
    'lubridate',
    'ncdf4',
    'PCICt',
    'pbapply'
)
version_pkgs = c( # Set requirements for package versions
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
# options(error = function() { flog.fatal(geterrmessage()) ; quit(runLast=FALSE)}) # Override R's default error handling

option_list = list( make_option(c("-n", "--nthreads"),
                                type="integer",
                                default=NA,
                                help="Number of threads to be used. Set to NA for automatic detection (all available threads in the node will be used) [default: %default]"),
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
                    make_option("--progress",
                                action="store_true",
                                help="Show a progress bar - this slightly decreases performance"),
                    make_option("--assume_monthly",
                                action="store_true",
                                help="Assume the input file has the correct monthly periodicity, and only use the time of the first timestep to define times"),
                    make_option("--compress",
                                action="store_true",
                                help="Activate netCDF compression (with deflate level 1) for the main variable"),
                    make_option("--debug",
                                action="store_true",
                                help="Print additional debug output. This flag is also useful if you want to check that the options were correctly understood"),
                    make_option("--dryrun",
                                action="store_true",
                                help="Perform a dry run, do not compute nor write anything to file")
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
    MISSING  PACKAGES: {missing}

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
        flog.fatal('Cannot write to log file, check permissions and path')
    }
    invisible(flog.appender(appender.tee(logfile)))
    if (file.exists(logfile)) flog.debug('Appending output to logfile %s', logfile)
}

dryrun = isTRUE(opt$dryrun)

flog.info(glue(' ### Starting {program_name} version {version} from {author} ({contact}) ### '))

progress = isTRUE(opt$progress)
if (progress) {
    flog.debug('Progress bar will be shown')
} else {
    flog.debug('Progress bar will not be shown')
}

assume_mon = isTRUE(opt$assume_monthly)
if (assume_mon) {
    flog.debug('Assuming the input file has the correct monthly periodicity, and only using the time of the first timestep to define times')
}

compress = isTRUE(opt$compress)
deflate_level = 1
if (compress) flog.debug('Deflate compression activated')

flog.debug("Parsing input arguments")

fn_in = arguments$args[1]
flog.debug("Input file: %s", fn_in)
if (!file.exists(fn_in)) flog.fatal('Input file does not exist')
if (file.access(fn_in, mode=4) == -1) flog.fatal('Input file cannot be read, check permissions and path')

fn_out = arguments$args[2]
flog.debug("Output file: %s", fn_out)
if (file.exists(fn_out)) flog.fatal('Output file already exists')
if (file.access(dirname(fn_in), mode=2) == -1) flog.fatal('Output file cannot be created, check permissions and path')

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
version_pkgs %>% names %>% sapply(function(x) {
    ver = p_ver(x)
    if ( ver < version_pkgs[x] ) flog.fatal('Version %s of package %s required, got %s. Maybe your (probably CRAN) version is behind the github one?', version_pkgs[x], x, as.character(ver))
}) %>% invisible
flog.debug("Loaded packages")

#============= READ INPUT DATA =============

flog.info("Reading metadata")
nc_in = nc_open(fn_in)

# Read time
time_var = 'time'
times = nc_in %>% ncvar_get(time_var)
time_units = nc_in %>% ncatt_get(time_var, 'units')
if (!time_units$hasatt) flog.fatal('Cannot find time units!')
flog.debug('Time units: %s', time_units$value)
time_units = strsplit(time_units$value, " ")[[1]]
flog.debug('Time unit parsed as:',  data.frame(WHAT = c('unit', 'since', 'date', 'time', 'timezone'), VALUE = c(time_units, rep(NA, 5 - length(time_units)))), capture=TRUE)

time_cal = nc_in %>% ncatt_get(time_var, 'calendar')
if (!time_cal$hasatt) {
    flog.info('No input calendar, assuming standard')
    time_cal = 'gregorian'
} else {
    time_cal = time_cal$value %>% tolower
    flog.debug('Input calendar: %s', time_cal)
}

pcict_calendars = c('365', '365_day', 'noleap', '360', '360_day', 'gregorian', 'proleptic_gregorian', 'standard')
if (time_cal %in% pcict_calendars) {
    # Use PCIct to deal with times
    time_offset_unit <- time_units[1]
    time_tz <- time_units[5]
    time_start <- strsplit(time_units[4], ":")[[1]]
    if (length(time_start) != 3 || time_start[1] > 24 || time_start[2] > 60 || time_start[3] > 60 || any(time_start < 0)) {
        flog.warn("%s is not a valid start time. Assuming 00:00:00", time_start)
        time_units[4] <- "00:00:00"
    }
    if (! time_tz %in% OlsonNames()) {
        flog.debug("%s is not a valid timezone. Assuming UTC", time_tz)
        time_tz <- "UTC"
    }
    time_start <- ymd_hms(paste(time_units[3], time_units[4]), tz=time_tz)
    # Find the correct lubridate time function based on the unit
    time_f <- switch(tolower(time_offset_unit),
        seconds=seconds, second=seconds, sec=seconds,
        minutes=minutes, minute=minutes, min=minutes,
        hours=hours,     hour=hours,     h=hours,
        days=days,       day=days,       d=days,
        months=months,   month=months,   m=months,
        years=years,     year=years,     yr=years,
        NA
    )
    if ( grepl("month", tolower(time_offset_unit)) && ( time_cal %in% c('365', '365_day', 'noleap', '360', '360_day') ) ) {
        second_offsets = switch(time_cal,
            '365' = (365/12) * 24 * 3600,
            '365_day' = (365/12) * 24 * 3600,
            'noleap' = (365/12) * 24 * 3600,
            '360' = 30 * 24 * 3600,
            '360_day' = 30 * 24 * 3600
        ) * floor(times)
    } else {
        second_offsets = as.numeric(time_f(floor(times))) # TODO this floor here is a hack, to work around files which have non-integer times. Works in most cases.
    }
    times = as.PCICt(time_start, cal=time_cal) + second_offsets
} else {
    flog.fatal('This program does not support calendar ""%s"', time_cal)
}

# Simplify dates, since we only care about months
times = suppressWarnings(times %>% as.POSIXct %>% round('day'))
if (assume_mon) {
    flog.debug('Assuming monthly data, since --assume_monthly was set')
    expected_times = times[1] + months(1:length(times)-1)
    times = expected_times
}
years = year(times)
r_years = years %>% rle
lr_years = r_years$lengths %>% length
check12 = all( r_years$lengths[c(-1, -lr_years)] == 12 )
if (!check12) flog.fatal('Some internal years were detected to not have exactly 12 timesteps. Missing any month?')
first_year = years[1]
last_year = years[length(years)]
if (any( unique(years) != (first_year:last_year) )) flog.fatal('Years appear to not be consecutive. Missing any years?')
flog.debug("First year in the file: %d", first_year)
flog.debug("Last  year in the file: %d", last_year)

# Function to get dimensions for a given variable
ncdim_get = function(nc, varid) {
    nc$var[[varid]]$dim %>% sapply(`[[`, 'name')
}

# Inspect input variable to read
nc_var_in = nc_in$var[[var_in]]
if (is.null(nc_var_in)) flog.fatal('Input file %s does not contain variable %s', fn_in, var_in)
nc_var_in_dims = nc_in %>% ncdim_get(var_in)

# Read data
flog.info("Reading data")
spi_in = nc_in %>% ncvar_get(var_in)

#============= COMPUTE =============

if (dryrun) {
    flog.info('Quitting, this was just a dry run')
    quit()
}

spi_spell = function(v, yrs, eventstart = -1, eventend = 0) {
    stopifnot( length(v) == length(yrs) )
    # Define what the labels mean, for clarity
    A = 1; B = 2; C = 3
    # Categorize input values:
    # (A) lower than eventstart, start an event
    # (B) higher than eventend, stop an event if 2 consecutives are found
    # (B) in-between, event status not affected
    # We treat any NA as Bs
    v[is.na(v)] <- (eventstart + eventend) / 2
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
if (progress) {
    pboptions('type' = 'timer')
} else {
    pboptions('type' = 'none')
}
spell_out = pbapply(spi_in,
    which(nc_var_in_dims != 'time'),
    spi_spell, cl = cluster,
    yrs = years, eventstart = estart, eventend = eend
)
if (nthreads > 1 ) stopCluster(cluster)
flog.info("Ended computation")

# Repermutate the output which has been swapped by apply
shift_vec = function(v) {
    len = length(v)
    return(c(v[2:len], v[1]))
}
spell_out = aperm(spell_out, shift_vec(1:length(dim(spell_out))))

#============= WRITE OUTPUT DATA =============

# Get dimensions and variables to be cloned from the source file
# coordinates and grid_mapping are special variables that act as pointers to other variables/dimensions
nc_coords = nc_in %>% ncatt_get(var_in, 'coordinates')
nc_grid_mapping = nc_in %>% ncatt_get(var_in, 'grid_mapping')
vars2copy = NULL
if (nc_coords$hasatt) {
    flog.debug('Found "coordinates" attribute for variable %s', var_in)
    vars2copy = nc_coords$value %>% strsplit(' ') %>% unlist
}
if (nc_grid_mapping$hasatt){
    flog.debug('Found "grid_mapping" attribute for variable %s', var_in)
    vars2copy = c(vars2copy, nc_grid_mapping$value %>% strsplit(' ') %>% unlist)
}
dims2copy = c(nc_var_in_dims, lapply(vars2copy, function(v) nc_in %>% ncdim_get(v)) %>% unlist) %>% unique
flog.debug("Variables  to be cloned: {paste(vars2copy, collapse=', ')}" %>% glue)
flog.debug("Dimensions to be cloned: {paste(dims2copy, collapse=', ')}" %>% glue)

# Creating output file
flog.info('Initialising output file')
# Remove time from output dimensions since we need to redefine out own
nc_var_out_dims = nc_var_in_dims[which(nc_var_in_dims != 'time')]
nc_var_out_time = ncdim_def(
    name = 'time',
    units = 'years since {first_year}-06-15 00:00:00' %>% glue,
    vals = unique(years) - first_year,
    unlim = TRUE,
    create_dimvar = TRUE,
    calendar = 'standard',
    longname = 'time'
)

nc_var_out = ncvar_def(
    var_out, '1', c(nc_in$dim[nc_var_out_dims], list(nc_var_out_time)),
    longname='Yearly spell count',
#     shuffle = compress,
    compression = ifelse(compress, deflate_level, NA)
)
# It is here necessary to reset chunksizes to NA, in some cases, possibly due to a bug in ncdf4
nc_in_vars2copy = lapply(nc_in$var[vars2copy], function(x) {x$chunksizes = NA; x})
nc_out = fn_out %>% nc_create(c(list(nc_var_out), nc_in_vars2copy), force_v4 = TRUE)

# Wrapper for ncatt_put to log in debug
ncatt_put = function(nc, varid, attname, attval, ...) {
    flog.debug('Setting %s attribute: %s = %s', ifelse(varid == 0, 'global', varid), attname, attval)
    ncdf4::ncatt_put(nc, varid, attname, attval, ...)
}

# Add time standard name
nc_out %>% ncatt_put('time', 'standard_name', 'time')

# Clone all the attributes from the relevant variables and dimensions in the source file, excluding time
for (v in c('global attributes', dims2copy[which(dims2copy != 'time')], vars2copy[which(vars2copy != 'time')])) {
    if (v == 'global attributes') v = 0
    if (v == 0) {
        flog.debug('Getting global attributes')
    } else {
        flog.debug('Getting attributes for variable %s', v)
    }

    atts = nc_in %>% ncatt_get(v)
    for (a in names(atts)) {
        nc_out %>% ncatt_put(v, a, atts[[a]])
    }
}

# Update netCDF history
nc_h = ncatt_get(nc_in, 0, 'history')
nc_new_h = '{date()}: yearly spell count calculated by {program_name} version {version} ({gh_url})' %>% glue
if (nc_h$hasatt) {
    nc_out %>% ncatt_put(0, 'history', paste(nc_new_h, nc_h$value, sep='\n'))
} else {
    nc_out %>% ncatt_put(0, 'history', nc_new_h)
}

# Fill in cloned variables
for (v in vars2copy) {
    flog.debug('Copying values for cloned variable %s', v)
    nc_out %>% ncvar_put(v, nc_in %>% ncvar_get(v))
}
flog.debug('Closing input file')
nc_close(nc_in)

# Clone coordinates and grid_mapping attributes, if present
if (nc_coords$hasatt) nc_out %>% ncatt_put(var_out, 'coordinates', nc_coords$value)
if (nc_grid_mapping$hasatt) nc_out %>% ncatt_put(var_out, 'grid_mapping', nc_grid_mapping$value)

# Fill SPI variable with data and attributes
flog.info('Filling output file')
nc_out %>% ncvar_put(var_out, spell_out)
nc_out %>% ncatt_put(var_out, 'event_start', estart)
nc_out %>% ncatt_put(var_out, 'event_end'  , eend)
nc_out %>% ncatt_put(var_out, 'program', '{program_name} version {version} ({gh_url})' %>% glue)
nc_out %>% ncatt_put(var_out, 'R_version', p_ver(R) %>% as.character)
nc_out %>% ncatt_put(var_out, 'ncdf4_pkg_version', p_ver(ncdf4) %>% as.character)

# Close
flog.debug('Closing output file')
nc_close(nc_out)
flog.info(' ### Goodbye! ### ')

