#!/home/jeff/linux_x86/bin/convsh

# 
# NAME:
#       Interpolate to 2.5x2.5
#
# PURPOSE:
#       Interpolate a given input netCDF file to a 2.5x2.5 regular grid (144x72),
#       using either bilinear or area-weighted interpolation.
#
# CALLING SEQUENCE:
#       interpolate_to_2.5x2.5.tcl -i inputfile -o outputfile -f fieldnum -[abm] -x tempfile
#
# ARGUMENTS AND SWITCHES:
#       -i : Defines the full path of the file(s) to be interpolated.
#       -o : Defines the full path of the output netCDF file 
#            (one output file only)
#       -f : The number of the field(s) in the input file(s) to be
#            interpolated (zero-based; default is field "0")
#       -b : Use bilinear interpolation
#       -a : Use area-weighted interpolation
#            (cannot be used with -m, but can be used with -x)
#       -x : Extrapolate over missing data before interpolation
#            (also sets -m).  Specify the path to a temporary file to
#            be created for the intermediate extrapolation step.
#       -m : Field contains missing data (assumed if -x is set)
#            (cannot be used with -a, but can be used with -b)
#       -h : Print help message
#
#
# RESTRICTIONS:
#       None
#
# SIDE EFFECTS:
#       Will overwrite any existing netCDF file with name [outputfile] or [tempfile]
#
# MODIFICATION HISTORY:
#       0.3 - Added help message and better comments; also fixed missing value handling
#             so that it actually works (NPK) (29/04/06 Euro)
#       0.2 - Added ability to extrapolate over missing data before interpolating (NPK) (22/04/06 Euro)
#       0.1 - Written by Nicholas Klingaman (21/04/06 Euro)
#

# Set version number
set version 0.3

# Set output format to netCDF
set outputformat netcdf

# Automatically work out input file type
set filetype 0

# Get command line arguments
set i false; set o false; set f false; set a false; set b false; \
    set m false; set x false
set bilinear false; set area_weighted false
set missing false; set extrapolate false

foreach arg $argv {
    switch -glob -- $arg {
	-i {set i true; set o false; set f false; set x false}
	-o {set i false; set o true; set f false; set x false}
	-f {set i false; set o false; set f true; set x false}
	-x {set i false; set o false; set f false; set x true}	
	-a - -b - -m - -h {set i false; set o false; set f false; set x false}
	-* {puts "unknown option $arg" ; set i false ; set o false; \
		set f false; set x false}
	default {
	    if {$i} {
		set infile [lappend infile $arg]
	    } elseif {$o} {
		set outfile [lappend outfile $arg]	
	    } elseif {$f} {
		set fieldlist [lappend field $arg]
	    } elseif {$x} {
		set extrapolate true
		set tempfile $arg
	    } else {
		puts "unknown option $arg"
	    }
	}
    }
}

set location [lsearch $argv "-h"]
if {$location >= 0} {
    puts " "
    puts "Interpolate to 2.5x2.5 (Convsh script) Version $version.  Written by Nicholas Klingaman."
    puts " "
    puts "Command-line options are as follows:"
    puts "-a : Use area-weighted interpolation (incompatible with -b and -m)"
    puts "-b : Use bilinear interpolation (incompatible with -b)"
    puts "-f \[field number(s)\] : Field codes in input netCDF file(s) to interpolate"
    puts "-h : Print this help message"
    puts "-i \[input file(s)\]  : Input netCDF file(s) to interpolate"
    puts "-m : Input field contain missing data that should not be interpolated"
    puts "     (missing_value must be defined in netCDF file)"
    puts "-o \[output file(s)\] : Output netCDF file (one file only)"
    puts "-x \[temporary file\] : Extrapolate over missing data to temporary file before interpolating"
    puts " "
    exit
}
set location [lsearch $argv "-b"]
if {$location >= 0} {
    set bilinear true
}
set location [lsearch $argv "-a"]
if {$location >= 0} {
    set area_weighted true
}
set location [lsearch $argv "-m"]
if {$location >= 0} {
    set missing true
}
set location [lsearch $argv "-x"]
if {$location >= 0} {
    set extrapolate true
}

if {! [info exists infile]} {
    puts "$argv0 : You must give at least one input filename (-i option)."
    puts $bilinear
    puts $area_weighted
    exit
}

if {$extrapolate && ! [info exists tempfile]} {
    puts "$argv0 : You specified the -x option, but didn't specify a "
    puts "$argv0 : location for the temporary file."
    exit
}

if {[info exists outfile]} {
    if {[llength $outfile] > 1} {
	set outfile [lindex $outfile 0]
	puts "$argv0 : You may only output to one file (-o option)."
    }
} else {
    puts "$argv0 : You must specify an output file (-o option)."
    exit
}

if {$bilinear && $area_weighted} {
    puts "$argv0 : You may not specify both bilinear and area-weighted interpolation."
    exit
} elseif {! $bilinear && ! $area_weighted } {
    puts "$argv0 : You must specify either bilinear (-b) or area-weighted (-a) interpolation."
    puts $bilinear
    puts $area_weighted
    exit
}

if {$area_weighted && $missing} {
    puts "$argv0 : You may not specify both area-weighted interpolation (-a) and missing data (-m)."
    puts "$argv0 : Use extrapolation over missing data (-x) instead."
    exit
}

# Read in each of the input files
foreach file $infile {
    readfile $filetype $file
}

if {$extrapolate} {
    # Extrapolate over missing data and save results to a temporary file before continuing
    puts "$argv0 : You have selected the -x option, so this program will now"
    puts "$argv0 : extrapolate over missing data before interpolating to 2.5x2.5."
    puts "$argv0 : A temporary file will be created at the location you specified."
    puts "$argv0 : You may delete this temporary file after the program finishes."
    setgridtype regular
    extrap $fieldlist 1
    writefile $outputformat $tempfile $fieldlist
    set infile $tempfile
}

set nsize 1000

# Set up 576x434 grid
#set nx 576
#set x0 0.0
#set dx 1.25
#set ny 434
#set y0 90.0
#set dy -0.4166666666

# Set up 2.5x2.5 grid
set nx 144
set x0 1.25
set dx 2.5
set ny 73
set y0 90
set dy -2.5

# Transform to a regular grid
setgridtype regular
# Define x-axis
setgrid 1 $nx $x0 $dx
# Define y-axis
setgrid 2 $ny $y0 $dy

if {$area_weighted} {
    setinterp area_weighted
} else {
    setinterp bilinear
    if {$missing} {
	setmiss true
    } else {
	setmiss false
    }
}

# Interpolate data from input grid to 2.5x2.5 output grid
interp_grid $fieldlist

# Write out all input fields to a single netCDF file
writefile $outputformat $outfile $fieldlist

