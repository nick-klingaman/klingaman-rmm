#!/bin/ksh

# Script to calculate the Wheeler and Hendon (WH04) RMM1 and RMM2 values 
# from netCDF input data.

# For HadGEM3-A AMIP2 run airxv
#run_id=airxv
#run_type=hadgem3a_morph3_final_n96_amip2
#input_directory=${HOME}/um_output3/hadgem3_monwg/airxv
#n_years=27
#year_range='1982-2008'
#ndays_per_year=360
#set -A nc_variables 'olr' 'u' 'u' 'precip'
#obs_start_year=7

# For HadGEM3-A at vn7.4 with blended FOAM/AMIP2 SSTs and no CMT
#run_id=nocmt_blendsst_vn74
#run_type=hadgem3a
#input_directory=${HOME}/um_output3/${run_type}_${run_id}
#n_years=6
#year_range='years1-6'
#ndays_per_year=360
#set -A nc_variables 'olr' 'u' 'u' 'precip'
#obs_start_year=7

# For single HadGEM3-A at vn7.4 with blended FOAM/AMIP2 SSTs and CMT
#run_id=xetug
#run_type=hadgem3a_foam_ctl_vn74_blendsst
#input_directory=$HOME/um_output3/xetug
#n_years=6
#year_range='i0-i5'
#ndays_per_year=360
#set -A nc_variables 'olr' 'u' 'u' 'precip'
#obs_start_year=7

# For single HadGEM3-A at vn7.8 with UKMO clim SST
run_id=xjhwb
run_type=hadgem3kpp_nrglobal_n216
input_directory=$HOME/um_output6/${run_id}
n_years=106
year_range='years1-106'
ndays_per_year=360
set -A nc_variables 'toa_outgoing_longwave_flux' 'eastward_wind' 'eastward_wind' 'PP_1_5216_vn708'
obs_start_year=7

# For NOAA OLR and ERA-Interim winds
#run_id=obs
#run_type=noaa_eraint
#input_directory=/home/ss901165/datasets/MJO_INDICES
#n_years=20
#year_range='1989-2008'
#ndays_per_year=365
#set -A nc_variables 'olr' 'U' 'U'
#obs_start_year=14

# Input file with observed values of RMM1 and RMM2 (for comparison to model data)
obs_rmm_indices_file=/home/ss901165/datasets/MJO_INDICES/MJO_rmm1_rmm2.jan-dec_dmeans.1975-2009.index_values.nc
#obs_start_year=3     # As an offset, e.g., to get the right years to compare to an AMIP run
obs_nyears=25         # Change this if the observations have fewer years than the model 
                      # (or if you want to compare different numbers of years in the model 
                      # and observations, for some other reason)
setup idl idl_guide5

# Input files with daily means of OLR, U850, U200
olr_dmean_file=${input_directory}/${run_type}.jan-dec_dmeans.${year_range}.olr.nc
u200_dmean_file=${input_directory}/${run_type}.jan-dec_dmeans.${year_range}.u200.nc
u850_dmean_file=${input_directory}/${run_type}.jan-dec_dmeans.${year_range}.u850.nc
precip_dmean_file=${input_directory}/${run_type}.jan-dec_dmeans.${year_range}.precip.nc

# Input files with daily climatologies of OLR, U850, U200
olr_clim_file=${input_directory}/${run_type}.jan-dec_dmean_clim.${year_range}.olr.nc
u200_clim_file=${input_directory}/${run_type}.jan-dec_dmean_clim.${year_range}.u200.nc
u850_clim_file=${input_directory}/${run_type}.jan-dec_dmean_clim.${year_range}.u850.nc
precip_clim_file=${input_directory}/${run_type}.jan-dec_dmean_clim.${year_range}.precip.nc

set -A variables 'olr' 'u200' 'u850' 'precip'
set -A field_multipliers '1' '1' '1' '86400'
set -A dmean_files ${olr_dmean_file} ${u200_dmean_file} ${u850_dmean_file} ${precip_dmean_file}
set -A clim_files ${olr_clim_file} ${u200_clim_file} ${u850_clim_file} ${precip_clim_file}
set -A all_input_files ${dmean_files[@]} ${clim_files[@]}
set -A dmean_field_numbers '3' '4' '4' '3'
set -A clim_field_numbers '3' '3' '3' '3'
set -A interpolation_flag '1' '1' '1' '1'

# Convert the input files to 2.5x2.5 resolution
interpolation_script=$HOME/src/convsh/interpolate_to_2.5x2.5_linux.tcl
set -A dmean_files_interpolated
set -A clim_files_interpolated
i=0
for this_file in ${dmean_files[@]} ; do
    output_file=${this_file%%.nc}.2.5x2.5.nc
    if [ ${interpolation_flag[i]} -eq 1 ] ; then
	command=${interpolation_script}' -i '${this_file}' -o '${output_file}' -a -f '${dmean_field_numbers[${i}]} ; echo $command ; $command
    fi
    command='ncrename -d t_1,t -v t_1,t '${output_file} ; echo $command ; $command
    set -A dmean_files_interpolated ${dmean_files_interpolated[@]} ${output_file}
    i=$((i+1))
done
i=0
for this_file in ${clim_files[@]} ; do    
    output_file=${this_file%%.nc}.2.5x2.5.nc
    if [ ${interpolation_flag[i]} -eq 1 ] ; then 
	command=${interpolation_script}' -i '${this_file}' -o '${output_file}' -a -f  '${clim_field_numbers[${i}]} ; echo $command ; $command
    fi
    command='ncrename -d t_1,t -v t_1,t '${output_file} ; echo $command ; $command
    set -A clim_files_interpolated ${clim_files_interpolated[@]} ${output_file}
    i=$((i+1))
done

# Run IDL procedure to remove daily climatology and first three harmonics from daily means
# Run IDL procedure to latitude average 15S-15N and remove mean of previous 120 days
i=0
set -A anomaly_files
 cat << EOF > ${input_directory}/idl.startup
 .compile remove_anncycle_three_harmonics.pro
EOF
for this_variable in ${variables[@]} ; do
    dmean_interpolated_file=${dmean_files_interpolated[i]}
    clim_interpolated_file=${clim_files_interpolated[i]}
    this_anomaly_file=${input_directory}/${run_type}.jan-dec_dmeans_anom-3harm.${year_range}.${this_variable}.2.5x2.5.nc
    this_latavg_anomaly_file=${input_directory}/${run_type}.jan-dec_dmeans_anom-3harm-120d_latavg-15S-15N.${year_range}.${this_variable}.2.5x2.5.nc
    cat <<EOF >> ${input_directory}/idl.startup
 remove_anncycle_three_harmonics,'${dmean_interpolated_file}','${clim_interpolated_file}','${nc_variables[i]}',${ndays_per_year},'${this_anomaly_file}'
 latitude_average_and_remove_120days,'${this_anomaly_file}','${nc_variables[i]}',${ndays_per_year},'${this_latavg_anomaly_file}',/ENSEMBLE
EOF
    set -A anomaly_files ${anomaly_files[@]} ${this_anomaly_file}
    set -A latavg_anomaly_files ${latavg_anomaly_files[@]} ${this_latavg_anomaly_file}
    i=$((i+1))
done
cat <<EOF >> ${input_directory}/idl.startup
exit
EOF
export IDL_STARTUP=${input_directory}/idl.startup
idl

# Run Fortran code to project onto observed WH04 EOFs
rmm_indices_file=${input_directory}/${run_type}.jan-dec_dmeans.${year_range}.rmm_indices.nc
src=${HOME}/src/mjo/projectRMM1RMM2_npk_netcdf.f90
exe=${HOME}/src/mjo/projectRMM1RMM2_npk_netcdf

# sed s/\\\//\\\\\\\\\\\//g
#       1235|123412341235|
#       
# 1 - Escapes the text file
# 2 - Escapes the shell
# 3 - Escapes sed
# 4 - Prints an actual forward slash
# 5 - Prints an actual backward slash

set -A latavg_anomaly_files_forsed
j=0
input_directory_forsed=`echo ${input_directory} | sed s/\\\//\\\\\\\\\\\//g`
set -A latavg_anomaly_files_symlinks 'olr.nc' 'u200.nc' 'u850.nc' 'precip.nc'
for this_file in ${latavg_anomaly_files[@]} ; do
    ln -sf ${this_file} ${input_directory}/${latavg_anomaly_files_symlinks[j]}
    j=$((j+1))
done
rmm_indices_file_symlink='rmm_indices.nc'
ln -sf ${rmm_indices_file} ${input_directory}/${rmm_indices_file_symlink}
src_forsed=`echo ${src} | sed s/\\\//\\\\\\\\\\\//g`

echo 'sed #1'
sed s/n_years=-1/n_years=${n_years}/1 ${src} > ${src}.${run_id}.1 || exit
echo 'sed #2'
sed s/n_days_per_year=-1/n_days_per_year=${ndays_per_year}/1 ${src}.${run_id}.1 > ${src}.${run_id}.2 || exit
sed s/n_days_removed=-1/n_days_removed=120/1 ${src}.${run_id}.2 > ${src}.${run_id}.1 || exit
sed s/olr_file/\'${latavg_anomaly_files_symlinks[0]}\'/1 ${src}.${run_id}.1 > ${src}.${run_id}.2 || exit
sed s/olr_name=olr_name/olr_name=\'${nc_variables[0]}\'/1 ${src}.${run_id}.2 > ${src}.${run_id}.1 || exit
sed s/u200_file/\'${latavg_anomaly_files_symlinks[1]}\'/1 ${src}.${run_id}.1 > ${src}.${run_id}.2 || exit
sed s/u200_name=u200_name/u200_name=\'${nc_variables[1]}\'/1 ${src}.${run_id}.2 > ${src}.${run_id}.1 || exit
sed s/u850_file/\'${latavg_anomaly_files_symlinks[2]}\'/1 ${src}.${run_id}.1 > ${src}.${run_id}.2 || exit
sed s/u850_name=u850_name/u850_name=\'${nc_variables[2]}\'/1 ${src}.${run_id}.2 > ${src}.${run_id}.1 || exit
sed s/output_file/\'${rmm_indices_file_symlink}\'/1 ${src}.${run_id}.1 > ${src}.${run_id}.2 || exit
sed s/directory=directory/directory=\'${input_directory_forsed}\'/1 ${src}.${run_id}.2 > ${src%%.f90}.${run_id}.f90 || exit
compile_script=${HOME}/src/mjo/projectRMM1RMM2_npk_netcdf.compile
sed s/${src_forsed}/${src_forsed%%.f90}.${run_id}.f90/1 ${compile_script} > ${compile_script}.${run_id}
echo ${compile_script}.${run_id} ; chmod u+x ${compile_script}.${run_id} ; ${compile_script}.${run_id}

${exe}

# Run IDL procedure to generate probability plots
setup idl_guide4
cat << EOF > ${input_directory}/idl.startup
 .compile FSC_COLOR
 .compile mjo_indices_wh04_probabilities.pro
 mjo_indices_wh04_probabilities,'${rmm_indices_file}','${run_id}',start_year=0,nyear=${n_years},obs_input_file='${obs_rmm_indices_file}',obs_start_year=${obs_start_year},obs_nyear=${obs_nyears}
 exit
EOF
idl
# Filter the gridded 2.5x2.5 anomalies to 20-80 day periods
set -A filtered_anomaly_files
cat << EOF > ${input_directory}/idl.startup
 .compile lanczos_filter_file.pro
EOF
j=0
for this_variable in ${variables[@]} ; do
    this_input_anomaly_file=${anomaly_files[j]}
    this_filtered_anomaly_file=${this_input_anomaly_file%%.${year_range}.${this_variable}.2.5x2.5.nc}_filter2080.${year_range}.${this_variable}.2.5x2.5.nc
    echo ${this_input_anomaly_file} ${this_filtered_anomaly_file}
    cat << EOF >> ${input_directory}/idl.startup
 lanczos_filter_file,'${this_input_anomaly_file}','${nc_variables[j]}',20,80,'${this_filtered_anomaly_file}',field_multiplier=${field_multipliers[j]},/CONTINUOUS
EOF
    set -A filtered_anomaly_files ${filtered_anomaly_files[@]} ${this_filtered_anomaly_file}
    j=$((j+1))
done
cat << EOF >> ${input_directory}/idl.startup
 exit
EOF
idl 

# Run IDL procedure to generate filtered variability plots
cat << EOF > ${input_directory}/idl.startup
 .compile FSC_COLOR.pro
 .compile mjo_indices_filtered_variance.pro
EOF
j=0
for this_variable in ${variables[@]} ; do
    command='ncrename -v '${nc_variables[j]}','${this_variable}' '${filtered_anomaly_files[j]} ; echo $command ; $command
    command='ncpdq -a time,member,latitude,longitude -O -o '${filtered_anomaly_files[j]}' '${filtered_anomaly_files[j]} ; echo $command ; $command
    cat << EOF >> ${input_directory}/idl.startup
 mjo_indices_filtered_variance,'${filtered_anomaly_files[j]}','${run_id}','${this_variable}',20,80,${n_years},${ndays_per_year},box=[-30,0,30,360]
EOF
    j=$((j+1))
done
cat << EOF >> ${input_directory}/idl.startup
 exit
EOF
idl

# Run Fortran code to compute space-time power spectra
src=${HOME}/src/mjo/space_time_powerspec.f90
exe=${HOME}/src/mjo/space_time_powerspec
compile_script=${HOME}/src/mjo/space_time_powerspec.ksh
compile_script_final=${compile_script%%.ksh}.${run_id}.ksh
sed1=${src%%.f90}.${run_id}.f90.1
sed2=${src%%.f90}.${run_id}.f90.2
src_final=${src%%.f90}.${run_id}.f90
echo 'sed #1'
src_forsed=`echo ${src} | sed s/\\\//\\\\\\\\\\\//g`
echo 'sed #2'
src_final_forsed=`echo ${src_final} | sed s/\\\//\\\\\\\\\\\//g`
set -A spacetime_powerspec_files
j=0
for this_variable in ${variables[@]} ; do
    this_powerspec_outfile=${anomaly_files[j]%%.${year_range}.${this_variable}.2.5x2.5.nc}_space-time_powerspec.${year_range}.${this_variable}.2.5x2.5.nc
    this_powerspec_outfile_symlink=${input_directory}/spacetime_spec_${this_variable}.nc
    this_powerspec_outfile_forsed=`basename ${this_powerspec_outfile_symlink}`
    ln -sf ${this_powerspec_outfile} ${this_powerspec_outfile_symlink}
    sed s/directory=directory/directory=\'${input_directory_forsed}\'/1 ${src} > ${sed1}
    sed s/input_file/\'${latavg_anomaly_files_symlinks[j]}\'/1 ${sed1} > ${sed2}
    sed s/output_file/\'${this_powerspec_outfile_forsed}\'/1 ${sed2} > ${sed1}    
    sed s/varname=variable_name/varname=\'${nc_variables[j]}\'/1 ${sed1} > ${sed2}
    sed s/ndays_per_year=ndays_per_year/ndays_per_year=${ndays_per_year}/1 ${sed2} > ${sed1}
    sed s/nyears=nyears/nyears=${n_years}/1 ${sed1} > ${sed2}
    sed s/field_multiplier=1/field_multiplier=${field_multipliers[j]}/1 ${sed2} > ${src_final}
    sed s/${src_forsed}/${src_final_forsed}/1 ${compile_script} > ${compile_script_final}    
    chmod u+x ${compile_script_final}
    ${compile_script_final}
    ${exe}
    set -A spacetime_powerspec_files ${spacetime_powerspec_files[@]} ${this_powerspec_outfile}
    j=$((j+1))
done

# Run IDL code to plot space-time power spectra
cat << EOF > ${input_directory}/idl.startup
 .compile mjo_indices_spacetime_power_spectra.pro
EOF
j=0
for this_variable in ${variables[@]} ; do
    cat << EOF >> ${input_directory}/idl.startup
 mjo_indices_spacetime_power_spectra,'${spacetime_powerspec_files[j]}','${run_id}','${this_variable}',[-10,0.005,10,0.1]
EOF
    j=$((j+1))
done
cat << EOF >> ${input_directory}/idl.startup
 exit
EOF
export IDL_STARTUP=${input_directory}/idl.startup
idl

# Computation of WK99 power spectra on OLR

echo ' -- '
echo ' -- Computation of WK99 spectra on OLR --'
echo ' -- '

# Step 1. Create segments and divide into symmetric/asymmetric components
src=${HOME}/src/mjo/create_segments.f90
exe=${HOME}/src/mjo/create_segments
compile_script=${HOME}/src/mjo/create_segments.ksh
compile_script_final=${compile_script%%.ksh}.${run_id}.ksh
sed1=${src%%.f90}.${run_id}.f90.1
sed2=${src%%.f90}.${run_id}.f90.2
src_final=${src%%.f90}.${run_id}.f90
echo 'sed #1'
src_forsed=`echo ${src} | sed s/\\\//\\\\\\\\\\\//g`
echo 'sed #2'
src_final_forsed=`echo ${src_final} | sed s/\\\//\\\\\\\\\\\//g`
asym_outfile=${anomaly_files[0]%%.${year_range}.olr.2.5x2.5.nc}_wk99_segments_asym.${year_range}.olr.2.5x2.5.nc
asym_symlink=${input_directory}/wk99_segments_asym.nc
ln -sf ${asym_outfile} ${asym_symlink}
sym_outfile=${anomaly_files[0]%%.${year_range}.olr.2.5x2.5.nc}_wk99_segments_sym.${year_range}.olr.2.5x2.5.nc
sym_symlink=${input_directory}/wk99_segments_sym.nc
ln -sf ${sym_outfile} ${sym_symlink}
ln -sf ${anomaly_files[0]} ${input_directory}/wk99_segments_input.nc
sed s/directory=directory/directory=\'${input_directory_forsed}\'/1 ${src} > ${sed1}
sed s/input_file/\'wk99_segments_input.nc\'/1 ${sed1} > ${sed2}
sed s/output_file_asym/\'wk99_segments_asym.nc\'/1 ${sed2} > ${sed1}
sed s/output_file_sym/\'wk99_segments_sym.nc\'/1 ${sed1} > ${sed2}
sed s/variable_name/\'${nc_variables[0]}\'/1 ${sed2} > ${sed1}
sed s/num_x/144/1 ${sed1} > ${sed2} 
sed s/num_y/73/1 ${sed2} > ${sed1}
sed s/ndpy/${ndays_per_year}/1 ${sed1} > ${sed2}
sed s/num_yr/${n_years}/1 ${sed2} > ${sed1}
sed s/num_s/96/1 ${sed1} > ${sed2}
sed s/sel_y/13/1 ${sed2} > ${sed1}
sed s/yoff/30/1 ${sed1} > ${sed2}
sed s/num_o/60/1 ${sed2} > ${src_final}
sed s/${src_forsed}/${src_final_forsed}/1 ${compile_script} > ${compile_script_final}    
chmod u+x ${compile_script_final}
${compile_script_final}
${exe}

# Step 2. Calculate power
# Step 2.1 Symmetric
src=${HOME}/src/mjo/wk99_power.f90
exe=${HOME}/src/mjo/wk99_power
compile_script=${HOME}/src/mjo/wk99_power.ksh
compile_script_final=${compile_script%%.ksh}.${run_id}.ksh
sed1=${src%%.f90}.${run_id}.f90.1
sed2=${src%%.f90}.${run_id}.f90.2
src_final=${src%%.f90}.${run_id}.f90
echo 'sed #1'
src_forsed=`echo ${src} | sed s/\\\//\\\\\\\\\\\//g`
echo 'sed #2'
src_final_forsed=`echo ${src_final} | sed s/\\\//\\\\\\\\\\\//g`
sym_power_outfile=${anomaly_files[0]%%.${year_range}.olr.2.5x2.5.nc}_wk99_power_sym.${year_range}.olr.2.5x2.5.nc
sym_power_symlink=${input_directory}/wk99_power_sym.nc
ln -sf ${sym_power_outfile} ${sym_power_symlink}
sed s/directory=directory/directory=\'${input_directory_forsed}\'/1 ${src} > ${sed1}
sed s/input_file/\'wk99_segments_sym.nc\'/1 ${sed1} > ${sed2}
sed s/output_file/\'wk99_power_sym.nc\'/1 ${sed2} > ${sed1}
sed s/variable_name/\'olr_sym\'/1 ${sed1} > ${sed2}
sed s/num_s/96/1 ${sed2} > ${sed1}
sed s/num_x/144/1 ${sed1} > ${sed2}
sed s/sel_y/13/g ${sed2} > ${src_final}
sed s/${src_forsed}/${src_final_forsed}/1 ${compile_script} > ${compile_script_final}    
chmod u+x ${compile_script_final}
${compile_script_final}
${exe}

# Step 2.2 Asymmetric
asym_power_outfile=${anomaly_files[0]%%.${year_range}.olr.2.5x2.5.nc}_wk99_power_asym.${year_range}.olr.2.5x2.5.nc
asym_power_symlink=${input_directory}/wk99_power_asym.nc
ln -sf ${asym_power_outfile} ${asym_power_symlink}
sed s/directory=directory/directory=\'${input_directory_forsed}\'/1 ${src} > ${sed1}
sed s/input_file/\'wk99_segments_asym.nc\'/1 ${sed1} > ${sed2}
sed s/output_file/\'wk99_power_asym.nc\'/1 ${sed2} > ${sed1}
sed s/variable_name/\'olr_asym\'/1 ${sed1} > ${sed2}
sed s/num_s/96/1 ${sed2} > ${sed1}
sed s/num_x/144/1 ${sed1} > ${sed2}
sed s/sel_y/13/g ${sed2} > ${src_final}
sed s/${src_forsed}/${src_final_forsed}/1 ${compile_script} > ${compile_script_final}    
chmod u+x ${compile_script_final}
${compile_script_final}
${exe}

# Step 2.3 Normalise and compute background
src=${HOME}/src/mjo/wk99_smoothing.f90
exe=${HOME}/src/mjo/wk99_smoothing
compile_script=${HOME}/src/mjo/wk99_smoothing.ksh
compile_script_final=${compile_script%%.ksh}.${run_id}.ksh
sed1=${src%%.f90}.${run_id}.f90.1
sed2=${src%%.f90}.${run_id}.f90.2
src_final=${src%%.f90}.${run_id}.f90
echo 'sed #1'
src_forsed=`echo ${src} | sed s/\\\//\\\\\\\\\\\//g`
echo 'sed #2'
src_final_forsed=`echo ${src_final} | sed s/\\\//\\\\\\\\\\\//g`
sym_power_norm_outfile=${anomaly_files[0]%%.${year_range}.olr.2.5x2.5.nc}_wk99_power_norm_sym.${year_range}.olr.2.5x2.5.nc
sym_power_norm_symlink=${input_directory}/wk99_power_sym_norm.nc
asym_power_norm_outfile=${anomaly_files[0]%%.${year_range}.olr.2.5x2.5.nc}_wk99_power_norm_asym.${year_range}.olr.2.5x2.5.nc
asym_power_norm_symlink=${input_directory}/wk99_power_asym_norm.nc
bkgrnd_power_norm_outfile=${anomaly_files[0]%%.${year_range}.olr.2.5x2.5.nc}_wk99_power_norm_bkgrnd.${year_range}.olr.2.5x2.5.nc
bkgrnd_power_norm_symlink=${input_directory}/wk99_power_bkgrnd_norm.nc
ln -sf ${sym_power_norm_outfile} ${sym_power_norm_symlink}
ln -sf ${asym_power_norm_outfile} ${asym_power_norm_symlink}
ln -sf ${bkgrnd_power_norm_outfile} ${bkgrnd_power_norm_symlink}
sed s/directory=directory/directory=\'${input_directory_forsed}\'/1 ${src} > ${sed1}
sed s/input_file_sym/\'wk99_power_sym.nc\'/1 ${sed1} > ${sed2}
sed s/input_file_asym/\'wk99_power_asym.nc\'/1 ${sed2} > ${sed1}
sed s/output_file_sym/\'wk99_power_sym_norm.nc\'/1 ${sed1} > ${sed2}
sed s/output_file_asym/\'wk99_power_asym_norm.nc\'/1 ${sed2} > ${sed1}
sed s/output_file_bkgrnd/\'wk99_power_bkgrnd_norm.nc\'/1 ${sed1} > ${sed2}
sed s/variable_name_sym/\'power\'/1 ${sed2} > ${sed1}
sed s/variable_name_asym/\'power\'/1 ${sed1} > ${sed2}
sed s/num_x/144/1 ${sed2} > ${sed1}
sed s/num_s/96/1 ${sed1} > ${src_final}
sed s/${src_forsed}/${src_final_forsed}/1 ${compile_script} > ${compile_script_final}    
chmod u+x ${compile_script_final}
${compile_script_final}
${exe}

#!/bin/ksh

f90=f90n
npt=${HOME}/src/nicks_subroutines.f90

${f90} -o /home/ss901165/src/mjo/create_segments ${npt} /home/ss901165/src/mjo/create_segments.xjhva.f90 -lfftpack

#!/bin/ksh

f90=f90n
npt=${HOME}/src/nicks_subroutines.f90

${f90} -o /home/ss901165/src/mjo/space_time_powerspec ${npt} /home/ss901165/src/mjo/space_time_powerspec.xjhva.f90 -lfftpack
#!/bin/ksh

f90=f90n
npt=${HOME}/src/nicks_subroutines.f90

${f90} -o /home/ss901165/src/mjo/wk99_power ${npt} /home/ss901165/src/mjo/wk99_power.xjhva.f90 -lfftpack

#!/bin/ksh

f90=f90n
npt=${HOME}/src/nicks_subroutines.f90

${f90} -o /home/ss901165/src/mjo/wk99_smoothing ${npt} /home/ss901165/src/mjo/wk99_smoothing.xjhva.f90 -lfftpack
