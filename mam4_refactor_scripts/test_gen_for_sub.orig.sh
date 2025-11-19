

main() {
    #-------------------------
    #User input Starts
    #-------------------------

    #subroutine name
    sub_name=seasalt_emisflx_calc

    #ins and in-outs variables
    #(a comma seprated string, like 'var1, var2, var3'- remember the quotes around the list of variables)
    ins='ncol, fi, ocnfrc, emis_scale, flx_type'

    #outs and in-outs variables
    #(a comma seprated string, same as above)
    outs='cflx'

    #module name
    module_name=seasalt_model

    #file path
    dir_path=components/eam/src/chemistry/modal_aero/

    #CPP directive to turn on file writing
    cpp_directive=YAML_DUST_SSLT_EM

    #-------------------------
    #USER INPUT ENDS
    #-------------------------

    #----------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------

    #ensure that script is run in the right directory (i.e. "e3sm_mam4_refactor") as
    #relative path code works only if script is run at the code root
    this_dir=${PWD##*/}

    if [ $this_dir != "e3sm_mam4_refactor" ]; then
        echo 'Please run this script at source code root in directory "e3sm_mam4_refactor"'
        exit 1
    fi

    #Yaml directory path
    yaml_dir_path=components/eam/src/chemistry/yaml

    #path to directory containing beg and end input stubs
    stub_dir=mam4_refactor_scripts

    #relative path
    relative_path=`realpath --relative-to="$dir_path" "$yaml_dir_path"`

    #create directory and files if not already created
    /bin/mkdir -p $yaml_dir_path/$module_name/f90_yaml

    #generate code in tmp files for the calls to add in the "beg" and "end" files
    ftmp_ins='tmp_in_calls.txt'
    ftmp_outs='tmp_out_calls.txt'
    create_code_for_in_outs "$ins" "$outs" $ftmp_ins $ftmp_outs

    #code to add
    echo '----------------------------'
    echo 'Code to add to the F90 file:'
    echo '----------------------------'

    echo '(If not already added, add the following line at module level at the top)'
    echo "#include \"$relative_path/common_files/common_uses.ymlf90\""


    #Subroutine beginning file
    newline
    create_file $module_name $sub_name $yaml_dir_path $relative_path $ftmp_ins "beginning"

    #subroutine end file
    newline
    create_file $module_name $sub_name $yaml_dir_path $relative_path $ftmp_outs "end"

    #remove temporary files
    rm $ftmp_ins $ftmp_outs
}

#---------------------
# Function Definitions
#---------------------

#Generate a newline
newline () { echo ''; }

create_code_for_in_outs () {
    newline
    newline
    echo 'Generating files with code to add to the include files ......'
    newline
    #call python script to generate code for the calls
    MY_PATH="$(dirname -- "${BASH_SOURCE[0]}")" # relative path
    MY_PATH="$(cd -- "$MY_PATH" && pwd)"        # absolute and normalized path
    /bin/python3 $MY_PATH/gen_input_write_vars.py  --finp $3 --fout $4 -i "$1" -o "$2"

    #check if the files exists that has intent-in, intent-inouts and intent-outs calls
    if [ ! -f $3 ]; then
        echo "File tmp_in_calls.txt doesn't exist, check gen_input_write_vars.py script, exiting"
        exit 1
    fi
    if [ ! -f $4 ]; then
        echo "File tmp_out_calls.txt doesn't exist, check gen_input_write_vars.py script, exiting"
        exit 1
    fi

}


create_file () {

    #create_file $module_name $sub_name $yaml_dir_path $relative_path $ftmp_ins/$ftmp_outs "beginning"/"end"
    echo "(Add the following line at the $6 of the subroutine)"

    #find beginning or end sub string
    if [ $6 == "beginning" ]; then
        sub_str=beg
    elif [ $6 == "end" ]; then
        sub_str=end
    else
        echo 'Invalid string (beginning or end)for function create_file'
        exit 1
    fi

    #full file name
    f_name=${sub_name}_${sub_str}_yml.f90

    #code to include
    echo "#include \"$4/$1/f90_yaml/$f_name\""

    #create this file if it doesn't exist
    f_path=$3/$1/f90_yaml/$f_name
    if test -f "$f_path"; then
        newline
        echo "[INFO ONLY]: $f_path already exists."
    else
        newline
        echo '[INFO ONLY]:Creating NEW file: '$f_path
        f_tmp=tmp.inp
        sed s/"!#ifdef YAML_CPP"/"#ifdef $cpp_directive"/g $stub_dir/${sub_str}.inp > $f_tmp
        sed -i s/"SUB_NAME"/"\'$2\'"/g $f_tmp
        sed -i -e "/WRITE_CALLS/r $5" -e "s///" $f_tmp
        cat $f_tmp>>$f_path
        /bin/rm $f_tmp
    fi
}

#call the main function
main
