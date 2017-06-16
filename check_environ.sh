#!/bin/bash
# specify commands to look for, and then try. The executable may be present, but crash at runtime (for missing shared libraries, for example).
#
#

# subfunction to capture output
check_execution() {
    retval=0
    errorfilename=".error.log"
    if [ -e $errorfilename ]
    then
        rm $errorfilename
    fi

    # run command
    output=$($1 $2 2>$errorfilename)  

    # check return status
    if [ $? == 0 ]
    then
        echo "[$1] COMMAND RAN SUCCESSFULLY."
    else
        echo "COMMAND FAILED: gave \`$output'" 
        if [ -e $errorfilename ]
        then
            echo "AND GAVE THIS ERROR: $(cat $errorfilename)"
        fi
        retval=1
    fi

    # clean up
    if [ -e $errorfilename ]
    then
        rm $errorfilename
    fi

    return $retval
}

# bowtie1

if hash bowtie >/dev/null
then
    check_execution bowtie --version
else
    echo "bowtie NOT FOUND."
fi
