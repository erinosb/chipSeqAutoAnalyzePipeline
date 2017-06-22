#!/bin/bash
# specify commands to look for, and then try. The executable may be present, but crash at runtime (for missing shared libraries, for example).
#
#

# ANSI escape codes to be used with echo -e
RED='\033[0;31m'
BLACK='\033[0;30m'

# emoji
ANGER='\U1F4A2'
HOPEFUL='\U1F91E'
THINKING='\U1F914'
SMILEY='\U1F600'
CRYING='\U1F62D'
FEARFUL='\U1F628'


# subfunction to do tests and capture std/err output
check_execution() {
    retval=0 # 0 if no error; resolves to "true" in if-clause
    errorfilename=".error.log"

    # remove error file if exists
    if [ -e $errorfilename ]
    then
        rm $errorfilename
    fi

    # test for command in path
    if hash $1 >/dev/null
    then 
        echo -e "[$1] FOUND in PATH."
    else
        echo -e "[$1] NOT FOUND in PATH."
    fi

    # run command
    output=$($1 $2 2>$errorfilename)  

    # check return status
    if [ $? == 0 ]
    then
        echo -e "[$1] COMMAND RAN SUCCESSFULLY."
    else
        # include stdout failure message, if exists
        output_msg="."
        if [ ${#output} != 0 ]; then output_msg=": gave \'$output'"; fi
        echo -e "[$1] COMMAND FAILED$output_msg" 

        # include stderr failure message, if exists
        if [ -e $errorfilename ]
        then
            echo -e "with ERROR:${RED} $(cat $errorfilename)${BLACK}"
        fi
        retval=1 # will resolve to "false" in an if-clause
    fi

    # clean up
    if [ -e $errorfilename ]
    then
        rm $errorfilename
    fi

    return $retval
}
chk() {
    ret=1
    echo -e "${HOPEFUL} ================ $PROG ============================================================== "
    if check_execution $1 $2
    then
        echo -e "$PROG PASSED. ${SMILEY}"
    else
        echo -e "$PROG FAILED ${FEARFUL}"
        ret=0
    fi
    return $ret
}

echo "Your PATH is:"
echo -e "\t$PATH" | sed 's/:/\n\t/g'

echo "Your LD_LIBRARY_PATH is:"
echo -e "\t$LD_LIBRARY_PATH" | sed 's/:/\n\t/g'

# -version
for PROG in java;
do
    chk $PROG -version
done
# --version
for PROG in bowtie bowtie2 macs2 tagdust bedtools samtools;
do
    chk $PROG --version
done


