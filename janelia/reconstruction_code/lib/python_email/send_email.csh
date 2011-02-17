#!/bin/csh
if(! $?EM_CODE_DIR) then
    set EM_CODE_DIR=/groups/chklovskii/chklovskiilab/em_reconstruction/code
endif

if($#argv < 6) then
    python -c "import sys; sys.path = ['$EM_CODE_DIR/lib/python_email'] + sys.path; import py_email; py_email.sendMail('$argv[1]','$argv[2]','$argv[3]','$argv[4]','$argv[5]');"
else
    set att="['$argv[6]'"
    if($#argv >6) then
        foreach a ($argv[7-])
            set att="$att ,'$a'"
        end
    endif
    set att="$att]"
    python -c "import sys; sys.path = ['$EM_CODE_DIR/lib/python_email'] + sys.path; import py_email; py_email.sendMailA('$argv[1]','$argv[2]','$argv[3]','$argv[4]','$argv[5]', $att);"
endif


