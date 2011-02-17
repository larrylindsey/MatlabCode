#!/bin/bash

if [ $# -le 6 ]; then
    python -c "import sys; sys.path = ["\'"$EM_CODE_DIR/lib/python_email"\'"] + sys.path; import py_email; py_email.sendMail("\'"$1"\',\'"$2"\',\'"$3"\',\'"$4"\',\'"$5"\'");"
fi


