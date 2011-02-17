#!/bin/csh
python -c "import sys; sys.path = ['./'] + sys.path; import py_email; py_email.sendMail('$argv[1]','$argv[2]','$argv[3]','$argv[4]','$argv[5]');"


