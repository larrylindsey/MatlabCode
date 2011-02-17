#! /bin/bash

"$@"

if [ $? -ne 0 ]; then
  quoted=$(echo -e "$@" | awk 'BEGIN {
  FS = "\n" } { printf "\"%s\" ", $1 }' | sed -e s#\"\"##)
  eval "$EM_CODE_DIR/lib/python_email/send_email.sh em.reconstruct janelia.em.reconstruct \"$LOGNAME@janelia.hhmi.org\" \"error: \"$(pwd)\" $1\" $quoted"
fi

