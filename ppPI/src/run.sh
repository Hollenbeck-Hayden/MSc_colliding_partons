#!/bin/bash

if [ ! $1 ]; then
    echo "Missing entry: usage ./runrep.sh NAMEEXP IREP PREP"
    exit
else
    NAMEEXP=$1
    a=( $2 )
    IREP=${a[0]}
    PREP=${a[1]}
fi

echo "$IREP $PREP" >> parallel_log.txt

./hadunp.x<<EOF
$NAMEEXP
$IREP
$PREP
EOF

exit 0
