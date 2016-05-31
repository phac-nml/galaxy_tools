#!/bin/bash
if [[ ! -d $INSTALL_DIR ]]
then
    echo "ERROR: Environment variable INSTALL_DIR not a directory!"
fi
export INSTALL_DIR=${INSTALL_DIR:-$PWD}
export INSTALL_DIR=`(cd "$INSTALL_DIR"; pwd)`
#============================================================
echo "Setting environment variables for spades version 3.8.0"
#============================================================
specifc_action_done=0
#------------------------------------------------------------
if [[ $specifc_action_done == 0 ]]
then
    echo "Non-platform-specific actions"
    export PATH=$INSTALL_DIR/bin:$PATH
fi
