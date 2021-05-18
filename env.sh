#!/usr/bin/env bash

COMMAND=$1

function set_up_apt_variables {
  echo "setting up variables for apt" &&
    export HDF5_USE_FILE_LOCKING="FALSE" &&
    export C_INCLUDE_PATH="/usr/local/include:/usr/include/hdf5/serial:/usr/include/suitesparse" &&
    export LIBRARY_PATH="/usr/lib/:/usr/local/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu/hdf5/serial"
}

function set_up_brew_variables {
  echo "setting up variables for Homebrew" &&
    export HDF5_USE_FILE_LOCKING="FALSE" &&
    export C_INCLUDE_PATH="/usr/local/include" &&
    export LIBRARY_PATH="/usr/local/lib"
}

function set_up_port_variables {
  echo "setting up variables for MacPort" &&
    export HDF5_USE_FILE_LOCKING="FALSE" &&
    export C_INCLUDE_PATH="/opt/local/include" &&
    export LIBRARY_PATH="/opt/local/lib"
}

if [[ $COMMAND == "brew" ]]; then
  set_up_brew_variables
elif [[ $COMMAND == "port" ]]; then
  set_up_port_variables
elif [[ $COMMAND == "apt" ]]; then
  set_up_apt_variables
else
  env
fi