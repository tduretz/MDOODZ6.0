#!/bin/bash

COMMAND=$1
OPTIONS=$2

function usage {
  cat <<EOF
CLI interface for automatic deployment of the MDOODZ6.0 project

Usage:
  ops.sh [install]

Examples:
  ./ops.sh install brew                         # install and sets up packages using Homebrew
  ./ops.sh install port                         # install and sets up packages using MacPort
  ./ops.sh install apt                          # install and sets up packages using apt (native Ubuntu)
EOF
}

function install_packages_with_apt {
    echo "installing packages using apt (Ubuntu)" &&
    sudo apt -y install git || echo "git is already installed" &&
    sudo apt -y install libsuitesparse-dev || echo "suitesparse is already installed" &&
    sudo apt -y install hdf5-tools || echo "hdf5 is already installed" &&
    sudo apt -y install gcc-7 || echo "gcc7 is already installed" &&
    ./env.sh apt
}

function install_packages_with_brew {
  echo "installing packages using Homebrew" &&
    brew install git || echo "git is already installed" &&
    brew install suite-sparse || echo "suitesparse is already installed" &&
    brew install hdf5 || echo "hdf5 is already installed" &&
    brew install gcc@7 || echo "gcc7 is already installed" &&
    ./env.sh brew
}

function install_packages_with_macport {
  echo "installing packages using MacPort" &&
    sudo port install git || echo "git is already installed" &&
    sudo port install gccX || echo "gccX is already installed" &&
    sudo port install suitesparse || echo "suitesparse is already installed" &&
    sudo port install hdf5 || echo "hdf5 is already installed" &&
    ./env.sh port
}

if [[ $COMMAND == "install" ]]; then
  if [[ $OPTIONS == "brew" ]]; then
    install_packages_with_brew
  elif [[ $OPTIONS == "port" ]]; then
    install_packages_with_macport
  elif [[ $OPTIONS == "apt" ]]; then
    install_packages_with_apt
  else
    usage
  fi
else
  usage
fi
