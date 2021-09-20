#! /bin/bash

### script to look for root installations, pythonpath etc and set to proper values needed for the project. 

## set the path of Theia directory
## This variable is needed in the CMakeLists.txt
export THEIAPATH=$PWD
export STORAGEDIR=$PWD/Executables/source/Executables
echo "Setting the THEIAPATH to "$THEIAPATH
echo "Setting the STORAGEDIR to "$STORAGEDIR



### Now prepare the build and install directory

BUILD_DIR="../build/"
if [ -d "${BUILD_DIR}" ]; then
  ### Take action if $DIR exists ###
  echo "Building files in ${BUILD_DIR}..."
else
  ###  Control will jump here if $DIR does NOT exists ###
  echo "${BUILD_DIR} directory not found. Making the directory first!"
  mkdir -p ${BUILD_DIR}
fi


RUN_DIR="../run/"
if [ -d "${RUN_DIR}" ]; then
  ### Take action if $DIR exists ###
  echo "Building files in ${RUN_DIR}..."
else
  ###  Control will jump here if $DIR does NOT exists ###
  echo "${RUN_DIR} directory not found. Making the directory first!"
  mkdir -p ${RUN_DIR}
fi


INSTALL_DIR="../install/"
if [ -d "${INSTALL_DIR}" ]; then
  ### Take action if $DIR exists ###
  echo "Installing files in ${INSTALL_DIR}..."
else
  ###  Control will jump here if $DIR does NOT exists ###
  echo "${INSTALL_DIR} directory not found. Making the directory first!"
  mkdir -p ${INSTALL_DIR}
fi


BIN_DIR="../install/bin/"
if [ -d "${BIN_DIR}" ]; then
  ### Take action if $DIR exists ###
  echo "Installing files in ${BIN_DIR}..."
else
  ###  Control will jump here if $DIR does NOT exists ###
  echo "${BIN_DIR} directory not found. Making the directory first!"
  mkdir -p ${BIN_DIR}
fi

export INSTALL_PATH=$PWD/../install/bin

### cmake and install
cd ../build
cmake ../Theia
cmake --build . --target install 

### insert the path of libraries
export PATH=$PATH:${INSTALL_PATH}

### go back to run directory
#cd ../run
