#!/bin/bash

#
# //  Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
# //  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
# //  https://www.leibniz-hki.de/en/applied-systems-biology.html
# //  HKI-Center for Systems Biology of Infection
# //  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
# //  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
# //
# //  This code is licensed under BSD 2-Clause
# //  See the LICENSE file provided with this code for the full license.
#

# use the same as in depency script!!
export CC=/usr/local/bin/gcc-7.3
export CXX=/usr/local/bin/g++-7.3

USER=$(whoami)

# BUILD_TYPE="Release"
BUILD_TYPE="Debug"

if [[ "$1" =~ ^(release|Release)$ ]]; then BUILD_TYPE="Release"
elif [[ "$1" =~ ^(debug|Debug)$ ]]; then BUILD_TYPE="Debug"
else
    echo "No build type given: setting build type to RELEASE"
    BUILD_TYPE="Release"
fi

if [[ "$BUILD_TYPE" =~ Release ]]
then
    echo "compiling Release"
    cd cmake-build-release
    echo $PWD
    rm -v CMakeCache.txt
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make clean
    make -j 20 #(zahl der threads)
    # mv -v ParameterEstimatorForDynamicModels ${HOME}/bin/parameterFitting_${HOSTNAME}
    # ldd ${HOME}/bin/parameterFitting_${HOSTNAME}
fi

if [[ "$BUILD_TYPE" =~ Debug ]]
then
    echo "compiling Debug"
    cd cmake-build-debug
    echo $PWD
    rm -v CMakeCache.txt
    cmake -DCMAKE_BUILD_TYPE=Debug ..
    make clean
    make -j 20 #(zahl der threads)
fi
