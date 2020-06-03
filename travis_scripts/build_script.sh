#!/usr/bin/env bash

# fwdpp is included as a submodule,
# so we need to initalize it here:
git submodule update --init --recursive

if [ "$USECONDA" == "1" ];
then
    export PATH="$HOME/miniconda/bin:$PATH"
    # if [ "$TRAVIS_OS_NAME" == "linux" ]; then export LD_LIBRARY_PATH=$HOME/miniconda/lib; fi;
    # if [ "$TRAVIS_OS_NAME" == "osx" ]; then export DYLD_FALLBACK_LIBRARY_PATH=$HOME/miniconda/lib; fi;
    # export CPPFLAGS="-I$HOME/miniconda/include $CPPFLAGS"
    # export LDFLAGS="-L$HOME/miniconda/lib $LDFLAGS"
    # if [ "$TRAVIS_OS_NAME" == "linux" ]; then CC=$CC CXX=$CXX python setup.py build_ext -i; fi
    # if [ "$TRAVIS_OS_NAME" == "linux" ];then CC=$CC CXX=$CXX python -m unittest discover -v tests; fi
    # # if [ "$TRAVIS_OS_NAME" == "linux" ];then cd doc && make doctest; fi
    # if [ "$TRAVIS_OS_NAME" == "osx" -a "$OSXGCC" == "0" ]; then python setup.py build_ext -i; fi
    # if [ "$TRAVIS_OS_NAME" == "osx" -a "$OSXGCC" == "1" ]; then CC=gcc CXX=g++ python setup.py build_ext -i --gcc; fi
    # if [ "$TRAVIS_OS_NAME" == "osx" -a "$OSXGCC" == "0" ]; then LDFLAGS='-stdlib=libc++ -mmacosx-version-min=10.7' CPPFLAGS="-stdlib=libc++ -mmacosx-version-min=10.7" python -m unittest discover -v tests; fi
    # if [ "$TRAVIS_OS_NAME" == "osx" -a "$OSXGCC" == "1" ]; then CC=gcc CXX=g++ python -m unittest discover -v tests; fi
    if [ "$CXXSTANDARD" == "17" ];
    then
        CC=$CC CXX=$CXX python setup.py build_ext -i --cpp17
    else
        CC=$CC CXX=$CXX python setup.py build_ext -i 
    fi
    if [ "$?" != "0" ];
    then
        exit 1
    fi
    CC=$CC CXX=$CXX python -W ignore -m unittest discover -v tests
    if [ "$?" != "0" ];
    then
        exit 1
    fi
else
    echo "compilers are $CC $CXX"
    echo "python version is $TRAVIS_PYTHON_VERSION"
    if [ "$CXXSTANDARD" == "17" ];
    then
        CC=$CC CXX=$CXX python setup.py build_ext -i --cpp17
    else
        CC=$CC CXX=$CXX python setup.py build_ext -i 
    fi
       export LD_LIBRARY_PATH=/usr/local/lib
    if [ "$?" != "0" ];
    then
        exit 1
    fi
    CC=$CC CXX=$CXX python -W ignore -m unittest discover -v tests
    if [ "$?" != "0" ];
    then
        exit 1
    fi
    # cd doc && make doctest
    # if [ "$?" != "0" ];
    # then
    #     exit 1
    # fi
fi
