pacman -Sy
pacman -Su --nocomfirm
pacman -S base-devel cmake git eigen3 python3 suitesparse tbb superlu libxt --noconfirm

cat <<EOF >config_matlab_container.txt

CMAKE_FLAGS="\
  -G \"Unix Makefiles\"\
  -DCMAKE_CXX_FLAGS=\"-O3 -std=c++20 -march=native -fno-common -ftemplate-backtrace-limit=0\"\
  -DCMAKE_BUILD_TYPE=Release\
  -DCMAKE_DISABLE_FIND_PACKAGE_MPI=TRUE\
  -DCMAKE_CXX_COMPILER=g++\
  -DCMAKE_C_COMPILER=gcc\
  -DCMAKE_FIND_LIBRARY_SUFFIXES=\".a .so\"\
  -DDUNE_REENABLE_ADD_TEST=TRUE \
  -DDUNE_PYTHON_VIRTUALENV_SETUP=TRUE \
  -DDUNE_PYTHON_ALLOW_GET_PIP=TRUE \
  -DMatlab_ROOT_DIR=/home/matlab/local_install/local_matlab \
"
MAKE_FLAGS="-j2"
EOF
