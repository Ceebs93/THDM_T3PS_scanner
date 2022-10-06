# .bash_profile

module load python/2.7.14
module load gcc/11.1.0
module load numpy
module load pandas

export FC=gfortran
export F77=gfortran
export CXX=g++
export CPP=cpp

# Get the aliases and functions
#if [ -f ~/.bashrc ]; then
#	. ~/.bashrc
#fi

# User specific environment and startup programs

PATH=$PATH:$HOME/bin:/local/bin:.
export PATH

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/local/software/gsl/2.6/include/gsl
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/local/software/gcc/11.1.0/bin/gcc


module load null
