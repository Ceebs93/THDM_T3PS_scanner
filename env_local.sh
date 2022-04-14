
# Load Python environment
#conda deactivate
module unload --all
module load conda
#conda activate THDM
source activate THDM

# Load libraries and compiler
module load gcc/11.1.0
module load gsl/2.6

# Load THDM and T3PS
#echo "Path is: ${PATH}"
echo "Scanner Dir: ${THDM_T3PS_SCANNER_DIR}"
~                                             
