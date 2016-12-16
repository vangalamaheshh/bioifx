# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

#alias python='/usr/local/python-2.7.6/bin/python'

alias python='/apps/python-2.7.2/bin/python'
#alias_python='/apps/python-2.7.6/bin/python'
export PYTHONPATH=$PYTHONPATH:/apps/python-2.7.2/lib/python2.7/site-packages/:/ifs/rcgroups/Apps/RNAseq_pipeline/RNAseqENV/lib/python2.7/site-packages/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64
export PATH=/zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/rna_seq/cfce_support_scripts/:/ifs/rcgroups/Apps/RNAseq_pipeline/RNAseqENV/bin/:$PATH


#alias Rscript='/usr/local/R-2.15.1/bin/Rscript'
alias R='/apps/R-2.15.1/bin/R'
alias Rscript='/apps/R-2.15.1/bin/Rscript'


