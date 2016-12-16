#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import string
import argparse
import subprocess

# -----------------------------------
# custom package
# -----------------------------------
import RNA_pipe

### global variable
from RNA_pipe.GlobalVar    import * 

### tool function
from RNA_pipe.Utility      import (sp,
                                   pdf_name,
                                   raise_error,
                                   wlog,
                                   readAnnotation,
                                   textformat
                                   )
### read and generate config file
from RNA_pipe.parse_config import (gen_conf,
                                   read_conf,
                                   make_conf
                                   )     


                                   
# -------------------
# main step
# -------------------
from RNA_pipe.step0_integrate_data   import step0_integrate_data
from RNA_pipe.step1_make_bam_star    import step1_make_bam_star
from RNA_pipe.step1_make_bam_tophat    import step1_make_bam_tophat
from RNA_pipe.step3_single_bam_QC    import step3_single_bam_QC
#from RNA_pipe.step4_summary_SingleQC import step4_summary_SingleQC
from RNA_pipe.step5_rpkm_table       import step5_rpkm_table
from RNA_pipe.step6_diffexp          import step6_diffexp
from RNA_pipe.stepFinal_summarize    import stepFinal_summarize
from RNA_pipe.step7_diff_AS          import step7_diff_AS
#from RNA_pipe.step01_filter_rRNA      import step01_remove_rRNA
#from RNA_pipe.step8_mutation         import step8_mutation

# ------------------------------------
# Misc functions
# ------------------------------------

    
#### read options

class ChiLinParser(argparse.ArgumentParser):
    '''
    Ex version for argparse , add raise error function .
    '''
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit()

def parse_args():
    '''
    Read parameter 
    '''
    description = ">.< RNAseq pipeline "
    parser = ChiLinParser(description = description)
    sub_parsers = parser.add_subparsers(help = "sub-command help", dest = "sub_command")

    ### generate config file
    template_parser = sub_parsers.add_parser("gen",  help = "generate a template of config file",
                                             description = "RNAseq config file generation")
    template_parser.add_argument("-s","--species", choices = ("hg19", "mm9"), required = True)
    template_parser.add_argument("-n","--name", dest="config_name",required = True)

    ### run config file
    pipe_parser = sub_parsers.add_parser("run", help = "run pipeline using a config file",
                                         description = "Run RNAseq pipeline using a config file")
    pipe_parser.add_argument("-c","--config", required = True,
                             help = "specify the config file to use", )
    pipe_parser.add_argument("-f","--force_overwrite",dest='fover',  default=False, action='store_true', 
                             help = "specify the config file to create output folder , this cmd will rm existing result if set True ~!! ", )
    pipe_parser.add_argument("-m","--mapping",dest='map',  default='star', 
                             help = "select mapping software: STAR or Tophat; default is STAR", ) 
    pipe_parser.add_argument("-as","--differential_alternative_splicing",dest='das',  default=False, action='store_true', 
                             help = "Using MATS call differential alternative splicing " )
    #  pipe_parser.add_argument("-mu","--mutaion",dest='mutation',  default=False, action='store_true', 
   #                          help = "Call mutation" )    
#    pipe_parser.add_argument("-r","--remove_rRNA",dest='ffpe',  default=False, action='store_true', 
#                             help = "remove rRNA" )                                                                                  
    ### simple mode
    simple_parser = sub_parsers.add_parser("simple", help = "run pipeline using simple mode",
                                         description = "(Run RNAseq pipeline using simple mode) Usage: RNA_QC.py -t a.fastq,b.fastq -c c_r1.fastq;c_r2.fastq,d.bam --treatlab t1,t2 --ctrllab c1,c3 -n yourname -s hg19 -f")
    simple_parser.add_argument("-t","--treat", dest = 'treat',required = True,
                             help = "treat data, comma(,) separate replicates , semicolon(;) separate pairend_fastq , format(fastq,bam) fixed by extension" )
    simple_parser.add_argument("-c","--ctrl",dest='ctrl',default = ' ',
                             help = "ctrl data, comma(,) separate replicates , semicolon(;) separate pairend_fastq , format(fastq,bam) fixed by extension" )
    simple_parser.add_argument("--treatlab", dest = 'treatlab',default = ' ',
                             help = "(no space inside)label of corresponded treat data, will use default if left blank or have different number with treat data, default : name of input file, separate by comma",)
    simple_parser.add_argument("--ctrllab", dest = 'ctrllab',default = ' ',
                             help = "(no space inside)label of corresponded ctrl data, will use default if left blank or have different number with ctrl data, default : name of input file, separate by comma",)
    simple_parser.add_argument("-n","--name", dest="name",required = True,help="name of you config file and output dir")
    simple_parser.add_argument("-s","--species",  choices = ("hg19", "mm9"), required = True,
                             help = "species ,choose from hg19 and mm9" )
    simple_parser.add_argument("-f","--force_overwrite",dest='fover',  default=False, action='store_true', 
                             help = "specify the config file to create output folder , this cmd will rm existing result if set True ~!! " )
    simple_parser.add_argument("-m","--mapping",dest='map',  default='star', 
                             help = "specify mapping software: STAR or Tophat; default is STAR", )    
    simple_parser.add_argument("-as","--alternative_splicing",dest='das',  default=False, action='store_true', 
                             help = "Using MATS call differential alternative splicing " )
    #simple_parser.add_argument("-mu","--mutaion",dest='mutation',  default=False, action='store_true', 
     #                        help = "Call mutation" ) 
#    simple_parser.add_argument("-r","--remove_rRNA",dest='ffpe',  default=False, action='store_true', 
#                             help = "remove rRNA" )                                                             
    args = parser.parse_args()

    if args.sub_command == "gen":
        gen_conf(args.species, args.config_name)
        sys.exit(0)

    if args.sub_command == "run":
        if os.path.isfile(args.config):
            return args
        else:
            print ('ERROR : -c input is not a config file\n')
            print pipe_parser.print_help()
            sys.exit()
    
    if args.sub_command == "simple":
        make_conf(args.treat,args.ctrl,args.treatlab,args.ctrllab,args.species,args.name,args.fover)
        args.config = args.name + '.conf'
        return args

# ------------------------------------
# Main function
# ------------------------------------

def main():
    args = parse_args()
    print args,args.map
    conf_dict = read_conf(args.config)
    ### read raw path of output dir
    conf_dict['General']['startdir'] = os.getcwd()+'/'
    if conf_dict['General']['outputdirectory'] == "":
        conf_dict['General']['outputdirectory'] = "./"
        current = 1
    else:
        current = 0
    conf_dict['General']['outputdirectory'] = sp('readlink -f %s'%(conf_dict['General']['outputdirectory']))[0].strip()+"/"
    ### creat output dir
    if current == 0:
		if args.fover :
			try:
			    if os.path.isdir(conf_dict['General']['outputdirectory']):
				    os.system('rm -r %s'%(conf_dict['General']['outputdirectory']))
			except:
				print 'can not create output dir , do you have write permission ?'
				sys.exit(1)
		if os.path.isdir(conf_dict['General']['outputdirectory']) or os.path.isfile(conf_dict['General']['outputdirectory']):
			print 'name of your output dir is already exist , pipeline exit to avoid overwrite file , you can add -f to force overwrite (see help)'
			sys.exit(1)
		else:
			os.system("mkdir %s"%(conf_dict['General']['outputdirectory']))
    else:
        print 'output dir is using default , output in current dir'
        
    os.chdir(conf_dict['General']['outputdirectory'])
    ## cp config file to output folder
    cmd = 'cp %s %s'%(conf_dict['General']['startdir']+args.config,conf_dict['General']['outputdirectory'])
    os.system(cmd)
    ### specify the main progress log file
    logfile = conf_dict['General']['outputdirectory'] + 'progress_log.txt'
    
    step0_integrate_data(conf_dict,logfile)
#    if args.ffpe:
#    	step01_remove_rRNA(conf_dict,logfile)
    print conf_dict
    if args.map.lower()=='star':
        step1_make_bam_star(conf_dict,logfile)
    else:   	
        step1_make_bam_tophat(conf_dict,logfile)
    step3_single_bam_QC(conf_dict,logfile)
    step5_rpkm_table(conf_dict,logfile)
    step6_diffexp(conf_dict,logfile)
    if args.das:
   		step7_diff_AS(conf_dict,logfile)
    stepFinal_summarize(conf_dict,logfile)
    cmd = 'chmod -R 777 %s'%conf_dict['General']['outputdirectory']
    os.system(cmd)
  #  if args.mutation:
  #  	step8_mutation(conf_dict,logfile)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(1)

