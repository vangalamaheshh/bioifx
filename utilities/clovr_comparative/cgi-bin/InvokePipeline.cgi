#!/usr/bin/perl
# Name: InvokePipeline.cgi
# @author: Mahesh Vangala
####################################

use strict;
use warnings;
use JSON::PP;
#use LWP::UserAgent;
#use HTTP::Request::Common;

my $refSeqFile = $ARGV[0];
my $mapFile = $ARGV[1];
my $pipeline = $ARGV[2];
my $optionalGenbankFile = $ARGV[3];
my $optionalMapFile = $ARGV[4];
my $pipelineName = $ARGV[5];
my $wrapper = "clovr_wrapper";
my $config_file = "/opt/clovr_pipelines/workflow/project_saved_templates/".$pipeline."/".$pipeline.".config";

sub getRefSeqInfo ($);
sub parseConfigFileIntoHash ($);

my $root;
my $refSeqInfo = getRefSeqInfo($refSeqFile);
$$root{'pipeline_config'} = parseConfigFileIntoHash($config_file);
$$root{'name'} = 'local';
$$root{'pipeline'} = $wrapper;
$$root{'pipeline_name'} = $pipelineName;
my $jsonText =  encode_json($root);
print $jsonText;
#print postURL();
exit(0);

sub postURL {
	my $ua = new LWP::UserAgent();
	#$ua->request(POST 'http://localhost/vappio/runPipeline_ws.py', $jsonText);
}

sub getRefSeqInfo ($) {
	my ($file) = @_;
	my $info = '';
	open(FH, "<$file") or die "Error in opening the file, $file, $!\n";
	while(my $line = <FH>) {
		$info .= $line;	
	}
	close FH;
	return $info;
}

sub parseConfigFileIntoHash ($) {
	my ($file) = @_;
	my $refHash = {};
	open(FH, "<$file") or die "Error in opening the file, $file, $!\n";
	my $line = <FH>;
	while($line) {
		chomp $line;
		if($line =~ /^\[(.+)\]/) {
			my $category = $1;
			$line = <FH>;
			while($line && $line !~ /^\[/) {
				chomp $line;
				if($line =~ /^(\S+?)=(.*)/) {
					if($1 eq 'INPUT_FILE') {
						$$refHash{$category."\.".$1}=$refSeqInfo;
					}
					elsif($1 eq 'ORGANISM_TO_DB_MAPPING') {
						$$refHash{$category."\.".$1}=$mapFile;
					}
					else {
						$$refHash{$category."\.".$1} = $2;
					}
				}
				$line = <FH>;
			}
		}
	}
	close FH;
	return $refHash;
}

