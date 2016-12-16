#!/usr/bin/perl
# Takes the wrapper pipeline name from the client
# runs another tasklet that runs JSON pipeline summary metric
# returns the task id to the client
# @author: Mahesh Vangala
###########################################################

use strict;
use warnings;
use CGI;
use JSON::PP;

my $PIPELINE_NAME = 'pipeline_name';
my $WRAPPER = 'wrapper';
my $WORKER = 'worker';
my $WRAPPER_STATE = 'wrapper_state';
my $WORKER_STATE = 'worker_state';

print "Content-type: text/html \n\n";

my $q = new CGI;
my $params = $q -> Vars;
#my $params = { $PIPELINE_NAME => shift @ARGV , $WRAPPER => shift @ARGV, $WORKER => shift @ARGV };
my ($wrapper_id, $wrapper_cluster, $worker_name) = run_conf_command( $$params{$PIPELINE_NAME} );
my ($worker_id, $worker_cluster, $dont_want_it) = run_conf_command( $worker_name );

my $result = {};

($$result{$WRAPPER_STATE}, $$result{$WRAPPER}) = get_state_and_task_id( run_JSON_tasklet( $wrapper_id, "local" ) ) if( $$params{$WRAPPER} );
($$result{$WORKER_STATE}, $$result{$WORKER}) = get_state_and_task_id( run_JSON_tasklet( $worker_id, $worker_cluster ) ) if( $$params{$WORKER} );

print encode_json( $result );

exit $?;

sub get_state_and_task_id {
	my ($arg) = @_;
	my ($state_line) = $$arg[0];
	my ($result_line) = $$arg[1];
	my ($state, $task_id);
	if( $state_line =~ /state:\s+(\w+)/i ) { $state = $1; }
	if( $result_line =~ /.+'(.+)\\n'/ ) { $task_id = $1; }
	return ($state, $task_id);
}

sub run_conf_command {

	my ($pipeline_name) = @_;
	my $command = "vp-describe-task --show-all --block --no-print-polling `vp-run-metrics -t \"get-pipeline-conf $pipeline_name\"`";

	open(COMMAND, "$command |") or die "Error in executing the command, $command, $!\n";

	my ($id, $cluster, $name);

	while(my $line = <COMMAND>) {
		if( $line =~ /PIPELINE_ID=(.+?)\\n/ ) {
			$id = $1;
		}
		if( $line =~ /CLUSTER_NAME=(.+?)\\n/ ) {
			$cluster = $1;
		}
		if( $line =~ /input\.PIPELINE_NAME=(.+?)\\n/ ) {
			$name = $1;
		}
	}

	close COMMAND;
	die "Error in executing the command, $command, $!\n" if( $? );
	return ($id, $cluster, $name);

}

sub run_JSON_tasklet {
	my ($id, $cluster) = @_;
	my $array = [];
	push @$array, `vp-describe-task --show-all --block --no-print-polling \`vp-run-metrics -t --name $cluster \"run-JSON-pipeline-summary $id\"\``;
	my $filter = [];
	foreach( @$array ) {
		if( $_ =~ /^(Task)|(Result)/i ) { push @$filter, $_; }
	}
	return $filter;
}

