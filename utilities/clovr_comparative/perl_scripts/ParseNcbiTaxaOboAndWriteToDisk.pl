#!/usr/bin/perl
# Name: ParseNcbiTaxaOboAndWriteToDisk.pl
# Desc: Gets the Ncbi taxonomy information into a tree data structure and writes to disk
#		Uses Storable module.
#		$$root{accession_id}{text} gives name of the current node visiting
#		$$root{accession_id}{children} gives a reference to another hash 
#			that has current node's children info
#		$$root{accession_id}{leaf} tells us whether the node being visited is a leaf or not
# NOTE: As we won't be using GO::Parser in the .cgi program that sends data to client, we
#		append the top root node id's into an array and save that along with $root hash into
#		a Final array and will store that array to disk using Storable
# Dependencies: GO::Parser, Storable
# @author: Mahesh Vangala
###################         LIBRARIES & PRAGMAS         #################

use strict;
use warnings;
use GO::Parser;
use Storable;

###################             CONSTATNS               #################

use constant node => 'node';
use constant id	=> 'id';
use constant children => 'children';
use constant text => 'text';
use constant leaf => 'leaf';
use constant leaf_count => 'leaf_count';
use constant true => 1;
use constant false => 0;

###################             GLOBAL VARIABLES        #################

my $root = {};
my $topNodes = [];
my $writeFile = 'NcbiOboTaxaDataStructure';

###################            SUBS DECLARATION         #################

sub populateNodeInfo ($);
sub getChildren ($);
sub getImmediateChildren ($);

###################            MAIN PROGRAM             #################

my $input_file = shift@ARGV || die "Usage: perl program.pl 'NCBI OBO TAXA FILE NAME'\n";

my $parser = new GO::Parser({handler=>'obj'});
$parser->parse($input_file);
my $graph = $parser->handler->graph;
foreach(@{$graph->find_roots}) {
	push @$topNodes, $_->acc;
	populateNodeInfo($_);
}

store([$root,$topNodes], $writeFile) or die "Error in writing to the file, $writeFile, $!\n";
print "Success\n";
exit(0);

##################             END OF MAIN               #################

## SUBS
# populates the node's name and it's children info
# also tells us whether it's a leaf or not by having LEAF boolean
# also populates the leaf count at every level
# recursive sub
sub populateNodeInfo ($) {
	my ($node) = @_;
	if(!$graph->is_leaf_node($node->acc)) {
		foreach(@{$graph->get_child_terms($node->acc)}) {
			$$root{$node->acc}{children}{$_->acc} = {};
			populateNodeInfo($_);
			$$root{$node->acc}{leaf_count} += $$root{$_->acc}{leaf_count};
		}
		$$root{$node->acc}{leaf} = false;
	}
	else { 
		$$root{$node->acc}{children} = {}; 
		$$root{$node->acc}{leaf} = true;
		$$root{$node->acc}{leaf_count} = 1;
	}
	$$root{$node->acc}{text} = $node->name;
}

# prints children by following a recursive call
# recursive sub
sub getChildren ($) {
	my ($node) = @_;
	while(my ($key,$value) = each %$node) {
		#print "ID: ",$key,"\t", "Name: ", $$root{$key}{text} || "Name not found","\n";
		getChildren($$root{$key}{children});
	}
}

# Prints Immediate children
# Not a recursive function
# just loops through all the keys in children hash and prints the name of each
sub getImmediateChildren ($) {
	my ($refHash) = @_;
	my $counter = 0;
	while(my ($key, $value) = each %$refHash) {
		#print "Name of child ", $counter++, " '$$root{$key}{text}' \t LEAF COUNT: '$$root{$key}{leaf_count}'\n";
	}	
}
