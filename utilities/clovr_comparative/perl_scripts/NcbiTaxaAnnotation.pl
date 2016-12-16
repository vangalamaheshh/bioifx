#!/usr/bin/perl
# Name: NcbiTaxaAnnotation.pl
# Desc: Reads total taxa data structure and the annotation data structure
#		filters the taxa data structure and keeps only those taxa that have
#		annotation information.
#		Leaf counts are programmed in such a way that at any level, the leaf
#		count tells us the number of refseqs present at that level
#		deletes those keys in the taxa data structure that don't have annotation
#		info. The final data struc is written to disk.
#		$$root{$id}{$name}
#		$$root{$id}{$annot}[[filepath, seq length, seq id], [...], []]
#		$$root{$id}{$children}{id of child1 ..}
#		The children hash also filtered by deleting those child id's with out 
#		annotation information
#		If at a level, no child has annotation then the leaf flag is set to true for that
#		level as long as there are some refseqs available at that level.
# Dependencies: Storable
# @author: Mahesh Vangala
###########################       LIBRARIES & PRAGMAS        ########################

use warnings;
use strict;
use Storable;

###########################           CONSTANTS              ########################

my $ROOT = '/';
my $NCBI_TAXA = 'NCBITaxon:';
my $TRUE = 1;
my $FALSE = 0;
my $ANNOT = 'ANNOTATION_INFO';
my $MAP = 'MAP_INFO';
my $NAME = 'text';
my $LEAF = 'leaf';
my $CHILDREN = 'children';
my $LEAF_COUNT = 'leaf_count';

############################      SUBS DECLARATION           ######################## 

sub generateAnnotationInfo ();
sub addAnnotationFields ($$);
sub generateFakeRootToSynchronizeTopNodes();
sub filterOutUnAnnotatedOnes ();
sub helperFilter ($);
sub initializeLeafCountsToZero();
sub makeLeafCountZero($);
sub addMapInfo ($);
sub helperAddMapInfo ($);

############################        MAIN PROGRAM             #########################

my ($taxa) = retrieve('../binary_files/NcbiOboTaxaDataStructure')
							|| die "Error retrieving the NcbiOboTaxaDataStructure\n";
							
my ($annot_info) = retrieve('../binary_files/NcbiAnnotationDataStructure')
							|| die "Error retrieving the NcbiAnnotationDataStruc\n";
						
my ($root) = $$taxa[0];
generateFakeRootToSynchronizeTopNodes ();
initializeLeafCountsToZero();						
generateAnnotationInfo ();
filterOutUnAnnotatedOnes ();
addMapInfo($root);
store($root, '../binary_files/NcbiTreeWidgetDataStructure') || 
	die "Error in writing the data structure to the disk\n";
print "Success:\n";
exit(0);

############################         END OF MAIN              ##########################

### SUBS

sub initializeLeafCountsToZero() {
	while(my ($key,$value) = each %{$$root{$ROOT}{$CHILDREN}}) {
		makeLeafCountZero($key);
	}
}

sub makeLeafCountZero ($) {
	my ($cur_node) = @_;
	unless($$root{$cur_node}{$LEAF}) {
		while(my ($key,$value) = each %{$$root{$cur_node}{$CHILDREN}}) {
			makeLeafCountZero($key);
		}
	}
	$$root{$cur_node}{$LEAF_COUNT} = 0;
}

sub helperFilter ($) {
	my ($cur_node) = @_;
	unless($$root{$cur_node}{$LEAF}) {
		while(my ($key,$value) = each %{$$root{$cur_node}{$CHILDREN}}) {
			my $ref_seq_count = helperFilter($key);
			$ref_seq_count ? $$root{$cur_node}{$LEAF_COUNT} += $ref_seq_count
						   : delete($$root{$cur_node}{$CHILDREN}{$key});	
		}
		if($$root{$cur_node}{$LEAF_COUNT} && $$root{$cur_node}{$ANNOT} &&
				$$root{$cur_node}{$LEAF_COUNT} == scalar@{$$root{$cur_node}{$ANNOT}}) {
			$$root{$cur_node}{$LEAF} = $TRUE;
		}
		
		my @children = keys(%{$$root{$cur_node}{$CHILDREN}});
		if($$root{$cur_node}{$LEAF_COUNT} && !(@children)) {
		    $$root{$cur_node}{$LEAF} = $TRUE;
		}
	}
		
	unless($$root{$cur_node}{$LEAF_COUNT}) {
		delete($$root{$cur_node});
	}
	
	return $$root{$cur_node}{$LEAF_COUNT} || 0;
}

sub filterOutUnAnnotatedOnes () {
	while(my ($key,$value) = each %{$$root{$ROOT}{$CHILDREN}}) {
		$$root{$key}{$LEAF_COUNT} = helperFilter($key);
		unless($$root{$key}{$LEAF_COUNT}) {
			delete($$root{$ROOT}{$CHILDREN}{$key});
		}
	}	
}

sub generateFakeRootToSynchronizeTopNodes () {
	foreach(@{$$taxa[1]}) {
		$$root{$ROOT}{$CHILDREN}->{$_} = {};
	}
}

sub generateAnnotationInfo () {
	while(my ($key,$value) = each %$annot_info) {
		if($key =~ /:(.+)/g) {
			my ($ncbi_key) = $NCBI_TAXA.$1;
			if($$root{$ncbi_key}) {
				addAnnotationFields($ncbi_key, $key);
			} 
			else {
				$$root{$ncbi_key}{$ANNOT} = [];
			}
		}
	}	
}

sub addAnnotationFields ($$) {
	my ($ncbi_key, $annot_key) = @_;
	$$root{$ncbi_key}{$ANNOT} = $$annot_info{$annot_key};
	$$root{$ncbi_key}{$LEAF_COUNT} = scalar@{$$annot_info{$annot_key}};	
}

sub addMapInfo ($) {	
	my ($root) = @_;
	while(my ($key,$value) = each %{$$root{$ROOT}{$CHILDREN}}) {
		helperAddMapInfo($key);
	}	
}

sub helperAddMapInfo ($) {
	my ($id) = @_;
	$$root{$id}{$MAP} = [];
	if($$root{$id}{$ANNOT} && scalar @{$$root{$id}{$ANNOT}} > 0) {
		my $temp = {};
		foreach my $refSeq(@{$$root{$id}{$ANNOT}}) {
			open(FH,"<$$refSeq[0]") or die "Error in reading the file, $$refSeq[0], $!\n";
			while(my $line = <FH>) {
				$line =~ s/^\s+//;
				$line =~ s/\s+$//;
				if($line =~ /^ORGANISM\s+(.+)/) {
					my @array = split(" ", $1);
					my $string = join(" ", (@array,$array[scalar@array-1]));
					$$temp{$string}++;
					last;
				}	
			}
			close FH;
		}		
		push @{$$root{$id}{$MAP}}, keys%$temp;
	}
	while(my ($key,$value) = each %{$$root{$id}{$CHILDREN}}) {
		helperAddMapInfo($key);
	}
}

