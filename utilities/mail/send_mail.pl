#!/usr/bin/perl

use strict;
use warnings;
use lib '/zfs/cores/mbcf/mbcf-storage/devel/umv/software/perl/lib';
use Mail::Sendmail;

my @files = @ARGV;
unless( @files ) {
	exit 0;
} 

foreach my $file( @files ) {
	my $content = get_content( $file );
	unlink $file; #time to delete the file
	my %mail = ( From => 'vangalamaheshh@gmail.com',
		To => 'zachary_herbert@dfci.harvard.edu',
		Bcc => 'maura_berkeley@dfci.harvard.edu, Leslie_Grimmett@dfci.harvard.edu, umam_vangala@dfci.harvard.edu, zherbert3@gmail.com', 'Leon_Sheynkman@dfci.harvard.edu', 'Shiprav_Gupta@dfci.harvard.edu',
		Subject => $$content{ 'subject' },
		Message => $$content{ 'body' } );
	sendmail( %mail ) or die $Mail::Sendmail::error;
	print "Mail successfully sent\n\n", $Mail::Sendmail::log; 
}


sub get_content {
	my( $file ) = @_;
	my $content = {};
	open( FH, "<$file" ) or die "Error opening the file, $file, $!\n";
	my $subject = <FH>;
	chomp $subject;
	my @body = <FH>;
	$$content{ 'subject' } = $subject;
	$$content{ 'body' } = join( "", @body );
	close FH or die "Error closing the file, $file, $!\n";
	return $content;
}
