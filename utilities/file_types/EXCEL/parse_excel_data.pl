#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Spreadsheet::ParseExcel;
use Spreadsheet::WriteExcel::Simple;

my $options = parse_options( );
my $parser = Spreadsheet::ParseExcel -> new( );
my $workbook = $parser -> parse( $$options{ 'input_excel_workbook' } );

if( not defined $workbook ) {
	die $parser -> error( ), "\n";
}

my $compound_worksheet = $workbook -> worksheet( 'Compound' );
my $peak_worksheet = $workbook -> worksheet( 'Peak' );
my $aa_info = get_aa_info( $compound_worksheet );
my $aa_vals = get_rest_info( $peak_worksheet, sort { $a <=> $b } keys %$aa_info );
write_data( $aa_info, $aa_vals, $$options{ 'output_excel_sheet' } );
exit $?;

sub parse_options {
	my $options = {};
	GetOptions( $options, 'input_excel_workbook|i=s', 'output_excel_sheet|o=s', 'help|h' );
	unless( $$options{ 'input_excel_workbook' } and $$options{ 'output_excel_sheet' } ) {
		die "Usage: $0 <--input_excel_workbook|-i> <--output_excel_sheet|-o>\n";
	}
	return $options;
}

sub get_aa_info {
	my( $worksheet ) = @_;
	my $info = {};
	my( $row_min, $row_max ) = $worksheet -> row_range( );
	my( $col_min, $col_max ) = $worksheet -> col_range( );
	my $total_amount = 0;
	#row_min which is '0' is a header so we start reading the excel file from row '1'
	foreach my $row( $row_min + 1 .. $row_max ) {
		next unless( $worksheet -> get_cell( $row, 12 ) );
		my( $name, $amount, $cal_compound ) = (	$worksheet -> get_cell( $row, 12 ) -> value( ),
                                     		       	$worksheet -> get_cell( $row, 13 ) -> value( ),
                                                	$worksheet -> get_cell( $row, 14 ) -> value( ) );
		$$info{ $row }{ 'name' } = $name;
		$$info{ $row }{ 'amount' } = $amount;
		$$info{ $row }{ 'comp' } = $cal_compound;
		$total_amount += $amount;
	}
	foreach my $row( sort { $a <=> $b } keys %$info ) {
		$$info{ $row }{ 'amount_percent' } = ( ( $$info{ $row }{ 'amount' } / $total_amount ) * 100 );
	}
	return $info;
}

sub get_rest_info {
	my( $worksheet, @rows ) = @_;
	my $info = {};
	my( $row_min, $row_max ) = $worksheet -> row_range( );
	my( $col_min, $col_max ) = $worksheet -> col_range( );

	my( $total_height, $total_area ) = ( 0, 0 );

	#we want to consider only those rows for whom there was name in the previous worksheet 'Compound'
	foreach my $row( @rows ) {
		$$info{ $row }{ 'meas_ret_time' } = $worksheet -> get_cell( $row, 10 ) -> value( );
		$$info{ $row }{ 'corr_exp_ret_time' } = $worksheet -> get_cell( $row, 11 ) -> value( );
		$$info{ $row }{ 'int_peak_type' } = $worksheet -> get_cell( $row, 12 ) ? $worksheet -> get_cell( $row, 12 ) -> value( )
											: 'BLANK';
		$$info{ $row }{ 'area' } = $worksheet -> get_cell( $row, 13 ) -> value( );
		$$info{ $row }{ 'height' } = $worksheet -> get_cell( $row, 14 ) -> value( );
		$$info{ $row }{ 'width' } = $worksheet -> get_cell( $row, 15 ) -> value( );
		$$info{ $row }{ 'symmetry' } = $worksheet -> get_cell( $row, 16 ) -> value( );
		$$info{ $row }{ 'cal_peak' } = $worksheet -> get_cell( $row, 17 ) -> value( );
		$total_height += $$info{ $row }{ 'height' };
		$total_area += $$info{ $row }{ 'area' };
	}
	foreach my $row( sort { $a <=> $b } keys %$info ) {
		$$info{ $row }{ 'height_percent' } = ( ( $$info{ $row }{ 'height' } / $total_height ) * 100 );
		$$info{ $row }{ 'area_percent' } = ( ( $$info{ $row }{ 'area' } / $total_area ) * 100 );
	}
	return $info;
}

sub write_data {
	my( $aa_info, $aa_vals, $file_name ) = @_;
	my $simple_sheet = Spreadsheet::WriteExcel::Simple -> new( );
	my @headings = qw( Name Time Height Height_Percentage Area Area_Percentage Quantity Quantity_Percentage );
	$simple_sheet -> write_bold_row( \@headings );
	foreach my $row( sort { $a <=> $b } keys %$aa_info ) {
		my @data = ( $$aa_info{ $row }{ 'name' }, $$aa_vals{ $row }{ 'meas_ret_time' }, $$aa_vals{ $row }{ 'height' },
				$$aa_vals{ $row }{ 'height_percent' }, $$aa_vals{ $row }{ 'area' }, 
				$$aa_vals{ $row }{ 'area_percent' }, $$aa_info{ $row }{ 'amount' }, $$aa_info{ $row }{ 'amount_percent' } );
		$simple_sheet -> write_row( \@data );
	}
	$simple_sheet -> save( $file_name );
}
