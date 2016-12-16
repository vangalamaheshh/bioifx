#!/usr/bin/env ruby

require 'bio'
require 'getoptlong'

def parse_options
	opts = GetoptLong.new( [ '--fasta_file', '-f', GetoptLong::REQUIRED_ARGUMENT ],
				[ '--id_file', '-d', GetoptLong::OPTIONAL_ARGUMENT ],
				[ '--id', '-i', GetoptLong::OPTIONAL_ARGUMENT ],
				[ '--reg_exp_to_fetch_id', '-r', GetoptLong::OPTIONAL_ARGUMENT ],
				[ '--help', '-h', GetoptLong::NO_ARGUMENT ] );
	fasta_file, id_file, usergiven_regex, id_list = nil, nil, nil, Array.new

	opts.each do | opt, arg |
		case opt
			when "--help"
				puts <<-EOF
	Usage: __FILE__ --fasta_file|-f <path to fasta file> --id_file|-i <file contiaining ids> 
			--reg_exp_to_fetch_id|-r <regexp to get the gene id> --help|-h <to print this help message>

		Options:

			--fasta_file|-f:
				File path to a multi or single fasta formatted file.

			--id_file|-d:
				File containing ids for whome the fasta seqeunces to be
				fetched from multifasta file given
			
			--id|-i:
				Id for whom the fasta sequence to be fetched. Multiple id's can be given in the
				format; -i id1 -i id2 ...
			
			--reg_exp_to_fetch_id|-r:
				Optional argument. Regexp to filter the gene id from the fasta description line.
				If not given, the < default regex: /(\w+)/ > is used.
						  --------------------------

			--help|-h:
				To see this help message.

		AUTHOR:
			Mahesh Vangala
			vangalamaheshh@gmail.com
			mvangala@som.umaryland.edu
				
				EOF
				exit 1

			when "--fasta_file"
				fasta_file = arg
			
			when "--id_file"
				id_file = arg

			when "--id"
				id_list << arg

			when "--reg_exp_to_fetch_id"
				usergiven_regex = arg
		end
	end

	unless fasta_file and ( id_file or id_list )
		puts "try --help"
		exit 1
	end
	
	unless usergiven_regex
		usergiven_regex = '(\w+)'
	end

	return fasta_file, id_file, usergiven_regex, id_list
end

def get_multi_fasta_objects( all_seqs_as_string )
	multi_obj = Bio::Probcons::DEFAULT_PARSER.new( all_seqs_as_string )
	array_of_fasta_obj = multi_obj.entries( )
	return array_of_fasta_obj
end

def concat_file_as_string( fasta_file )
	file_obj = File.new( fasta_file, "r" );
	file_as_string = String.new
	while line = file_obj.gets
		file_as_string.concat( line );
	end
	file_obj.close
	return file_as_string
end

def get_ids_in_hash( id_file, id_list )
	ids_info = Hash.new
	unless id_file
		id_list.each do | id |
			ids_info[ id ] = false
		end
		return ids_info
	end
	file_obj = File.new( id_file, "r" )
	while id = file_obj.gets
		id.chomp!
		unless ids_info.has_key?( id )
			ids_info[ id ] = false
		end
	end
	file_obj.close
	return ids_info
end

def print_matched_fasta_seqs( array_of_fasta_objects, ids_info, usergiven_regex )
	array_of_fasta_objects.each do | fasta_obj |
		definition = fasta_obj.definition( )
		match_data = /#{ usergiven_regex }/.match( definition )
		acc_id = match_data[ 1 ]
		if ids_info.has_key?( acc_id )
			$stdout.puts fasta_obj.entry( )
			ids_info[ acc_id ] = true
		end
	end	
	return ids_info
end

def print_unmatched_ids( ids_info )
	ids_info.each do |key,value|
		unless value
			$stderr.puts key
		end
	end
end

fasta_file, id_file, usergiven_regex, id_list = parse_options( )
array_of_fasta_objects = get_multi_fasta_objects( concat_file_as_string( fasta_file ) )
ids_info = get_ids_in_hash( id_file, id_list )
ids_info = print_matched_fasta_seqs( array_of_fasta_objects, ids_info, usergiven_regex )
print_unmatched_ids( ids_info )

exit 0
