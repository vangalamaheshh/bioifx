#!/bin/bash

#-------------------------------------------------------------------
#  Dependencies: sffinfo should be in the path (Newbler suite)
#  Input: Takes a list of absolute path of sff files
#  Output: Generates both fasta and quality files for each sff file
#  @author: Mahesh Vangala
#           vangalamaheshh@gmail.com
#
#  Usage: ./sff2fasta.bash $(find $PWD -type f -name "*sff")
#--------------------------------------------------------------------

for file in "$@"; do
  echo -n "Processing $file ... "; >&2
  file_without_ext=$(echo "$file" | ruby -e ' print $stdin.gets.split( "." )[0..-2].join(".")')
  sffinfo -seq -notrim "$file" > "$file_without_ext.notrim.fa"
  sffinfo -qual -notrim "$file" > "$file_without_ext.notrim.qual"
  sffinfo -seq "$file" > "$file_without_ext.trim.fa"
  sffinfo -qual "$file" > "$file_without_ext.trim.qual"
  echo "Done" >&2
done

exit $?

#---------------------------------------------------------------------
