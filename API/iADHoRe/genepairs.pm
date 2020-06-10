package iADHoRe::genepairs;
#Created by Cedric Simillion on Sat Jun 26 12:03:21 CEST 2004

use strict;
use iADHoRe::settings;

my $settings=iADHoRe::settings->get_settings;

#===============================================================================
# Package methods
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub read_pairs
# Creates an access object to a hash containing all gene pairs.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub read_pairs
 {
  my %pairs;

  if ( !scalar( keys %pairs ) ) #Read the gene pairs if they haven't been read yet.
   {
    open (IN,$settings->blast_table) || die "Could not open BLAST table file ".settings->blast_table.": $!\n";
    while (<IN>)
     {
      chomp;
      my ($ga,$gb)= /^(\S*?)[+-]?\t(\S*?)[+-]?$/;
      ($ga eq $gb) && next;
      $pairs{$ga}{$gb}=-1;
      $pairs{$gb}{$ga}=-1;
     } #while (<IN>)
    close IN; 
   } #if ( scalar( keys %pairs ) )
   
  bless \%pairs, $_[0];
  return\%pairs;
   
 } #sub read_pairs


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub pairs_of
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub pairs_of
 {
  my ($pairs,$gene)=@_;
  exists( $$pairs{$gene} ) ? return( keys %{$$pairs{$gene}} ) : return;
 } #sub pairs_of

1;
