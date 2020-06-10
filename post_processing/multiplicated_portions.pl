#!/usr/bin/perl -w
#This script was written by Cedric Simillion on Tue Jun 21 14:23:00 CEST 2005

=head1 Usage

 multiplicated_portions.pl settings_file

=head1 Description

This program analyses the results from an i-ADHoRe run and prints a tab-delimited table to
STDOUT that gives for each genome the fraction covered by multiplicons of different
multiplication levels as well as the overall duplicated fraction. This table lists for
each gene list first the total number of positions followed by the number of positions in
each multiplication level as well as the fraction in percent.

=head1 Author

This script was written by Cedric Simillion.

=cut

use strict;
use lib "/home/sepro/i-adhore2.1/API";
use iADHoRe;

my ($data, $max_level);
my (%listsizes,%segment_level,%level_portion);

#Read the dataset
$data=iADHoRe->read_dataset;
$max_level=0;

#===============================================================================
# Get the sizes of each list and prepare a datastructure that stores the
# multiplication level of each position
#===============================================================================
foreach my $list ( $data->genelists )
 {
  my ($genome,$listname);
  $genome = $list->genome;
  $listname = $list->listname;
  $listsizes{$genome}{$listname} = $list->unmapped_size;
  $listsizes{$genome}{"total"} += $listsizes{$genome}{$listname};
  foreach my $i ( 0 .. $listsizes{$genome}{$listname} - 1 )
   {
    $segment_level{$genome}{$listname}[$i] = 1;
   } #foreach my $i ( 0 .. $listsizes{$genome}{$listname} - 1 )
 } #foreach my $list ( $data->genelists )

#===============================================================================
# Compute the multiplication level of each position
#===============================================================================
foreach my $multiplicon ( $data->multiplicons )
 {
  my %level_per_genome;
  
  $multiplicon->is_redundant && next;
  
  foreach my $segment ( @{$multiplicon->profile->segments} )
   {
    $level_per_genome{ $segment->genome }++;
   } #foreach my $segment ( @{$multiplicon->profile->segments} )
  
  foreach my $segment ( @{$multiplicon->profile->segments} )
   {
    my ($first_gene,$last_gene,$genome,$listname,$level);
    $first_gene = $segment->first->gene;
    $last_gene = $segment->last->gene;
    $genome = $segment->genome;
    $listname = $segment->listname;
    $level = $level_per_genome{$genome};
    foreach my $i ( $first_gene->coordinate .. $last_gene->coordinate )
     {
      ( $level > $segment_level{$genome}{$listname}[$i] ) && ( $segment_level{$genome}{$listname}[$i] = $level );
     } #foreach my $i ( $first_gene->coordinate .. $last_gene->coordinate )
   } #foreach my $segment ( @{$multiplicon->profile->segments} )
 } #foreach my $multiplicon ( $data->multiplicons )

#===============================================================================
# Compute the number of positions of each multiplication level
#===============================================================================
foreach my $genome (keys %segment_level)
 {
  foreach my $listname ( keys %{$segment_level{$genome}} )
   {
    foreach my $level ( @{$segment_level{$genome}{$listname}} )
     {
      $level_portion{$genome}{$listname}{$level}++;
      ( $level > $max_level ) && ( $max_level = $level );
     } #foreach my $level ( @{$segment_level{$genome}{$listname}} )
   } #foreach my $listname ( keys %{$segment_level{$genome}} )
 } #foreach my $genome (keys %segment_level)

#===============================================================================
# Output the results
#===============================================================================
foreach my $genome ( sort(keys %level_portion) )
 {
  my $portion_duplicated;
  print STDOUT "$genome\n";
  print STDOUT "list\tlist size\t2";
  foreach ( 3 .. $max_level ) { print STDOUT "\t\t$_" }
  print STDOUT "\n";
  foreach my $listname ( sort( keys %{$level_portion{$genome}} ) )
   {
    my $size = $listsizes{$genome}{$listname};
    print STDOUT "$listname\t$size";
    foreach my $level ( 2 .. $max_level )
     {
      my $portion = defined( $level_portion{$genome}{$listname}{$level} ) ? $level_portion{$genome}{$listname}{$level} : 0;
      print STDOUT "\t$portion\t".($portion/$size*100)."%";
      $portion_duplicated+=$portion;
     } #foreach my $level ( 2 .. $max_level )
    print STDOUT "\n";
   } #foreach my $listname ( keys %{$level_portion{$listname}} )
  print STDOUT ($portion_duplicated/$listsizes{$genome}{"total"}*100) ."% of the genome of $genome is duplicated.\n\n";
 } #foreach my $genome (keys %level_portion)
