#!/usr/bin/perl -w
#This script was written by Cedric Simillion on Thu Jun 23 17:36:26 CEST 2005


=head1 Usage

 $ colinear_portions.pl <inifile>
 
<inifile> is a i-ADHoRe settings file.

=head1 Description

In a multi-genome dataset, lists for each genome wich portion is colinear to any of the
other genomes.

=head1 Author

This script was written by Cedric Simillion.

=cut

use strict;
use lib "/home/sepro/i-adhore2.3/API";
use iADHoRe;

my ($dataset);
my (@genomes);
my (%colinear_positions,%listsize);

$dataset = iADHoRe->read_dataset;

@genomes = sort( $dataset->genomes );
( scalar(@genomes) == 1 ) && die "Error: Only one genome in the dataset!\n";

foreach my $genelist ( $dataset->genelists )
 {
  $listsize{ $genelist->genome }{ $genelist->listname }=$genelist->unmapped_size;
  $listsize{ $genelist->genome }{"total"}+=$genelist->unmapped_size
 } #foreach my $genelist ( $dataset->genelists )

foreach my $multiplicon ( $dataset->multiplicons )
 {
  my @genomes_in_multiplicon;
  my %segments;
  
  $multiplicon->is_redundant && next;

  foreach my $segment ( @{$multiplicon->profile->segments} )
   {
    push( @{$segments{ $segment->genome }}, $segment );
   } #foreach my $segment ( {@{$multiplicon->profile->segments} )
  @genomes_in_multiplicon = keys %segments;
  
  foreach my $ref_genome ( @genomes_in_multiplicon )
   {
    foreach my $other_genome ( @genomes_in_multiplicon )
     {
      ($ref_genome eq $other_genome) && next;
      foreach my $ref_segment ( @{$segments{$ref_genome}} )
       {
        foreach my $i ( $ref_segment->first->gene->coordinate .. $ref_segment->last->gene->coordinate )
	 {
	  $colinear_positions{ $ref_genome }{ $ref_segment->listname }{ $other_genome }{$i}=1;
	 } #foreach my $i ( $ref_segment->first->gene->coordinate .. $ref_segment->last->gene->coordinate )
       } #foreach my $ref_segment ( @{$segments{$ref_genome}} )
     } #foreach my $other_genome ( @genomes_in_multiplicon )
   } #foreach my $ref_genome ( @genomes_in_multiplicon )
   
 } #foreach my $multiplicon ( $dataset->multiplicons )

foreach my $ref_genome ( @genomes )
 {
  my @other_genomes;
  my %total_colinearity;
  
  foreach ( @genomes ) { ($_ eq $ref_genome) || push( @other_genomes, $_ ) }
  print STDOUT "$ref_genome\n";
  print STDOUT "list\tsize\t".join("\t\t", @other_genomes)."\n";
  foreach my $listname ( sort( keys %{$colinear_positions{$ref_genome}} ) )
   {
    my $size = $listsize{$ref_genome}{$listname};
    print STDOUT "$listname\t$size";
    foreach my $other_genome( @other_genomes )
     {
      my $portion = scalar( keys %{$colinear_positions{$ref_genome}{$listname}{$other_genome}} );
      print STDOUT "\t$portion\t".($portion/$size*100)."%";
      $total_colinearity{$other_genome}+=$portion;
     } #foreah my $other_genome( @other_genomes )
    print STDOUT "\n"; 
   } #foreach my $listname ( keys %{$listsize{$ref_genome}} )
  print STDOUT "\n";
  foreach my $other_genome( @other_genomes ) 
   {
    my $fraction = $total_colinearity{$other_genome}/$listsize{$ref_genome}{"total"}*100;
    print STDOUT "$fraction% of the genome of $ref_genome is colinear to $other_genome.\n";
   } #foreah my $other_genome( @other_genomes ) 
  print STDOUT "\n"; 
 } #my $ref_genome ( @genomes )
