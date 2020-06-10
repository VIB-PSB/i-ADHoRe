#!/usr/bin/perl -w
#This script was written by Cedric Simillion on Fri Jul  8 14:44:31 CEST 2005

=head1 Usage

  $ generate_overview.pl <inifile>

<inifile> is a i-ADHoRe settings file.

=head1 Description

Creates in the output directory a folder called 'html'. This folder contains a set of HTML
files that show the profile of each multiplicon detected and also lists the multiplicons
that have been found using this profile. The first file in this folder is linked to as
'index.html'

=head1 Author

This script was written by Cedric Simillion.

=cut

use strict;
use lib "/home/sepro/i-adhore2.3/API";
use iADHoRe;
use iADHoRe::ghm;

#===============================================================================
# Preparation
#===============================================================================
my ($dataset,$output_path,$html_header,$first_id);
my (@top_multiplicons);
my (%children_multiplicons,%child_ids,%prev,%next,%up);

print STDERR "Reading dataset...";
$dataset = iADHoRe->read_dataset;
print STDERR "\n";
$dataset->get_gene_pairs;

$output_path = $dataset->get_settings->output_path."html/";
(-d $output_path) || ( mkdir( $output_path, 0777 ) || die "Could not create $output_path directory: $!\n" );

#===============================================================================
# Main
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Map child multiplicons to parent multiplicons
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
foreach my $multiplicon ( $dataset->multiplicons )
 {
  if ( $multiplicon->level > 2 )
   {
    my $parent_id=$multiplicon->x_object->multiplicon->id;
    push( @{$children_multiplicons{$parent_id}}, $multiplicon );
    $up{ $multiplicon->id } = $parent_id;
    } #if ( $multiplicon->level > 2 ) 
  else 
   {
    push( @top_multiplicons, $multiplicon );
   } #else 
 } #foreach my $multiplicon ( $dataset->multiplicons )

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sort multiplicons and create datastructures to quickly retrieve
# parent, child and sibling multiplicons
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
foreach my $set ( \@top_multiplicons, values %children_multiplicons )
 {
  ( scalar( @$set ) > 1 ) || next;
  @$set = sort by_multiplicon_size ( @$set );
  foreach my $i ( 1 .. $#{$set} )
   {
    my ($first_id,$second_id);
    $first_id = $$set[$i-1]->id;
    $second_id = $$set[$i]->id;
    $next{$first_id} = $second_id;
    $prev{$second_id} = $first_id;
   } #foreach my $i ( 1 .. $#{$set} ) 
 } #my $set ( \@top_multiplicons, values %children_multiplicons )
$first_id = $top_multiplicons[0]->id;

foreach my $id (keys  %children_multiplicons )
 {
  @{$child_ids{$id}} = map { $_->id } ( @{$children_multiplicons{$id}} )
 } #foreach my $id (keys  %children_multiplicons )

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create HMTL files
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print STDERR "Generating HTML files";
foreach my $multiplicon ( $dataset->multiplicons )
 {
  $multiplicon->is_redundant && next;
  my ($id,$family_matrix,$family_colors,$matrix,$img);
  $id = $multiplicon->id;

  open( OUT, ">$output_path$id.html" ) || die "Could not create file: $!\n";
  print OUT "<HTML>\n <HEAD>\n  <TITLE>multiplicon $id</TITLE>\n </HEAD>\n <BODY>\n";

  print OUT " <h3>General information</h3>\n <TABLE  border='1'>\n";
  print OUT "  <TR>\n   <TH align='right'>ID:</TH>\n   <TH align='center'>".$id."</TH>\n  </TR>\n";
  print OUT "  <TR>\n   <TH align='right'>Level:</TH>\n   <TD align='left'>".$multiplicon->level."</TD>\n  </TR>\n";
  print OUT "  <TR>\n   <TH align='right'>Number of anchorpoints:</TH>\n   <TD align='left'>"
    .$multiplicon->number_of_anchorpoints."</TD>\n  </TR>\n";

  if ( exists $up{$id} )
   {
    my $temp_id = $id;
    print OUT "  <TR>\n   <TH align='right'>Detection history:</TH>\n";
    print OUT "   <TD align='left'>\n    <UL>\n";
    while ( exists $up{$temp_id} )
     {
      print OUT "     <LI><A HREF='$up{$temp_id}.html'>$up{$temp_id}</A></LI>\n";
      $temp_id = $up{$temp_id};
     } #while ( exists $up{$temp_id} )
    print OUT "    </UL>\n   </TD>\n  </TR>\n";
   } #if ( exists $up{$id} )

  print OUT " </TABLE>\n <HR>\n";


  ($family_matrix,$family_colors) = $multiplicon->profile->matrix_single_linkage;
  $img = $multiplicon->profile->init_image( $family_colors, 1 );
  $matrix = $multiplicon->profile->matrix;
  
  print OUT " <h3>Alignment plot</h3>\n <IMG SRC = '../alignment_plots/$id.png'>\n <hr>\n";
  
  print OUT " <h3>Alignment table</h3>\n <TABLE border='1'>\n";
  foreach my $i ( 0 .. $#{$matrix} )
   {
    print OUT "  <TR>\n";
    print OUT "   <TH>".${$multiplicon->profile->segments}[$i]->genome." "
     .${$multiplicon->profile->segments}[$i]->listname."</TH>\n";
    foreach my $j ( 0 .. $#{$$matrix[$i]} )
     {
      if ( defined $$matrix[$i][$j] )
       {
        my ($family,$color_index,$color_code);
	$family = defined( $$family_matrix[$i][$j] ) ? $$family_matrix[$i][$j] : "white";
	$color_index = $$family_colors{$family};
	$color_code = &rgb2hex( $img->rgb($color_index) );
	print OUT "   <TD bgcolor ='$color_code'>".$$matrix[$i][$j]->gene->ID."</TD>\n";
       } #if ( defined $$matrix[$i][$j] )
      else 
       {
        print OUT "   <TD></TD>\n";
       } #else 
     } #foreach my $j ( 0 .. $#{$$matrix[$i]} )
    print OUT "  </TR>\n"; 
   } #foreach my $i ( 0 .. $#{$matrix} )
  print OUT " </TABLE>\n <hr>\n";

  if ( exists $child_ids{$id} )
   {
    print OUT " <TABLE border='1'>\n  <TR><TH>Multiplicons detected with this profile</TH></TR>\n";
    foreach my $child_id ( @{$child_ids{$id}} )
     {
      print OUT "   <TR><TD><A HREF='$child_id.html'>$child_id</A></TD></TR>\n";
     } #foreach my $child_id ( @{$child_ids{$id}} )
    print OUT " </TABLE>\n <hr>\n";
   } #if ( exists $child_ids{$id} )

  exists( $prev{$id} ) && print OUT " <A HREF='$prev{$id}.html'>previous</a><br>\n";
  exists( $next{$id} ) && print OUT " <A HREF='$next{$id}.html'>next</a><br>\n";
  exists( $up{$id} ) && print OUT " <A HREF='$up{$id}.html'>up</a><br>\n";

  print OUT " </BODY>\n</HTML>\n";
  close OUT;
  print STDERR ".";
 } #foreach my $multiplicon ( $dataset->multiplicons )

system( "ln -s $first_id.html ${output_path}index.html" );
print STDERR "\n";

#===============================================================================
# Subroutines
#===============================================================================
sub rgb2hex
 {
  my $hex='#';
  foreach ( @_ )
   {
    $hex.=sprintf("%02X", $_);
   } #foreach ( @_ )
  return $hex; 
 } #sub rgb2hex

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub by_multiplicon_size
# A sorting routine to sort multiplicons first by their number of 
# anchorpoints, then by DPD distance they span.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub by_multiplicon_size
 {
  ( $b->number_of_anchorpoints <=> $a->number_of_anchorpoints )
   ||
  (
   iADHoRe::ghm->dpd( $b->begin_x, $b->begin_y, $b->end_x, $b->end_y )
    <=>
   iADHoRe::ghm->dpd( $a->begin_x, $a->begin_y, $a->end_x, $a->end_y )
  ) 
 } #sub by_multiplicon_size
