#!/usr/bin/perl -w
#This script was written by Cedric Simillion on Wed Jun 29 10:45:51 CEST 2005

=head1 Usage
 
 multiplicon_plot.pl settings_file [ multiplicons_ids ]
 
=head1 Description
 
Makes diagrams of multiplicons detected with i-ADHoRe in .eps format. The settings_file is
the same as used with i-ADHoRe.pl If no additional command line arguments are given, then
a .eps file will be created for each multiplicon in the multiplicons table. If one or more
multiplicon IDs (from the ID field of the multiplicons table) are given, then diagrams
will be generated for these multiplicons only. All .eps files will be generated in a 
subdirectory called "multiplicon_plots" located in the output directory specified in the
settings file.

=head1 Requirements

This scripts relies on the dupliviz.pl script that also comes with the i-ADHoRe
package. Its path must be set in the $dupliviz_path variable in this script.
 
=head1 Author

This script was written by Cedric Simillion.

=cut

use strict;
use lib "/home/sepro/i-adhore2.1/API";
use iADHoRe;

my ($dataset,$output_path,$plot_selection,$dupliviz_path);
my (%plot);

#Change this variable to the correct path of dupliviz.pl
$dupliviz_path = "/nas/biocomp/projects/adhore/dupliviz.pl"

print STDERR "Reading dataset...";
$dataset = iADHoRe->read_dataset;
$output_path = iADHoRe->get_settings->output_path;

if ( scalar( @ARGV ) > 1 ) 
 {
  $plot_selection=-1;
  shift( @ARGV );
  %plot = map { ( $_, -1 ) } ( @ARGV );
 } #if ( scalar( @ARGV ) > 1 )
else 
 {
  $plot_selection=0;
 } #else

unless (-e "${output_path}multiplicon_plots")
 {
  mkdir( "${output_path}multiplicon_plots", 0777 ) || die "Could not create output directory: $!\n";
 } #unless (-e "${output_path}alignmentplots")


print STDERR "\nGenerating images";
foreach my $multiplicon ( $dataset->multiplicons )
 {
  ( $plot_selection && !exists( $plot{ $multiplicon->id } ) ) && next;

  my ($s,$m,$id);
  my (@anchors2plot);
  my (%seg_num,%seg_enc,%in_multi);

  $id=$multiplicon->id;
  $s=1;
  open (OUT,">".$output_path."multiplicon_plots/$id");

  foreach my $segment ( @{$multiplicon->profile->segments} )
   {
    my ($list,$listname,$list_id,$begin,$end);

    $listname=$segment->genome.$segment->listname;
    $seg_enc{$listname}++;
    $list_id=$listname.sprintf("%02d",$seg_enc{$listname});

    print OUT "segment $s\nchromosome $list_id\n\n#geneid\tstart\tstop\torientation\n";

    $seg_num{$list_id}=$s;
    $begin = $segment->first->gene->coordinate;
    $end = $segment->last->gene->coordinate;
    $list = $segment->first->gene->list;

    foreach my $c ( $begin .. $end )
     {
      my $gene = ${$list->elements}[$c]->gene;
      $in_multi{$gene->ID}=$list_id;
      print OUT join( "\t", $gene->ID, ($c-$begin)*20, ($c-$begin)*20+15, $gene->orientation )."\n";
     } #foreach my $c ($$segment{"begin"}..$$segment{"end"})

    $s++;
    print OUT "\n";
   } #foreach my $segment (@{$multi_segments{$meta}})

  print OUT "generelations\n#genex\tsegmentx\tgeney\tsegmenty\n";
  $m=$multiplicon;

  while ( defined $m )
   {
    foreach my $anchor ( $m->get_anchorpoints )
     {
      ( $in_multi{$anchor->gene_x->ID} && $in_multi{$anchor->gene_y->ID} ) && push(@anchors2plot,$anchor);
     } #foreach my $anchor (@{$anchorpoints{$meta}})
    $m = ( ref( $m->x_object ) eq 'iADHoRe::profile' ) ? $m->x_object->multiplicon : undef;
   } #while ($meta_clusters{$m}{"parent"})

  foreach my $anchor (@anchors2plot) 
   {
    my ($list_id_x,$list_id_y);
    $list_id_x = $in_multi{$anchor->gene_x->ID};
    $list_id_y = $in_multi{$anchor->gene_y->ID};
    print OUT join("\t", $anchor->gene_x->ID, $seg_num{$list_id_x}, $anchor->gene_y->ID, $seg_num{$list_id_y} )."\n";
   } #foreach my $anchor (@{$anchorpoints{$meta}})

  close OUT;
  `$dupliviz_path ${output_path}multiplicon_plots/$id`;
  unlink ("${output_path}multiplicon_plots/$id");
  print STDERR ".";

 } #foreach my $multiplicon ( $dataset->multiplicons )

print STDERR "\nDone!\n";
