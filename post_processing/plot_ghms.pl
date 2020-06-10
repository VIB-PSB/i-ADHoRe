#!/usr/bin/perl -w
#This script was written by Cedric Simillion on Tue Aug 23 10:31:24 CEST 2005

=head1 Description

Generates plots of GHMs between a selected pair of gene lists in a dataset or between all
lists in a dataset. The plots are in .png format and written in a subdirectory of the
output path (as listed in the <settings_file>) called "ghm_plots". Each file is named
<genome1>.<list1>__<genome2>.<list2>.png.

On the plot, bright red dots represent points of the positive orientation class (i.e. both
genes involved have the same orientation) and dark red dots are negative points (genes
have opposite orientations). Note that the GHMs plotted are based on the *remapped* gene
lists.

=head1 Usage

  > plot_ghms.pl <settings_file> [ <genome1> <list1> <genome2> <list2> ]

=over

=item <settings_file>

Mandatory. The settings file used to generate the dataset.

=item <genome1> <list1> <genome2> <list2>

Optional. The genome and list names of two gene lists between which GHMs should be
plotted. 

=head1 !!!WARNING!!!

If no genelists are specified, all GHMs between all gene lists will be plotted possibly
resulting in a very large number of files. If your dataset contains x genelists, x*(x+1)/2
GHM plots will be created. You want to avoid this if you want to remain friends with your
sysadmin. :-)

=back

=head1 Author

This script was written by Cedric Simillion.

=cut

use strict;
use lib "/home/sepro/i-adhore2.1/API";
use iADHoRe;
use iADHoRe::ghm;

my ($dataset,$output_path,$genome1,$list1,$genome2,$list2);
my (@genelists);

$dataset = iADHoRe->read_dataset;
$output_path = $dataset->get_settings->output_path;

unless (-e "${output_path}ghm_plots")
 {
  mkdir( "${output_path}ghm_plots", 0777 ) || die "Could not create output directory: $!\n";
 } #unless (-e "${output_path}alignmentplots")

$dataset->get_gene_pairs;

($genome1,$list1,$genome2,$list2) = @ARGV[1 .. 4];

@genelists = defined( $genome1 ) ?
 grep {
	( ($_->genome eq $genome1) && ($_->listname eq $list1) )
       ||
	( ($_->genome eq $genome2) && ($_->listname eq $list2) )
       } ( $dataset->genelists ) :
$dataset->genelists;

foreach my $i ( 0 .. $#genelists )
 {
  foreach my $j ( $i .. $#genelists )
   {
    my ($ghm,$filename);
    $filename = "${output_path}ghm_plots/".
     $genelists[$i]->genome.".".$genelists[$i]->listname."__".$genelists[$j]->genome.".".$genelists[$j]->listname.".png";
    $ghm=iADHoRe::ghm->new( $genelists[$i], $genelists[$j] );
    $ghm->plot_ghm( $filename );
    print STDERR ".";
   } #foreach my $j ( $i .. $#genelists )
  print STDERR "\n"; 
 } #foreach my $i ( 0 .. $#genelists )
