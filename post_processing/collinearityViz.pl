#!/usr/bin/perl -w


=head1 Description

This script generates a diagram showing a lineair representation of the collinear blocks.
The regions of the reference genome that are collinear with the 'mapping' genome are
colored according to the chromosome of the mapping genome they are collinear to. Above the
colored regions a number is plotted that refers to the multiplicon ID. The output file is
called <reference genome>_<mapping genome>.svg and is generated as svg format in the
output directory as specified in the settings file.

This script should only be used when doing a level_2_only analysis. 
  
=head1 Usage

   collinearityViz.pl <settings_file> <name reference genome> <name mapping genome>
 
The names of the reference and the mapping genome should be identical to those used in the
settings file.
 
=head1 Requirements
 
This program requires the Perl SVG library, which can be downloaded from CPAN at:
 
   http://search.cpan.org/CPAN/authors/id/R/RO/RONAN/SVG-2.33.tar.gz
 
=head1 Author

This script was written by Lieven Sterck
  
=cut

BEGIN {
	#Check if there are arguments given, otherwise print usage-message
	&usage("\n\t!!!!!!not enough parameters!!!!!!!!!!\n") if(scalar(@ARGV)<3);

	#------------------------------------------------------------
	sub usage( $ ){
    print "$_[0]\n";
    system("pod2text $0");
    exit(1);
  	}
	#------------------------------------------------------------
	
}

#===============================================================================
# Initialisation
#===============================================================================
use strict;
use lib "API";
use lib "/home/sepro/i-adhore2.3/API";
use iADHoRe;
use SVG;
use POSIX;

my $dataset = '';
my $output_path = '';
my $level2 = '';

my $file = $ARGV[0];
my $genomeA = $ARGV[1];
my $genomeB = $ARGV[2];

#=========================================================================================

print STDERR "Reading dataset...";
$dataset = iADHoRe->read_dataset;
$output_path = iADHoRe->get_settings->output_path;
$level2 = $dataset->get_settings->level_2_only;
print STDERR " done\n";

if ($level2 == 0) { die "\n\tThis is not a level_2_only i-ADHoRe run!!\n"; }

unless (-e "${output_path}CollinearityViz") {
  mkdir( "${output_path}CollinearityViz", 0777 ) || die "Could not create output directory: $!\n";
}

my %genomes;
my @genelists = $dataset->genelists;
foreach my $list (@genelists){
	$genomes{$list->genome}{$list->listname}{'list'} = $list;
	
}

#create number of different colors
$genomes{$genomeB} = &colors ($genomes{$genomeB});

my $max_list = 1;
foreach my $key (keys %{$genomes{$genomeA}}) {
	if ($genomes{$genomeA}{$key}{'list'}->unmapped_size > $max_list) {
		$max_list = $genomes{$genomeA}{$key}{'list'}->unmapped_size;
	}
}

my $width= $max_list + 150;
my $height= 100 * (scalar(keys %{$genomes{$genomeA}})+2);

#build SVG
my $svg = SVG->new(width=>$width ."px", height=>$height ."px");

my @sort_keys = sort (keys %{$genomes{$genomeA}});	
for (my $x=0; $x < scalar(@sort_keys); $x++){
	my $list = $genomes{$genomeA}{$sort_keys[$x]}{'list'};
	$svg->rect(x=>125, y=>(($x+1)*100), height=>10, width=>$list->unmapped_size,rx=>2,ry=>2, style=> {'fill'=>'black'} );
	$svg->text(x=>10, y=>(($x+1)*100))->cdata($sort_keys[$x]);
	$genomes{$genomeA}{$sort_keys[$x]}{'Y_coord'} = (($x+1)*100);
}

foreach my $multiplicon ( $dataset->multiplicons ){
	if ($multiplicon->level > 2 || $multiplicon->x_object->genome eq $multiplicon->y_list->genome) {next;}
	 	my $segment_X = ${$multiplicon->profile->segments}[0];
		my $segment_Y = ${$multiplicon->profile->segments}[1];
		if ($multiplicon->x_object->genome eq $genomeA) {
			$svg->rect(x=>(125 + $segment_X->first->gene->coordinate), y=>($genomes{$genomeA}{$segment_X->listname}{'Y_coord'} - 11), height=>'10', width=>((25 + $segment_X->last->gene->coordinate)-(25 + $segment_X->first->gene->coordinate)),rx=>2,ry=>2,style=> {'fill'=>$genomes{$genomeB}{$segment_Y->listname}{'color'}});
			$svg->text(x=>(125 + $segment_X->first->gene->coordinate) -5, y=>($genomes{$genomeA}{$segment_X->listname}{'Y_coord'} -15))->cdata($multiplicon->id);
		} elsif ($multiplicon->y_list->genome eq $genomeA) {
			$svg->rect(x=>(125 + $segment_Y->first->gene->coordinate), y=>($genomes{$genomeA}{$segment_Y->listname}{'Y_coord'} - 11), height=>'10', width=>((25 + $segment_Y->last->gene->coordinate)-(25 + $segment_Y->first->gene->coordinate)),rx=>2,ry=>2,style=> {'fill'=>$genomes{$genomeB}{$segment_X->listname}{'color'}});
			$svg->text(x=>(125 + $segment_Y->first->gene->coordinate) -5, y=>($genomes{$genomeA}{$segment_Y->listname}{'Y_coord'} -15))->cdata($multiplicon->id);
		}
}

$svg->rect(id=>'scale' ,x=>25, y=>(100 * (scalar(keys %{$genomes{$genomeA}})+1)), height=>'10', width=>100,rx=>2,ry=>2, style=> {'fill'=>'black'});
$svg->text(id=>'scale_name', x=>25, y=>((100 * (scalar(keys %{$genomes{$genomeA}})+1))-10))->cdata('Scale');
$svg->text(id=>'scale_length', x=>130, y=>((100 * (scalar(keys %{$genomes{$genomeA}})+1))+5))->cdata('(100 genes)');

my $t =0;
foreach my $key (sort (keys %{$genomes{$genomeB}})) {
	$svg->rect(id=>"scale_". $key ,x=>(300 + (150 * $t)), y=>(100 * (scalar(keys %{$genomes{$genomeA}})+1)), height=>'10', width=>'100',rx=>2,ry=>2, style=> {'fill'=>$genomes{$genomeB}{$key}{'color'}});
	$svg->text(id=>"color_". $key, x=>(300 + (150 * $t)), y=>((100 * (scalar(keys %{$genomes{$genomeA}})+1))-10))->cdata($key);
	$t++;
}

my $svg_XML = $svg->xmlify();
my $out_file = $genomeA . "_" . $genomeB;
open (SVG_OUT, ">${output_path}CollinearityViz/${out_file}.svg") ||  die "Could not write svg !\n";
 print SVG_OUT $svg_XML;
close (SVG_OUT);

exit;
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sub routines
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub random ( $ $ )
	{
		my $max = $_[0];
		my $test = $_[1];
		my $NR = '';
		do {
			$NR = int(rand ($max));
		} until ($test !~ m/\,$NR\,/);
		
	return ($NR);
	}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub colors ( $ )
	{
		my $genome = $_[0];
		my $number_of_colors = (scalar (keys %$genome));		
		my (@rgb)=(255,255,255);
		my @colors;
		my $step=int( 255 / ( ceil( ($number_of_colors+1)**(1/3) ) - 1 ) );
		foreach my $aantal (1..$number_of_colors) {
  			my $i=0;
  			$rgb[$i]-=$step;
  			until ($rgb[$i]>=0){
   			$rgb[$i]=255;
    			$rgb[$i+1]-=$step;
    			$i++;
   		} 
   		#acces @rgb here
			if ($rgb[0] eq '255' && $rgb[1] eq '255' && $rgb[2] eq '255') {
				redo;
			} elsif ($rgb[0] eq '0' && $rgb[1] eq '0' && $rgb[2] eq '0') {
				redo;
			} else {
				push (@colors, [@rgb]);
			}
		}
		foreach my $key (keys %$genome){
			my $col = shift @colors;
			$col = join (",", @$col);
			my $color_code = "rgb(" . $col .")";
			$$genome{$key}{'color'} = $color_code;
		}
		
		return ($genome);
}



