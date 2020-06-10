#!/usr/bin/perl -w

=head1 Description

Makes a diagram of the duplication level for all regions of a genome in SVG format. The
height of the stack reflects the duplication level. An empty bar indicates that the
identified block is a hidden block. A block is considered hidden if it does not have more
than <minimum_number_anchorpoints> (cfr settings file) with the reference segment, meaning
that it can only be found by using profiles.

The settings_file is the same as used with i-ADHoRe. The diagrams are written to a
subdirectory named "dupliStacks" in the output-path specified in the settings_file and are
named dupliStack_<genome>.svg 

A lightgrey color denotes regions that are found in 2 copies, darkgrey regions are found
in 3 or 4 copies, colored regions are found in more than 4 copies in the genome, regions
present in the same multiplicon have the same color. Black regions denotes genomic
sequence that is not duplicated. On the top of the figure a mark is drawn every 500 genes.

=head1 Usage

  dupliStacks.pl <settings_file>

=head1 Requirements
 
This program requires the Perl GD library, which can be downloaded from CPAN at:
 
 http://search.cpan.org/CPAN/authors/id/L/LD/LDS/GD-2.12.tar.gz
 
=head1 Author 

This script was written by Lieven Sterck
  
=cut

BEGIN {
	#Check if there are arguments given, otherwise print usage-message
	&usage("\n\t!!!!!!not enough parameters!!!!!!!!!!\n") if(scalar(@ARGV) < 1);

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
use lib "/home/sepro/i-adhore2.3/API";
use iADHoRe;
use SVG;

my @colors;
my ($m,$i,$white,$black,$blue,$light_gray,$dark_gray,$multi_id);

my $dataset = '';
my $output_path = '';
my $anchorpoints = '';
my @colors_svg;

#=============================================================

print STDERR "Reading dataset...";
$dataset = iADHoRe->read_dataset;
$output_path = iADHoRe->get_settings->output_path;
$anchorpoints = $dataset->get_settings->anchorpoints;
print STDERR " done\n";

my %genomes;
my @genelists = $dataset->genelists;
foreach my $list (@genelists){
	$genomes{$list->genome}{$list->listname}{'list'} = $list;
	
}

unless (-e "${output_path}dupliStacks") {
  mkdir( "${output_path}dupliStacks", 0777 ) || die "Could not create output directory: $!\n";
}

foreach my $genome_key (sort keys %genomes) {
	
	my $img;
	my $max_list = 0;
	my $i;
	
	foreach my $key (keys %{$genomes{$genome_key}}) {
		if ($genomes{$genome_key}{$key}{'list'}->unmapped_size > $max_list) {
			$max_list = $genomes{$genome_key}{$key}{'list'}->unmapped_size;
		}
	}

	print STDERR "\nGenerating figure for $genome_key";
	my $height= 200 * (scalar(keys %{$genomes{$genome_key}})+1);
	my $svg = SVG->new(width=>$max_list , height=>$height, 'xmlns:xlink'=>"http://www.w3.org/1999/xlink", xmlns=>"http://www.w3.org/2000/svg", 'xmlns:ev'=>"http://www.w3.org/2001/xml-events" );

	#create number of different colors
	foreach my $i (0..(scalar $dataset->multiplicons)+1) {
		my ($r,$g,$b);
		$r=int ( rand (256) );
		$g=int ( rand (256) );
		$b=int ( rand (256) );
		if ($r == 255 && $g == 255 && $b == 255) {redo;} #skip white color
		if ($r == 0 && $g == 0 && $b == 0) {redo;} # skip black color
		
		push (@colors_svg, "rgb($r,$g,$b)");
	}
	
	my @sort_keys = sort (keys %{$genomes{$genome_key}});	
	for (my $x=0; $x < scalar(@sort_keys); $x++){
		my $list = $genomes{$genome_key}{$sort_keys[$x]}{'list'};
		$svg->rect(x=>0,y=>(200*($x+1))-10 ,height=>10, width=>($list->unmapped_size)-1, style=>{'fill'=>'black'});
	}
	
	foreach $i (1..int($max_list/500)) {
		$svg->line(x1=>$i*500, y1=>0, x2=>$i*500,y2=>15, style=>{'stroke'=>'black', 'stroke-width'=>'5'});
		$svg->line(x1=>$i*500, y1=>$height, x2=>$i*500,y2=>$height-15, style=>{'stroke'=>'black', 'stroke-width'=>'5'});
	}
	
	foreach my $multiplicon ( $dataset->multiplicons ){
	        $multiplicon->is_redundant && next;
		if ($multiplicon->level == 2 && $multiplicon->x_object->genome ne $multiplicon->y_list->genome) {next;}
		
		my $level=$multiplicon->level;
		($multiplicon->profile->size >= 9) || next; #multiplicons smaller than 9 genes are skipped
		
		my %blocks;
		my $s = scalar @{$multiplicon->profile->segments};
		my $Mu_tmp = $multiplicon;
		my $q =0;
		foreach my $segment ( reverse @{$multiplicon->profile->segments} ) {
			$blocks{$s}{'listname'} = $segment->listname;
			$blocks{$s}{'begin'} = $segment->first->gene->coordinate;
			$blocks{$s}{'end'} = $segment->last->gene->coordinate;
			$blocks{$s}{'list_ref'} = $segment->first->gene->list;
			$blocks{$s}{'segm_id'} = $segment->id;
			$blocks{$s}{'genome'} = $segment->genome;
			$blocks{$s}{'size'} = $segment->size;
			
			if ($segment->genome eq $genome_key) { $q++; }
			
			if ($s != 1){
				$blocks{$s}{'multiplicon'} = $Mu_tmp;
				
				if ($s > 2) {
					$Mu_tmp = $Mu_tmp->x_object->multiplicon;
				}
				
			} 
			
			%{$blocks{$s}{'hidden'}} = ();
			%{$blocks{$s}{'ghost'}} = ();
			
			$s--;
		
		}
		$blocks{'1'}{$genome_key} = $q;
			
		my %APs;
		if ($multiplicon->level > 2 ) { 
			for (my $x = 2; $x <= $level; $x++) {
				foreach my $anchor ( $blocks{$x}{'multiplicon'}->get_anchorpoints ) { 
					my %segment;
					$segment{'x'}=&which_segment( $multiplicon,$anchor->gene_x->list->genome,$anchor->gene_x->list->listname,$anchor->gene_x->coordinate );
					$segment{'y'}=&which_segment( $multiplicon,$anchor->gene_y->list->genome,$anchor->gene_y->list->listname,$anchor->gene_y->coordinate );
					if ($segment{'x'} ne '' && $segment{'y'} ne '') {
						$APs{$segment{'x'}."_".$segment{'y'}} += 1;
						$APs{$segment{'y'}."_".$segment{'x'}} += 1;
					}
				}
			}
			for (my $x = 1; $x <= $level; $x++) {
				my $base_id = $blocks{$x}{'segm_id'};
				for (my $y = 1; $y <= $level; $y++) {
					if ($x == $y) { next;}
					my $id = $blocks{$y}{'segm_id'};
					my $hidden = 1;
					if (($APs{$base_id."_".$id}) && $APs{$base_id."_".$id} >= $anchorpoints) {
						$hidden = 0;
					}
					$blocks{$x}{'hidden'}{$y} = $hidden;
				}
			}
		}
		
			

		#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# Draw figure
		#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		
		#select colors
		my $dupli_color;
		if ($blocks{'1'}{$genome_key} <=2){
			$dupli_color="darkgrey";
		} elsif ($blocks{'1'}{$genome_key} <=4){
			$dupli_color="gray";
		} else {
			$dupli_color=pop(@colors_svg);
		}
	
	
		#draw figure
		for (my $x = 1; $x <= $level; $x++) {
	
			if ($blocks{$x}{'genome'} ne $genome_key) {next;}
		
			my $chr;
			for (my $y=0;$y < scalar @sort_keys;$y++) {
				if ($sort_keys[$y] eq $blocks{$x}{'listname'}) {
					$chr = $y+1;
					$y = scalar @sort_keys;
				}
			}
			if ($blocks{'1'}{$genome_key} > 1){
				$svg->rect(x=>$blocks{$x}{'begin'}, y=>($chr*200)-10, height=>10, width=>($blocks{$x}{'end'} - $blocks{$x}{'begin'}), style=>{'fill'=>$dupli_color});
			}
			my @gen_sort = keys %blocks;
			if ($level > 2) {
				my @tmp_sort = sort { $a <=> $b} (keys %blocks);
				splice (@tmp_sort,$x-1,1);
				@gen_sort = sort { ($blocks{$x}{'hidden'}{$a} <=> $blocks{$x}{'hidden'}{$b}) } @tmp_sort;
				
			}
			
			my $d=1;
			for (my $w=0;$w<scalar(@gen_sort);$w++) {
				my $j = $gen_sort[$w];
				if ($x == $j) { next;}
				if ($blocks{$j}{'genome'} eq $genome_key) { 
					my ($x1,$y1,$x2,$y2,$h,$w);
					$x1=$blocks{$x}{'begin'};
					$y1=(($chr*200)-(($d)*15))-10;
					$h = 10;
					$w = $blocks{$x}{'end'} - $blocks{$x}{'begin'};
					if (($blocks{$x}{'hidden'}{$j}) ) {
						$svg->rect(x=>$x1+1.5,y=>$y1+1.5, height=>$h-3, width=>$w-3, style=>{'fill'=>'white', 'stroke'=>$dupli_color, 'stroke-width'=>'3'});
					} else {
						$svg->rect(x=>$x1,y=>$y1, height=>$h, width=>$w, style=>{'fill'=>$dupli_color});
					}
					$d++;
				} 
			} #for (my $w=0;$w<scalar(@gen_sort);$w++)
		} #for (my $x = 1; $x <= $level; $x++)
	} #foreach my $multiplicon ( $dataset->multiplicons )
	
	print STDOUT "\nWriting .svg file...";
	my $svg_XML = $svg->xmlify();
	open (SVG_OUT, ">$output_path/dupliStacks/dupliStack_$genome_key.svg");
	print SVG_OUT $svg_XML;
	close (SVG_OUT);
	print STDOUT " done!\n";

	
} #foreach my $genome_key (sort keys %genomes)

exit;

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#for each anchorpoint check on which segments the gene is present and then 
#return the name of that segment
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub which_segment
 {
	my ($mu,$genome,$list,$coord)=@_;
	
	foreach my $segment ( @{$mu->profile->segments} ){
		if ( ($segment->genome eq $genome) &&
				($segment->listname eq $list) &&
				($segment->first->gene->coordinate <= $coord) &&
				($segment->last->gene->coordinate >= $coord)){
		
			return( $segment->id );
	
		}
	}
 }
	

