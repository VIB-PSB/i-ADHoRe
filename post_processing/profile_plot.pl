#!/usr/bin/perl -w

=head1 Usage

  profile_plot.pl <iADHoRe settings_file> -real T/F [-prot <list prot coordinates> -TE <list TE coordinates>] [-ID "multiplicons_ids" ]

=head1 Description

This script takes the output of i-ADHoRe and generates a diagram in SVG format every
multiplicon in the dataset. All segments in the block are drawn in their real genomic
configuration. Normal protein coding genes have olive color and transposons have a red
color. When the command line option -real T (see below) is used then the output files are
named <multipicon_id>_real.svg, otherwise they are named <multipicon_id>.svg. All files
will be generated in a subdirectory (profile_plot/) in the output directory as specified
in the settings file. 

=head1 Command line options

=over

=item -real

Obligatory. Indicates if you want to draw an image reflecting the real situation (-real T,
genes and intergenics have the actual spacing like they are present on the chromosomes )
or just a schematic representation (-real F, all genes/intergenics have equal lengths).

=item -prot

if -real T is used then this is obligatory. You must specify a file containing the
positions of the proteins. This must be a tab-delimited file, containing the following
columns:

=over

=item gene-name

The identifier of the gene. Must be identical to the identifier of the gene used in the
gene lists file.

=item Orientation

A '+' or '-' sign, indicating the orientation of the gene

=item Start

The starting nucleotide position, in reference to the + strand, of the gene.

=item Stop

The end nucleotide position, in reference to the + strand, of the gene.

=back

Note that the column names should not be included in the file.

=item -TE

Optional, but only useful when -real T is used. You can specify a file containing the
positions of the transposable elements. The regions indicated in this file will be
inserted into each segment of the multiplicon. Again, this must be a tab-delimited file
with the following columns:

=over

=item Chromo-name

Name of the chromosome. This must be the same as the name of the gene lists listed in the
i-ADHoRe configuration file.

=item Start

The starting nucleotide position, in reference to the + strand, of the transposable
element.

=item Stop

The end nucleotide position, in reference to the + strand, of the transposable element.

=item Orientation

A '+' or '-' sign, indicating the orientation of the transposable element.

=item TE-name

Name of the transposable element.

=back

Again, column names should not be included in the file.

=item -ID

Optional. If this option is ommitted, then a file will be created for each multiplicon in
the multiplicons table. If one or more multiplicon IDs (from the ID field of the
multiplicons table) are given, then only these multiplicons are processed.

=head1 Requirements
 
This program requires the Perl SVG library, which can be downloaded from CPAN at:
 
 http://search.cpan.org/CPAN/authors/id/R/RO/RONAN/SVG-2.33.tar.gz
 
=head1 Author

This script was written by Lieven Sterck
  
=cut

BEGIN {
	#Check if there are arguments given, otherwise print usage-message
	&usage("\n\t!!!!!!not enough parameters!!!!!!!!!!\n") if(scalar(@ARGV) < 2);

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
use lib "/home/sepro/i-adhore2.1/API";
use iADHoRe;
use SVG;

my $dataset = '';
my $output_path = '';
my $plot_selection = '';
my $flip = '';
my (%prot,%TE,%plot);

#=============================================================

print STDERR "Reading dataset...";
$dataset = iADHoRe->read_dataset;
$output_path = iADHoRe->get_settings->output_path;
print STDERR " done\n";

shift @ARGV;
my %params = @ARGV;

if (uc ($params{"-real"}) eq "T") {
	#check if prot file can be opened
	unless ($params{"-prot"}) {
		print "\n\t!!!Please give in a file containing the protein positions!!!\n\n";
		system("pod2text $0");
		exit 1;
	}
	unless (open(PROT_LIST,$params{"-prot"})) {
		print "\n\tCannot open file ". $params{"-prot"} ."!!\n\n";
		exit;
	}
	while (<PROT_LIST>) {
		$_ =~ m/^\>?(\S+)\t([+|-])\t(\d+)\t(\d+)$/;
		$prot{$1}{'strand'} = $2;
		$prot{$1}{'start'} = $3;
		$prot{$1}{'stop'} = $4;
	}
	close (PROT_LIST);

	#check if TE file is given and can be opened
	if (defined $params{"-TE"}){
		unless (open(TE_LIST,$params{"-TE"})) {
		print "\n\tCannot open file ". $params{"-TE"} ."!!\n\n";
		exit;
		}
		while (<TE_LIST>) {
			$_ =~ m/^(\S+)\t(\d+)\t(\d+)\t([+|-])\t(\S+)$/;
			$TE{$1}{$2}{'start'} = $2;
			$TE{$1}{$2}{'stop'} = $3;
			$TE{$1}{$2}{'strand'} = $4;
			$TE{$1}{$2}{'name'} = $5;
		}
		close (TE_LIST);
	}
}

if ( exists $params{"-ID"} ) {
	my @IDs = split (" ",$params{"-ID"});
	$plot_selection=-1;
	%plot = map { ( $_, -1 ) } ( @IDs );
} else {
  $plot_selection=0;
}

unless (-e "${output_path}profile_plot") {
  mkdir( "${output_path}profile_plot", 0777 ) || die "Could not create output directory: $!\n";
}

my @file = ();

foreach my $multiplicon ( $dataset->multiplicons ){
        $multiplicon->is_redundant && next;
	if ($plot_selection && !exists( $plot{ $multiplicon->id })) {next;} 
	my ($id,$level,%blocks);
	$id = $multiplicon->id;
	$level = $multiplicon->level;

	my $output_name;
	if (uc ($params{"-real"}) eq "T") {
		$output_name = "${id}_real.txt";
	} else {
		$output_name = "${id}.txt";
	}
	
	#check if textual output file exists, if so then go directly to generating the SVG-images
	if (!-e "${output_path}profile_plot/$output_name" || -z "${output_path}profile_plot/$output_name") {
		
		open (OUT,"> ${output_path}profile_plot/$output_name");
	
		my $s = scalar @{$multiplicon->profile->segments};
		my $Mu_tmp = $multiplicon;
		my %segm_list = ();
	
		foreach my $segment ( reverse @{$multiplicon->profile->segments} ) {
			$blocks{$s}{'listname'} = $segment->listname;
			$blocks{$s}{'begin'} = $segment->first->gene->coordinate;
			$blocks{$s}{'end'} = $segment->last->gene->coordinate;
			$blocks{$s}{'list_ref'} = $segment->first->gene->list;
			$blocks{$s}{'segm_len'} = $prot{${$blocks{$s}{'list_ref'}->elements}[$blocks{$s}{'end'}]->gene->ID}{'stop'} - $prot{${$blocks{$s}{'list_ref'}->elements}[$blocks{$s}{'begin'}]->gene->ID}{'start'};
		
			if ($s != 1){
				$blocks{$s}{'flip'} = &flip ($Mu_tmp->get_baseclusters);
				$blocks{$s}{'multiplicon'} = $Mu_tmp;
			
				if ($s > 2) {
					$Mu_tmp = $Mu_tmp->x_object->multiplicon;
				}
			} else {
				$blocks{$s}{'flip'} = 0;
			}
		
			$segm_list{$segment->id} = $s;
			$s--;
		
		}
	
		for (my $x = 1; $x <= $level; $x++) {
			my $segm_begin = $prot{${$blocks{$x}{'list_ref'}->elements}[$blocks{$x}{'begin'}]->gene->ID}{'start'};
			my $diff_begin = $prot{${$blocks{$x}{'list_ref'}->elements}[$blocks{$x}{'begin'}]->gene->ID}{'start'} -500;
			my $segm_end = $prot{${$blocks{$x}{'list_ref'}->elements}[$blocks{$x}{'end'}]->gene->ID}{'stop'};
			
			my @TEs;
			if (defined $params{"-TE"} && ((uc ($params{"-real"})) eq "T")) {
				@TEs = &get_TEs ($blocks{$x}{'listname'},$segm_begin,$segm_end,) ;
				my $tot_len =0;
				for (my $q = 0; $q < scalar @TEs; $q++) {
					$tot_len += ($TE{$blocks{$x}{'listname'}}{$TEs[$q]}{'stop'} - $TE{$blocks{$x}{'listname'}}{$TEs[$q]}{'start'});
				}
				print OUT "id\tlistname\t#TEs\tlength_segm\tmasked_length\n";
				print OUT $id ."\t". $blocks{$x}{'listname'} ."\t". scalar @TEs ."\t". $blocks{$x}{'segm_len'} ."\t". $tot_len ."\n";
				print OUT "Segment ". $x ."\t". $blocks{$x}{'listname'} ."\t". ($blocks{$x}{'segm_len'}+1000)."\n";
			
			} elsif (uc ($params{"-real"}) ne "T") {
				print OUT "id\tlistname\tlength_segm\n";
				print OUT $id ."\t". $blocks{$x}{'listname'} ."\t". (((($blocks{$x}{'end'} - $blocks{$x}{'begin'})+1)*2000)+500) ."\n";
				print OUT "Segment ". $x ."\t". $blocks{$x}{'listname'} ."\t". (((($blocks{$x}{'end'} - $blocks{$x}{'begin'})+1)*2000)+500)."\n";
			
			} else {                                         
				print OUT "id\tlistname\tlength_segm\n";
				print OUT $id ."\t". $blocks{$x}{'listname'} ."\t". $blocks{$x}{'segm_len'} ."\n";
				print OUT "Segment ". $x ."\t". $blocks{$x}{'listname'} ."\t". ($blocks{$x}{'segm_len'}+1000)."\n";
			}
		
			
			if ($blocks{$x}{'flip'}) {
				my $offset = $segm_end;
				my $prev_gene_start = '';
				my $n = 1;
				foreach my $c (reverse( $blocks{$x}{'begin'} .. $blocks{$x}{'end'} )) {
					my $gene = ${$blocks{$x}{'list_ref'}->elements}[$c]->gene;
					my $or = $gene->orientation;
					$or = ( $or eq '+' ) ? '-' : '+';
						
					if (uc ($params{"-real"}) eq "T") {
						if (exists $params{"-TE"}) {
							if (scalar @TEs != 0 && $prev_gene_start ne '') {
								for (my $y = (scalar @TEs )-1; $y >= 0; $y--) {
									if ( $TEs[$y] < $prev_gene_start && $TEs[$y] > $prot{$gene->ID}{'stop'} ) {
										my $TE_or = $TE{$blocks{$x}{'listname'}}{$TEs[$y]}{'strand'};
										$TE_or = ( $or eq '+' ) ? '-' : '+';
										print OUT join ( "\t", $TE{$blocks{$x}{'listname'}}{$TEs[$y]}{'name'},(500 + ($offset - $TE{$blocks{$x}{'listname'}}{$TEs[$y]}{'stop'})),(500 + ($offset - $TE{$blocks{$x}{'listname'}}{$TEs[$y]}{'stop'})) + ($TE{$blocks{$x}{'listname'}}{$TEs[$y]}{'stop'} - $TE{$blocks{$x}{'listname'}}{$TEs[$y]}{'start'}), $TE_or , "TE" ) ."\n";
										pop @TEs;
									}
								}
							}
						}
						print OUT join( "\t", $gene->ID, (500 + ($offset - $prot{$gene->ID}{'stop'})), (500 + ($offset - $prot{$gene->ID}{'stop'})) + ($prot{$gene->ID}{'stop'} - $prot{$gene->ID}{'start'}), $or , "protein" ) ."\n";
						$prev_gene_start = $prot{$gene->ID}{'start'};
						
					} else {
						print OUT join( "\t", $gene->ID, ($n*500 + 1500*($n-1)), ($n*500 + 1500*$n), $or , "protein" ) ."\n";
						$n++;
					}
				}
				print OUT "\n";
			} else {
				my $prev_gene_stop = '';
				my $n = 1;
				foreach my $c ( $blocks{$x}{'begin'} .. $blocks{$x}{'end'} ) {
					my $gene = ${$blocks{$x}{'list_ref'}->elements}[$c]->gene;
					
					if (uc ($params{"-real"}) eq "T") {
						if (exists $params{"-TE"}) {
							if (scalar @TEs != 0 && $prev_gene_stop ne '') {
								my $aant = 0;
								for (my $y=0; $y < scalar @TEs; $y++) {
									if ( $TEs[$y] > $prev_gene_stop && $TEs[$y] < $prot{$gene->ID}{'start'} ) {
										print OUT join ( "\t", $TE{$blocks{$x}{'listname'}}{$TEs[$y]}{'name'},($TE{$blocks{$x}{'listname'}}{$TEs[$y]}{'start'} - $diff_begin),($TE{$blocks{$x}{'listname'}}{$TEs[$y]}{'stop'} - $diff_begin), $TE{$blocks{$x}{'listname'}}{$TEs[$y]}{'strand'}, "TE" ) ."\n";
										shift @TEs;
									}
								}
							}
						}
						
						print OUT join( "\t", $gene->ID, ($prot{$gene->ID}{'start'} - $diff_begin),($prot{$gene->ID}{'stop'} - $diff_begin), $gene->orientation, "protein" ) ."\n";
						$prev_gene_stop = $prot{$gene->ID}{'stop'};
					} else {
						print OUT join( "\t", $gene->ID, ($n*500 + 1500*($n-1)), ($n*500 + 1500*$n), $gene->orientation , "protein" ) ."\n";
						$n++;
					}
				}
				print OUT "\n";
			}
		}
		
		print OUT "gene relations:\n";
		print OUT "segment\tAPx-id\tsegment\tAPy-id\n";
		for (my $x = 2; $x <= $level; $x++) {
			foreach my $anchor ( $blocks{$x}{'multiplicon'}->get_anchorpoints ) { 
  				my %segment;
				$segment{'x'}=&which_segment( $multiplicon,$anchor->gene_x->list->genome,$anchor->gene_x->list->listname,$anchor->gene_x->coordinate, %segm_list );
      		$segment{'y'}=&which_segment( $multiplicon,$anchor->gene_y->list->genome,$anchor->gene_y->list->listname,$anchor->gene_y->coordinate, %segm_list );
				if ($segment{'x'} ne '' && $segment{'y'} ne '') {
					print OUT join("\t", $segment{'x'}, $anchor->gene_x->ID, $segment{'y'}, $anchor->gene_y->ID )."\n";
   			}
			}
		}
		%blocks = ();
	
		close (OUT);
	} else {
		print STDERR "File ". $output_name ." exists! will use existing one!\n";
	}
	
	my $cursegment;
	my $gene_relat;
	my %lengths;
	my %segmentnames;
	my %txtF;
	my $maxL = 1;
	my $maxTxt = 1;
	
	print STDERR ".";
	
	open (IN, "< ${output_path}profile_plot/$output_name") || die  "Couldn't open file ${output_path}multiplicon_lineair_plots_highMu/". $output_name ."\n";;
	# read textual output file with info for the image
	while (<IN>){
		my $line=$_;
		chomp($line);
		if($line =~ /^Segment\s(\d+)\t(\S+)\t(\d+)$/){
 			$cursegment = $1;
			$segmentnames{$cursegment}=$2;
   		$lengths{$1}{'length'}=$3;
			if ($3 > $maxL){
				$maxL = $3;
			}
			if (length $2 > $maxTxt){
				$maxTxt = length $2;
			}
			$line = <IN>;
		}
		if ($cursegment) {
			chomp $line;
			my @vals= split(/\t/,$line);
			if (scalar(@vals) == 5){
				$txtF{$cursegment}{$vals[0]}{'start'} = $vals[1];
				$txtF{$cursegment}{$vals[0]}{'stop'} = $vals[2];
				$txtF{$cursegment}{$vals[0]}{'str'} = $vals[3];
				$txtF{$cursegment}{$vals[0]}{'type'} = $vals[4];
			} else {
				undef $cursegment;
			}
		}
		if($line =~ /^gene\srelations:$/){
			$gene_relat=1;
			$line = <IN>;
		}
		if ($gene_relat){
			chomp $line;
			if ($line =~ m/^\d+\t\S+\t\d+\t\S+$/) {
				$txtF{'gene_relations'}{$line} = 1;	
			} elsif ($line !~ m/^\S+\t\S+\t\S+\t\S+$/) {
				undef $gene_relat;
			}
		}
	}
	close (IN);
        unlink( "${output_path}profile_plot/$output_name" ) || print STDERR "Could not delete $output_name: $!\n";

	my $margin = 10;
	my $width= ($maxL /100) + (3 * $margin) + ($maxTxt * 10);
	my $height=( scalar(keys (%txtF)) ) * 200;
	my $poly = '';
	my $line = '';
	my $text = '';

	#rescale all coordinates
	foreach my $key (keys (%txtF)) {
		if ($key ne "gene_relations"){
		
			$lengths{$key}{'length'} = $lengths{$key}{'length'} /100;
			$lengths{$key}{'offset'} = ((($maxL /100) - ($lengths{$key}{'length'}))/2) + (($maxTxt * 10) + 5);
		
			foreach my $key2 (keys (%{$txtF{$key}})) {
				$txtF{$key}{$key2}{'start'} = ($txtF{$key}{$key2}{'start'} /100) + $margin +$lengths{$key}{'offset'};
				$txtF{$key}{$key2}{'stop'} = ($txtF{$key}{$key2}{'stop'} /100) + $margin +$lengths{$key}{'offset'};
			}
		} else {
			my @keys_tmp = ( sort keys %{$txtF{'gene_relations'}});
			for ( my $k=0; $k < scalar @keys_tmp; $k++) {
				my @tmp = split (/\t/,$keys_tmp[$k]);
				my @redun = grep (/^$tmp[2]\t$tmp[3]\t/,(keys %{$txtF{'gene_relations'}})); 
				if ($redun[0]) {
					my @tmp1 = split (/\t/,$redun[0]);
					delete $txtF{'gene_relations'}{$tmp[0] ."\t". $tmp[1] ."\t". $tmp1[2] ."\t". $tmp1[3] };
				}
			}
		}
	}

	#make the svg-image
	my $svg = SVG->new(width=>$width , height=>$height, 'xmlns:xlink'=>"http://www.w3.org/1999/xlink", xmlns=>"http://www.w3.org/2000/svg", 'xmlns:ev'=>"http://www.w3.org/2001/xml-events" );

	#first draw out all gene relations so that they are in the background
	foreach my $APxAPy (sort keys %{$txtF{'gene_relations'}}){
		my @tmp = split (/\t/, $APxAPy);
			
		($svg) = &draw_relation (
			$poly, $svg,
			[ $txtF{$tmp[0]}{$tmp[1]}{'start'}, $txtF{$tmp[0]}{$tmp[1]}{'stop'}, $tmp[0], $txtF{$tmp[0]}{$tmp[1]}{'str'}, $tmp[1] ],
			[ $txtF{$tmp[2]}{$tmp[3]}{'start'}, $txtF{$tmp[2]}{$tmp[3]}{'stop'}, $tmp[2], $txtF{$tmp[2]}{$tmp[3]}{'str'}, $tmp[3] ] );		
	}

	#then draw the gene
	foreach my $key (sort (keys (%txtF))) {
		if ($key ne "gene_relations"){
			my @sorted = sort { $txtF{$key}{$a}{'start'} <=> $txtF{$key}{$b}{'start'} } keys (%{$txtF{$key}}); 
			my $y1 = ($key * 200);
		
			$line = $svg->line(id=>"segment".$key, x1=>($margin + $lengths{$key}{'offset'}), y1=>($y1-7), x2=>($lengths{$key}{'length'} +$margin +$lengths{$key}{'offset'}), y2=>($y1 -7), style=> {'stroke'=>'black','stroke-width'=>'5','stroke-linecap'=>'round'});
			$text = $svg->text(x=>(($margin + $lengths{$key}{'offset'})-(($maxTxt * 10) + 5)), y=>($y1-3), textLength=>($maxTxt * 10), style =>{'font-size'=>'10pt','font-weight'=>'bold'})->cdata($segmentnames{$key});		
			
			for (my $x=0; $x < scalar @sorted; $x++){
				($svg) = &draw_gene (
						$poly, $svg,
						[ $txtF{$key}{$sorted[$x]}{'start'},$txtF{$key}{$sorted[$x]}{'stop'} ],[ $y1,7 ],
						$sorted[$x], $key, $txtF{$key}{$sorted[$x]}{'str'}, $txtF{$key}{$sorted[$x]}{'type'} );
			
				my $X = ($txtF{$key}{$sorted[$x]}{'start'} + ($txtF{$key}{$sorted[$x]}{'stop'} - $txtF{$key}{$sorted[$x]}{'start'})/2);
				my $Y = ($y1-16);
				my $Gname = $sorted[$x];
				$Gname =~ s/$segmentnames{$key}_//;
				if ($txtF{$key}{$sorted[$x]}{'type'} eq 'protein') {
					$text = $svg->text(x=>$X, y=>$Y, transform=>"rotate(295,$X,$Y)", style => {'font-size'=>'6pt'})->cdata($Gname);		
				}
			}
		}	
	}

	my $svg_XML = $svg->xmlify();
	$output_name =~ s/\.txt$//;
	open (SVG_OUT, "> ${output_path}profile_plot/${output_name}.svg");
		print SVG_OUT $svg_XML;
	close (SVG_OUT);

}
print STDERR "\nDone!\n";


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub routines
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#select from the complete TE file only the TEs that are present on the 
#current segment
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_TEs 
	{
		my ($listname, $begin, $end, $diff) = @_;
		my @sub_TEs = ();
		
		foreach my $key (sort {$a <=> $b} (keys %{$TE{$listname}})) {
			if ($key > $begin && $key < $end) {
				push (@sub_TEs, $key);
			}
		}
		
		return (@sub_TEs);
	}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# check if the baseclusters from iADHoRe were twisted, if so return $flp
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub flip
	{
		my @baseclusters = @_;
		my $min = 0;
		my $plus = 0;
		my $flp = '';
		foreach my $baseclust (@baseclusters) {
			if ($baseclust->orientation eq '-') {
				$min += $baseclust->number_of_anchorpoints;
			} else {
				$plus += $baseclust->number_of_anchorpoints;
			}
		}
		if ($plus >= $min) {
			$flp=0;
 		} else {
			$flp=1;
		#	print STDERR "* FLIPPING Y segment\n";
		}
	
		return ($flp);
	}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#for each anchorpoint check on which segments the gene is present and then 
#return the name of that segment
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub which_segment
 {
  my ($mu,$genome,$list,$coord,%segm_lst)=@_;
  
  foreach my $segment ( @{$mu->profile->segments} ){
    if ( ($segment->genome eq $genome) &&
	 		($segment->listname eq $list) &&
	 		($segment->first->gene->coordinate <= $coord) &&
	 		($segment->last->gene->coordinate >= $coord)){
      	
			return( $segm_lst{$segment->id} );
     
	  } 
   } 
  
 }
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#first draw the relations between the anchorpoints
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub draw_relation
{
	my ($p,$s,$APx,$APy) = @_;
	
	if ( ${$APx}[3] eq ${$APy}[3] ){
		my $xv;
		if ( ${$APx}[3] eq '+'){
			$xv = [ ${$APx}[0], ${$APx}[1]-2, ${$APy}[1]-2, ${$APy}[0]];
		} else {
			$xv = [ ${$APx}[0]+2, ${$APx}[1], ${$APy}[1], ${$APy}[0]+2];
		} 
			
		my $yv = [ (${$APx}[2] *200), (${$APx}[2] *200), (${$APy}[2] *200)-14, (${$APy}[2] *200)-14 ];
		my $points = $s->get_path(x=>$xv, y=>$yv, -type=>'polygon');
		
		$p = $s->polygon(id=>${$APx}[4] ."_". ${$APy}[4], %$points, style=>{'fill'=>'gray' , 'fill-opacity'=>0.3,'stroke'=>'grey'});
		
	} else {
		my $xv;
		if ( ${$APx}[3] eq '+'){
			$xv = [ ${$APx}[0], ${$APx}[1]-2, ${$APy}[0]+2, ${$APy}[1]];
		} else {
			$xv = [ ${$APx}[0]+2, ${$APx}[1], ${$APy}[0], ${$APy}[1]-2];
		} 
			
		my $yv = [ (${$APx}[2] *200), (${$APx}[2] *200), (${$APy}[2] *200)-14, (${$APy}[2] *200)-14 ];
		my $points = $s->get_path(x=>$xv, y=>$yv, -type=>'polygon');
		
		$p = $s->polygon(id=>${$APx}[4] ."_". ${$APy}[4], %$points, style=>{'fill'=>'grey' , 'fill-opacity'=>0.3,'stroke'=>'grey'});
		
	}
	
	return ($s);
	
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#draw each gene (and transposon) on the segment
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub draw_gene
{
	my ($p,$s,$xcoords,$ycoords,$name,$k,$str,$type) = @_;
	
	if ( (${$xcoords}[1] - ${$xcoords}[0]) >3 ){
	
		if ($str eq '-') {
			my $xv = [ ${$xcoords}[0],${$xcoords}[0]+2,${$xcoords}[1],${$xcoords}[1],${$xcoords}[0]+2 ];
			my $yv = [ ${$ycoords}[0]-${$ycoords}[1],${$ycoords}[0],${$ycoords}[0],${$ycoords}[0]-(2* ${$ycoords}[1]),${$ycoords}[0]-(2* ${$ycoords}[1]) ];
			my $points = $s->get_path(x=>$xv, y=>$yv, -type=>'polygon');
		
			if ($type eq 'protein'){
				$p = $s->polygon(id=>$k ."_". $name, %$points, style=>{'fill'=>'olive','stroke'=>'black'});
			} elsif ($type eq 'TE'){
				$p = $s->polygon(id=>$k ."_". $name, %$points, style=>{'fill'=>'red','stroke'=>'darkred'});
			} else {
				$p = $s->polygon(id=>$k ."_". $name, %$points, style=>{'fill'=>'black','stroke'=>'darkred'});
			}

		} else {
			my $xv = [ ${$xcoords}[0],${$xcoords}[1]-2,${$xcoords}[1],${$xcoords}[1]-2,${$xcoords}[0] ];
			my $yv = [ ${$ycoords}[0],${$ycoords}[0],${$ycoords}[0]-${$ycoords}[1],${$ycoords}[0]-(2* ${$ycoords}[1]),${$ycoords}[0]-(2* ${$ycoords}[1]) ];
			my $points = $s->get_path(x=>$xv, y=>$yv, -type=>'polygon');
		
			if ($type eq 'protein'){
				$p = $s->polygon(id=>$k ."_". $name, %$points, style=>{'fill'=>'olive','stroke'=>'black'});
			} elsif ($type eq 'TE'){
				$p = $s->polygon(id=>$k ."_". $name, %$points, style=>{'fill'=>'red','stroke'=>'darkred'});
			} else {
				$p = $s->polygon(id=>$k ."_". $name, %$points, style=>{'fill'=>'black','stroke'=>'darkred'});
			}
		}
	} else {
		if ($str eq '-') {
			my $xv = [ ${$xcoords}[0],${$xcoords}[1],${$xcoords}[1] ];
			my $yv = [ ${$ycoords}[0]-${$ycoords}[1],${$ycoords}[0],${$ycoords}[0]-(2* ${$ycoords}[1]) ];
			my $points = $s->get_path(x=>$xv, y=>$yv, -type=>'polygon');
		
			if ($type eq 'protein'){
				$p = $s->polygon(id=>$k ."_". $name, %$points, style=>{'fill'=>'olive','stroke'=>'black'});
			} elsif ($type eq 'TE'){
				$p = $s->polygon(id=>$k ."_". $name, %$points, style=>{'fill'=>'red','stroke'=>'darkred'});
			} else {
				$p = $s->polygon(id=>$k ."_". $name, %$points, style=>{'fill'=>'black','stroke'=>'darkred'});
			}

		} else {
			my $xv = [ ${$xcoords}[0],${$xcoords}[1],${$xcoords}[0] ];
			my $yv = [ ${$ycoords}[0],${$ycoords}[0]-${$ycoords}[1],${$ycoords}[0]-(2* ${$ycoords}[1]) ];
			my $points = $s->get_path(x=>$xv, y=>$yv, -type=>'polygon');
		
			if ($type eq 'protein'){
				$p = $s->polygon(id=>$k ."_". $name, %$points, style=>{'fill'=>'olive','stroke'=>'black'});
			} elsif ($type eq 'TE'){
				$p = $s->polygon(id=>$k ."_". $name, %$points, style=>{'fill'=>'red','stroke'=>'darkred'});
			} else {
				$p = $s->polygon(id=>$k ."_". $name, %$points, style=>{'fill'=>'black','stroke'=>'darkred'});
			}
		}
		
	}
	
	return ($s);
	
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
