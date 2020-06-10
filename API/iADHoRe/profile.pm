package iADHoRe::profile;
#Created by Cedric Simillion on Tue Dec 28 13:46:07 CET 2004

use strict;
use iADHoRe::settings;
use iADHoRe::record_structure;
use iADHoRe::nw_path;
use iADHoRe::ghm;

our @ISA=qw(iADHoRe::record_structure);
our $settings=iADHoRe::settings->get_settings;

our $structure=iADHoRe::profile->create_structure(qw(multiplicon valid segments unaligned_x_object unaligned_y_list));

#===============================================================================
# Package methods
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub create
# Creates a new profile from a multiplicon
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub create
 {
  my ($pkg, $multiplicon,$x_object,$y_list,$profile,$segments);
  ($pkg, $multiplicon)=@_;  

  $x_object=$multiplicon->x_object;
  $y_list=$multiplicon->y_list;

  $profile=$pkg->new( $multiplicon, -1, [] );
  $profile->unaligned_x_object( $x_object->extract_segment( $multiplicon->begin_x, $multiplicon->end_x, $profile ) );
  $profile->unaligned_y_list( $y_list->extract_segment( $multiplicon->begin_y, $multiplicon->end_y, $profile ));
  $profile->apply_y_list_permutation;  
  
  $profile->check_masking;
  
  $profile->valid && $profile->align;
  $profile->valid && $profile->check_alignment;

  return($profile); 
 } #sub create

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub align
# Aligns the segments in a profile
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub align
 {
  my ($profile,$path,$first,$last,$segments);
  $profile=$_[0];
  $path=iADHoRe::nw_path->get_alignment_path($profile->unaligned_x_object, $profile->unaligned_y_list);
  
  if ( !$path->valid )
   {
    $profile->valid( 0 );
    print STDERR " unalignable multiplicon!\n";
    return;
   } #if ( !$path->valid )
  
  $first=0;
  $last=0;

  while ( my ($x,$y) = $path->next_offset )
   {

    if ($x == $y)
     {
      $first+=$x;
      $last+=$x;
     } #if ($x == $y)
    else 
     {
      my ($insert_object,$number_of_gaps);

      if ( $x > $y )
       {
        $insert_object=$profile->unaligned_y_list; #Introduce gaps in the y_list
	$number_of_gaps=$x-$y;
	$last+=$y;
       } #if ( $x > $y )
      else # $y > $x
       {
        $insert_object=$profile->unaligned_x_object; #Introduce gaps in the x_object;
	$last+=$x;
	$number_of_gaps=$y-$x;
       } #else

      $insert_object->introduce_gaps($first,$last,$number_of_gaps); 

      $last+=$number_of_gaps;
      $first=$last;  
     } #else 

   } #while ( my ($x,$y) = $path->next_offset )
  
  #Put the now aligned x_object and y_list together in the segments array
  $segments=$profile->segments;
  #For multiplicons with level 2, x_object is a single list, for higher levels it is an array of lists.
  ( $profile->multiplicon->level > 2 )
   ? push( @$segments, @{$profile->unaligned_x_object->segments} )
   : push( @$segments, $profile->unaligned_x_object );
  push( @$segments, $profile->unaligned_y_list );
  
 } #sub align

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub apply_y_list_permutation
# Permutates the unaligned_y_list object of a profile according to the
# orientation and was_twisted attribute of each basecluster of its
# multiplicon
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub apply_y_list_permutation
 {
  my ($profile,$offset,$baseclusters,$multiplicon,$y_list_inverted);
  my (@twisted_clusters);

  $profile=$_[0];
  $multiplicon=$profile->multiplicon;
  $offset=$multiplicon->begin_y;
  $baseclusters=$multiplicon->baseclusters;

  #First, check if the entire y_list needs to be reversed.
  if (
       ( ( $$baseclusters[0]->orientation eq '-' ) && ( !$$baseclusters[0]->was_twisted ) )
      ||
       ( ( $$baseclusters[0]->orientation eq '+' ) && ( $$baseclusters[0]->was_twisted ) )
     )
   {
    $profile->unaligned_y_list->invert_section( 0, $multiplicon->end_y - $offset );
    $y_list_inverted= -1;
   } #if ( ... )
  else 
   {
    $y_list_inverted= 0;
   } #else 
  
  #Make sure that the outermost baseclusters that are twisted first.
  foreach my $cluster ( @$baseclusters) 
   {
    $cluster->was_twisted || next;
    ( ( $cluster->begin_y == $offset ) || ( $cluster->end_y == $multiplicon->end_y ) )
     ? unshift( @twisted_clusters, $cluster )
     : push( @twisted_clusters, $cluster ); 
   } #foreach my $cluster ( grep { $_->was_twisted } (@$baseclusters) )
   
  #Next, reverse these sections that correspond to the twisted baseclusters
  foreach my $twisted_cluster ( @twisted_clusters )
   {
    my ($begin,$end);
    $begin=$twisted_cluster->begin_y - $offset;
    $end=$twisted_cluster->end_y - $offset;
    
    #If the entire y_list was inverted, the y-coordinates of the twisted cluster need to be mirrored first.
    if ( $y_list_inverted )
     {
      my $highest_y_coordinate = $profile->unaligned_y_list->size - 1;
      $begin = $highest_y_coordinate - $begin;
      $end = $highest_y_coordinate - $end;
      ($begin, $end) = ($end, $begin);
     } #if ( $y_list_inverted )

    $profile->unaligned_y_list->invert_section( $begin, $end );

   } #foreach my $twisted_cluster ( grep { $_->was_twisted } (@$baseclusters) )

 } #sub apply_y_list_permutation

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub size
# Returns the length of the profile
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub size
 {
  my ($profile,$segments);
  $profile=$_[0];
  $segments=$profile->segments;
  return( $$segments[0]->size );
 } #sub size

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub matrix
# Returns an array of the remapped lists of each segment in the
# profile.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub matrix
 {
  my ($profile,$segments);
  my (@matrix);
  $profile=$_[0];
  $segments=$profile->segments;
  
  foreach my $segment ( @$segments )
   {
    push( @matrix, $segment->remapped );
   } #foreach ( @$segments )
  
  return( \@matrix );
  
 } #sub matrix

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub extract_segment
# Creates a profile object of which all segments are subsegments of
# the calling profile object. Calls the
# iADHoRe::genelist::extract_segment mothod to achieve this.
# Arguments:
# - first position of the subsegment to be extracted
# - last position of the subsegment to be extracted 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub extract_segment
 {
  my ($profile,$subprofile,$profile_segments,$subprofile_segments,$begin,$end);
  ($profile,$begin,$end)=@_;
  $subprofile=iADHoRe::profile->new($profile->multiplicon, -1, []);
  
  $profile_segments=$profile->segments;
  $subprofile_segments=$subprofile->segments;
  
  foreach my $segment ( @$profile_segments )
   {
    push( @$subprofile_segments, $segment->extract_segment( $begin, $end, $profile ) );
   } #foreach my $segment ( @$profile_segments )
  
  return( $subprofile ); 
 } #sub extract_segment

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub introduce_gaps
# Introduces between two positions of all segments in a profile
# a number of gaps. Calls the iADHoRe::genelist::introduce_gaps method
# to achieve this.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub introduce_gaps
 {
  my ($profile, $segments,$begin,$end,$number_of_gaps);
  ($profile,$begin,$end,$number_of_gaps)=@_;
  $segments=$profile->segments;

  foreach my $segment ( @$segments )
   {
    $segment->introduce_gaps( $begin,$end,$number_of_gaps );
   } #foreach my $segment ( @$segments )

 } #sub introduce_gaps

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub check_masking
# Checks if the last segment added wasn't completely masked. If so,
# the profile's valid flag will be set to false. If not, the entire
# segment will be completely unmasked.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub check_masking
 {
  my ($profile,$elements,$all_masked);
  $profile=$_[0];
  $elements=$profile->unaligned_y_list->remapped;

  unless ( $settings->level_2_only )
   {
    $all_masked=-1;
    foreach my $element ( @$elements )
     {
      if (!$element->masked)
       {
	$all_masked=0;
	last;
       } #if (!$element->masked)
     } #foreach my $element ( @$elements )    
   } #unless ( $settings->level_2_only )
  else 
   {
    $all_masked=0;
   } #else 
   
  if ( $all_masked ) 
   {
    $profile->valid( 0 );
    print STDERR " masking check failed.\n";
   } #if ( $all_masked )
  else  
   {
    foreach my $element ( @$elements )
     {
      $element->masked( 0 );
     } #foreach my $element ( @$elements )
   } #else  

 } #sub check_masking

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub check_alignment
# Checks the alignment of a profile. This means checking if each
# segment has at least a minimum number of homologous genes on other
# segments.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub check_alignment
 {
  my ($profile,$segments,$max_nr_of_pairs);
  my (@in,@most_connected_pair);
  my (%pair_identities,%out,%multi_identities);
  
  $profile=$_[0];
  $segments=$profile->segments;
  
  $max_nr_of_pairs=0;

  foreach my $i ( 0 .. $#{$segments} -1 )
   {
    foreach my $j ( $i+1 .. $#{$segments} )
     {
      $pair_identities{$i}{$j} = iADHoRe::ghm->new( $$segments[$i], $$segments[$j] );
      $pair_identities{$j}{$i} = $pair_identities{$i}{$j};

      if ( $pair_identities{$i}{$j}->number_of_points > $max_nr_of_pairs )
       {
        @most_connected_pair=($i,$j);
	$max_nr_of_pairs=$pair_identities{$i}{$j}->number_of_points;
       } #if ( $pair_identities{$i}{$j} > $max_nr_of_pairs )
     } #foreach my $j ( 0 .. $#{$segments} )
   } #foreach my $i ( 0 .. $#{$segments} )

  if ( $max_nr_of_pairs < $settings->anchorpoints )
   {
    $profile->valid( 0 );
    print STDERR " alignment check failed.\n";
    return;
   } #if ( $max_nr_of_pairs < $settings->anchorpoints )
  
  @in=@most_connected_pair;
  foreach ( 0 .. $#{$segments} )
   {
    ( ($_ == $in[0]) || ($_ == $in[1]) ) || ($out{$_}=-1);
   } #foreach ( 0 .. $#{$segments} )

  foreach my $out_segment (keys %out)
   {
    $multi_identities{$out_segment}
     = iADHoRe::ghm->merge( $pair_identities{$in[0]}{$out_segment}, $pair_identities{$in[1]}{$out_segment} );
   } #foreach my $out_segment (keys %out)
  
  until( scalar( @in ) == scalar( @$segments ) )
   {
    my ($most_connected_segment);
    $max_nr_of_pairs=0;

    foreach my $out_segment (keys %out)
     {
      if ( $multi_identities{$out_segment}->number_of_y_hits > $max_nr_of_pairs ) 
       {
        $most_connected_segment = $out_segment;
	$max_nr_of_pairs = $multi_identities{$out_segment}->number_of_y_hits;
       } #if ($nr_of_pairs > $max_nr_of_pairs) 
     } #foreach my $out_segment (keys %out)

    if ( $max_nr_of_pairs >= $settings->anchorpoints ) 
     {
      push( @in, $most_connected_segment );
      delete $out{$most_connected_segment};
      
      foreach my $out_segment (keys %out)
       {
        $multi_identities{$out_segment}
	 = iADHoRe::ghm->merge( $multi_identities{$out_segment}, $pair_identities{$most_connected_segment}{$out_segment} );
       } #foreach my $out_segment (keys %out)
      
     } #if ( $max_nr_of_pairs >= $settings->anchorpoints )
    else 
     {
      $profile->valid( 0 );
      print STDERR " alignment check failed.\n";
      return;
     } #else 

   } #until( scalar( @in ) == scalar( @$segments ) )
  
 } #sub check_alignment

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub plot_profile
# Returns a GD image object containing an alignplot of the profile
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub plot_profile
 {
  my ($profile,$family_matrix,$family_colors,$image,$size,$sp,$segments,$img,$file);
  ($profile, $size) = @_;
  
  $file = defined( $_[2] ) ? $_[2] : "profile.png";
  
  ($family_matrix,$family_colors) = $profile->matrix_single_linkage;
  $img = $profile->init_image( $family_colors, $size );
  $sp = int($size/20);
  $segments = $profile->segments;
  
  foreach my $s ( 0 .. $#{$segments} )
   {
    foreach my $i ( 0 .. $#{$$segments[$s]->remapped} )
     {
      my ($id, $color);
      defined( ${$$segments[$s]->remapped}[$i] ) || next;
      $id = defined( $$family_matrix[$s][$i] ) ? $$family_matrix[$s][$i] : "white";
      $color = $$family_colors{$id};
      $img->filledRectangle($size+$i*$size+$sp, $s*$size+$sp, $size+$i*$size+($size-1)-$sp, $s*$size+($size-1)-$sp, $color);
      ( $id eq 'white' ) && $img->rectangle( $size+$i*$size+$sp, $s*$size+$sp, $size+$i*$size+($size-1)-$sp,
                                             $s*$size+($size-1)-$sp, $$family_colors{"black"} );
     } #foreach my $i ( 0 .. $#{$$segments[$s]->remapped} )
   } #foreach my $s ( 0 .. $#{$segments} )
  
  open ( OUT, ">$file" ) || die "Cannot write file $file: $!\n";
  binmode OUT;
  print OUT $img->png;
  close OUT;

  return( "profile plot saved in $file at".localtime );
  
 } #sub plot_profile

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub matrix_single_linkage
# Performs single linkage clustering on all genes in a profile.
# Returns:
# - a reference to a 2D array where each cells contains a scalar
#   reference to the family ID of the gene on the corresponding
#   position in the matrix
# - a reference to a hash where the keys are the IDs of the families
#   created. The values are all true.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub matrix_single_linkage
 {
  my ($profile,$matrix,$family_count,$size);
  my @family_matrix;
  my (%fam_colors,%family_of_gene,%family);
  
  $profile = $_[0];
  $matrix=$profile->matrix;
  $size = $profile->size;
  $family_count = 0;

  #Perform single linkage clustering on all genes in the profile. Each family will then be given a different color. 
  foreach my $i ( 0 .. $#{$matrix} - 1 )
   {
    foreach my $j ( $i + 1 .. $#{$matrix} )
     {
      foreach my $i_index ( 0 .. $size - 1 )
       {
        defined( $$matrix[$i][$i_index] ) || next;
        foreach my $j_index ( 0 .. $size - 1 )
	 {
	  defined( $$matrix[$j][$j_index] ) || next;
	  $$matrix[$i][$i_index]->is_pair_with( $$matrix[$j][$j_index] ) || next;
	  my ($gene_i, $gene_j);
	  $gene_i = $$matrix[$i][$i_index]->gene->ID;
	  $gene_j = $$matrix[$j][$j_index]->gene->ID;

          if ( !exists( $family_of_gene{$gene_i} ) && !exists( $family_of_gene{$gene_j} ) )
	   {
	    my $family_id = $family_count;
	    $family_count++;
	    $family_of_gene{$gene_i} = $family_id;
	    $family_of_gene{$gene_j} = $family_id;
	    $family{$family_id} = [ $gene_i, $gene_j ];
	   } #if ( !exists( $family_of_gene{$gene_i} ) && !exists( $family_of_gene{$gene_j} ) )
          elsif ( exists( $family_of_gene{$gene_i} ) && !exists( $family_of_gene{$gene_j} ) )
	   {
	    $family_of_gene{$gene_j} = $family_of_gene{$gene_i};
	    push( @{$family{$family_of_gene{$gene_i}}}, $gene_j );
	   } #elsif ( exists( $family_of_gene{$gene_i} ) && !exists( $family_of_gene{$gene_j} ) )
          elsif ( !exists( $family_of_gene{$gene_i} ) && exists( $family_of_gene{$gene_j} ) )
	   {
	    $family_of_gene{$gene_i} = $family_of_gene{$gene_j};
	    push( @{$family{$family_of_gene{$gene_j}}}, $gene_i );
	   } #elsif ( !exists( $family_of_gene{$gene_i} ) && exists( $family_of_gene{$gene_j} ) )
	  elsif ( $family_of_gene{$gene_i} != $family_of_gene{$gene_j} ) 
	   {
	    my $family_i = $family_of_gene{$gene_i};
	    my $family_j = $family_of_gene{$gene_j};
	    foreach my $gene ( @{$family{$family_j}} )
	     {
	      $family_of_gene{$gene} = $family_i;
	     } #foreach my $gene ( @{$family{$family_j}} )
	    push( @{$family{$family_i}}, @{$family{$family_j}} );
	    delete $family{$family_j}; 
	   } #elsif ( $family_of_gene{$gene_i} != $family_of_gene{$gene_j} ) 
	 } #foreach my $j_index ( 0 .. $size - 1 )
       } #foreach my $i_index ( 0 .. $size - 1 )
     } #foreach my $j ( $i + 1 .. $#{$matrix} )
   } #foreach my $i ( 0 .. $#{$matrix} )
  
  foreach my $i ( 0 .. $#{$matrix}  )
   {
    foreach my $i_index ( 0 .. $size - 1 )
     {
      defined( $$matrix[$i][$i_index] ) || next;
      my $gene = $$matrix[$i][$i_index]->gene->ID;
      if ( exists $family_of_gene{$gene} )
       {
        $family_matrix[$i][$i_index] = $family_of_gene{$gene};
	$fam_colors{$family_of_gene{$gene}} = -1;
       } #if ( exists $family_of_gene{$gene} )
     } #foreach my $i_index ( 0 .. $size - 1 )
   } #foreach my $i ( 0 .. $#{$matrix} )
  
  return( \@family_matrix, \%fam_colors );
  
 } #sub matrix_single_linkage

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub init_image
# Takes as input the hash output by matrix_single_linkage as well
# as the size of each cell in the plot to initialize a GD object
# to draw an alignment plot of a profile.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub init_image
 {
  use POSIX;
  use GD;
  my ($img,$number_of_segments,$longest_size,$fam_colors,$number_of_colors,$white,$black,$step,$size,$sp,$profile);
  my (@rgb)=(255,255,255);
  ($profile,$fam_colors,$size)=@_;
  $longest_size = $profile->size;
  $number_of_segments = scalar( @{$profile->segments} );
  $number_of_colors=scalar (keys %$fam_colors);
  $sp = int( $size/ 20 );
  $step=int( 255 / ( ceil( ($number_of_colors+1)**(1/3) ) - 1 ) );
  $img=new GD::Image($size*($longest_size+1),$size*$number_of_segments);
  $white=$img->colorAllocate(255,255,255);
  $black=$img->colorAllocate(0,0,0);
  foreach my $fam (keys %$fam_colors)
   {
    my $i=0;
    $rgb[$i]-=$step;
    until ($rgb[$i]>=0)
     {
      $rgb[$i]=255;
      $rgb[$i+1]-=$step;
      $i++;
     } #until ($rgb[$i]>=0)
    $$fam_colors{$fam}=$img->colorAllocate(@rgb); 
   } #foreach my $fam (keys %$fam_colors)
  $$fam_colors{"white"}=$white;
  
  foreach (0..$number_of_segments-1)
   {
    $img->filledRectangle(0,(int($size/2)-$sp)+$size*$_,$size*($longest_size+1),(int($size/2)+$sp)+$size*$_,$black);
   } #foreach (1..$number_of_segments)
  $$fam_colors{"black"}=$black;
  return ($img);

 } #sub init_image

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_segments
# Returns a list of all the segments of a profile
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_segments
 {
  return ( @{$_[0]->segments} );
 } #sub get_segments

#===============================================================================
# Accessor methods
#===============================================================================
 
sub multiplicon
 {
  my ($profile,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'multiplicon'}[ $$profile ] = $set_value );
  return ( $$structure{'multiplicon'}[ $$profile ] );
 } #sub multiplicon
 
#------------------------------------------------------------
 
sub valid
 {
  my ($profile,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'valid'}[ $$profile ] = $set_value );
  return ( $$structure{'valid'}[ $$profile ] );
 } #sub valid
 
#------------------------------------------------------------
 
sub segments
 {
  my ($profile,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'segments'}[ $$profile ] = $set_value );
  return ( $$structure{'segments'}[ $$profile ] );
 } #sub segments
 
#------------------------------------------------------------
 
sub unaligned_x_object
 {
  my ($profile,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'unaligned_x_object'}[ $$profile ] = $set_value );
  return ( $$structure{'unaligned_x_object'}[ $$profile ] );
 } #sub unaligned_x_object
 
#------------------------------------------------------------
 
sub unaligned_y_list
 {
  my ($profile,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'unaligned_y_list'}[ $$profile ] = $set_value );
  return ( $$structure{'unaligned_y_list'}[ $$profile ] );
 } #sub unaligned_y_list
 
1;
