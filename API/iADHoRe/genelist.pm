package iADHoRe::genelist;
#Created by Cedric Simillion on Fri Jun 25 16:35:47 CEST 2004

use strict;
use iADHoRe::gene;
use iADHoRe::list_element;
use iADHoRe::settings;
use iADHoRe::record_structure;
use iADHoRe::ghm;

our @ISA=qw(iADHoRe::record_structure);

our $structure=iADHoRe::genelist->create_structure(qw(listname genome elements remapped profile id));

my $settings=iADHoRe::settings->get_settings;

#===============================================================================
# Package methods
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub create_from_listfile
# Creates a genelist object from a listfile object by reading in the
# listfile.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub create_from_listfile
 {
  my ($pkg,$genelist,$listfile,$coord,$elements);
  ($pkg,$listfile)=@_;
  
  $genelist=$pkg->new($listfile->listname,$listfile->genome,[],[]);
  
  $elements=$genelist->elements;
  
  $coord=0;
  open (IN, $listfile->file ) || die("Could not open gene list file ".$listfile->file.": $!\n");
  while (<IN>)
   {
    my ($gene);
    chomp;
    s/([+-])$//;
    $gene=iADHoRe::gene->new($_,$listfile->genome,$genelist,$coord,$1,undef,0,0,undef,0,{});
    push @$elements,iADHoRe::list_element->new($gene,$1,0);
    $coord++;
   } #while (<IN>)
  close IN;
  
  return $genelist;
 } #sub create_from_listfile

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub extract_segment
# Creates a new gene list by extracting a segment from another gene
# list. The original gene list object, start and end coordinates are
# passed as arguments.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub extract_segment
 {
  my ($parent_list,$begin,$end,$new_segment,$new_elements,$parent_elements,$profile);
  ($parent_list,$begin,$end,$profile)=@_;
  
  $parent_elements=$parent_list->remapped;
  $new_segment=iADHoRe::genelist->new( $parent_list->listname, $parent_list->genome, [], [], $profile );
  $new_elements=$new_segment->remapped;

  foreach my $i ( $begin .. $end )
   {
    push( @$new_elements,
     ( ref($$parent_elements[$i]) eq 'iADHoRe::list_element' ) ? $$parent_elements[$i]->duplicate : undef );
   } #foreach my $i ( $begin .. $end )
  
  return( $new_segment );
 } #sub extract_segment

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub elementlist
# Returns a list of all elements (gene objects) of the genelist
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub elementlist
 {
  my $genelist=$_[0];
  my $element_list=$genelist->elements;
  return( @{$element_list} )
 } #sub elementlist

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub remapped_list
# Returns the remapped list of gene objects of the genelist
# You need to do tandem remapping first before using this method
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub remapped_list
 {
  my $genelist=$_[0];
  my $remapped_list=$genelist->remapped;
  return( @{$remapped_list} )
 } #sub remapped_list

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub remap_tandems
# Performs remapping of tandem repeated genes in a gene list object.
# Genes in a tandem array are remapped onto the first gene of that
# array which is marked as tandem_representative.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub remap_tandems
 {
  my ($genelist,$element_list,$gap,$substract,$remapped_list);
  $genelist=$_[0];
  $element_list=$genelist->elements;
  $remapped_list=$genelist->remapped;
  
  $gap=int( $settings->gap_size / 2 );
  
  #Detect tandem duplicated genes and remap then onto their representatives
  for ( my $i=scalar( @{$element_list} )-1; $i > 0 ; $i-- )
   {
    $$element_list[$i]->has_pairs || next;
    for ( my $j=$i-1; ( $j >= $i-($gap+1) ) and ( $j >= 0 ); $j-- )
     {
      $$element_list[$j]->has_pairs || next;
      if ( $$element_list[$i]->is_indirect_pair_with( $$element_list[$j] ) ) 
       {
        $$element_list[$i]->remap_to( $$element_list[$j] );
	last;
       } #if ( $$element_list[$i]->is_indirect_pair_with( $$element_list[$j] ) ) 
     } #for ( my $j=$i-1; ( $j>=$i-($gap+1) ) and ( $j >= 0 ); $j-- )
   } #for ( my $i=scalar( @{$element_list} )-1; $i > 0 ; $i-- )
  
  #Calculate remapped coordinate for each gene
  $substract=0;
  foreach my $element ( @{$element_list} )
   {
    $element->gene->remapped && ($substract++);
    $element->gene->remapped || push (@$remapped_list, $element); #The remapped list only contains the genes that are not remapped.
    my $coordinate=$element->gene->coordinate;
    $element->gene->remapped_coordinate( $coordinate-$substract );
   } #foreach my $gene ( @{$element_list} )
   
 } #sub remap_tandems

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub invert_section
# Inverses the order of a section of a gene list. First and last
# coordinate of this section are passed as arguments.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub invert_section
 {
  my ($segment,$begin,$end,$elements);
  ($segment,$begin,$end)=@_;
  $elements=$segment->remapped;

  foreach my $element ( @$elements[ $begin .. $end ] ) 
   {
    $element->invert_orientation
   } #foreach my $element ( @$elements[ $begin .. $end ] ) 

  @$elements[ $begin .. $end ] = reverse( @$elements[ $begin .. $end ] );
 } #sub invert_section

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub introduce_gaps
# Introduces between two positions in a segment a number of gaps. The
# list elements in between these positions will be evenly dispersed.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub introduce_gaps
 {
  my ($segment,$begin,$end,$number_of_gaps,$elements,$number_of_elements,$avg_gaps,$rest_gaps,$gaps_inserted);
  ($segment,$begin,$end,$number_of_gaps)=@_;
  
  $elements=$segment->remapped;
  
  $number_of_elements=($end-$begin)-1;
  $gaps_inserted=0;
  
  if ( $number_of_elements > 0 )
   {
    $avg_gaps=int( $number_of_gaps / $number_of_elements ); #The average number of gaps to be inserted per position
    $rest_gaps= $number_of_gaps % $number_of_elements; #The remaining number of gaps to be inserted
   } #if ( $number_of_elements > 0 )
  else 
   {
    $rest_gaps=$number_of_gaps;
    $avg_gaps=0;
   } #else 
  
  for (my $n=0; ($n < $number_of_elements) && ($avg_gaps > 0); $n++)
   {
    my ($position,$element);
    my (@substitute);

    $position=($begin+1)+($n * $avg_gaps) + $n;
    $element=$$elements[$position];
    
    #Substitute the element at position $position with a list consisting of $avg_gaps gaps and the element
    #itself embedded in the middle
    $#substitute=$avg_gaps-1; #Creates a list of $avg_gaps gaps (undef values) in @substitute
    splice( @substitute, int( $avg_gaps / 2 ), 0, $element ); #Insert the element in the middle of @substitute
    splice( @$elements,  $position, 1, @substitute );
    $gaps_inserted+=$avg_gaps
   } #for (my $n=0; ($n < $number_of_elements) && ($avg_gaps > 0); $n++)
  
  #If any gaps left: first introduce one at the first position
  if ($rest_gaps > 0 )
   {
    splice( @$elements, $begin+1 , 0, undef );
    $rest_gaps--;
    $gaps_inserted++;
   } #if ($rest_gaps > 0 )
  #Then at the last position 
  if ($rest_gaps > 0 )
   {
    splice( @$elements, $end + $gaps_inserted , 0, undef );
    $rest_gaps--;
    $gaps_inserted++;
   } #if ($rest_gaps > 0 )
  #Divide remaining gaps evenly across the section  
  for (my $n=0; $n < $rest_gaps; $n++)   
   {
    splice( @$elements, $begin+1 + ($n * ($number_of_elements + $gaps_inserted)/$rest_gaps) + $n , 0, undef );
   } #for (my $n=0; $n < $rest_gaps; $n++)   

 } #sub introduce_gaps

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub size
# Returns the number of elements (both gaps and list elements) of the
# *remapped* list.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub size
 {
  my ($list,$remapped);
  ($list)=$_[0];
  $remapped=$list->remapped;
  return( scalar( @$remapped ) );
 } #sub size

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub unmapped_size
# Returns the number of elements of the *unmapped* list. Does not
# apply to profile segments.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub unmapped_size
 {
  my ($list,$elements);
  ($list)=$_[0];
  $elements=$list->elements;
  return( scalar( @$elements ) );
 } #sub unmapped_size

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub mask
# Masks a section of the genelist. Arguments:
# - begin: first position to be masked
# - end: last position to be masked
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub mask
 {
  my ($list,$remapped,$begin,$end);
  ($list,$begin,$end)=@_;
  $remapped=$list->remapped;
  
  foreach my $i ($begin .. $end)
   {
    $$remapped[$i]->masked( -1 );
   } #foreach my $i ($begin .. $end)

 } #sub mask

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub number_of_masked_elements
# Returns the number of masked elements in a genelist
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub number_of_masked_elements
 {
  my ($list,$number_of_masked_elements);
  $list=$_[0];
  $number_of_masked_elements=0;

  foreach my $element ( @{$list->remapped} )
   {
    $element->masked && $number_of_masked_elements++;
   } #foreach my $element ( @{$ghm->y_list->remapped} )

  return( $number_of_masked_elements );

 } #sub number_of_masked_elements

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub first
# Returns the list_element that has the lowest coordinate on the 
# original gene lists.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub first
 {
  my ($genelist,$first,$lowest_coordinate);
  $genelist=$_[0];
  
  foreach my $element ( @{$genelist->remapped} )
   {
    if ( defined( $element ) )
     {
      $first = $element;
      $lowest_coordinate = $element->gene->remapped_coordinate;
      last;
     } #if ( defined( $element ) )
   } #foreach my $element ( @{$genelist->remapped} )
  
  foreach my $element ( @{$genelist->remapped} )
   {
    defined( $element ) || next;
    if ( $element->gene->remapped_coordinate < $lowest_coordinate )
     {
      $first = $element;
      $lowest_coordinate = $element->gene->remapped_coordinate;
     } #if ( $element->gene->remapped_coordinate < $lowest_coordinate )
   } #foreach my $element ( @{$genelist->remapped} )
  
  return( $first );
  
 } #sub first

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub last
# Returns the list_element that has the highest coordinate on the 
# original gene lists.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub last
 {
  my ($genelist,$last,$highest_coordinate);
  $genelist=$_[0];
  
  $highest_coordinate = -1;
  
  foreach my $element ( @{$genelist->remapped} )
   {
    defined( $element ) || next;
    if ( $element->gene->remapped_coordinate > $highest_coordinate )
     {
      $last = $element;
      $highest_coordinate = $element->gene->remapped_coordinate;
     } #if ( $element->gene->remapped_coordinate > $highest_coordinate )
   } #foreach my $element ( @{$genelist->remapped} )
  
  return( $last );
  
 } #sub last

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub unmapped_coordinates
# Returns the unmapped begin and end coordinates of a profile segment.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub unmapped_coordinates
 {
  my $genelist=$_[0];
  
  return( $genelist->first->gene->coordinate, $genelist->last->gene->coordinate );
 } #sub unmapped_coordinates

#===============================================================================
# Accessor methods
#===============================================================================
 
sub listname
 {
  my ($genelist,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'listname'}[ $$genelist ] = $set_value );
  return ( $$structure{'listname'}[ $$genelist ] );
 } #sub listname
 
#------------------------------------------------------------
 
sub genome
 {
  my ($genelist,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'genome'}[ $$genelist ] = $set_value );
  return ( $$structure{'genome'}[ $$genelist ] );
 } #sub genome
 
#------------------------------------------------------------
 
sub elements
 {
  my ($genelist,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'elements'}[ $$genelist ] = $set_value );
  return ( $$structure{'elements'}[ $$genelist ] );
 } #sub elements
 
#------------------------------------------------------------
 
sub remapped
 {
  my ($genelist,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'remapped'}[ $$genelist ] = $set_value );
  return ( $$structure{'remapped'}[ $$genelist ] );
 } #sub remapped
 
#------------------------------------------------------------
 
sub profile
 {
  my ($genelist,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'profile'}[ $$genelist ] = $set_value );
  return ( $$structure{'profile'}[ $$genelist ] );
 } #sub profile
 
#------------------------------------------------------------
 
sub id
 {
  my ($genelist,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'id'}[ $$genelist ] = $set_value );
  return ( $$structure{'id'}[ $$genelist ] );
 } #sub id

1;
