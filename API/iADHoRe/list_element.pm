package iADHoRe::list_element;
#Created by Cedric Simillion on Wed Jan  5 21:42:02 GMT 2005

use strict;
use iADHoRe::record_structure;

our @ISA=qw(iADHoRe::record_structure);

our $structure=iADHoRe::list_element->create_structure(qw(gene orientation masked id));

#===============================================================================
# Package methods
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub is_pair_with
# Returns true if the gene object of a list_element forms a gene pair
# with the gene of another list element that is supplied as an
# argument.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub is_pair_with
 {
  my ($element,$other_element)=@_;

  return( $element->gene->is_pair_with( $other_element->gene ) );
 } #sub is_pair_with

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub is_indirect_pair_with
# Returns true if a gene object of a list element forms an indirect
# gene pair with the gene object of another list element that is 
# supplied as an argument. An indirect gene pair is formed when two 
# genes A and B pair with each other directly or if both pair with the
# same gene C. Used in tandem remapping.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub is_indirect_pair_with
 {
  my ($element,$other_element)=@_;

  return( $element->gene->is_indirect_pair_with( $other_element->gene ) );
 } #sub is_pair_with

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub has_pairs
# Returns true if a gene object forms pairs with other genes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub has_pairs
 {
  my $element=$_[0];
  
  return( $element->gene->has_pairs );
 } #sub has_pairs

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub matching_positions
# Returns for a gene all coordinates of matching genes in a list.
# The list is passed as a reference.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub matching_positions
 {
  my ($element,$list,$pairs);
  my (@positions);
  my (%in_list);
  
  ($element,$list)=@_;
  
  foreach my $i ( 0..$#{$list} )
   {
    ( ref( $$list[$i] ) eq 'iADHoRe::list_element' ) || next;
    $$list[$i]->masked && next;
    $in_list{ ${$$list[$i]->gene} }=$i;
   } #foreach my $i ( 0..$#{$list} )
   
  $pairs=$element->gene->pair_with; 
  
  foreach my $match (keys %$pairs)
   {
    exists( $in_list{$match} ) && push( @positions, $in_list{$match} );
   } #foreach my $match (keys %$pairs)
  
  return ( @positions );
  
 } #sub matching_positions

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub remap_to
# Remaps a gene object 'a' onto a second one 'b' that is supplied as
# an argument.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub remap_to
 {
  my ($gene_a,$gene_b,$element_a,$element_b,$gene_a_pairs,$gene_b_pairs);
  ($element_a,$element_b)=@_;
  $gene_a=$element_a->gene;
  $gene_b=$element_b->gene;

  $gene_a->is_tandem(-1);
  $gene_a->is_tandem_representative(0);
  $gene_a->tandem_representative($gene_b);
  $gene_a->remapped(-1);
  
  $gene_b->is_tandem(-1);
  $gene_b->is_tandem_representative(-1);
  $gene_b->tandem_representative($gene_b);
  $gene_b->remapped(0);
  
 } #sub remap_to

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub invert_orientation
# Inverts the orientation of a list element object;
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub invert_orientation
 {
  my($element,$new_orientation);
  $element=$_[0];
  
  $new_orientation= ($element->orientation eq '+') ? '-' : '+';
  
  $element->orientation( $new_orientation );
 } #sub invert_orientation

#===============================================================================
# Accessor methods
#===============================================================================
 
sub gene
 {
  my ($list_element,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'gene'}[ $$list_element ] = $set_value );
  return ( $$structure{'gene'}[ $$list_element ] );
 } #sub gene
 
#------------------------------------------------------------
 
sub orientation
 {
  my ($list_element,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'orientation'}[ $$list_element ] = $set_value );
  return ( $$structure{'orientation'}[ $$list_element ] );
 } #sub orientation
 
#------------------------------------------------------------
 
sub masked
 {
  my ($list_element,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'masked'}[ $$list_element ] = $set_value );
  return ( $$structure{'masked'}[ $$list_element ] );
 } #sub masked
 
#------------------------------------------------------------
 
sub id
 {
  my ($list_element,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'id'}[ $$list_element ] = $set_value );
  return ( $$structure{'id'}[ $$list_element ] );
 } #sub id

1;
