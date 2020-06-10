package iADHoRe::gene;
#Created by Cedric Simillion on Fri Jun 25 16:30:05 CEST 2004

use strict;
use iADHoRe::record_structure;

our @ISA=qw(iADHoRe::record_structure);

our $structure=iADHoRe::gene->create_structure(qw(ID genome list coordinate orientation remapped_coordinate is_tandem 
                                   is_tandem_representative tandem_representative remapped pair_with));

#===============================================================================
# Package methods
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub add_pairs
# Stores for a gene object all genes with which it forms a gene
# pair. Takes a list of gene objects as argument.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub add_pairs
 {
  my ($gene,$pair_with);
  my (@pair_genes);
  
  ($gene,@pair_genes)=@_;
  $pair_with=$gene->pair_with;
  
  foreach (@pair_genes)
   {
    defined($_) || next; #Sometimes genenames are present in the BLAST table file that are not in a genelist. Ignore these.
    my $pair_nr=$$_;
    $$pair_with{$pair_nr}=-1;
   } #foreach my $pair (@pair_genes)

 } #sub add_pairs

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub is_pair_with
# Returns true if a gene object forms a gene pair with another gene
# object that is supplied as an argument
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub is_pair_with
 {
  my ($gene,$other_gene_nr,$pairs);
  $gene=$_[0];
  $other_gene_nr=${$_[1]};
  
  $pairs=$gene->pair_with;
  
  ( exists( $$pairs{$other_gene_nr} ) && $$pairs{$other_gene_nr} ) ? ( return -1) : ( return 0 );
  
 } #sub is_pair_with

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub is_indirect_pair_with
# Returns true if a gene object forms an indirect gene pair with
# another gene object that is supplied as an argument. An indirect
# gene pair is formed when two genes A and B pair with each other
# directly or if both pair with the same gene C. Used in tandem
# remapping.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub is_indirect_pair_with
 {
  my ($gene_a,$gene_b,$gene_a_pairs,$gene_b_pairs);
  ($gene_a,$gene_b)=@_;
  
  if ( $gene_a->is_pair_with($gene_b) )
   {
    return -1;
   } #if ( $gene_a->is_pair_with($gene_b) )
  else 
   {
    $gene_a_pairs=$gene_a->pair_with;
    $gene_b_pairs=$gene_b->pair_with;
    foreach my $gene_nr (keys %$gene_a_pairs)
     {
      if ( exists( $$gene_b_pairs{$gene_nr} ) && $$gene_b_pairs{$gene_nr} )
       {
        return -1;
	last;
       } #if ( exists( $$gene_b_pairs{$gene_nr} ) && $$gene_b_pairs{$gene_nr} )
     } #foreach my $gene_nr (keys %$gene_a_pairs)
   } #else
  
  return 0;
     
 } #sub is_indirect_pair_with

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub has_pairs
# Returns true if a gene object forms pairs with other genes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub has_pairs
 {
  my ($gene,$pairs);
  $gene=$_[0];

  $pairs=$gene->pair_with;
  scalar( keys %$pairs ) ? (return -1) : (return 0);

 } #sub has_pairs

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub output_table
# Returns a reference to a 2D table (an array of arrays) containing
# all the ouput data for all objects in this class.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub output_table
 {
  my $table;
  my @fields=qw(ID genome list coordinate orientation remapped_coordinate is_tandem is_tandem_representative tandem_representative remapped);
  $table = iADHoRe::gene->dump_table(@fields);

  foreach my $row ( @$table )
   {
    shift( @$row );
    $$row[2]= $$row[2]->listname;
    $$row[8]= defined($$row[8]) ? $$row[8]->ID : "";
   } #foreach my $row ( @$table )
  $fields[0] = 'id'; #The id field name must be in lower case in the output file. 
  unshift( @$table, \@fields );
  
  return($table); 
 } #sub output_table

#===============================================================================
# Accessor methods
#===============================================================================
 
sub ID
 {
  my ($gene,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'ID'}[ $$gene ] = $set_value );
  return ( $$structure{'ID'}[ $$gene ] );
 } #sub ID
 
#------------------------------------------------------------
 
sub genome
 {
  my ($gene,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'genome'}[ $$gene ] = $set_value );
  return ( $$structure{'genome'}[ $$gene ] );
 } #sub genome
 
#------------------------------------------------------------
 
sub list
 {
  my ($gene,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'list'}[ $$gene ] = $set_value );
  return ( $$structure{'list'}[ $$gene ] );
 } #sub list
 
#------------------------------------------------------------
 
sub coordinate
 {
  my ($gene,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'coordinate'}[ $$gene ] = $set_value );
  return ( $$structure{'coordinate'}[ $$gene ] );
 } #sub coordinate
 
#------------------------------------------------------------
 
sub orientation
 {
  my ($gene,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'orientation'}[ $$gene ] = $set_value );
  return ( $$structure{'orientation'}[ $$gene ] );
 } #sub orientation
 
#------------------------------------------------------------
 
sub remapped_coordinate
 {
  my ($gene,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'remapped_coordinate'}[ $$gene ] = $set_value );
  return ( $$structure{'remapped_coordinate'}[ $$gene ] );
 } #sub remapped_coordinate
 
#------------------------------------------------------------
 
sub is_tandem
 {
  my ($gene,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'is_tandem'}[ $$gene ] = $set_value );
  return ( $$structure{'is_tandem'}[ $$gene ] );
 } #sub is_tandem
 
#------------------------------------------------------------
 
sub is_tandem_representative
 {
  my ($gene,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'is_tandem_representative'}[ $$gene ] = $set_value );
  return ( $$structure{'is_tandem_representative'}[ $$gene ] );
 } #sub is_tandem_representative
 
#------------------------------------------------------------
 
sub tandem_representative
 {
  my ($gene,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'tandem_representative'}[ $$gene ] = $set_value );
  return ( $$structure{'tandem_representative'}[ $$gene ] );
 } #sub tandem_representative
 
#------------------------------------------------------------
 
sub remapped
 {
  my ($gene,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'remapped'}[ $$gene ] = $set_value );
  return ( $$structure{'remapped'}[ $$gene ] );
 } #sub remapped
 
#------------------------------------------------------------
 
sub pair_with
 {
  my ($gene,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'pair_with'}[ $$gene ] = $set_value );
  return ( $$structure{'pair_with'}[ $$gene ] );
 } #sub pair_with
 

1;

