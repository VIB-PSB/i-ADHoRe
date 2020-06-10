package iADHoRe::record_structure;
#Created by Cedric Simillion on Fri Jun 25 14:01:25 CEST 2004

use strict;

our (%structure,%free);

our $AUTOLOAD;

#===============================================================================
# Package methods
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub create_structure
# creates a new data structure;
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub create_structure
 {
  my ($structure_name);
  my (@attributes);
  
  $structure_name=shift @_;
  @attributes=@_;
  
  #create an empty array for each attribute
  foreach my $attr (@attributes)
   {
    $structure{$structure_name}{$attr}=[];
   } #foreach my $attr (@attributes)
  
  #Store the order in which the attributes were specified
  $structure{$structure_name}{"attributes"}=[ @attributes ];
  
  #Initialize the free array for the datastructure that contains the list of free slots in the datastructure;
  $free{$structure_name}=[ 0 ];
   
  bless $structure{$structure_name}, 'iADHoRe::record_structure';
  return $structure{$structure_name};
  
 } #sub create_structure

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub base_upon
# Creates a new data structure with all the attributes of an existing
# one. The name of the original data structure is passed as argument,
# as wel as any additional attributes.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub base_upon
 {
  my ($original,$structure_name);
  my (@attributes);
  
  $structure_name=shift @_;
  $original= shift @_;
  
  @attributes= @{$structure{$original}{"attributes"}};
  push(@attributes, @_);
  
  $structure_name->create_structure( @attributes );
  
 } #sub base_upon
 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub new
# Adds a new record to an existing data structure
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub new
 {
  my ($record_nr,$structure_name);
  my (@attributes,@values);
  
  $structure_name=shift @_;
  @attributes=@{$structure{$structure_name}{"attributes"}};
  @values=@_;
  
  $record_nr=shift @{$free{$structure_name}};
  foreach my $i (0..$#attributes)
   {
    if ( exists $values[$i] )
     {
      $structure{ $structure_name }{ $attributes[$i] }[ $record_nr ]=$values[$i];
     } #if (exists($values[$i])
    else 
     {
      $structure{ $structure_name }{ $attributes[$i] }[ $record_nr ]=undef;
     } #else 
   } #foreach my $i (0..$#attributes)

  scalar( @{$free{$structure_name}} ) || ( @{$free{$structure_name}} = scalar( @{ $structure{$structure_name}{$attributes[0]} } ) );
  
  bless \$record_nr, $structure_name;
  
 } #sub new

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub duplicate
# Creates a copy of an existing object;
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub duplicate
 {
  my ($object,$class,$duplicate);
  my (@attribute_names,@new_attributes);

  $object=$_[0];
  $class= ref( $object );
  @attribute_names=@{ $structure{$class}{"attributes"} };
  
  foreach my $attribute ( @attribute_names )  
   {
    push( @new_attributes, $structure{$class}{$attribute}[ $$object ] );
   } #foreach my $attribute ( @attribute_names )  
   
  $duplicate=$class->new( @new_attributes );
  
  return( $duplicate ); 
 } #sub duplicate

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub dump_table
# Returns a 2D-array containing selected attributes of all elements
# in a datastructure. The first column is the index number, all other
# columns are the selected attributes, in the specified order.
# Arguments: a list of attributes to be returned
# Returns: A reference an array of arrays.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub dump_table
 {
  my ($class,$number_of_rows);
  my (@attributes, @table);
  ($class,@attributes)=@_;
  
  $number_of_rows=scalar( @{ $structure{$class}{$attributes[0]} } );
  foreach my $i ( 0 .. $number_of_rows-1)
   {
    my (@row);
    defined( $structure{$class}{$attributes[0]}[$i] ) || next;
    push( @row, $i );
    foreach ( @attributes )
     {
      push( @row, $structure{$class}{$_}[$i] );
     } #foreach ( @attributes )
    push(@table, \@row );
   } #foreach my $i ( 0 .. $number_of_rows-1)
  
  return(\@table); 
 } #sub dump_table

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_index
# Returns the index of an object in the datastructure, which is a
# simple matter of de-referencing.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_index
 {
  my $object=$_[0];
  return( $$object );
 } #sub get_index

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub dump_object
# Returns a hash containing all attributes of a current object.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub dump_object
 {
  my ($structure_name,$object);
  my (%dump);
  
  $object=$_[0];
  $structure_name=ref $object;
  
  foreach my $attr ( @{$structure{$structure_name}{"attributes"}} )
   {
    $dump{$attr} = $structure{$structure_name}{$attr}[ $$object ];
   } #foreach my $attr ( @{$structure{$structure_name}{"attributes"}} )

  return (\%dump);

 } #sub dump_object

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub DESTROY
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub DESTROY
 {
  my ($structure_name,$object,$lastfree);
  
  $object=$_[0];
  $structure_name=ref $object;
  $lastfree=$free{$structure_name}[-1];
  
  #When an object is destroyed manually, and gets subsequently out of scope, the DESTROY method is invoked
  #twice. Therefore, we first need to check if the value of $$object is defined.
  if (defined $$object )
   {
    foreach my $attr ( @{$structure{$structure_name}{"attributes"}} )
     {
      undef $structure{$structure_name}{$attr}[ $$object ];
     } #foreach my $attr ( @{$structure{$structure_name}{"attributes"}} )

    push (@{$free{$structure_name}},$$object);
    undef $$object;
   } #if ( $last_free != $$object )
        
 } #sub DESTROY

1;
