package iADHoRe::multiplicon;
#Created by Cedric Simillion on Tue Aug 31 13:04:36 CEST 2004

use strict;
use iADHoRe::record_structure;
use iADHoRe::cluster;
use iADHoRe::profile;

our @ISA=qw(iADHoRe::cluster iADHoRe::record_structure);

our $structure=iADHoRe::multiplicon->base_upon("iADHoRe::cluster", qw(level baseclusters number_of_baseclusters profile id is_redundant) );

#===============================================================================
# Package methods
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub add_baseclusters
# Adds a list of baseclusters to a multiplicon
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub add_baseclusters
 {
  my ($multiplicon,$basecluster,$basecluster_array);
  my (@basecluster_list);
  ($multiplicon, @basecluster_list)=@_;
  
  $basecluster_array = $multiplicon->baseclusters;

  #Add each basecluster of the list to the multiplicon
  foreach my $basecluster ( @basecluster_list )
   {
    push( @$basecluster_array, $basecluster );
    $basecluster->multiplicon( $multiplicon );
    $multiplicon->number_of_baseclusters( $multiplicon->number_of_baseclusters + 1 );
    $multiplicon->number_of_anchorpoints( $multiplicon->number_of_anchorpoints + $basecluster->number_of_anchorpoints );

    #Tell the anchorpoints of the basecluster that they're now also part of a multiplicon;
    foreach my $anchorpoint ( $basecluster->get_anchorpoints )
     {
      $anchorpoint->multiplicon( $multiplicon );
     } #foreach my $anchorpoint ( $basecluster->get_anchorpoints )
     
   } #foreach my $basecluster ( @basecluster_list )
   
 } #sub add_baseclusters

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_anchorpoints
# Returns a list of all anchorpoints of all baseclusters in a
# multiplicon.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_anchorpoints
 {
  my ($multiplicon,$baseclusters);
  my (@anchorpoints);
  $multiplicon=$_[0];
  $baseclusters=$multiplicon->baseclusters;
  
  foreach my $cluster ( @$baseclusters )
   {
    push(@anchorpoints, $cluster->get_anchorpoints );
   } #foreach my $baseclusters ( @$baseclusters )
   
  return (@anchorpoints);
  
 } #sub get_anchorpoints

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_baseclusters
# Returns a list of all baseclusters in a multiplicon.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_baseclusters
 {
  my ($multiplicon,$baseclusters);
  $multiplicon=$_[0];
  $baseclusters=$multiplicon->baseclusters;
  
  return (@$baseclusters);
  
 } #sub get_baseclusters

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub empty
# Empties a multiplicon: i.e. the array of baseclusters is emptied.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub empty
 {
  my $multiplicon=$_[0];

  $multiplicon->baseclusters( [] );
  $multiplicon->number_of_baseclusters( 0 );

 } #sub empty

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub create_profile
# Creates a profile from the segments in the multiplicon
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub create_profile
 {
  my ($multiplicon)=$_[0];
  $multiplicon->profile( iADHoRe::profile->create( $multiplicon ) );
 } #sub create_profile


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub mask_segments
# Mask the segments in a multiplicon. If the multiplicon has level 2,
# both the x_object and the y_segment will be masked. Otherwise only
# the y_segment will be masked.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub mask_segments
 {
  my ($multiplicon)=$_[0];
  
  ($multiplicon->level == 2) && $multiplicon->x_object->mask( $multiplicon->begin_x , $multiplicon->end_x );
  $multiplicon->y_list->mask( $multiplicon->begin_y , $multiplicon->end_y );
  
 } #sub mask_segments

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub unmapped_coordinates
# Returns the unmapped coordinates of a cluster.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub unmapped_coordinates
 {
  my $multiplicon=$_[0];
  my @coordinates;
  
  if ( $multiplicon->level == 2 )
   {
    foreach my $x ( $multiplicon->begin_x, $multiplicon->end_x )
     {
       push( @coordinates, ${$multiplicon->x_object->remapped}[$x]->gene->coordinate );
     } #foreach my $x ( $multiplicon->begin_x, $multiplicon->end_x )
   } #if ( $multiplicon->level == 2 )
  else 
   {
    push( @coordinates, $multiplicon->begin_x, $multiplicon->end_x )
   } #else 

  foreach my $y ( $multiplicon->begin_y, $multiplicon->end_y )
   {
    push( @coordinates, ${$multiplicon->y_list->remapped}[$y]->gene->coordinate );
   } #foreach my $coordinate ( $multiplicon->begin_x, $multiplicon->end_x, $multiplicon->begin_y, $multiplicon->end_y )

  return( @coordinates );
   
 } #sub unmapped_coordinates

#===============================================================================
# Accessor methods
#===============================================================================
                                                                                                                              
sub level
 {
  my ($multiplicon,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'level'}[ $$multiplicon ] = $set_value );
  return ( $$structure{'level'}[ $$multiplicon ] );
 } #sub level
                                                                                                                              
#------------------------------------------------------------
                                                                                                                              
sub baseclusters
 {
  my ($multiplicon,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'baseclusters'}[ $$multiplicon ] = $set_value );
  return ( $$structure{'baseclusters'}[ $$multiplicon ] );
 } #sub baseclusters
                                                                                                                              
#------------------------------------------------------------
                                                                                                                              
sub number_of_baseclusters
 {
  my ($multiplicon,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'number_of_baseclusters'}[ $$multiplicon ] = $set_value );
  return ( $$structure{'number_of_baseclusters'}[ $$multiplicon ] );
 } #sub number_of_baseclusters
                                                                                                                              
#------------------------------------------------------------
                                                                                                                              
sub profile
 {
  my ($multiplicon,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'profile'}[ $$multiplicon ] = $set_value );
  return ( $$structure{'profile'}[ $$multiplicon ] );
 } #sub profile
                                                                                                                              
#------------------------------------------------------------
                                                                                                                              
sub id
 {
  my ($multiplicon,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'id'}[ $$multiplicon ] = $set_value );
  return ( $$structure{'id'}[ $$multiplicon ] );
 } #sub id

#------------------------------------------------------------

sub is_redundant
 {
  my ($multiplicon,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'is_redundant'}[ $$multiplicon ] = $set_value );
  return( $$structure{'is_redundant'}[ $$multiplicon ] );
 } #sub is_redundant

1;
