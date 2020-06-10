package iADHoRe::cluster;
#Created by Cedric Simillion on Sun Jul 11 17:27:01 CEST 2004

use strict;
use iADHoRe::record_structure;
use iADHoRe::anchorpoint;

our @ISA=qw(iADHoRe::record_structure);

iADHoRe::cluster->create_structure(qw(x_object y_list begin_x end_x begin_y end_y number_of_anchorpoints));

#===============================================================================
# Package methods
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub lowest_x
# Returns for all anchorpoints of a basecluster the lowest
# x-coordinate
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub lowest_x
 {
  my ($cluster,$lowest_x);
  my (@anchorpoints);

  $cluster=$_[0];
  @anchorpoints=$cluster->get_anchorpoints;
  
  $lowest_x=$anchorpoints[0]->x;
  shift (@anchorpoints);
  
  foreach my $ap (@anchorpoints)  
   {
    ($ap->x < $lowest_x) && ($lowest_x=$ap->x);
   } #foreach my $ap (@anchorpoints)
   
  return $lowest_x; 
 } #sub lowest_x

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub highest_x
# Returns for all anchorpoints of a basecluster the highest
# x-coordinate
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub highest_x
 {
  my ($cluster,$highest_x);
  my (@anchorpoints);

  $cluster=$_[0];
  @anchorpoints=$cluster->get_anchorpoints;
  
  $highest_x=$anchorpoints[0]->x;
  shift (@anchorpoints);
  
  foreach my $ap (@anchorpoints)  
   {
    ($ap->x > $highest_x) && ($highest_x=$ap->x);
   } #foreach my $ap (@anchorpoints)
   
  return $highest_x; 
 } #sub highest_x

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub lowest_y
# Returns for all anchorpoints of a basecluster the lowest
# y-coordinate
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub lowest_y
 {
  my ($cluster,$lowest_y);
  my (@anchorpoints);

  $cluster=$_[0];
  @anchorpoints=$cluster->get_anchorpoints;
  
  $lowest_y=$anchorpoints[0]->y;
  shift (@anchorpoints);
  
  foreach my $ap (@anchorpoints)  
   {
    ($ap->y < $lowest_y) && ($lowest_y=$ap->y);
   } #foreach my $ap (@anchorpoints)
   
  return $lowest_y; 
 } #sub lowest_y

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub highest_y
# Returns for all anchorpoints of a basecluster the highest
# y-coordinate
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub highest_y
 {
  my ($cluster,$highest_y);
  my (@anchorpoints);

  $cluster=$_[0];
  @anchorpoints=$cluster->get_anchorpoints;
  
  $highest_y=$anchorpoints[0]->y;
  shift (@anchorpoints);
  
  foreach my $ap (@anchorpoints)  
   {
    ($ap->y > $highest_y) && ($highest_y=$ap->y);
   } #foreach my $ap (@anchorpoints)
   
  return $highest_y; 
 } #sub highest_y

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub set_bounds
# Sets the begin_x, end_x, begin_y and end_y attributes of a cluster
# object.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub set_bounds
 {
  my $cluster=$_[0];
  
  $cluster->begin_x( $cluster->lowest_x );
  $cluster->end_x( $cluster->highest_x );
  $cluster->begin_y( $cluster->lowest_y );
  $cluster->end_y( $cluster->highest_y );
  
 } #sub set_bounds

#===============================================================================
# Accessor methods
#===============================================================================
 
sub x_object
 {
  my ($cluster,$set_value,$structure,$package);
  ($cluster,$set_value)=@_;
  $package = ref( $cluster );
  $structure = $iADHoRe::record_structure::structure{$package};
  
  defined( $set_value ) && ( $$structure{'x_object'}[ $$cluster ] = $set_value );
  return ( $$structure{'x_object'}[ $$cluster ] );
 } #sub x_object
 
#------------------------------------------------------------
 
sub y_list
 {
  my ($cluster,$set_value,$structure,$package);
  ($cluster,$set_value)=@_;
  $package = ref( $cluster );
  $structure = $iADHoRe::record_structure::structure{$package};
  
  defined( $set_value ) && ( $$structure{'y_list'}[ $$cluster ] = $set_value );
  return ( $$structure{'y_list'}[ $$cluster ] );
 } #sub y_list
 
#------------------------------------------------------------
 
sub begin_x
 {
  my ($cluster,$set_value,$structure,$package);
  ($cluster,$set_value)=@_;
  $package = ref( $cluster );
  $structure = $iADHoRe::record_structure::structure{$package};
  
  defined( $set_value ) && ( $$structure{'begin_x'}[ $$cluster ] = $set_value );
  return ( $$structure{'begin_x'}[ $$cluster ] );
 } #sub begin_x
 
#------------------------------------------------------------
 
sub end_x
 {
  my ($cluster,$set_value,$structure,$package);
  ($cluster,$set_value)=@_;
  $package = ref( $cluster );
  $structure = $iADHoRe::record_structure::structure{$package};
  
  defined( $set_value ) && ( $$structure{'end_x'}[ $$cluster ] = $set_value );
  return ( $$structure{'end_x'}[ $$cluster ] );
 } #sub end_x
 
#------------------------------------------------------------
 
sub begin_y
 {
  my ($cluster,$set_value,$structure,$package);
  ($cluster,$set_value)=@_;
  $package = ref( $cluster );
  $structure = $iADHoRe::record_structure::structure{$package};
  
  defined( $set_value ) && ( $$structure{'begin_y'}[ $$cluster ] = $set_value );
  return ( $$structure{'begin_y'}[ $$cluster ] );
 } #sub begin_y
 
#------------------------------------------------------------
 
sub end_y
 {
  my ($cluster,$set_value,$structure,$package);
  ($cluster,$set_value)=@_;
  $package = ref( $cluster );
  $structure = $iADHoRe::record_structure::structure{$package};
  
  defined( $set_value ) && ( $$structure{'end_y'}[ $$cluster ] = $set_value );
  return ( $$structure{'end_y'}[ $$cluster ] );
 } #sub end_y
 
#------------------------------------------------------------
 
sub number_of_anchorpoints
 {
  my ($cluster,$set_value,$structure,$package);
  ($cluster,$set_value)=@_;
  $package = ref( $cluster );
  $structure = $iADHoRe::record_structure::structure{$package};
  
  defined( $set_value ) && ( $$structure{'number_of_anchorpoints'}[ $$cluster ] = $set_value );
  return ( $$structure{'number_of_anchorpoints'}[ $$cluster ] );
 } #sub number_of_anchorpoints
 

1;
