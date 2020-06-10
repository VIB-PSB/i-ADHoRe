package iADHoRe::listfile;
#Created by Cedric Simillion on Fri Jun 25 13:13:37 CEST 2004

use strict;
use iADHoRe::record_structure;

our @ISA=qw(iADHoRe::record_structure);

our $structure=iADHoRe::listfile->create_structure(qw(genome listname file));

#===============================================================================
# Accessor methods
#===============================================================================
 
sub genome
 {
  my ($listfile,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'genome'}[ $$listfile ] = $set_value );
  return ( $$structure{'genome'}[ $$listfile ] );
 } #sub genome
 
#------------------------------------------------------------
 
sub listname
 {
  my ($listfile,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'listname'}[ $$listfile ] = $set_value );
  return ( $$structure{'listname'}[ $$listfile ] );
 } #sub listname
 
#------------------------------------------------------------
 
sub file
 {
  my ($listfile,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'file'}[ $$listfile ] = $set_value );
  return ( $$structure{'file'}[ $$listfile ] );
 } #sub file
 

1;
