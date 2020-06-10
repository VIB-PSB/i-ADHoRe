package iADHoRe::anchorpoint;
#Created by Cedric Simillion on Mon Jul  5 19:55:59 CEST 2004

use strict;
use iADHoRe::record_structure;

our @ISA=qw(iADHoRe::record_structure);

our $structure=iADHoRe::anchorpoint->create_structure( qw( gene_x gene_y x y is_real_anchorpoint basecluster multiplicon id) );

#===============================================================================
# Accessor methods
#===============================================================================
 
sub gene_x
 {
  my ($anchorpoint,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'gene_x'}[ $$anchorpoint ] = $set_value );
  return ( $$structure{'gene_x'}[ $$anchorpoint ] );
 } #sub gene_x
 
#------------------------------------------------------------
 
sub gene_y
 {
  my ($anchorpoint,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'gene_y'}[ $$anchorpoint ] = $set_value );
  return ( $$structure{'gene_y'}[ $$anchorpoint ] );
 } #sub gene_y
 
#------------------------------------------------------------
 
sub x
 {
  my ($anchorpoint,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'x'}[ $$anchorpoint ] = $set_value );
  return ( $$structure{'x'}[ $$anchorpoint ] );
 } #sub x
 
#------------------------------------------------------------
 
sub y
 {
  my ($anchorpoint,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'y'}[ $$anchorpoint ] = $set_value );
  return ( $$structure{'y'}[ $$anchorpoint ] );
 } #sub y
 
#------------------------------------------------------------
 
sub is_real_anchorpoint
 {
  my ($anchorpoint,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'is_real_anchorpoint'}[ $$anchorpoint ] = $set_value );
  return ( $$structure{'is_real_anchorpoint'}[ $$anchorpoint ] );
 } #sub is_real_anchorpoint
 
#------------------------------------------------------------
 
sub basecluster
 {
  my ($anchorpoint,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'basecluster'}[ $$anchorpoint ] = $set_value );
  return ( $$structure{'basecluster'}[ $$anchorpoint ] );
 } #sub basecluster
 
#------------------------------------------------------------
 
sub multiplicon
 {
  my ($anchorpoint,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'multiplicon'}[ $$anchorpoint ] = $set_value );
  return ( $$structure{'multiplicon'}[ $$anchorpoint ] );
 } #sub multiplicon
#------------------------------------------------------------
 
sub id
 {
  my ($anchorpoint,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'id'}[ $$anchorpoint ] = $set_value );
  return ( $$structure{'id'}[ $$anchorpoint ] );
 } #sub id
 

1;

