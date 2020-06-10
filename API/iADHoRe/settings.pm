package iADHoRe::settings;
#Created by Cedric Simillion on Fri Jun 25 09:57:43 CEST 2004

use strict;
use iADHoRe::listfile;

our $AUTOLOAD;
our %settings;

#===============================================================================
# Package methods
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_settings (the settings object creator)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_settings
 {
  scalar( keys %settings ) || &read_settings;
  my $settings_object= \%settings;
  bless $settings_object, $_[0];
  return $settings_object;
 } #sub get_settings


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub read_settings
# Reads and parses the settings file into %settings;
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub read_settings
 {
  my $settings;
  
  #-----------------------------------------------------------
  # Read settings file (passed as first parameter to the main
  # program)
  #-----------------------------------------------------------
  open (IN,$ARGV[0]) || die "Could not open settings file $_[0]: $!\n";
  $settings=join("",<IN>);
  close IN;

  #-----------------------------------------------------------
  # Parse listfile names for each genome
  #-----------------------------------------------------------
  while ($settings=~ s/genome=((?:.|\n)*?)\n(\w+=)/$2/)
   {
    my ($genome,$list_string);
    my (@lists);
    $list_string=$1;
    ($genome,@lists)=split(/\n/,$list_string);
    $genome=~ s/ +//g;
    foreach (@lists)
     {
      /(\w+)\s+(\S+)/;
      $settings{"genomes"}{$genome}{$1}=$2;
      (-e $2) || die "Could not locate chromosome list file $2.\n";
     } #foreach (@lists)
   } #while ($settings=~ s/genome=(.|\n)*?\n(\w+=)/$2/)

  #-----------------------------------------------------------
  # Parse BLAST file name
  #-----------------------------------------------------------
  $settings=~ /blast_table= *(\S+)/;
  $settings{"blast_table"}=$1;
  (-f $settings{"blast_table"}) || (die "Could not find blast_table file ".$settings{"blast_table"}.".\n");

  #-----------------------------------------------------------
  # Parse output directory (if none given, current directory
  # is taken)
  #-----------------------------------------------------------
  if ($settings=~ /output_path= *(\S+)/)
   {
    $settings{"output_path"}=$1;
    unless (-d $1)
     {
      mkdir ($1,0777) || die "Could not create output directory $1: $!\n";
     } #unless (-d $1)
    $settings{"output_path"}=~ s/([^\/])$/$1\//;     
   } #if ($settings=~ /output_path= *(\S+)/)
  else
   {
    $settings{"output_path"}='';
   } #else

  #-----------------------------------------------------------
  # Parse algorithm parameters
  #-----------------------------------------------------------
  ($settings=~ /gap_size= *(\d+)/) ? ($settings{"gap_size"}=$1) : (die "Gap size not specified!");
  ($settings=~ /cluster_gap= *(\d+)/) ? ($settings{"cluster_gap"}=$1) : (die "Cluster gap size not specified!");
  ($settings=~ /q_value= *(\d*\.?\d+)/) ? ($settings{"q_value"}=$1) : (die "Q value not specified!");
  ($settings=~ /prob_cutoff= *(\d*\.?\d+)/) ? ($settings{"prob_cutoff"}=$1) : (die "Probability cutoff not specified!");
  ($settings=~ /anchor_?points= *(\d+)/) ? ($settings{"anchorpoints"}=$1)
   : (die "Minimum number of anchor points per segment size not specified!");

  #-----------------------------------------------------------
  # Other settings
  #-----------------------------------------------------------
  ($settings=~ /level_2_only= *true/) ? ($settings{"level_2_only"}=-1) : ($settings{"level_2_only"}=0);
  ($settings=~ /database_name= *(\S+)/) ? ($settings{"database_name"}=$1) : ($settings{"database_name"}='');
  ($settings=~ /database_host= *(\S+)/) ? ($settings{"database_host"}=$1) : ($settings{"database_host"}='');
  ($settings=~ /database_login= *(\S+)/) ? ($settings{"database_login"}=$1) : ($settings{"database_login"}='');
  
 } #sub read_settings

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub AUTOLOAD
# Function that returns an attribute of a settings object when this
# attribute is not explicitly defined.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub AUTOLOAD
 {
  my ($object,$attribute);
  $object=$_[0];
  $attribute=$AUTOLOAD;
  $attribute=~ s/^.*:://;
  
  if (exists $$object{$attribute})
   {
    return $$object{$attribute};
   } #if (exists $$object{$attribute})
  else 
   {
    die "The attribute $attribute does not exist for ".ref($object)." objects";
   } #else 
 } #sub AUTOLOAD

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub listfiles
# Returns all listfiles described in the settings files as an array
# of listfile objects.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub listfiles
 {
  my $object=$_[0];
  my (@listfiles);
  
  foreach my $genome (keys %{$$object{"genomes"}})
   {
    foreach my $list ( keys %{$$object{"genomes"}{$genome}} )
     {
      push (@listfiles, iADHoRe::listfile->new( $genome, $list, $$object{"genomes"}{$genome}{$list} ) )
     } #foreach my $list ( keys %{$$object{"genomes"}{$genome}} )
   } #foreach my $genome (keys %{$object}{"genomes"})
  
  return (@listfiles); 
 } #sub listfiles

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub gap_sizes
# Returns 10 or fewer integer values ranging from 3 to $max_dist, 
# ascending exponentially.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub gap_sizes
 {
  my ($step,$i,$value,$settings);
  my (@gap);
  
  $settings=$_[0];
  $step=( log($settings->gap_size) / log(3)-1 ) / 9;
  for ($i=0;$i<10;$i++)
   {
    $value=sprintf("%.0f",3**(1+($i*$step)));
    if($i==0)
     {
      push (@gap,$value);
     } #if($i==0)
    elsif ($value!=$gap[-1]) 
     {
      push (@gap,$value);
     } #elsif ($value!=$gap[-1]) 
   } #for ($i=0;$i<10;$i++) 
  return (@gap); 
 } #gap_sizes

1;
