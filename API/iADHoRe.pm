package iADHoRe;
#Created by Cedric Simillion on Fri Jun 25 23:46:18 CEST 2004

use strict;
use lib "/Users/cesim/code/i-adhore/API";
use iADHoRe::settings;
use iADHoRe::genelist;
use iADHoRe::genepairs;
use iADHoRe::ghm;

my $settings=iADHoRe::settings->get_settings;

#===============================================================================
# Package methods
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub new_dataset
# Creates a new iADHoRe dataset. A dataset consists of an array of
# gene lists, an array of multiplicons that need to be evaluated
# and an array of already evaluated multiplicons.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub new_dataset
 {
  my ($pkg)=$_[0];
  my (%dataset);
  
  foreach my $listfile ($settings->listfiles)
   {
    push (@{$dataset{"genelists"}},iADHoRe::genelist->create_from_listfile($listfile));
   } #foreach my $listfile ($settings->listfiles)
  
  $dataset{"multiplicons2evaluate"}=[];
  $dataset{"evaluated_multiplicons"}=[];
   
  bless \%dataset, $pkg;
  
  return \%dataset; 
  
 } #sub new_dataset

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub genelists
# Returns an array of all genelists objects of the iADHoRe
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub genelists
 {
  return (sort { ($a->genome cmp $b->genome) || ($a->listname cmp $b->listname)  } ( @{$_[0]{"genelists"}} ) )
 } #sub genelists

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_gene_pairs
# Adds for every gene of every genelist in the iADHoRe all homologous
# genes as read from the BLAST table.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_gene_pairs
 {
  my ($dataset,$pairs);
  my (%gene_name2object);

  print STDERR "Mapping gene pairs...";
  
  $dataset=$_[0]; 
  $pairs=iADHoRe::genepairs->read_pairs;
  
  foreach my $list ( $dataset->genelists )
   {
    foreach my $element ( $list->elementlist )
     {
      my $name=$element->gene->ID;
      $gene_name2object{$name}=$element->gene;
     } #foreach my $gene ( $list->elements )
   } #foreach my $list ( $dataset->genelists )
  
  foreach my $list ( $dataset->genelists )
   {
    foreach my $element ( $list->elementlist )
     {
      my $name=$element->gene->ID;
      my @pair_names;
      
      @pair_names=$pairs->pairs_of($name);
      $element->gene->add_pairs( @gene_name2object{@pair_names} )
      
     } #foreach my $gene ( $list->elements )
   } #foreach my $list ( $dataset->genelists )

  print STDERR "\n";

 } #sub get_gene_pairs

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub remap_tandems
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub remap_tandems
 {
  print STDERR "Remapping tandem duplicates";
  foreach my $genelist ($_[0]->genelists)
   {
    $genelist->remap_tandems;
    print STDERR ".";
   } #foreach ($_[0]->genelists)
  print STDERR "\n"; 
 } #sub remap_tandems

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub level1_adhore
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub level1_adhore
 {
  my ($dataset);
  my (@genelists,@multiplicons);
  
  print STDERR "Detecting level 2 multiplicons\n";
  
  $dataset=$_[0];
  @genelists=$dataset->genelists;
  
  foreach my $x (0..$#genelists)
   {

    foreach my $y ($x..$#genelists)
     {
      my $ghm=iADHoRe::ghm->new( $genelists[$x], $genelists[$y] );
      push( @multiplicons, $ghm->adhore ); 
      print STDERR ".";
     } #foreach my $y ($x..$#genelists)

    print STDERR "\n"; 

   } #foreach my $x (0..$#genelists)

  @multiplicons = sort by_multiplicon_size ( @multiplicons );
  $dataset->add_multiplicons( @multiplicons );

 } #sub level1_adhore

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub profile_detection
# Detects higher level multiplicons by creating profiles from each
# multiplicon in the array "multiplicons2evaluate" of an iADHoRe
# dataset and searching all genelists in the "genelist" array.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub profile_detection
 {
  my $dataset=$_[0];
  print STDERR "\nHigher level detection\n";
  while (my $multiplicon = $dataset->next_multiplicon)
   {
    my (@new_multiplicons);
    print STDERR $dataset->multiplicons_left." multiplicons to evaluate - ";
    print STDERR "evaluating level ".$multiplicon->level." multiplicon...";
    $multiplicon->create_profile;

    #If a profile is not valid, destroy it and the segments referencing to it.
    if (!$multiplicon->profile->valid)
     {
      foreach my $segment ( @{$multiplicon->profile->segments} )
       {
        $segment->DESTROY;
       } #foreach my $segment ( @{$multiplicon->profile->segments} )
      $multiplicon->DESTROY; 
      next;
     } #if (!$multiplicon->profile->valid)

    $multiplicon->mask_segments;

    unless( $settings->level_2_only )
     {
      @new_multiplicons = $dataset->profile_search( $multiplicon->profile );
      @new_multiplicons = sort by_multiplicon_size ( @new_multiplicons );
      print STDERR " ".scalar( @new_multiplicons ). " new multiplicons found.";
      $dataset->add_multiplicons( @new_multiplicons );
     } #unless( $settings->level_2_only )

    print STDERR "\n";
    $dataset->store_multiplicon( $multiplicon );
   } #while (my $multiplicon = $dataset->next_multiplicon)
  
 } #sub profile_detection

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub profile_search
# Takes a profile object as input and searches all gene lists for
# matching segments. These segments are added and returned as new
# multiplicons
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub profile_search
 {
  my ($dataset,$profile)=@_;
  my (@multiplicons);
  
  foreach my $genelist ( $dataset->genelists )
   {
    ($genelist->size - $genelist->number_of_masked_elements < $settings->anchorpoints) && next;
    my $ghm=iADHoRe::ghm->new( $profile, $genelist );
    push( @multiplicons, $ghm->adhore );
   } #foreach my $genelist ( $dataset->genelists )
  
  return( @multiplicons );
   
 } #sub profile_search

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub add_multiplicons
# Adds new multiplicons to the "multiplicons2evaluate" array of a
# dataset.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub add_multiplicons
 {
  my ($dataset)=shift(@_);
  push( @{ $$dataset{"multiplicons2evaluate"} }, @_ );
 } #sub add_multiplicons

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub store_multiplicon
# Stores a multiplicon to the "evaluated_multiplicons" array of a
# dataset.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub store_multiplicon
 {
  my ($dataset)=shift(@_);
  push( @{ $$dataset{"evaluated_multiplicons"} }, $_[0] );
 } #sub store_multiplicon

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub next_multiplicon
# Returns the next multiplicon to be evaluated. This multiplicon
# is then removed from the "multiplicons2evaluate" array.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub next_multiplicon
 {
  my ($dataset)=shift(@_);
  return( pop( @{ $$dataset{"multiplicons2evaluate"} } ) );
 } #sub next_multiplicon

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub multiplicons_left
# Returns the number of multiplicons left in the dataset.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub multiplicons_left
 {
  my ($dataset)=shift(@_);
  #Add 1 because the next multiplicon has already been popped of the array
  return( scalar( @{ $$dataset{"multiplicons2evaluate"} } ) + 1 ); 
 } #sub multiplicons_left

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub output
# Outputs the processed dataset to the different files in the
# designated output directory.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub output
 {
  my ($dataset)=shift(@_);
  my (%tables);

  print STDERR "\n\nGenerating output files...";
  #We can get the genes table complete from the iADHoRe::gene::output_table class method.
  $tables{"genes"}=iADHoRe::gene->output_table;

  #Add header rows to the other output tables
  push ( @{$tables{"multiplicons"}},
   [qw(id genome_x list_x parent genome_y list_y level number_of_anchorpoints profile_length begin_x end_x begin_y end_y)]);
  push ( @{$tables{"baseclusters"}}, [qw(id multiplicon number_of_anchorpoints orientation was_twisted random_probability)]);
  push ( @{$tables{"anchorpoints"}}, [qw(id multiplicon basecluster gene_x gene_y coord_x coord_y is_real_anchorpoint)]);
  push ( @{$tables{"segments"}}, [qw(id multiplicon genome list first last order)]);
  push ( @{$tables{"list_elements"}}, [qw(id segment gene position orientation)]);
  
  #Go through all stored multiplicons
  foreach my $multiplicon ( @{ $$dataset{"evaluated_multiplicons"} } )
   {
    #Store the multiplicon data
    my @record=($multiplicon->get_index);
    ($multiplicon->level == 2)
     ? push( @record, $multiplicon->x_object->genome, $multiplicon->x_object->listname, "" )
     : push( @record, "", "", $multiplicon->x_object->multiplicon->get_index );
    push( @record, $multiplicon->y_list->genome, $multiplicon->y_list->listname);
    push( @record, $multiplicon->level, $multiplicon->number_of_anchorpoints, $multiplicon->profile->size );
    push( @record, $multiplicon->unmapped_coordinates );
    
    push(@{$tables{"multiplicons"}}, \@record );
    
    #Store data of the baseclusters of the multiplicon
    foreach my $basecluster ( $multiplicon->get_baseclusters )
     {
      push ( @{$tables{"baseclusters"}}, [ $basecluster->get_index, $multiplicon->get_index,
       $basecluster->number_of_anchorpoints, $basecluster->orientation, $basecluster->was_twisted,
       $basecluster->random_probability ] );
     } #foreach my $basecluster ( $multiplicon->get_baseclusters )

    #Store data of the anchorpoints
    foreach my $anchorpoint ( $multiplicon->get_anchorpoints ) 
     {
      push( @{$tables{"anchorpoints"}}, [ $anchorpoint->get_index, $multiplicon->get_index,
       $anchorpoint->basecluster->get_index, $anchorpoint->gene_x->ID, $anchorpoint->gene_y->ID, $anchorpoint->x, 
       $anchorpoint->y, $anchorpoint->is_real_anchorpoint ] );
     } #foreach my $anchorpoint ( $multiplicon->get_anchorpoints )

    #Store the segments of the multiplicon's profile
    my $order=0;
    foreach my $segment ( @{$multiplicon->profile->segments} )  
     {
      my($first, $last, $i);

      foreach ( @{$segment->remapped} )
       {
        if ( ref($_) eq 'iADHoRe::list_element' )
	 {
	  $first = $_->gene;
	  last;
	 } #if ( ref($_) eq 'iADHoRe::list_element' )
       } #foreach ( @{$segment->remapped} )

      $last=$first;
      $i=-1;

      #Store all list elements of each segment
      foreach my $element ( @{$segment->remapped} )
       {
        $i++;
        ( ref( $element ) eq 'iADHoRe::list_element' ) || next; #Skip gaps
        push( @{$tables{"list_elements"}}, [$element->get_index, $segment->get_index, $element->gene->ID, $i, 
	                                   $element->orientation] );
	#Determine the first and last genes of the segment on the original genelist				   
	($element->gene->coordinate < $first->coordinate) && ($first=$element->gene);
	($element->gene->coordinate > $last->coordinate) && ($last=$element->gene);				   
       } #foreach my $element ( @{$segment->remapped} )
      
      push( @{$tables{"segments"}}, [$segment->get_index, $multiplicon->get_index, $segment->genome, $segment->listname,
                                    $first->ID, $last->ID, $order] );
      $order++;
     } #foreach my $segment ( @{$multiplicon->profile->segments} )  

   } #foreach my $multiplicon ( @{ $$dataset{"evaluated_multiplicons"} } )
  
  #Finally, write out all tables as tab-delimited files.
  print STDERR "\nWriting output files";
  foreach my $table (keys %tables)
   {
    open ( OUT, ">".$settings->output_path."$table.txt" ) || die "Could not write output file $table.txt: $!\n";
    foreach ( @{$tables{$table}} )
     {
      print OUT join("\t", @$_ )."\n";
     } #foreach ( @{$tables{$table}} )
    close OUT;
    print STDERR ".";
   } #foreach my $table (keys %tables)
  
 } #sub output

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub by_multiplicon_size
# A sorting routine to sort multiplicons first by their number of 
# anchorpoints, then by DPD distance they span.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub by_multiplicon_size
 {
  ( $a->number_of_anchorpoints <=> $b->number_of_anchorpoints )
   ||
  (
   iADHoRe::ghm->dpd( $b->begin_x, $b->begin_y, $b->end_x, $b->end_y )
    <=>
   iADHoRe::ghm->dpd( $a->begin_x, $a->begin_y, $a->end_x, $a->end_y )
  ) 
 } #sub by_multiplicon_size

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub read_dataset
# Reads in the output from a finished i-ADHoRe-run and stores all the
# multiplicons in the evaluated_multiplicons array. These can then be
# accessed for postprocessing
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub read_dataset
 {
  my ($pkg)=$_[0];
  my (@fields);
  my (%dataset,%genelists,%genes, %tables, %keyfield, %profiles);
  
  $dataset{"multiplicons2evaluate"}=[];
  $dataset{"evaluated_multiplicons"}=[];

  #Read in the genes and gene lists.
  open( IN, $settings->output_path."genes.txt" ) || die "Could not open genes table: $!\n";
  @fields = split(/\t/, <IN> );
  chomp @fields;
  while( <IN> )
   {
    my ($gene,$element,$remapped,$unmapped);
    my %line;
    chomp;
    @line{ @fields }= split(/\t/);
    if ( !exists( $genelists{ $line{'genome'} }{ $line{'list'} } ) )
     {
      $genelists{ $line{'genome'} }{ $line{'list'} } = iADHoRe::genelist->new( @line{'list', 'genome'}, [], [], undef );
      push( @{$dataset{"genelists"}}, $genelists{ $line{'genome'} }{ $line{'list'} } );
     } #if ( !exists( $genelists{ $line{'genome'} }{ $line{'list'} } ) )
    $gene = iADHoRe::gene->new( @line{qw(id genome)}, $genelists{ $line{'genome'} }{ $line{'list'} },
     @line{qw(coordinate orientation remapped_coordinate is_tandem is_tandem_representative tandem_representative remapped)}, {} );
    $genes{ $line{'id'} }= $gene;
    $element = iADHoRe::list_element->new( $gene, $line{'orientation'}, 0 );
    $remapped = $genelists{ $line{'genome'} }{ $line{'list'} }->remapped;
    $unmapped = $genelists{ $line{'genome'} }{ $line{'list'} }->elements;
    $gene->remapped || ($$remapped[ $line{'remapped_coordinate'} ] = $element);
    $$unmapped[ $line{'coordinate'} ] = $element;
   } #while( <IN> )
  close IN;
  
  #Read in the other tables
  %keyfield = ( 'multiplicons' => 'id', baseclusters => 'multiplicon', 'anchorpoints' => 'basecluster',
                'segments' => 'multiplicon', 'list_elements' => 'segment' );
  foreach my $table ( qw(multiplicons baseclusters anchorpoints segments list_elements) )
   {
    my $key = $keyfield{$table};
    my %table;

    open( IN, $settings->output_path.$table.".txt" ) || die "Could not open $table table: $!\n";
    @fields = split(/\t/, <IN> );
    chomp @fields;
    while( <IN>) 
     {
      my %line;
      chomp;
      @line{ @fields }= split(/\t/);
      defined( $line{$key} ) || next;
      ( $table eq 'multiplicons' ) ? ( $table{ $line{$key} } = \%line ) : push( @{$table{ $line{$key} }}, \%line );
     } #while( <IN>) 
    close IN;
    
    $tables{$table} = \%table;
   } #foreach my $table ( qw(multiplicons baseclusters anchorpoints segments list_elements) )
  
  #Create all multiplicon objects and all objects associated to multiplicons
  foreach my $record ( sort { $$a{'level'} <=> $$b{'level'} } ( values %{$tables{'multiplicons'}} ) )
   {
    my ($multiplicon,$x_object,$y_list,$begin_x,$begin_y,$end_x,$end_y,$id,$is_redundant);
    
    $id = $$record{'id'};
    $is_redundant = $$record{'is_redundant'};
    
    if ( $$record{'level'} == 2 )
     {
      $x_object = $genelists{ $$record{'genome_x'} }{ $$record{'list_x'} };
      $begin_x = ${$x_object->elements}[ $$record{'begin_x'} ]->gene->remapped_coordinate;
      $end_x = ${$x_object->elements}[ $$record{'end_x'} ]->gene->remapped_coordinate;
     } #if ( $$record{'level'} == 2 )
    else 
     {
      $x_object = $profiles{ $$record{'parent'} };
      $begin_x = $$record{'begin_x'};
      $end_x = $$record{'end_x'};
     } #else 

    $y_list = $genelists{ $$record{'genome_y'} }{ $$record{'list_y'} };
    $begin_y = ${$y_list->elements}[ $$record{'begin_y'} ]->gene->remapped_coordinate;
    $end_y = ${$y_list->elements}[ $$record{'end_y'} ]->gene->remapped_coordinate;
    
    $multiplicon = iADHoRe::multiplicon->new( $x_object, $y_list, $begin_x, $end_x, $begin_y, $end_y, 
     $$record{'number_of_anchorpoints'}, $$record{'level'}, [], 0, undef, $id, $is_redundant );
     
    #Create the basecluster objects for the multiplicon
    foreach my $cluster_data  ( @{$tables{'baseclusters'}{$id}} )
     {
      my $basecluster = iADHoRe::basecluster->new( $x_object, $y_list, undef, undef, undef, undef,
       $$cluster_data{'number_of_anchorpoints'}, [], @{$cluster_data}{'random_probability', 'orientation'}, $multiplicon, undef,
       undef, undef, undef, undef, undef, undef, undef, undef, $$cluster_data{'was_twisted'}, $$cluster_data{'id'} );
       
      foreach my $anchorpoint_data ( @{$tables{'anchorpoints'}{ $$cluster_data{'id'} }} ) 
       {
        my $anchorpoint = iADHoRe::anchorpoint->new( $genes{ $$anchorpoint_data{'gene_x'} },
	 $genes{ $$anchorpoint_data{'gene_y'} }, @{$anchorpoint_data}{qw(coord_x coord_y is_real_anchorpoint)}, $basecluster,
	 $multiplicon, $$anchorpoint_data{'id'} );
	push( @{$basecluster->anchorpoints}, $anchorpoint );
       } #foreach my $anchorpoint_data ( @{$tables{'anchorpoints'}{ $$cluster_data{'id'} }} )
      
      $basecluster->begin_x( $basecluster->lowest_x );
      $basecluster->end_x( $basecluster->highest_x );
      $basecluster->begin_y( $basecluster->lowest_y );
      $basecluster->end_y( $basecluster->highest_y );
      
      push( @{$multiplicon->baseclusters}, $basecluster );
       
     } #foreach my $cluster_data  ( @{$tables{'baseclusters'}{$id}} )
    
    #Create the profile
    my $profile =  iADHoRe::profile->new( $multiplicon, -1, [] );
    foreach my $segment_data ( sort { $$a{'order'} <=> $$b{'order'} } ( @{$tables{'segments'}{$id}} ) )
     {
      my ( $segment, $remapped );
      $segment = iADHoRe::genelist->new( @{$segment_data}{'list', 'genome'}, [], [], $profile, $$segment_data{'id'} );
      $remapped = $segment->remapped;
      foreach my $element_data ( @{ $tables{'list_elements'}{ $$segment_data{'id'} } } )
       {
        $$remapped[ $$element_data{'position'} ]
	 = iADHoRe::list_element->new( $genes{ $$element_data{'gene'} }, $$element_data{'orientation'}, 0, $$element_data{'id'} );
       } #foreach my $element_data ( @{ $tables{'list_elements'}{ $$segment_data{'id'} } } )
      push( @{$profile->segments}, $segment );
     } #foreach my $segment_data ( sort { $$a{'order'} <=> $$b{'order'} } ( @{$tables{'segments'}{$id}} ) )
     
    $profiles{$id}= $profile;
    $multiplicon->profile( $profile );
    
    $multiplicon->number_of_baseclusters( scalar( @{$multiplicon->baseclusters} ) );
    
    push( @{$dataset{'evaluated_multiplicons'}}, $multiplicon );
   } #foreach my $record ( sort { $$a{'level'} <=> $$b{'level'} } ( values %{$table{'multiplicons'}} ) )
  
  bless \%dataset, $pkg;
  return \%dataset; 
  
 } #sub read_dataset

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub multiplicons
# Returns a list of all evaluated multiplicons.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub multiplicons
 {
  my $dataset=$_[0];
  
  return( @{$$dataset{'evaluated_multiplicons'}} );
 } #sub multiplicons

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub anchorpoints
# Returns a list of all anchorpoints in the dataset.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub anchorpoints
 {
  my $dataset=$_[0];
  my (@anchorpoints);
  
  foreach my $multiplicon ( $dataset->multiplicons )
   {
    push( @anchorpoints, $multiplicon->get_anchorpoints );
   } #foreach my $multiplicon ( $dataset->multiplicons )
  
  return( @anchorpoints ); 
 } #sub anchorpoints

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub genomes
# Returns a list of all genomes in the dataset.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub genomes
 {
  my $dataset=$_[0];
  my %genomes;
  
  foreach my $list ( $dataset->genelists )
   {
    $genomes{ $list->genome }=-1;
   } #foreach my $list ( $dataset->genelists )
  
  return( sort( keys %genomes ) );
  
 } #sub genomes

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_settings
# Returns the get_settings object.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_settings
 {
  return $settings;
 } #sub get_settings

1;
