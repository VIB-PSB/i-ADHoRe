package iADHoRe::ghm;
#Created by Cedric Simillion on Tue Jun 29 14:09:26 CEST 2004

use strict;
use iADHoRe::basecluster;
use iADHoRe::multiplicon;
use iADHoRe::settings;
use iADHoRe::anchorpoint;

my $settings=iADHoRe::settings->get_settings;
our $AUTOLOAD;

#===============================================================================
# Package methods
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub new
# Builds a GHM from 2 objects, given as arguments. The first object is
# called the x_object and is either a profile object or a genelist
# object. The second is always a genelist object and is called y_list. 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub new
 {
  my ($x_object,$y_list,$pkg,$ghm);
  ($pkg,$x_object,$y_list)=@_;
  
  $ghm={ "x_object" => $x_object, 
         "y_list" => $y_list,
	 "+", => {}, #A hash containing the positive points in the matrix 
	 "-", => {}, #A hash containing the negative points in the matrix
	 "number_of_points" => 0, #The number of points in the matrix
	 "area" => undef,
         "baseclusters" => { '+' => [], '-' => [] }, #A hash containing arrays to store baseclusters
	 "multiplicons" => [] };
	 
  bless $ghm,$pkg;

  $$ghm{"size_x"}=$x_object->size;
  $$ghm{"size_y"}=$y_list->size;

  #Fill the matrix
  if ( ref($x_object) eq "iADHoRe::genelist" )
   {
    $$ghm{"lists_identical"} =
     ( ($x_object->listname eq $y_list->listname) && ($x_object->genome eq $y_list->genome) && !defined($x_object->profile) ) ? -1 : 0;
    $$ghm{"level"}=2;
    $ghm->build_matrix( $x_object, $y_list, 1 );
   } #if ( $x_object->isa("genelist") )
  elsif( ref($x_object) eq "iADHoRe::profile" )
   {
    my $segments=$x_object->segments;
    foreach my $i ( 0 .. $#{$segments} )
     {
      $ghm->build_matrix( $$segments[$i], $y_list, $i+1 );
     } #foreach my $i ( 0 .. $#{$segments} )
    $$ghm{"level"}=$x_object->multiplicon->level+1;
    $$ghm{"lists_identical"} = 0
   } #elsif( ref($x_object) eq "iADHoRe::profile" )

  return $ghm;

 } #sub new

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub merge
# Merges a list of GHMs together. The x-objects of all GHMs must have
# the same size and the y_lists must be identical.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub merge
 {
  my ($package,$new_ghm,$x_index);
  my (@ghms);

  ($package,@ghms)=@_;
  
  $new_ghm={ "x_object" => [ map { $_->x_object } (@ghms) ],
       "y_list" => $ghms[0]->y_list,
       "+", => {}, #A hash containing the positive points in the matrix 
       "-", => {}, #A hash containing the negative points in the matrix
       "number_of_points" => 0, #The number of points in the matrix
       "level" => 0,
       "area" => undef };

  bless $new_ghm, $package;

  $x_index=0;
  
  foreach my $ghm ( @ghms )
   {
    foreach my $or ( '+', '-' )
     {
      foreach my $x ( keys %{$$ghm{$or}} )
       {
        foreach my $y ( keys %{$$ghm{$or}{$x}} )
	 {
	  exists( $$new_ghm{$or}{$x}{$y} ) || ( $$new_ghm{"number_of_points"}++ );
	  $$new_ghm{$or}{$x}{$y} = $$ghm{$or}{$x}{$y} + $x_index;
	 } #foreach my $y ( keys %{$$ghm{$or}{$x}} )
       } #foreach my $x ( keys %{$$ghm{$or}} )
     } #foreach my $or ( '+', '-' )
    
    $x_index += $ghm->level; 
    $$new_ghm{"level"}+= $ghm->level-1;
   } #foreach $ghm ( $ghm_a, $ghm_b )
  
  $$new_ghm{"level"}++;
  
  return $new_ghm;  

 } #sub merge

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub number_of_y_hits
# Returns the number of elements in the y_list that have a homolog
# with the x_object.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub number_of_y_hits
 {
  my $ghm = $_[0];
  
  unless( exists $$ghm{"number_of_y_hits"} )
   {
    my %y_hit;
    foreach my $or ( '+', '-' )
     {
      foreach my $x ( keys %{$$ghm{$or}} )
       {
        foreach my $y ( keys %{$$ghm{$or}{$x}} )
	 {
	  $y_hit{$y}=-1;
	 } #foreach my $y ( keys %{$$ghm{$or}{$x}} )
       } #foreach my $x ( keys %{$$ghm{$or}} )
     } #foreach my $or ( '+', '-' )
    $$ghm{"number_of_y_hits"} = scalar( keys %y_hit ); 
   } #unless( exists $$ghm{"number_of_y_hits"} )
  
  return $$ghm{"number_of_y_hits"};
  
 } #sub number_of_y_hits

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub build_matrix
# Fills in the actual coordinates of a GHM that is being created with
# the new method. The arguments passed are:
# - the x_list (the gene list that will mapped on the x-axis). This
#   can be an individual gene list or a segment part of a profile. To
#   construct a GHM between a profile and a gene list, a superposition
#   must be made of the GHMs between each segment in the profile and
#   the gene list.
# - the y_list (the gene list that will mapped on the y-axis)
# - the x_index. When the x_list passed is part of a profile, this
#   indicates which number this is. This should be 1 for the first
#   segment and 2,3,4,... for all additional segments. For single
#   x_lists, this should simply be 1.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub build_matrix
 {
  my ($ghm,$x_list,$y_list,$x_elements,$y_elements,$ghm_area,$x_index);
  ($ghm,$x_list,$y_list,$x_index)=@_;
  
  $x_elements=$x_list->remapped;
  $y_elements=$y_list->remapped;
  
  foreach my $x ( 0..$#{$x_elements} )
   {

    ( ref( $$x_elements[$x] ) eq 'iADHoRe::list_element') || next; #Check if the current position is not a gap.
    $$x_elements[$x]->has_pairs || next;
    
    my $max_y= $$ghm{"lists_identical"} ? ($x-1) : $#{$y_elements};

    foreach my $y ( $$x_elements[$x]->matching_positions( $y_elements ) )
     {
      ($y > $max_y) && next;
      my $or;
      ( $$x_elements[$x]->orientation eq $$y_elements[$y]->orientation ) ? ( $or='+' ) : ( $or='-');
      $$ghm{$or}{$x}{$y}=$x_index;
      $$ghm{"number_of_points"}++; #!!! This will count some anchor points from profiles multiple times. Leave as is???
     } #foreach my $y ( $$x_elements[$x]->matching_positions( $y_elements ) )

   } #foreach my $x ( 0..$#{$x_elements} )  
  
 } #sub build_matrix

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub adhore
# runs the basic ADHoRe algorithm on a ghm object. Returns a list of
# detected multiplicons
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub adhore
 {
  my ($ghm)=$_[0];

  ( $ghm->number_of_points < $settings->anchorpoints ) && return();
  
  $ghm->prepare_for_statistical_validation;
  
  foreach my $gap ( $settings->gap_sizes )
   {
    foreach my $orientation ( '+', '-' )
     {
      $ghm->seed_base_clusters( $gap, $orientation );
      $ghm->enrich_clusters( $gap, $orientation, $orientation );
      $ghm->join_clusters( $gap, $orientation );
     } #foreach my $orientation ( '+', '-' )
   } #foreach my $gap ( settings->gap_sizes )
  
  $ghm->enrich_clusters( $settings->gap_size, '+', '-' );
  $ghm->enrich_clusters( $settings->gap_size, '-', '+' );
  
  $ghm->filter_baseclusters;
  
  $ghm->cluster_clusters;
  
  return( $ghm->get_multiplicons );
  
 } #sub adhore

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub prepare_for_statistical_validation
# Calculates a number of properties of the GHM used in the statistical
# validation of ADHoRe clusters.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub prepare_for_statistical_validation
 {
  my ($ghm,$number_of_masked_rows);
  $ghm=$_[0];

  #Get the number of masked rows
  $number_of_masked_rows= ( $ghm->level > 2 ) ? $ghm->y_list->number_of_masked_elements : 0;

  #Calculate the area of the ghm
  $ghm->lists_identical
   ? $ghm->area( $ghm->size_x * ($ghm->size_x - 1) / 2 )
   : $ghm->area( $ghm->size_x * ($ghm->size_y - $number_of_masked_rows) );

 } #sub prepare_for_statistical_validation

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub seed_base_clusters
# Creates initial base clusters from a GHM object, with the gap size
# and orientation class ('+' or '-') supplied as argument. If the
# detected clusters contain 3 or more anchor points and their
# r_squared is equal to or higher then the q_value specified by the
# user, the cluster is retained. Note that checking if a cluster
# contains the user specified minimum number of anchor points is done
# at a later stage.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub seed_base_clusters
 {
  my ($ghm,$gap,$orientation,$sign);

  ($ghm,$gap,$orientation)=@_;
  
  ($orientation eq "-") ? ($sign = -1) : ($sign=1);
  
  #Try to build clusters starting from each point in the current orientation class of the ghm
  foreach my $start_x ( sort { $a <=> $b} (keys %{$$ghm{$orientation}}) )
   {
    foreach my $start_y ( sort { ($a <=> $b) * $sign } (keys %{$$ghm{$orientation}{$start_x}}) )
     {
      my ($cluster, $ref_x, $ref_y);
      
      #Initiate a base cluster with the current starting point ($start_x, $start_y) as first anchor point
      $cluster=iADHoRe::basecluster->new
        ( $ghm->x_object, $ghm->y_list, undef, undef, undef, undef, 0, [], undef, $orientation, undef, undef, undef, undef, undef,
	  undef, undef, undef, undef, undef, 0 );
      $cluster->add_anchorpoint( $start_x, $start_y );
      
      #The starting point is the first reference point ($ref_x, $ref_y)
      ($ref_x, $ref_y) = ($start_x, $start_y);
      
      #Look for a next anchor point that lies within the designated gap distance from the reference point
      while ( defined $ref_x )
       {
        
	my ($closest_distance, $closest_x, $closest_y);
	$closest_distance=$gap + 1;
	
	foreach my $x ( $ref_x+1 .. $ref_x+$gap )
	 {

          exists( $$ghm{$orientation}{$x} ) || next;

          for ( my $y = $ref_y + $sign; abs($ref_y - $y) <= $gap; $y+= $sign )
	   {

            exists( $$ghm{$orientation}{$x}{$y} ) || next;

	    my $distance = iADHoRe::ghm->dpd( $ref_x, $ref_y, $x, $y );
            #If more than one point lies within a distance $gap of the reference point, the closest one is chosen.
	    if ( $distance < $closest_distance )
	     {
	      $closest_distance=$distance;
	      $closest_x=$x;
	      $closest_y=$y
	     } #if ( $distance < $closest_distance )

	   } #for ( my $y = $ref_y + $sign; $y+= $sign; $y == $ref_y + ($gap*$sign) )

	 } #foreach my $x ( $ref_x+1 .. $ref_x+$gap )

        #If a new anchor point is found, it is added to the cluster and designated as the next reference point
        if ( defined $closest_x )
	 {
	  $cluster->add_anchorpoint( $closest_x, $closest_y );
	  ($ref_x,$ref_y) = ($closest_x, $closest_y);
	 } #if ( defined $closest_x )
        else
	 {
	  $ref_x = undef;
	 } #else
	 
       } #while ( defined $ref_x )

      #If a cluster contains 3 or more anchor points (for now) and its r_squared value is higher than or equal to the
      #q_value specified by the user, it is retained.
      if ( ( $cluster->number_of_anchorpoints >= 3 ) && ( $cluster->r_squared >= $settings->q_value ) )
       {
        $ghm->add_cluster( $cluster );
	$cluster->confidence_interval;
       } #if ( $cluster->number_of_anchorpoints >= 3 )
      else 
       {
        #Since the anchorpoints objects of a cluster object also point to the cluster object itself,
	#the reference count to the cluster object will not drop to 0 once $cluster is out of scope.
	#Therefore, we must destroy the $cluster object manually.
        $cluster->DESTROY;
       } #else 
     } #foreach my $y ( sort { $a <=> $b } (keys %{$$ghm{$orientation}{$x}}) )
   } #foreach my $x ( sort { $a <=> $b } (keys %{$$ghm{$orientation}}) )
  
 } #sub seed_base_clusters

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub enrich_clusters
# Adds anchor points of a given orientation class to clusterset of a
# given orientation class, using a specified gap size.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub enrich_clusters
 {
  my ($ghm,$cluster_orientation,$ghm_orientation,$gap,$point_found);
  my (@clusters_2_evaluate);

  ($ghm,$gap,$cluster_orientation,$ghm_orientation)=@_;
  @clusters_2_evaluate=$ghm->get_baseclusters($cluster_orientation);

  do
   {
    my @keep_4_next_iteration;
   
    $point_found=0;
    
    #Check for every point...
    foreach my $x ( sort { $a <=> $b } (keys %{$$ghm{$ghm_orientation}}) )
     {
      foreach my $y ( sort { $a <=> $b } (keys %{$$ghm{$ghm_orientation}{$x}}) )
       {
        my ($closest_cluster,$closest_distance);
	$closest_distance = $gap+1;
	
        #and cluster in the GHM if...
	foreach my $cluster (@clusters_2_evaluate)
	 {
	  $cluster->confidence_interval;
	  my $distance=$cluster->distance_to_point($x,$y);
	  if ( 
	       ( $distance < $closest_distance ) #the point lies within the designated distance of the cluster, ...
	      &&
	       ( $cluster->r_squared($x,$y) >= $settings->q_value ) #the r_squared is high enough when the point is added...
	      &&
	       ( $cluster->in_interval($x,$y) ) #and if the point lies within the confidence interval of the cluster.
	     )
	   {
	    $closest_distance=$distance;
	    $closest_cluster=$cluster;
	   } #if ...  
	 } #foreach my $cluster (@clusters_2_evaluate)

        if ( defined $closest_cluster ) #If such a point and cluster are found, add the point to the cluster...
	 {
	  $closest_cluster->add_anchorpoint($x,$y);
	  $closest_cluster->confidence_interval;
	  delete $$ghm{$ghm_orientation}{$x}{$y};
	  push (@keep_4_next_iteration, $closest_cluster); #and keep the cluster for the next iteration.
	  $point_found=-1;
	 } #if ( defined $closest_cluster)

       } #foreach my $y ( sort { $a <=> $b } (keys %{$$ghm{$ghm_orientation}{$x}}) )
     } #foreach my $x ( sort { $a <=> $b } (keys %{$$ghm{$ghm_orientation}}) )

    @clusters_2_evaluate = @keep_4_next_iteration;

   } while ($point_found);

 } #sub enrich_clusters

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub join_clusters
# Merges pairs of clusters of the same orientation class that are
# within a given distance from each other by adding all anchorpoints
# of one cluster to the other and deleting the first cluster. The
# distance and orientation class are passed as arguments
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub join_clusters
 {
  my ($ghm,$gap,$orientation,$join_flag,$baseclusters);
  ($ghm,$gap,$orientation)=@_;
  
  $baseclusters=$$ghm{"baseclusters"}{$orientation};
  
  do
   {
    my ($i,$j)=$ghm->get_closest_clusters($gap,$orientation);
    $join_flag=0;
    
    if (defined $i) 
     {
      $$baseclusters[$i]->merge_with( $$baseclusters[$j] );
      splice (@$baseclusters, $j,1);
      $join_flag=-1;
     } #if (defined $closest_i) 

   } while ($join_flag);
  
 } #sub join_clusters

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_closest_clusters
# For all pairs of clusters that are within a specified gap distance,
# the function returns the pair of clusters that are the closest to
# each other. If only one orientation is passed as argument, only the
# corresponding orientation class is searched. If two orientations are
# given, the pairs evaluated will contain one cluster from each class.
# The values returned are the indices of each cluster in their
# respective orientation classes. If one orientation class was given,
# both indices come from the same orientation class, if not the first
# index comes from the first orientation class specified, the second
# index from the second class.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_closest_clusters
 {
  my ($ghm,$gap,$ref_orientation,$comp_orientation,$clusters_i,$clusters_j,$sameset,$closest_i,$closest_j,
      $closest_distance);
  ($ghm,$gap,$ref_orientation,$comp_orientation)=@_;
  
  $clusters_i = $$ghm{"baseclusters"}{$ref_orientation};

  #Figure it if two different orientation classes were given. If $comp_orientation was not given,
  #the second set of clusters to be searched will be same as the first.
  if ( defined($comp_orientation) && ($comp_orientation ne $ref_orientation) )
   {
    $clusters_j = $$ghm{"baseclusters"}{$comp_orientation};
    $sameset = 0;
   } #if ( defined $comp_orientation )
  else
   {
    $clusters_j = $$ghm{"baseclusters"}{$ref_orientation};
    $sameset = -1;
   } #else

  $closest_distance=$gap+1;

  #Detect the closest clusters, if any.
  for ( my $i=0; $i < $#{$clusters_i}; $i++ )
   {
    for ( my $j = $sameset ? ( $i+1 ) : 0; $j <= $#{$clusters_j}; $j++ )
     {
      #Skip this iteration if cluster i and cluster j are assigned to the same multiplicon
      ( defined( $$clusters_i[$i]->multiplicon ) && defined( $$clusters_j[$j]->multiplicon )  )
       &&
      ( $$clusters_i[$i]->multiplicon == $$clusters_j[$j]->multiplicon )
       &&
      next; 

      my $distance = $$clusters_i[$i]->distance_to_cluster( $$clusters_j[$j] );

      #If the distance is zero, check if one of the clusters lies completely in the confidence interval of the other
      ($distance == 0)
       &&
      ( $$clusters_i[$i]->overlapping_interval( $$clusters_j[$j] ) || ( $distance = $closest_distance + 1 ) );

      if (
	   ( $distance < $closest_distance )
	    &&
	   ( $$clusters_i[$i]->r_squared( $$clusters_j[$j] ) >= $settings->q_value )
	 )
       {
	$closest_i=$i;
	$closest_j=$j;
	$closest_distance=$distance;
       } #if ( ... )

     } #for ( my $j = $samesets ? ( $i+1 ) : 0; $j <= $#{$clusters_j}; $j++ )
   } #for ( my $i=0; $i < scalar(@$baseclusters)-1; $i++ )

  return ( $closest_i, $closest_j );

 } #sub get_closest_clusters

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub filter_baseclusters
# Performs statistical filtering on all baseclusters of a ghm object.
# Clusters with a probability higher then the cutoff specified in the
# settings file are removed. The value of the probability is stored
# for the remaining clusters.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub filter_baseclusters
 {
  my ($ghm,$cutoff,$baseclusters);
  
  $ghm=$_[0];
  
  #Transform the cutoff so it can be compared with the value returned by the calculate_probability method
  $cutoff = log( 1- $settings->prob_cutoff );
  
  $baseclusters=$ghm->baseclusters;
  
  foreach my $or ( '+', '-' )
   {
    for ( my $i=0; $i <= $#{ $$baseclusters{$or} }; $i++ )
     {
      my $log_one_minus_prob =
       $$baseclusters{$or}[$i]->calculate_probability( $ghm->area, $ghm->number_of_points, $ghm->level );
     
      if ( $log_one_minus_prob  < $cutoff )
       {
        my $removed_cluster=splice( @{ $$baseclusters{$or} }, $i, 1 );
#	$removed_cluster->DESTROY;
	$i--;
       } #if ( $log_one_minus_prob  < $cutoff )
      else 
       {
        #If the cluster is retained (i.e. its probability is below the cutoff), store its probability
        $$baseclusters{$or}[$i]->random_probability( 1 - exp( $log_one_minus_prob ) );
       } #else 
       
     } #for ( my $i=0; $i <= $#{ $$baseclusters{$or} }; $i++ )
   } #foreach my $or ( '+', '-' )
  
 } #sub filter_clusters

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub cluster_clusters
# Performs metaclustering, i.e., grouping baseclusters together. This
# is done by first trying to merge clusters from one orientation class
# with clusters from the opposite class by 'twisting' them (mirroring
# around their central y-value) and try if they can be merged together
# just as is done in the join_clusters method. Next, clusters of the
# same class are also grouped together if they are close enough to
# each other.
# These groups of baseclusters are called metaclusters and are stored.
# in multiplicon objects. Thus, this function generates the
# multiplicon objects. The terms metacluster and multiplicon are used
# interchangebly
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub cluster_clusters
 {
  my ($ghm,$join_flag,$cluster_gap);
  
  $ghm=$_[0];
  $cluster_gap=$settings->cluster_gap;
  
  do
   {
    $join_flag=0;
    
    foreach ( ['+', '-'], ['+','+'], ['-', '+'], ['-', '-'] )
     {
      my ($x, $y, $i, $j,$was_twisted);

      ($x, $y)= @$_; #x and y are the two orientation classes in this iteration

      ($x ne $y) && $ghm->twist_clusters($y); #Twist clusters from the opposite orientation class
      ($x ne $y) ? ($was_twisted=-1) : ($was_twisted=0); #Set the flag if clusters from class y were twisted
      ($i,$j)=$ghm->get_closest_clusters($cluster_gap,$x,$y);  #Get any clusters that could be grouped
      ($x ne $y) && $ghm->twist_clusters($y); #Untwist clusters from the opposite orientation class
      
      if (defined $i) #Check if there are two clusters to be grouped.
       {
        my ($cluster_i, $cluster_j);
	$cluster_i = $$ghm{"baseclusters"}{$x}[$i];
	$cluster_j = $$ghm{"baseclusters"}{$y}[$j];
	
	( $was_twisted && $cluster_i->was_twisted ) && next;
	
	$cluster_j->was_twisted( $was_twisted );
        $join_flag=-1;
	
        #Check if both i and j are already assigned to a metacluster
        if ( defined( $cluster_i->multiplicon ) && defined( $cluster_j->multiplicon ) ) 
	 {
          #If so, join these two metaclusters:
	  my $old_multiplicon_of_j= $cluster_j->multiplicon;
	  #First, add all clusters of the multiplicon of cluster_j to the multiplicon of cluster_i
          $cluster_i->multiplicon->add_baseclusters( $cluster_j->multiplicon->get_baseclusters );
	  #Next, empty the old multiplicon of cluster_j
	  $old_multiplicon_of_j->empty;
	 } #if ( defined( $metacluster_of{$i} ) && defined( $metacluster_of{$j} ) )

        #If not, check if only one cluster has been asigned to a metacluster
	elsif ( defined( $cluster_i->multiplicon ) || defined( $cluster_j->multiplicon ) ) 
	 {
          #Make sure cluster i is the one with its metacluster defined
	  defined( $cluster_i->multiplicon ) || ( ($cluster_i,$cluster_j) = ($cluster_j,$cluster_i) );
          #Add cluster j to the metacluster of cluster i
	  $cluster_i->multiplicon->add_baseclusters( $cluster_j );
	 } #elsif ( defined( $metacluster_of{$i} ) || defined( $metacluster_of{$j} ) )

        #If both clusters are not assigned to a metacluster yet, create a new one containing both clusters
	else
	 {
	  my $multiplicon = iADHoRe::multiplicon->new
	   ( $ghm->x_object, $ghm->y_list, undef, undef, undef, undef, 0, $ghm->level, [], 0 );
	  $ghm->add_cluster( $multiplicon );
	  $multiplicon->add_baseclusters( $cluster_i, $cluster_j );
	 } #else
	 
       } #if (defined $i)

     } #foreach ( ['+', '-'], ['+','+'], ['-', '+'], ['-', '-'] )

   } while ($join_flag);

  #Check the configuration of any obtained metaclusters
  $ghm->check_metacluster_configuration;
 
  #Assign all baseclusters that were not assigned to a metacluster to their own multiplicon
  foreach my $basecluster ( $ghm->get_baseclusters )
   {
    ref( $basecluster->multiplicon ) && next;
    my $multiplicon = iADHoRe::multiplicon->new
     ( $ghm->x_object, $ghm->y_list, undef, undef, undef, undef, 0, $$ghm{"level"}, [], 0 );
    $basecluster->was_twisted( 0 );
    $multiplicon->add_baseclusters( $basecluster );
    $ghm->add_cluster( $multiplicon );
   } #foreach my $baseclusters ( $ghm->get_baseclusters )
  
  #Remove empty multiplicons from the ghm
  for (my $m=0; $m < scalar( @{ $$ghm{"multiplicons"} } ); $m++)
   {
    if (
         ( $$ghm{"multiplicons"}[$m]->number_of_baseclusters == 0 )
	||
	 ( $$ghm{"multiplicons"}[$m]->number_of_anchorpoints < $settings->anchorpoints )
       )
     {
      splice( @{ $$ghm{"multiplicons"} }, $m, 1 );
      $m--;
     } #if ( ... )
   } #for (my $m=0; $m < scalar( @{ $$ghm{"multiplicons"} } ); $m++)
   
 } #sub cluster_clusters

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub check_metacluster_configuration
# This checks for all metaclusters (i.e. multiplicons) if the
# configuration of the individual basecluster still allows the detect
# a possible alignment path later on.
# This is done by checking if there are any two clusters of the same
# orientation that have overlapping coordinates but are not located
# within each other's intervals. 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub check_metacluster_configuration
 {
  my ($ghm);
  $ghm = $_[0];
  
  foreach my $multiplicon ( @{$ghm->multiplicons} )
   {
    my ( $baseclusters, $configuration_ok );
    $baseclusters = $multiplicon->baseclusters;
    $configuration_ok = -1;
    
    foreach my $i ( 0 .. $#{$baseclusters} - 1 )
     {
      my $cluster_i = $$baseclusters[$i];
      foreach my $j ( $i + 1 .. $#{$baseclusters} )
       {
        my $cluster_j = $$baseclusters[$j];
	
        ( $cluster_i->orientation eq $cluster_j->orientation ) || next;
	if ( $cluster_i->overlapping_coordinates( $cluster_j ) && !$cluster_i->overlapping_interval( $cluster_j ) )
	 {
	  $configuration_ok = 0;
	  last;
	 } #if ( $cluster_i->overlapping_coordinates( $cluster_j ) && !$cluster_i->overlapping_interval( $cluster_j ) )
       } #foreach my $j ( $i + 1 .. $#{$baseclusters} )
      $configuration_ok || last; 
     } #foreach my $i ( 0 .. $#{$baseclusters} - 1 )

    if( !$configuration_ok )
     {
      foreach my $cluster ( @$baseclusters )
       {
        $cluster->multiplicon( 0 );
       } #foreach my $cluster ( @$baseclusters )
      $multiplicon->empty; 
     } #if( !$configuration_ok )

   } #foreach my $multiplicon ( $ghm->multiplicons )
  
 } #sub check_metacluster_configuration

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub twist_clusters
# Mirrors all clusters of a given orientation class around their
# central y-value. !!!This changes the y-values of the individual
# anchorpoints!!!
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub twist_clusters
 {
  my ($ghm,$orientation)=@_;
  foreach my $cluster ( $ghm->get_baseclusters($orientation) )
   {
    $cluster->twist_cluster;
   } #foreach $cluster ( $ghm->get_baseclusters($orientation) )
 } #sub twist_clusters

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub add_cluster
# Adds a basecluster or multiplicon object passed as argument to a ghm
# object
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub add_cluster
 {
  my ($ghm,$cluster)=@_;

  if (ref $cluster eq 'iADHoRe::basecluster')
   {
    my $orientation = $cluster->orientation;

    push (@{$$ghm{"baseclusters"}{$orientation}}, $cluster);

    #Points in the matrix that are added to a cluster are removed from the matrix
    foreach my $anchorpoint ($cluster->get_anchorpoints)
     {
      my ($x,$y)=($anchorpoint->x, $anchorpoint->y);
      delete $$ghm{$orientation}{$x}{$y};
     } #foreach my $anchorpoint ($cluster->get_anchorpoints)

   } #if (ref $cluster eq 'basecluster')
  elsif (ref $cluster eq 'iADHoRe::multiplicon')
   {
    push (@{$$ghm{"multiplicons"}}, $cluster);
   } #if (ref $cluster eq 'base_cluster') 
   
 } #sub add_cluster

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_baseclusters
# If no orientation ('+' or '-') is given as argument, returns a list
# of all baseclusters of both orientation classes. Otherwise, a list
# of baseclusters of the specified orientation class is returned.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_baseclusters
 {
  my ($ghm,$orientation)=@_;
  
  if ( defined $orientation )
   {
    return ( @{ $$ghm{"baseclusters"}{$orientation} } )
   } #if ( defined $orientation )
  else 
   {
    return ( @{ $$ghm{"baseclusters"}{'+'} }, @{ $$ghm{"baseclusters"}{'-'} } )
   } #else 
 } #sub get_baseclusters

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_multiplicons
# Returns a list of all multiplicons (metaclusters) detected in the
# in a GHM. Before the list is returned, some additional data of the
# multiplicon, its baseclusters and its anchorpoints need to be filled
# in.
# Note: still needs to be updated to able to handle GHMs where the
# x-object is a profile.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_multiplicons
 {
  my ($ghm,$x_elements,$y_elements);

  $ghm=$_[0];
  $x_elements = ( ref( $ghm->x_object ) eq "iADHoRe::genelist" ) ? $ghm->x_object->remapped : $ghm->x_object->matrix;
  $y_elements = $ghm->y_list->remapped;

  foreach my $multiplicon ( @{$$ghm{"multiplicons"}} )
   {
    $multiplicon->set_bounds;

    foreach my $basecluster ( $multiplicon->get_baseclusters )
     {
      my $anchorpoint_list= $basecluster->anchorpoints;
      my (@false_anchorpoints);
      $basecluster->set_bounds;
      
      foreach my $anchorpoint ( @$anchorpoint_list )
       {
        my ($x,$y) = ( $anchorpoint->x, $anchorpoint->y );
        if ( ref( $ghm->x_object ) eq "iADHoRe::genelist" )
	 {
	  #Process the anchor point if the GHM is a 'simple' (genelist vs genelist) GHM
	  $anchorpoint->gene_x( $$x_elements[$x]->gene );
	  $anchorpoint->gene_y( $$y_elements[$y]->gene );
	 } #if ( ref( $ghm->x_object ) eq "iADHoRe::genelist" )
	else
	 {
	  #Deal with 'multiple' (profile vs genelist) GHMs. In such a GHM, an anchorpoint can match with 
	  #multiple positions in a column. If this is the case, we create a set of 'false' (real_anchorpoint
	  #flag set to false) anchorpoints for each additional match.
	  my $first_match_encountered=0;
	  foreach my $list ( @$x_elements )
	   {
	    ( ref( $$list[$x] ) eq 'iADHoRe::list_element' ) || next;
	    if ( $$list[$x]->is_pair_with( $$y_elements[$y] ) )
	     {
	      if ( $first_match_encountered )
	       {
	        my $false_anchorpoint=iADHoRe::anchorpoint->new
		 ( $$list[$x]->gene, $$y_elements[$y]->gene, $x, $y, 0, $basecluster, $multiplicon );
		push( @false_anchorpoints, $false_anchorpoint );   
	       } #if ( $first_match_encountered )
	      else
	       {
	        $anchorpoint->gene_x( $$list[$x]->gene );
		$anchorpoint->gene_y( $$y_elements[$y]->gene );
	        $first_match_encountered=-1;
	       } #else
	     } #if ( $$list[$x]->is_pair_with( $$y_elements[$y] ) )
	   } #foreach my $list ( @$elements )
	  $first_match_encountered || die "Anchorpoint without matching x-gene!!!"; 
	 } #else 
       } #foreach my $anchorpoint ( $anchorpoint_list )
      
      push( @$anchorpoint_list, @false_anchorpoints );
      
     } #foreach my $basecluster ( $multiplicon->get_baseclusters )

   } #foreach my $multiplicon ( @{$$ghm{"multiplicons"}} )

  return( @{$$ghm{"multiplicons"}} );
   
 } #sub get_multiplicons

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub identity
# Returns true if for any orientation class there is a value in the
# cell at a given set of coordinates x,y in a GHM.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub identity
 {
  my ($ghm,$x,$y,$identity);
  ($ghm,$x,$y)=@_;
  
  foreach my $or ("+","-")
   {
    if (exists $$ghm{$or}{$x})
     {
      if (exists $$ghm{$or}{$x}{$y})
       {
        $identity= $$ghm{$or}{$x}{$y};
	last;
       } #if (exists $$ghm{$or}{$x}{$y})
     } #if (exists $$ghm{$or}{$x})
   } #foreach my $or ("+","-")
  
  return( $identity );
   
 } #sub identity

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub dpd
# Returns the discrete pseudo distance between two points.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub dpd
 {
  shift(@_);
  my ($x1,$y1,$x2,$y2)=@_;
  
  ( abs( $x1-$x2 ) < abs( $y1-$y2 ) ) && ( ($x1,$x2,$y1,$y2) = ($y1,$y2,$x1,$x2) );
  
  return ( 2 * abs( $x1-$x2 ) - abs( $y1-$y2 ) );
  
 } #sub dpd

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub plot_ghm
# Creates a .png visualisation of a GHM. For debugging purposes.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub plot_ghm
 {
  use GD;
  my ($ghm,$img,$filename);
  my (%color);
  
  $ghm=$_[0];
  $filename = defined( $_[1] ) ? $_[1] : "ghm.png";
  
  
  $img = new GD::Image( $ghm->size_x, $ghm->size_y );
  $color{"black"}=$img->colorAllocate(0,0,0);
  $color{'+'}{'1'}=$img->colorAllocate(255,0,0);
  $color{'-'}{'1'}=$img->colorAllocate(127,0,0);
  $color{'+'}{'2'}=$img->colorAllocate(255,153,0);
  $color{'-'}{'2'}=$img->colorAllocate(127,76,0);
  $color{'+'}{'3'}=$img->colorAllocate(203,255,0);
  $color{'-'}{'3'}=$img->colorAllocate(101,127,0);
  $color{'+'}{'4'}=$img->colorAllocate(50,255,0);
  $color{'-'}{'4'}=$img->colorAllocate(25,127,0);
  $color{'+'}{'5'}=$img->colorAllocate(0,255,102);
  $color{'-'}{'5'}=$img->colorAllocate(0,127,51);
  $color{'+'}{'6'}=$img->colorAllocate(0,255,255);
  $color{'-'}{'6'}=$img->colorAllocate(0,127,127);
  $color{'+'}{'7'}=$img->colorAllocate(0,101,255);
  $color{'-'}{'7'}=$img->colorAllocate(0,50,127);
  $color{'+'}{'8'}=$img->colorAllocate(50,0,255);
  $color{'-'}{'8'}=$img->colorAllocate(25,0,127);
  $color{'+'}{'9'}=$img->colorAllocate(204,0,255);
  $color{'-'}{'9'}=$img->colorAllocate(102,0,127);
  $color{'+'}{'10'}=$img->colorAllocate(255,0,153);
  $color{'-'}{'10'}=$img->colorAllocate(127,0,76);
  $color{"white"}=$img->colorAllocate(255,255,255);
  $color{"gray+"}=$img->colorAllocate(170,170,170);
  $color{"gray"}=$img->colorAllocate(128,128,128);
  $color{"gray-"}=$img->colorAllocate(85,85,85);

  #Mark masked segments
  if ( $ghm->level > 2 )
   {
    my $y_elements=$ghm->y_list->remapped;
    foreach my $y (0 .. $#{$y_elements} )
     {
      $$y_elements[$y]->masked && $img->line(0,$y,$ghm->size_x-1,$y,$color{"gray-"});
     } #foreach my $y (@$y_elements)
    } #if ( $ghm->level > 2 )

  #Draw a rectangle around each multiplicon (metacluster) in the ghm
  foreach my $multiplicon ( @{$$ghm{"multiplicons"}} )
   {
    $img->rectangle( $multiplicon->lowest_x, $multiplicon->lowest_y, $multiplicon->highest_x, $multiplicon->highest_y,
                     $color{"gray-"} );
   } #foreach my $multiplicon ( @{$$ghm{"multiplicons"}} )

  foreach my $or ("+", "-")
   {
   
    #Plot each basecluster in the ghm
    foreach my $cluster ( @{$$ghm{"baseclusters"}{$or}} )
     {
      
      #The best-fit line
      $img->line( $cluster->x_end1, $cluster->y_end1, $cluster->x_end2, $cluster->y_end2, $color{"gray"} );
      
      #The confidence interval
      foreach my $x ( $cluster->x_end1 - $settings->gap_size .. $cluster->x_end2 + $settings->gap_size )
       {
        my ($upper,$lower)=$cluster->interval_bounds($x);
	$img->setPixel( $x, $upper, $color{"gray+"} );
	$img->setPixel( $x, $lower, $color{"gray+"} );
       } #foreach my $x ( $cluster->x_end1 - $settings->gap_size .. $cluster->x_end2 + $settings->gap_size )
      
      #And the individual anchorpoints
      foreach my $anchorpoint ( $cluster->get_anchorpoints )
       {
        $img->setPixel( $anchorpoint->x, $anchorpoint->y, $color{"white"} );
       } #foreach my $anchorpoint ( $cluster->get_anchorpoints )
       
     } #foreach my $cluster ( @{$$ghm{"baseclusters"}{$or}} )
     
    #Plot all singletons.
    foreach my $x ( keys %{$$ghm{$or}} )
     {
      foreach my $y ( keys %{$$ghm{$or}{$x}} )
       {
        $img->setPixel($x,$y,$color{$or}{ $$ghm{$or}{$x}{$y} } );
       } #foreach my $y ( keys %{$$ghm{$or}{$x}} )
     } #foreach my $x ( keys %{$$ghm{$or}} )
     
   } #foreach my $or ("+", "-")

  open (PNG,">$filename");
  binmode PNG;
  print PNG $img->png;
  close PNG;
  
  return ("$filename written at ".localtime);
  
 } #sub plot_ghm

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub AUTOLOAD
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub AUTOLOAD
 {
  my($object, $setvalue, $attribute);

  $object=$_[0];
  $setvalue=$_[1];
  ($attribute)= ($AUTOLOAD=~ /::([^:]+)$/);
  
  ($attribute eq 'DESTROY') && return;
  
  exists( $$object{$attribute} ) || die "Attribute $attribute not defined for ".ref($object).".\n";

  defined($setvalue) && ( $$object{$attribute}=$setvalue );
  return ( $$object{$attribute} );
  
 } #sub AUTOLOAD

1;
