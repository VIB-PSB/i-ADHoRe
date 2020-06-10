package iADHoRe::basecluster;
#Created by Cedric Simillion on Mon Jul  5 20:04:02 CEST 2004

use strict;
use Statistics::Distributions;
use iADHoRe::record_structure;
use iADHoRe::cluster;
use iADHoRe::ghm;

our @ISA=qw(iADHoRe::cluster iADHoRe::record_structure);

our $structure=iADHoRe::basecluster->base_upon("iADHoRe::cluster",qw(anchorpoints random_probability orientation multiplicon a b avg_x 
                                                      var_x variance x_end1 y_end1 x_end2 y_end2 was_twisted id) );

#===============================================================================
# Package methods
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub add_anchorpoint
# Adds an anchor point to a basecluster object. The x and y
# coordinates of the anchorpoint are passed as arguments
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub add_anchorpoint
 {
  my ($cluster,$x,$y,$anchorpoint,$anchorpoint_array,$number_of_anchorpoints);
  ($cluster,$x,$y)=@_;
  
  $anchorpoint=iADHoRe::anchorpoint->new(undef, undef, $x, $y, -1, $cluster, undef );
  $anchorpoint_array=$cluster->anchorpoints;
  
  push( @$anchorpoint_array, $anchorpoint );
  
  $number_of_anchorpoints = $cluster->number_of_anchorpoints;
  $cluster->number_of_anchorpoints( $number_of_anchorpoints+1 );
  $cluster->b(0); #Set the coefficient of the best fit line to 0 to indicate that the statistics need to be updated
      
 } #sub add_anchorpoint

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_anchorpoints
# Returns a list of all anchorpoint objects of a cluster
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_anchorpoints
 {
  my ($cluster,$anchor_point_array);
  
  $cluster=$_[0];
  $anchor_point_array=$cluster->anchorpoints;
  
  return ( @$anchor_point_array );
 } #sub get_anchorpoints

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub r_squared
# Calculates the squared Pearson value for all x and y-coordinates
# of all anchor points in a cluster when no additional argument is
# passed. If additional x and y coordinates or another cluster are
# passed as arguments then the r_squared is calculated for the
# cluster with the additional coordinates included or with the
# anchorpoints of the other cluster included.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub r_squared
 {
  my ($cluster,$n,$sum_x,$sum_y,$sum_xy,$sum_x2,$sum_y2,$r);
  my (@values);
  
  $cluster=shift(@_);
  $n=0;
  
  #Add the x and y values in pairs to the @values array
  foreach my $anchorpoint ($cluster->get_anchorpoints)
   {
    push( @values, $anchorpoint->x);
    push( @values, $anchorpoint->y);
   } #foreach my $anchorpoint ($cluster->get_anchorpoints)
  

  #If another cluster is given as argument, add its x and y values as well to the @values array
  if ( defined($_[0]) &&  ( ref($_[0]) eq 'iADHoRe::basecluster' ) )
   {
    my $other_cluster=$_[0];
    foreach my $anchorpoint ($other_cluster->get_anchorpoints)
     {
      push( @values, $anchorpoint->x);
      push( @values, $anchorpoint->y);
     } #foreach my $anchorpoint ($other_cluster->get_anchorpoints)
   } #if ( defined($_[0]) &&  ( ref($_[0]) eq 'basecluster' ) )
  else #If not, add the contents of @_ to @values
   {
    push( @values, @_ );
   } #else
  
  while ( scalar(@values) )
   {
    my ($x,$y);
    $x=shift(@values);
    $y=shift(@values);
    $sum_x+=$x; $sum_y+=$y;
    $sum_xy+=$x*$y;
    $sum_x2+=$x**2;
    $sum_y2+=$y**2;
    $n++;
   } #foreach ($cluster->get_anchorpoints)

  if ( ( ($sum_x2-$sum_x**2/$n)*($sum_y2-$sum_y**2/$n) ) == 0 )
   {
    return -1
   } #if ((($sum_x2-$sum_x**2/$n)*($sum_y2-$sum_y**2/$n))==0)
  else
   {
    $r= ( $sum_xy - ($sum_x*$sum_y/$n) ) / sqrt( ($sum_x2-$sum_x**2/$n) * ($sum_y2-$sum_y**2/$n) );
    return ($r**2)
   } #else
   
 } #sub r_squared

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub in_interval
# Checks if a given point in a ghm (x and y coordinates passed as
# arguments) lies within the confidence interval of a cluster.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub in_interval
 {
  my ($cluster,$x,$y,$c_up,$c_down);
  ($cluster,$x,$y)=@_;

  ($c_up,$c_down)=$cluster->interval_bounds($x);

  ( ($y>=$c_down) and ($y<=$c_up) ) ? (return -1) : (return 0);

 } #sub in_interval

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub interval_bounds
# Returns for a given x-coordinate of a cluster the upper and lower
# y-coordinates of the confidence interval
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub interval_bounds
 {
  my ($cluster,$x,$n,$a,$b,$avg_x,$var_x,$variance,$sdev_y_est,$c_up,$c_down,$lambda);
  ($cluster,$x)=@_;

  $cluster->confidence_interval;
  $n=$cluster->number_of_anchorpoints;
  $a=$cluster->a;
  $b=$cluster->b;
  $avg_x=$cluster->avg_x;
  $var_x=$cluster->var_x;
  $variance=$cluster->variance;

  $lambda=Statistics::Distributions::tdistr($n,(0.01)/2); #Student-t distribution 99.9% confidence interval
  $sdev_y_est=sqrt( $variance * ( 1 + ( (1/$n) * (1+(($x-$avg_x)**2)/$var_x) ) ) );
  $c_up= ($a+$b*$x) + ($lambda*$sdev_y_est); #upper limit of the confidence interval
  $c_down=($a+$b*$x) - ($lambda*$sdev_y_est); #lower limit of the confidence interval
  
  return ($c_up, $c_down);
  
 } #sub interval_bounds

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub confidence_interval
# Updates, if necessary, all parameters describing the confidence
# interval of a cluster.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub confidence_interval
 {
  my ($cluster);
  $cluster=$_[0];
  
  if ( !$cluster->b )
   {
    my ($avg_x,$n,$a,$b,$y_est,$variance,$var_x,$i,$x_end1,$y_end1,$x_end2,$y_end2);

    $cluster->regression;
    $n=$cluster->number_of_anchorpoints; #number of points in cluster
    $avg_x=$cluster->avg_x;
    $a=$cluster->a;
    $b=$cluster->b;
    
    foreach my $anchorpoint ($cluster->get_anchorpoints)
     {
      my ($x,$y);
      $x=$anchorpoint->x;
      $y=$anchorpoint->y;
      $y_est=$a+($b*$x);
      $variance+= ($y-$y_est)**2;
      $var_x+=($x-$avg_x)**2;
     } #foreach my $anchorpoint ($cluster->get_anchorpoints)

    $variance=$variance/($n-2); #variance on deviation of y from the best-fit line
    $var_x=$var_x/$n; #variance on x
    $x_end1=$cluster->lowest_x;
    $x_end2=$cluster->highest_x;
    $y_end1=$x_end1*$b+$a;
    $y_end2=$x_end2*$b+$a;  

    $cluster->a($a);
    $cluster->b($b);
    $cluster->avg_x($avg_x);
    $cluster->var_x($var_x);
    $cluster->variance($variance);
    $cluster->x_end1($x_end1);
    $cluster->y_end1($y_end1);
    $cluster->x_end2($x_end2);
    $cluster->y_end2($y_end2);

   } #if ( !$cluster->b )

 } #sub confidence_interval

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub regression
# Calculates the best-fit line through the anchorpoints of a
# basecluster.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub regression
 {
  my ($cluster,$n,$sum_x,$sum_y,$sum_xy,$sum_x2,$x_avg,$y_avg,$a,$b);

  $cluster=$_[0];
  $n=$cluster->number_of_anchorpoints;

  ($n < 3) && return( -1 );

  foreach my $anchorpoint ($cluster->get_anchorpoints)
   {
    my ($x,$y);
    $x=$anchorpoint->x;
    $y=$anchorpoint->y;
    $sum_x+=$x; $sum_y+=$y;
    $sum_xy+=$x*$y;
    $sum_x2+=$x**2;
   } #foreach my $anchorpoint ($cluster->get_anchorpoints)
   
  $x_avg=$sum_x/$n;
  $y_avg=$sum_y/$n;

  $b=($sum_xy-($n*$x_avg*$y_avg))/($sum_x2-($n*$x_avg**2));
  $a=$y_avg-($b*$x_avg);    

  $cluster->a($a);
  $cluster->b($b);
  $cluster->avg_x($x_avg);

 } #sub regression

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub distance_to_point
# Calculates the distance from a cluster to a point in a ghm; x and
# y-coordinates of the point are given as arguments.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub distance_to_point
 {
  my ($cluster,$x,$y)=@_;
  
  if ( ( $x >= $cluster->x_end1 ) && ( $x <= $cluster->x_end2 ) )
   {
    return 0;
   } #if ( ( $x >= $cluster->x_end1 ) && ( $x <= $cluster->x_end2 ) )
  elsif ( $x < $cluster->x_end1 ) 
   {
    return iADHoRe::ghm->dpd( $x, $y, $cluster->x_end1, $cluster->y_end1 );
   } #elsif ( $x < $cluster->x_end1 )
  else 
   {
    return iADHoRe::ghm->dpd( $x, $y, $cluster->x_end2, $cluster->y_end2 );
   } #else 
  
 } #sub distance_to_point

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub distance_to_cluster
# Returns the distance from one basecluster to another in a GHM.
# Returns 0 if the x-coordinates or the y-coordinates of both clusters
# overlap; returns the dpd distance between the most adjacent
# endpoints in all other cases.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub distance_to_cluster
 {
  my ($cluster_a,$cluster_b,$lowest_distance);
  my (@distances);
  my (%a,%b);
  
  ($cluster_a,$cluster_b)=@_;
  
  $cluster_a->confidence_interval;
  $cluster_b->confidence_interval;
  
  $a{"x"}{"1"}=$cluster_a->x_end1;
  $a{"y"}{"1"}=$cluster_a->y_end1;
  $a{"x"}{"2"}=$cluster_a->x_end2;
  $a{"y"}{"2"}=$cluster_a->y_end2;
  $b{"x"}{"1"}=$cluster_b->x_end1;
  $b{"y"}{"1"}=$cluster_b->y_end1;
  $b{"x"}{"2"}=$cluster_b->x_end2;
  $b{"y"}{"2"}=$cluster_b->y_end2;

  #Check first if both the x and y coordinates overlap, IOW, if both clusters overlap
  if ( $cluster_a->overlapping_coordinates( $cluster_b ) == 2 )
   {
    $lowest_distance=0; #If so, the distance between the clusters is 0.
   } #$cluster_a->overlapping_coordinates( $cluster_b ) == 2
  else #If not, the lowest of all 4 possible distances between end points is returned.
   {      

    my ($first,$last)=(1,2);
    foreach ( 0,1 )
     {

      push( @distances, iADHoRe::ghm->dpd( $a{"x"}{$first}, $a{"y"}{$first}, $b{"x"}{$first}, $b{"y"}{$first} ) );
      push( @distances, iADHoRe::ghm->dpd( $a{"x"}{$first}, $a{"y"}{$first}, $b{"x"}{$last}, $b{"y"}{$last} ) );

      ($first,$last)=($last,$first);
     } #foreach ( 0,1 )

     $lowest_distance=shift @distances;
     foreach my $d ( @distances)
      {
       ( $d < $lowest_distance ) && ( $lowest_distance = $d );
      } #foreach my $d ( @distances)

   } #else

  return $lowest_distance;

 } #sub distance_to_cluster

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub overlapping_interval
# Returns true if all anchorpoints of one cluster lie within the
# confidence interval of another cluster.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub overlapping_interval
 {
  my ($cluster_a,$cluster_b,$overlapping);
  ($cluster_a,$cluster_b)=@_;
  
  $cluster_a->confidence_interval;
  $cluster_b->confidence_interval;
  
  $overlapping=-1;
  
  foreach ( 0, 1 )
   {

    #First, check if all anchorpoints of cluster_a lie within the confidence interval of cluster_b
    foreach my $anchorpoint ( $cluster_a->get_anchorpoints )
     {
      if ( !$cluster_b->in_interval( $anchorpoint->x, $anchorpoint->y ) )
       {
        $overlapping=0;
	last;
       } #if ( !$cluster_b->in_interval( $anchorpoint->x, $anchorpoint->y ) )
     } #foreach my $anchorpoint ( $cluster->get_anchorpoints )
    
    if ($overlapping) #If so, quit the loop
     {
      last;
     } #if ($overlapping)
    elsif ( $_ == 0 )  #If not, swap cluster_a and cluster_b and try again
     {
      ($cluster_a,$cluster_b) = ($cluster_b,$cluster_a);
      $overlapping=-1;
     } #else 
     
   } #foreach ( 0, 1 )

  return $overlapping;
  
 } #sub overlapping_interval

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub overlapping_coordinates
# Returns 1 if either the x- or the y-coordinates of two clusters
# overlap, 2 if both the x and y coordinates overlap and 0 if there is
# no overlap at all.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub overlapping_coordinates
 {
  my ($cluster_a,$cluster_b,$overlap_status);
  my (%a,%b);

  ($cluster_a,$cluster_b)=@_;
  $overlap_status = 0;

  $a{"x1"}=$cluster_a->x_end1;
  $a{"y1"}=$cluster_a->y_end1;
  $a{"x2"}=$cluster_a->x_end2;
  $a{"y2"}=$cluster_a->y_end2;
  $b{"x1"}=$cluster_b->x_end1;
  $b{"y1"}=$cluster_b->y_end1;
  $b{"x2"}=$cluster_b->x_end2;
  $b{"y2"}=$cluster_b->y_end2;

  ($cluster_a->orientation eq '-') && ( @a{"y1", "y2"} =  @a{"y2", "y1"} );
  ($cluster_b->orientation eq '-') && ( @b{"y1", "y2"} =  @b{"y2", "y1"} );

  ( $a{"x1"} < $b{"x2"} ) && ( $a{"x2"} > $b{"x1"} ) && $overlap_status++;
  ( $a{"y1"} < $b{"y2"} ) && ( $a{"y2"} > $b{"y1"} ) && $overlap_status++;

  return( $overlap_status );

 } #sub overlapping_coordinates

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub merge_with
# Adds all anchorpoints of a basecluster object passed as argument to
# the basecluster object the method is invoked on.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub merge_with
 {
  my ($cluster_a,$cluster_b)=@_;
  
  foreach my $anchorpoint ( $cluster_b->get_anchorpoints )
   {
    $cluster_a->add_anchorpoint( $anchorpoint->x, $anchorpoint->y );
   } #foreach my $anchorpoint ( $cluster_b->get_anchorpoints )
  
  $cluster_a->b(0); #Reset the cluster confidence interval statistics.
  $cluster_a->confidence_interval;
   
 } #sub merge_with

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub calculate_probability
# Calculates for a basecluster the probability to be generated by
# chance. The area and number of points in the ghm must be passed
# as arguments. The value returned is actually log ( 1- probability ).
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub calculate_probability
 {
  use POSIX;
  use Math::BigFloat;
  
  my ($cluster,$area,$points,$density,$n,$probability,$cells,$level,$anchorpoint_array);
  
  ($cluster, $area, $points, $level)=@_;
  
  $n= $cluster->number_of_anchorpoints;
  $points = $points - $n;
  $density = $points / $area;
  $anchorpoint_array=$cluster->anchorpoints;
  
  @$anchorpoint_array = sort { $a->x <=> $b->x } ( @$anchorpoint_array );
  
  $probability = Math::BigFloat->bone; #The initial probability of the cluster is set to 1

  for (my $i=1; $i < $n; $i++)
   {
    my ($c,$p_i,$distance);
    
    $distance = iADHoRe::ghm->dpd( $$anchorpoint_array[$i-1]->x, $$anchorpoint_array[$i-1]->y, 
                                   $$anchorpoint_array[$i]->x,   $$anchorpoint_array[$i]->y );

    #Remove any $anchorpoint_array that have the same coordinates (resulting in zero distance)
    if ( $distance == 0 )
     {
      my $removed_anchorpoint=splice( @$anchorpoint_array, $i, 1);
      $removed_anchorpoint->DESTROY;
      $n--;
      $cluster->number_of_anchorpoints( $n );
      $i--;
      next;
     } #if ( $distance == 0 )

    #Calculate the probability for the current anchorpoint
    $c=ceil($distance**2/2);
    $cells+=$c;
    $p_i=$c*$density*(1-$density)**($c-1);

    #Adjust the probability for higher-level multiplicons.
    $p_i=1-$p_i;
    $p_i=$p_i**($level - 1);
    $p_i=1-$p_i;

    $probability->bmul($p_i);
    
   } #for (my $i=1; $i < $n; $i++)

  $probability=1-$probability;
  $probability->blog;
  $probability->bmul($points);
  $probability->accuracy(10);  

  return $probability->bsstr;

 } #sub calculate_probability

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub twist_cluster
# Mirrors a cluster round its central y-value. All y-values of the
# anchor points are thus changed.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub twist_cluster
 {
  my ($cluster,$max_y,$min_y);

  $cluster=$_[0];
  
  $max_y=$cluster->highest_y;
  $min_y=$cluster->lowest_y;
  
  foreach my $anchorpoint ( $cluster->get_anchorpoints )
   {
    $anchorpoint->y( $max_y + $min_y - $anchorpoint->y );
   } #foreach my $anchorpoint ( $cluster->get_anchorpoints )

  $cluster->b(0);
  $cluster->confidence_interval;

 } #sub twist_cluster

#===============================================================================
# Accessor methods
#===============================================================================
 
sub anchorpoints
 {
  my ($basecluster,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'anchorpoints'}[ $$basecluster ] = $set_value );
  return ( $$structure{'anchorpoints'}[ $$basecluster ] );
 } #sub anchorpoints
 
#------------------------------------------------------------
 
sub random_probability
 {
  my ($basecluster,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'random_probability'}[ $$basecluster ] = $set_value );
  return ( $$structure{'random_probability'}[ $$basecluster ] );
 } #sub random_probability
 
#------------------------------------------------------------
 
sub orientation
 {
  my ($basecluster,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'orientation'}[ $$basecluster ] = $set_value );
  return ( $$structure{'orientation'}[ $$basecluster ] );
 } #sub orientation
 
#------------------------------------------------------------
 
sub multiplicon
 {
  my ($basecluster,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'multiplicon'}[ $$basecluster ] = $set_value );
  return ( $$structure{'multiplicon'}[ $$basecluster ] );
 } #sub multiplicon
 
#------------------------------------------------------------
 
sub a
 {
  my ($basecluster,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'a'}[ $$basecluster ] = $set_value );
  return ( $$structure{'a'}[ $$basecluster ] );
 } #sub a
 
#------------------------------------------------------------
 
sub b
 {
  my ($basecluster,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'b'}[ $$basecluster ] = $set_value );
  return ( $$structure{'b'}[ $$basecluster ] );
 } #sub b
 
#------------------------------------------------------------
 
sub avg_x
 {
  my ($basecluster,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'avg_x'}[ $$basecluster ] = $set_value );
  return ( $$structure{'avg_x'}[ $$basecluster ] );
 } #sub avg_x
 
#------------------------------------------------------------
 
sub var_x
 {
  my ($basecluster,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'var_x'}[ $$basecluster ] = $set_value );
  return ( $$structure{'var_x'}[ $$basecluster ] );
 } #sub var_x
 
#------------------------------------------------------------
 
sub variance
 {
  my ($basecluster,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'variance'}[ $$basecluster ] = $set_value );
  return ( $$structure{'variance'}[ $$basecluster ] );
 } #sub variance
 
#------------------------------------------------------------
 
sub x_end1
 {
  my ($basecluster,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'x_end1'}[ $$basecluster ] = $set_value );
  return ( $$structure{'x_end1'}[ $$basecluster ] );
 } #sub x_end1
 
#------------------------------------------------------------
 
sub y_end1
 {
  my ($basecluster,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'y_end1'}[ $$basecluster ] = $set_value );
  return ( $$structure{'y_end1'}[ $$basecluster ] );
 } #sub y_end1
 
#------------------------------------------------------------
 
sub x_end2
 {
  my ($basecluster,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'x_end2'}[ $$basecluster ] = $set_value );
  return ( $$structure{'x_end2'}[ $$basecluster ] );
 } #sub x_end2
 
#------------------------------------------------------------
 
sub y_end2
 {
  my ($basecluster,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'y_end2'}[ $$basecluster ] = $set_value );
  return ( $$structure{'y_end2'}[ $$basecluster ] );
 } #sub y_end2
 
#------------------------------------------------------------
 
sub was_twisted
 {
  my ($basecluster,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'was_twisted'}[ $$basecluster ] = $set_value );
  return ( $$structure{'was_twisted'}[ $$basecluster ] );
 } #sub was_twisted

#------------------------------------------------------------
 
sub id
 {
  my ($basecluster,$set_value)=@_;
  defined( $set_value ) && ( $$structure{'id'}[ $$basecluster ] = $set_value );
  return ( $$structure{'id'}[ $$basecluster ] );
 } #sub id

1;
