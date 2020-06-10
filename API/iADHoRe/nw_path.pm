package iADHoRe::nw_path;
#Created by Cedric Simillion on Tue Dec 28 22:04:09 CET 2004

use strict;

#===============================================================================
# Package methods
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_alignment_path
# - solve $min_x
# - work around extra row/column (all identities are offset +1,+1)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_alignment_path
 {
  my($pkg,$x_object,$y_list,$identities,$gap,$mismatch,$match,$min_x,$max_x,$max_y,$alignment_score,$back_x,$back_y,
     $prev_match_x);
  my(@matrix,@path);
  my(%backtrack);
  
  ($pkg,$x_object,$y_list)=@_;
  
  #We create a new GHM to get the identities between to objects to be aligned
  $identities=iADHoRe::ghm->new($x_object,$y_list);

  $gap=-2;
  $mismatch=-1;
  $match=100;

  $max_x=$identities->size_x;
  $max_y=$identities->size_y;

  foreach my $x ( 1 .. $max_x )  
   {
    if ( $identities->identity($x-1,0) )
     {
      $min_x=$x;
      last;
     } #if ( $identities->identity($x-1,0) )
   } #foreach my $x ( 1 .. $max_x ) 
    
  #Initialize first row of the matrix with zeroes
  foreach my $x (0..$max_x) 
   {
    $matrix[$x][0]=0;
   } #foreach my $x (0..$max_x)
  #Initialize first column of the matrix with zeroes
  foreach my $y (0..$max_y) 
   {
    $matrix[0][$y]=0;
   } #foreach my $y (0..$max_y)
   
  #Fill in the rest of the matrix
  foreach my $x (1..$max_x)
   {
    foreach my $y (1..$max_y)
     {
      my ($highest_score);
      my (@neighbours,@penalties,@scores);

      @neighbours= ( {"x" => $x-1, "y" => $y-1} , {"x" => $x, "y" => $y-1}, {"x" => $x-1, "y" => $y} );
      @penalties= ( $identities->identity($x-1,$y-1) ? $match : $mismatch ,$gap,$gap);
      foreach (0..2)      
       {
        $scores[$_]= $matrix[ $neighbours[$_]{"x"} ][ $neighbours[$_]{"y"} ] + $penalties[$_];
       } #foreach (0..2)
       
      ($highest_score)=sort {$b <=> $a} (@scores);
      $matrix[$x][$y]=$highest_score;
      
      #Find out which neighbouring cells had the highest score
      foreach (0..2) 
       {
        if ($scores[$_] == $highest_score)
	 {
	  $backtrack{$x}{$y}=$neighbours[$_]; #If different options, we prefer first up and left, then up, then left.
	  last;
	 } #if ($scores[$_] == $highest_score)
       } #foreach (0..2)
       
     } #foreach my $y (1..$max_y)
   } #foreach my $x (1..$max_x)
  
  #Determine which cell on the bottom row will be the starting point for backtracking. This must be an identity position.
  $alignment_score=$matrix[1][$max_y];
  foreach my $x (2..$max_x)
   {
    if ( ($matrix[$x][$max_y] > $alignment_score) && $identities->identity($x-1, $max_y-1) )
     {
      $alignment_score=$matrix[$x][$max_y];
      $back_x=$x;
     } #if ($matrix[$x][$max_y] > $aligment_score)
   } #foreach my $x (1..$max_x)
  

  #Check if the alignment path is valid
  if ( !defined( $back_x ) )
   {
    my $path = [];
    bless $path, $pkg;
    return( $path );
   } #if ( !defined( $back_x ) )

  #Perform backtracking.
  $back_y=$max_y;
  $prev_match_x=0;
  do
   {
    #If the current position is on a match position add it to the path provided that the previous position added
    #to the path doesn't have the same x-coordinate
    if ( $identities->identity($back_x-1,$back_y-1) && ($prev_match_x != $back_x) )
     {
       push(@path, [ $back_x-1 , $back_y-1 ] );
       $prev_match_x=$back_x;
     } #if ( $is_match{$back_x}{$back_y} && ($prev_match_x != $back_x) )
    ($back_x,$back_y)=@{$backtrack{$back_x}{$back_y}}{"x","y"};
   } until ( ($back_y==0) || ($back_x<$min_x) );
  
  #If the alignment path does not run to the top of the matrix, append the top left most identity position to it.
  if ($back_y!=0)
   {
    ($prev_match_x==$min_x) && pop(@path);
    push(@path, [ $min_x-1, 0 ] );
   } #if ($back_y!=0)
  
  ( ( $path[-1][0] == 0 ) && ( $path[-1][1] == 0 ) ) || push( @path, [0, 0] ); #Make sure the path begins with (0,0)
  #Temporary hack to remove points on the path that have y=0 or x=0
  until ( ( scalar( @path ) < 2 ) || ( ( $path[-2][0] > 0 ) && ( $path[-2][1] > 0 ) ) )
   {
    splice( @path, -2, 1);
   } #( ( scalar( @path ) < 2 ) || ( ( $path[-2][0] > 0 ) && ( $path[-2][1] > 0 ) ) )
  
  ( scalar( @path ) < 2 ) && ( @path = () );
  
  @path=reverse( @path );
  bless \@path, $pkg;

 } #get_alignment_path

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub next_offset
# Returns for an alignment path the relative  offset
# (i.e. the difference) between the current and next (x,y)
# coordinates.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub next_offset
 {
  my $path=$_[0];

  if ( scalar( @$path ) > 1 )
   {
    my $x=$$path[1][0] - $$path[0][0];
    my $y=$$path[1][1] - $$path[0][1];
    shift( @$path );
    return( $x,$y);
   } #if ( scalar( @$path ) > 1 )
  else 
   {
    return();
   } #else 
   
 } #sub next_offset

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub valid
# Returns true if an alignment could be found, false if not.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub valid
 {
  my $path = $_[0];
  return scalar( @$path ) ? -1 : 0;
 } #sub valid

1;
