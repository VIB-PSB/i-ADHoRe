#!/usr/bin/perl

=head1 Description

This script serves no purpose by itself but it called by the multiplicon_plot.pl script.

=cut

my @segments;
my %chrominfo;
my %relations;
my %cursegment;
my $curseg;
my $currel;
my $numsegments;
my $inpfile=$ARGV[0];
my $outputfile=$inpfile.".eps";

srand(time);

####################################################################
##  PARSING INPUT FILE         
####################################################################

undef %cursegment;
$currel=0;
$numsegments=0;
open(INP, "<$inpfile") || die  "Couldn't open file $inpfile\n";
while(<INP>){
  my $line=$_;
  chomp($line);
  $_=$line;
  if(/^segment\s(\d+)/){
    $curmode=1;
    $curseg=$1;
    #if cursegment exists, save it in the hash and start new one
    #if(defined %cursegment){
    #  push @segments, %cursegment;
      $numsegments++;
    #  undef %cursegment;
    #}
  }
  elsif(/^generelations/){
    $curmode=2;
  }
  elsif((not (/^#/)) and (length($line) > 4)){
    chomp($line);
    if($curmode==1){
      if(/^chromosome\s(.+)/){
        $chrominfo{$curseg}=$1;
      }
      else{
        my @vals=split/\t/,$line;
        if (scalar(@vals) == 4){
          #parse into geneid, start, stop and orientation
          #push(@{${%{$segments[$curseg-1]}}{$vals[0]}},$vals[1]);
          #push(@{${%{$segments[$curseg-1]}}{$vals[0]}},$vals[2]);
          #push(@{${%{$segments[$curseg-1]}}{$vals[0]}},$vals[3]);
          @{${%{$segments[$curseg-1]}}{$vals[0]}}=($vals[1],$vals[2],substr($vals[3],0,1));
        }
        else {die "Incorrect format in line : \n$line\n";}
      }
    }
    elsif($curmode==2){
      chomp($line);
      my @vals=split/\t/,$line;
      if (scalar(@vals) == 4){
        #parse into genex segmentx geney segmenty
        $_=$vals[3];/^(\d+)/;$vals[3]=$1;
        @{$relations{$currel}}=($vals[0],$vals[1],$vals[2],$vals[3]);
        $currel++;
      }
      else {die "Incorrect format in line : \n$line\n";}
    }
    else {die "Inconsistent mode : curmode = $curmode\n";}
  }
}
close(INP);
#if(defined %cursegment){
#  push @segments, %cursegment;
#  $numsegments++;
#  undef %cursegment;
#}

#print "$numsegments segments and $currel relations\n";
#for($iter=0; $iter < $numsegments; $iter++){
#  print "Segment $iter ".($segments[$iter])."\n";
#  foreach $key1 (keys %{$segments[$iter]}){
#    print "$key1 : ".(${%{$segments[$iter]}}{$key1})." ";
#    #foreach my $val (@{${%{$segments[$iter]}}{$key1}}){
#    # print "$val ";
#    #}
#    print "\n";
#  }
#}

####################################################################
##  DISTRIBUTING SEGMENTS ON CIRClE
####################################################################

my @lengths;
#determine the length of each segment
my $min;
my $max;
for(my $iter=0;$iter<$numsegments;$iter++){
  #determine the minimum of all startpostions
  my @allkeys=keys %{$segments[$iter]};
  $firstkey = $allkeys[0];
  #print "first key = ".($firstkey)."\n";
  $min=${@{${%{$segments[$iter]}}{$firstkey}}}[1];
  $max=${@{${%{$segments[$iter]}}{$firstkey}}}[0];
  #print "min = $min max=$max\n";
  foreach my $key (keys %{$segments[$iter]}){
    my $start=${@{${%{$segments[$iter]}}{$key}}}[0];
    my $stop=${@{${%{$segments[$iter]}}{$key}}}[1];
    if ($start < $min){$min=$start;}
    if ($stop > $max){$max=$stop;}
    #print "curstart = $start curstop = $stop min = $min max = $max\n";
  }
  $lengths[$iter]=$max-$min+1;
  #print "Len $iter = ".($lengths[$iter])."\n";
  #Now add all relative start en stop sites of the genes
  foreach my $key (keys %{$segments[$iter]}){
    my $start=${@{${%{$segments[$iter]}}{$key}}}[0];
    my $stop=${@{${%{$segments[$iter]}}{$key}}}[1];
    push(@{${%{$segments[$iter]}}{$key}},$start-$min);
    push(@{${%{$segments[$iter]}}{$key}},$stop-$min);
  }
}

# Drawing parameters
my $mx=297.5, $my=421;		#centre of the circle
my $radius=200;			#radius of the circle     
my $gapdeg=10;			#angle for gaps
my $triangdeg=0.5;		#angle for triangle part of arrows
my $PI=3.14159265358979323846264338327950288419716939937510 ;
my $r1=(2*$PI*$radius)*(1-($gapdeg/360));	#distance left for segments
my $arrowhight=5;			#5 pts high arrows
my $fontsize=8;			#standard fontsize for text

my %start_angles;
my %rad_distances;
my %end_angles;

my $degleft = 360 - ($numsegments*$gapdeg);
my $totlen=0;
for(my $iter=0;$iter<$numsegments;$iter++){
  $totlen=$totlen + $lengths[$iter];
}
for(my $iter=0;$iter<$numsegments;$iter++){
  if($iter==0){
    $start_angles{$iter}=90;
  }
  else{
    $start_angles{$iter}=$end_angles{$iter-1}-$gapdeg;
    #if($start_angles{$iter} < 0){$start_angles{$iter} +=360;}
  }
  $rad_distances{$iter}=($lengths[$iter] / $totlen)*$r1;
  $end_angles{$iter}=$start_angles{$iter}-($degleft*$lengths[$iter] / $totlen);
  #if($end_angles{$iter} < 0){$end_angles{$iter}+=360;}
}    

open(OUT,">$outputfile");
print OUT "%!PS-Adobe-3.0 EPSF-1.2
%%Title: $outputfile
%%LanguageLevel: 1
%%Creator: dupliviz.pl (C) Yvan Saeys
%%CreationDate: Wed Feb 27 10:55:04 2002
%%For: yvsae
%%DocumentMedia: A4 595.27559 841.88976 0 ( ) ( )
%%Orientation: Portrait
%%BoundingBox: 0 0 595 841
%%EndComments
%%BeginProlog
%%EndProlog
0 setlinewidth\n";

#Visualize the segments
print OUT "\n% Segments\n";
for(my $iter=0;$iter<$numsegments;$iter++){
  my $start=$start_angles{$iter};
  my $stop=$end_angles{$iter};
  #if($stop < $start){
  #  my $tmp=$stop;
  #  $stop=$start;
  #  $start=$tmp;
  #}
  #print "$start $stop\n";
  print OUT "% Segment $iter\n";
  print OUT "newpath\n";
  print OUT "$mx $my $radius $start $stop arcn\n";
  print OUT "stroke\n";
}


####################################################################
##  DETERMINING GENES ON SEGMENTS
####################################################################

for(my $iter=0;$iter<$numsegments;$iter++){
  foreach my $key (keys %{$segments[$iter]}){
    my $start=${@{${%{$segments[$iter]}}{$key}}}[3];
    my $stop=${@{${%{$segments[$iter]}}{$key}}}[4];
    push(@{${%{$segments[$iter]}}{$key}},$start_angles{$iter}-($start*$degleft/$totlen));
    push(@{${%{$segments[$iter]}}{$key}},$start_angles{$iter}-($stop*$degleft/$totlen));
  }
}

sub drawgene{
  (my $sa, my $ea, my $or, my $cl)= @_;  #start_angle, end_angle, orientation, color
  if(abs($ea-$sa) < $triangdeg){
    my $phi=0;
    my $p1x=0;my $p1y=$radius+$arrowhight;
    my $theta=abs($ea-$sa);
    if($theta < 90){$phi=90-$theta;}else{$phi=$theta-90;}
    #print "phi=$phi\n";
    $phi=$phi*$PI/180;
    my $p3x=$radius*cos($phi);my $p3y=$radius*sin($phi);
    print STDOUT "$sa $ea\n";
    print "phi=$phi p3x = $p3x  p3y = $p3y\n cos(phi)=".(cos($phi))."sin(phi)=".(sin($phi))."\n";
    my $radiusmin=$radius-$arrowhight;
    my $radiusplus=$radius+$arrowhight;
    print OUT "gsave\n";
    print OUT "$mx $my translate\n";
    print OUT (-(90-$sa))." rotate\n";
    print OUT "newpath\n";
    my @vals=split/ /,$cl;
    #print STDOUT ($vals[0])." ".($vals[1])." ".($vals[2])." setrgbcolor\n";
    print OUT ($vals[0])." ".($vals[1])." ".($vals[2])." setrgbcolor\n";
    print OUT "$p1x $p1y moveto\n";
    print OUT "0 0 $radiusplus 90 ".(90-abs($ea-$sa))." arcn\n";
    print OUT "$p3x $p3y lineto\n";
    print OUT "0 0 $radiusmin ".(90-abs($ea-$sa))." 90 arc\n";
    print OUT "closepath\nfill\ngrestore\n";
  }
  elsif($or eq "+"){
    my $phi=0;
    my $p1x=0;my $p1y=$radius+$arrowhight;
    my $theta=abs($ea-$sa);
    if($theta < 90){$phi=90-$theta;}else{$phi=$theta-90;}
    #print "phi=$phi\n";
    $phi=$phi*$PI/180;
    my $p3x=$radius*cos($phi);my $p3y=$radius*sin($phi);
    print STDOUT "$sa $ea\n";
    print "phi=$phi p3x = $p3x  p3y = $p3y\n cos(phi)=".(cos($phi))."sin(phi)=".(sin($phi))."\n";
    my $radiusmin=$radius-$arrowhight;
    my $radiusplus=$radius+$arrowhight;
    print OUT "gsave\n";
    print OUT "$mx $my translate\n";
    print OUT (-(90-$sa))." rotate\n";
    print OUT "newpath\n";
    my @vals=split/ /,$cl;
    #print STDOUT ($vals[0])." ".($vals[1])." ".($vals[2])." setrgbcolor\n";
    print OUT ($vals[0])." ".($vals[1])." ".($vals[2])." setrgbcolor\n";
    print OUT "$p1x $p1y moveto\n";
    print OUT "0 0 $radiusplus 90 ".(90-abs($ea-$sa)+$triangdeg)." arcn\n";
    print OUT "$p3x $p3y lineto\n";
    print OUT "0 0 $radiusmin ".(90-abs($ea-$sa)+$triangdeg)." 90 arc\n";
    print OUT "closepath\nfill\ngrestore\n";
  }    
  else{
    my $phi=0;
    my $p1x=0;my $p1y=$radius+$arrowhight;
    my $theta=abs($ea-$sa);
    if($theta < 90){$phi=90-$theta;}else{$phi=$theta-90;}
    #print "phi=$phi\n";
    $phi=$phi*$PI/180;
    my $p3x=-($radius*cos($phi));my $p3y=$radius*sin($phi);
    print STDOUT "$sa $ea\n";
    print "phi=$phi p3x = $p3x  p3y = $p3y\n cos(phi)=".(cos($phi))."sin(phi)=".(sin($phi))."\n";
    my $radiusmin=$radius-$arrowhight;
    my $radiusplus=$radius+$arrowhight;
    print OUT "gsave\n";
    print OUT "$mx $my translate\n";
    print OUT (-(90-$ea))." rotate\n";
    print OUT "newpath\n";
    my @vals=split/ /,$cl;
    #print STDOUT ($vals[0])." ".($vals[1])." ".($vals[2])." setrgbcolor\n";
    print OUT ($vals[0])." ".($vals[1])." ".($vals[2])." setrgbcolor\n";
    print OUT "$p1x $p1y moveto\n";
    print OUT "0 0 $radiusplus 90 ".(90+abs($ea-$sa)-$triangdeg)." arc\n";
    print OUT "$p3x $p3y lineto\n";
    print OUT "0 0 $radiusmin ".(90+abs($ea-$sa)-$triangdeg)." 90 arcn\n";
    print OUT "closepath\nfill\ngrestore\n";
  }    
    
}

my @colors=('1.0 0.0 0.0','0.0 1.0 0.0','0.0 0.0 1.0','0.8 0.0 0.0','0.0 0.8 0.0','0.0 0.0 0.8','0.5 0.0 0.0','0.0 0.5 0.0','0.0 0.0 0.5','0.5 0.5 0.5','0.2 0.4 0.6');
my @greycolors=('0.8 0.8 0.8','0.82 0.82 0.82', '0.84 0.84 0.84', '0.86 0.86 0.86', '0.88 0.88 0.88', '0.90 0.90 0.90');

sub drawrelation{
  (my $gene1_sa, my $gene1_ea, my $gene2_sa, my $gene2_ea, my $o1, my $o2, my $cl)= @_;  
  print "g1sa=$gene1_sa g1ea=$gene1_ea g2sa=$gene2_sa g2ea=$gene2_ea $o1 $o2\n";
  #gene1 start_angle, gene1 end_angle,gene2 start_angle, gene2 end_angle,g1-g2 orientation, color
  #construct the exact border points to connect them
  print OUT "gsave\n$mx $my translate\n";
  my $s1x=$radius*cos($gene1_sa*$PI/180);  my $s2x=$radius*cos($gene2_sa*$PI/180);
  my $s1y=$radius*sin($gene1_sa*$PI/180);  my $s2y=$radius*sin($gene2_sa*$PI/180);
  my $e1x=$radius*cos($gene1_ea*$PI/180);  my $e2x=$radius*cos($gene2_ea*$PI/180);
  my $e1y=$radius*sin($gene1_ea*$PI/180);  my $e2y=$radius*sin($gene2_ea*$PI/180);
  if($o1 ne $o2){
    print OUT "newpath\n";
    my @vals=split/ /,$cl;
    print OUT ($vals[0])." ".($vals[1])." ".($vals[2])." setrgbcolor\n";
    print OUT "$s1x $s1y moveto\n";
    print OUT "0 0 $radius $gene1_sa $gene1_ea arcn\n";
    print OUT "$e2x $e2y lineto\n";
    print OUT "0 0 $radius $gene2_ea $gene2_sa arc\n";
    print OUT "closepath\nfill\ngrestore\n";
  }
  else{
    print OUT "newpath\n";
    my @vals=split/ /,$cl;
    print OUT ($vals[0])." ".($vals[1])." ".($vals[2])." setrgbcolor\n";
    print OUT "$s1x $s1y moveto\n";
    print OUT "0 0 $radius $gene1_sa $gene1_ea arcn\n";
    print OUT "$s2x $s2y lineto\n";
    print OUT "0 0 $radius $gene2_sa $gene2_ea arcn\n";
    print OUT "closepath\nfill\ngrestore\n";
  }
}    

sub drawrelation2{
  (my $gene1_sa, my $gene1_ea, my $gene2_sa, my $gene2_ea, my $o1, my $o2, my $cl)= @_;  
  print "g1sa=$gene1_sa g1ea=$gene1_ea g2sa=$gene2_sa g2ea=$gene2_ea $o1 $o2\n";
  #gene1 start_angle, gene1 end_angle,gene2 start_angle, gene2 end_angle,g1-g2 orientation, color
  #construct the exact border points to connect them
  my $h1=abs($gene1_ea-$gene1_sa);my $h2=abs($gene2_ea-$gene2_sa);
  my $theta1=$gene1_sa-($h1/2);
  my $theta2=$gene2_sa-($h2/2);
  print OUT "gsave\n$mx $my translate\n";
  my $s1x=$radius*cos($theta1*$PI/180); my $s1y=$radius*sin($theta1*$PI/180);
  my $s2x=$radius*cos($theta2*$PI/180); my $s2y=$radius*sin($theta2*$PI/180);
  if($o1 ne $o2){
    print OUT "newpath\n";
    print OUT "0.25 0.25 0.25 setrgbcolor\n";
    print OUT "$s1x $s1y moveto\n";
    print OUT "$s2x $s2y lineto\n";
    print OUT "closepath\nstroke\ngrestore\n";
  }
  else{
    print OUT "newpath\n";
    print OUT "0.0 0.0 0.0 setrgbcolor\n";
    print OUT "$s1x $s1y moveto\n";
    print OUT "$s2x $s2y lineto\n";
    print OUT "closepath\nstroke\ngrestore\n";
  }
}

sub draw_straight_text{
  (my $angle, my $text, my $size)= @_;  
  my $neg=0;
  #calculate startpoint
  if ($angle < -90){
    $angle+=180;
    $neg=1;
  }
  print "anlge=$angle\n";
  my $p1x=($radius+(2*$arrowhight))*cos($angle*$PI/180);
  my $p1y=($radius+(2*$arrowhight))*sin($angle*$PI/180);
  print OUT "/Times-Roman findfont $fontsize scalefont setfont\n";
  print OUT "gsave\n";
  print OUT "$mx $my translate\n";
  print OUT ($angle)." rotate\n";
  print OUT "newpath\n";
  if(!$neg){
    print OUT ($radius+(2*$arrowhight))." "."-2 moveto\n";
  }
  else{
    print OUT (-($radius+(length($text)*$fontsize*2/3)))." "."-2 moveto\n";
  }
  print OUT "($text) show stroke\n";
  print OUT "grestore\n";
}

sub draw_horizontal_text{
  (my $angle, my $text, my $size)= @_;  
  #calculate startpoint
  print OUT "gsave\n";
  print OUT ($mx-($size/2))." ".($my-($size/2))." translate\n";
  my $p1x=(($radius+(2*$arrowhight)+$size)*cos($angle*$PI/180));
  my $p1y=(($radius+(2*$arrowhight)+$size)*sin($angle*$PI/180));
  print OUT "/Times-Roman findfont $fontsize scalefont setfont\n";
  print OUT "newpath\n";
  print OUT "$p1x $p1y moveto\n";
  print OUT "($text) show stroke\ngrestore\n";
}

print OUT "\n% Relations\n";
foreach my $key (keys %relations){
  my $g1=$relations{$key}[0];
  my $s1=${@{$relations{$key}}}[1];
  my $g2=${@{$relations{$key}}}[2];
  my $s2=${@{$relations{$key}}}[3];
  my $gene1_sa=${@{${%{$segments[$s1-1]}}{$g1}}}[5];
  my $gene1_ea=${@{${%{$segments[$s1-1]}}{$g1}}}[6];
  my $gene2_sa=${@{${%{$segments[$s2-1]}}{$g2}}}[5];
  my $gene2_ea=${@{${%{$segments[$s2-1]}}{$g2}}}[6];
  my $o1=${@{${%{$segments[$s1-1]}}{$g1}}}[2];
  my $o2=${@{${%{$segments[$s2-1]}}{$g2}}}[2];
  #print "g1=[$g1] g2=[$g2] s1=[$s1] s2=[$s2] g1sa=[$gene1_sa] g1ea=[$gene1_ea] g2sa=[$gene2_sa] g2ea=[$gene2_ea] [$o1] [$o2]\n";
  print OUT "%Relation : $g1 $g2\n";
  drawrelation2($gene1_sa,$gene1_ea,$gene2_sa,$gene2_ea,$o1,$o2,$greycolors[rand(5)]);
}

print OUT "\n% Genes\n";
#drawgene(90,60,'+','0.5 1.0 0.5');
for(my $iter=0;$iter<$numsegments;$iter++){
#for(my $iter=0;$iter<1;$iter++){
  $cl=rand(10);
  foreach my $key (keys %{$segments[$iter]}){
    #print "key $key\n";
    my $start=${@{${%{$segments[$iter]}}{$key}}}[5];
    my $stop=${@{${%{$segments[$iter]}}{$key}}}[6];
    my $or=${@{${%{$segments[$iter]}}{$key}}}[2];
    #my $cl='0.5 1.0 0.5';
    print OUT "%Gene $key on segment $iter\n";
    drawgene($start,$stop,$or,$colors[$iter%scalar(@colors)]);
  }
}

print OUT "\n% Genes start and stop positions\n";
for(my $iter=0;$iter<$numsegments;$iter++){
  print OUT "%Segment $iter\n";
  my $minangle3; my $maxangle3;
  my $mintext;my $maxtext;
  my $counter=0;
  foreach my $key (keys %{$segments[$iter]}){
     $angle1=${@{${%{$segments[$iter]}}{$key}}}[5];
     $angle2=${@{${%{$segments[$iter]}}{$key}}}[6];
     $angle3=$angle2+(abs($angle2-$angle1))*2/3;
     #$text=${@{${%{$segments[$iter]}}{$key}}}[0];
     $text=$key;
     if ($counter ==0){
       $minangle3=$angle3;
       $maxangle3=$angle3;
       $mintext=$text;
       $maxtext=$text;
     }
     #draw_straight_text($angle3,$text,$fontsize);
     if ($angle3 < $minangle3){
       $minangle3 = $angle3;
       $mintext=$text;
     }
     if ($angle3 > $maxangle3){
       $maxangle3 = $angle3;
       $maxtext=$text;
     }
     $counter++;
  }
  draw_straight_text($minangle3,$mintext,$fontsize);
  draw_straight_text($maxangle3,$maxtext,$fontsize);
}
print OUT "%%EOF\n";
close(OUT);
