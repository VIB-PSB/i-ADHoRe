package iADHoRe::central;
#Created by Cedric Simillion on Wed May 18 18:22:42 CEST 2005

use strict;
use iADHoRe::settings;
use iADHoRe::agent;
use DBI;

our ($scheme,$settings);
our @ISA=qw(iADHoRe::agent);

$settings=iADHoRe::settings->get_settings;

#Prepare the database scheme
$scheme="CREATE TABLE `_agents` (
  `id` smallint(5) unsigned NOT NULL auto_increment,
  `host` varchar(100) NOT NULL default '',
  `status` enum('idle','busy','terminated') NOT NULL default 'idle',
  PRIMARY KEY  (`id`),
  KEY `status` (`status`)
);
CREATE TABLE `_new_anchorpoints` (
  `id` int(10) unsigned NOT NULL auto_increment,
  `multiplicon` int(10) unsigned NOT NULL default '0',
  `basecluster` int(10) unsigned NOT NULL default '0',
  `gene_x` varchar(80) NOT NULL default '',
  `gene_y` varchar(80) NOT NULL default '',
  `coord_x` mediumint(8) unsigned NOT NULL default '0',
  `coord_y` mediumint(8) unsigned NOT NULL default '0',
  `is_real_anchorpoint` tinyint(4) NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `multiplicon` (`multiplicon`,`basecluster`,`gene_x`,`gene_y`,`coord_x`,`coord_y`,`is_real_anchorpoint`)
);
CREATE TABLE `_new_baseclusters` (
  `id` int(10) unsigned NOT NULL auto_increment,
  `multiplicon` int(10) unsigned NOT NULL default '0',
  `number_of_anchorpoints` mediumint(8) unsigned NOT NULL default '0',
  `orientation` enum('+','-') NOT NULL default '+',
  `was_twisted` tinyint(4) NOT NULL default '0',
  `random_probability` double unsigned NOT NULL default '0',
  `begin_x` mediumint(8) unsigned NOT NULL default '0',
  `end_x` mediumint(8) unsigned NOT NULL default '0',
  `begin_y` mediumint(8) unsigned NOT NULL default '0',
  `end_y` mediumint(8) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `multiplicon` (`multiplicon`,`number_of_anchorpoints`,`orientation`,`was_twisted`,`random_probability`)
);
CREATE TABLE `_new_multiplicons` (
  `id` int(10) unsigned NOT NULL auto_increment,
  `genome_x` varchar(80) NULL default '',
  `list_x` varchar(80) NULL default '',
  `parent` int(10) unsigned NULL default '0',
  `genome_y` varchar(80) NOT NULL default '',
  `list_y` varchar(80) NOT NULL default '',
  `level` smallint(5) unsigned NOT NULL default '0',
  `number_of_anchorpoints` mediumint(8) unsigned NOT NULL default '0',
  `begin_x` mediumint(8) unsigned NOT NULL default '0',
  `end_x` mediumint(8) unsigned NOT NULL default '0',
  `begin_y` mediumint(8) unsigned NOT NULL default '0',
  `end_y` mediumint(8) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `genome_x` (`genome_x`,`list_x`,`parent`,`genome_y`,`list_y`,`level`,`number_of_anchorpoints`)
);
CREATE TABLE `_tasks` (
  `id` bigint(20) unsigned NOT NULL auto_increment,
  `x_object` int(10) unsigned NOT NULL default '0',
  `x_object_type` enum('genelist','profile') NOT NULL default 'genelist',
  `y_list` int(10) unsigned NOT NULL default '0',
  `status` enum('queued','running') NOT NULL default 'queued',
  `agent`smallint(5) unsigned NULL default NULL,
  PRIMARY KEY  (`id`),
  KEY `status` (`status`)
);
CREATE TABLE `anchorpoints` (
  `id` int(10) unsigned NOT NULL auto_increment,
  `multiplicon` int(10) unsigned NOT NULL default '0',
  `basecluster` int(10) unsigned NOT NULL default '0',
  `gene_x` varchar(80) NOT NULL default '',
  `gene_y` varchar(80) NOT NULL default '',
  `coord_x` mediumint(8) unsigned NOT NULL default '0',
  `coord_y` mediumint(8) unsigned NOT NULL default '0',
  `is_real_anchorpoint` tinyint(4) NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `multiplicon` (`multiplicon`,`basecluster`,`gene_x`,`gene_y`,`coord_x`,`coord_y`,`is_real_anchorpoint`)
);
CREATE TABLE `baseclusters` (
  `id` int(10) unsigned NOT NULL auto_increment,
  `multiplicon` int(10) unsigned NOT NULL default '0',
  `number_of_anchorpoints` mediumint(8) unsigned NOT NULL default '0',
  `orientation` enum('+','-') NOT NULL default '+',
  `was_twisted` tinyint(4) NOT NULL default '0',
  `random_probability` double unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `multiplicon` (`multiplicon`,`number_of_anchorpoints`,`orientation`,`was_twisted`,`random_probability`)
);
CREATE TABLE `genes` (
  `id` varchar(80) NOT NULL default '',
  `genome` varchar(80) NOT NULL default '',
  `list` varchar(80) NOT NULL default '',
  `coordinate` mediumint(8) unsigned NOT NULL default '0',
  `orientation` enum('+','-') NOT NULL default '+',
  `remapped_coordinate` mediumint(8) unsigned NOT NULL default '0',
  `is_tandem` tinyint(4) NOT NULL default '0',
  `is_tandem_representative` tinyint(4) NOT NULL default '0',
  `tandem_representative` varchar(80) NULL default '',
  `remapped` tinyint(4) NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `genome` (`genome`,`list`,`coordinate`,`orientation`,`remapped_coordinate`,`is_tandem`,`is_tandem_representative`,`tandem_representative`,`remapped`)
);
CREATE TABLE `list_elements` (
  `id` bigint(20) unsigned NOT NULL auto_increment,
  `segment` int(10) unsigned NOT NULL default '0',
  `gene` varchar(80) NOT NULL default '',
  `position` mediumint(8) unsigned NOT NULL default '0',
  `orientation` enum('+','-') NOT NULL default '+',
  `masked` tinyint(4) NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `segment` (`segment`,`gene`,`position`),
  KEY `masked` (`masked`)
);
CREATE TABLE `multiplicons` (
  `id` int(10) unsigned NOT NULL auto_increment,
  `genome_x` varchar(80) NULL default '',
  `list_x` varchar(80) NULL default '',
  `parent` int(10) unsigned NULL default '0',
  `genome_y` varchar(80) NOT NULL default '',
  `list_y` varchar(80) NOT NULL default '',
  `level` smallint(5) unsigned NOT NULL default '0',
  `number_of_anchorpoints` mediumint(8) unsigned NOT NULL default '0',
  `profile_length` mediumint(8) unsigned NOT NULL default '0',
  `begin_x` mediumint(8) unsigned NOT NULL default '0',
  `end_x` mediumint(8) unsigned NOT NULL default '0',
  `begin_y` mediumint(8) unsigned NOT NULL default '0',
  `end_y` mediumint(8) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `genome_x` (`genome_x`,`list_x`,`parent`,`genome_y`,`list_y`,`level`,`number_of_anchorpoints`)
);
CREATE TABLE `segments` (
  `id` int(10) unsigned NOT NULL auto_increment,
  `multiplicon` int(10) unsigned NULL default NULL,
  `genome` varchar(80) NOT NULL default '',
  `list` varchar(80) NOT NULL default '',
  `first` varchar(80) NOT NULL default '',
  `last` varchar(80) NOT NULL default '',
  `order` mediumint(8) unsigned NULL default NULL,
  PRIMARY KEY  (`id`),
  KEY `multiplicon` (`multiplicon`,`genome`,`list`,`order`)
);";

#===============================================================================
# Package methods
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub setup
# Sets up the central by taking an iADHoRe::dataset object as input
# and setting up a database containing the dataset
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub setup
 {
  my ($package,$dbase,$dbase_name,$host,$login,$password,$dataset,$central);
  ($package,$dataset,$password)=@_;
  
  print STDERR "Setting up the central...\n";
  $host = $settings->database_host;
  $host || die "No database host specified!\n";
  $dbase_name = $settings->database_name;
  $dbase_name || die "No database name specified!\n";
  $login = $settings->database_login ? $settings->database_login : $ENV{'USER'};
  
  #Set up the database
  $dbase=DBI->connect("dbi:mysql::$host",$login,$password)
   || die "Could not connect to database server $host as $login: $!\n";
  $dbase->{PrintError}=1;
  $dbase->{RaiseError}=1;
  
  foreach my $statement ("CREATE DATABASE $dbase_name", "USE $dbase_name", split(/;/,$scheme) )
   {
    $dbase->do( $statement );
   } #foreach my $statement ("CREATE DATABASE $dbase_name", "USE $dbase_name", split(/;/,$scheme) )
  
  #Create the $central handle object
  $central = { 'connection' => $dbase,       #Contains the connection to the database
               'dataset' => $dataset,        #The iADHoRe dataset (genes + genelists);
	       'genelist_ids' => {},         #A hash table to retrieve IDs of the primary genelists in the database.
	       'genelists' => {},            #A hash table to retrieve gene lists objects by their genome and list name.
	       'gene_hash' => {},            #A hash table to retrieve gene objects by their database identifiers.
	       'parent' => {},               #Stores the parent multiplicon of each level
	       'parent_id' => {},            #Stores the database id of the parent multiplicon of each level
	       'current_level' => undef,     #The current search level
	       'store_multiplicon' => undef, # -+
	       'store_basecluster' => undef, #  |
	       'store_anchorpoint' => undef, #  |
	       'store_segment' => undef,     #  +--> often-used SQL statement handles
	       'store_list_element' =>undef, #  |
	       'store_gene' => undef,        #  |
	       'mask_elements' => undef      # -+ 
	       }; 
  bless $central, $package;
  
  #Store the dataset
  $central->prepare_queries;
  print STDERR "\nStoring the dataset...";
  $central->store_dataset;
  $central->map_gene_ids;
  print STDERR "\n";  
  
  return $central;
  
 } #sub setup

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub prepare_queries
# Prepares a set of often-used statement handles
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub prepare_queries
 {
  my ($central,$dbase,$sql);
  $central=$_[0];
  $dbase=$central->connection;
  
  $sql="INSERT INTO multiplicons (genome_x, list_x, parent, genome_y, list_y, level, number_of_anchorpoints, profile_length,
         begin_x, end_x, begin_y, end_y)
	VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ? ,?)";
  $central->store_multiplicon( $dbase->prepare( $sql ) );

  $sql="INSERT INTO baseclusters (multiplicon, number_of_anchorpoints, orientation, was_twisted, random_probability)
        VALUES (?, ?, ?, ?, ?)";
  $central->store_basecluster( $dbase->prepare( $sql ) );

  $sql="INSERT INTO anchorpoints (multiplicon, basecluster, gene_x, gene_y, coord_x, coord_y, is_real_anchorpoint)
        VALUES (?, ?, ?, ?, ?, ?, ?)";  
  $central->store_anchorpoint( $dbase->prepare( $sql ) );
  
  $sql="INSERT INTO segments (`multiplicon`, `genome`, `list`, `first`, `last`, `order`) VALUES (?, ?, ?, ?, ?, ?)";
  $central->store_segment( $dbase->prepare( $sql ) );
  
  $sql="INSERT INTO list_elements (segment, gene, position, orientation, masked) VALUES (?, ?, ?, ?, ?)";
  $central->store_list_element( $dbase->prepare( $sql ) );
  
  $sql="INSERT INTO genes VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
  $central->store_gene( $dbase->prepare( $sql ) );
  
  $sql="UPDATE list_elements SET masked = -1 WHERE segment = ? AND position >= ? AND position <= ?";
  $central->mask_elements( $dbase->prepare( $sql ) );
  
 } #sub prepare_queries

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub store_dataset
# Stores the genes and genelists in a iADHoRe::dataset in the
# database.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub store_dataset
 {
  my ($central,$dataset);
  $central=$_[0];
  $dataset=$central->dataset;
  
  foreach my $genelist ( $dataset->genelists )
   {
    my ($position,$segment_id);
   
    foreach my $gene ( map {$_->gene} ( @{$genelist->elements} ) )
     {
      $central->store_gene->execute( $gene->ID, $gene->genome, $gene->list->listname, $gene->coordinate, $gene->orientation,
        $gene->remapped_coordinate, $gene->is_tandem, $gene->is_tandem_representative,
	defined($gene->tandem_representative ) ? $gene->tandem_representative->ID : undef, $gene->remapped );
     } #foreach my $gene ( map {$_->gene} ( @{$genelist->elements} ) )

    $central->store_segment->execute( undef, $genelist->genome, $genelist->listname, ${$genelist->remapped}[0]->gene->ID,
      ${$genelist->remapped}[-1]->gene->ID, undef );
    ($segment_id)=$central->connection->selectrow_array("SELECT last_insert_id()");
    $position=0;
    foreach my $list_element ( @{$genelist->remapped} ) 
     {
      $central->store_list_element->execute
        ( $segment_id, $list_element->gene->ID, $position, $list_element->orientation, $list_element->masked );
      $position++;	
     } #foreach my $list_element ( @{$genelist->remapped} )
    
    ${$central->genelist_ids}{ $genelist->genome }{ $genelist->listname } = $segment_id;
    ${$central->genelists}{ $genelist->genome }{ $genelist->listname } = $genelist;
   } #foreach my $genelist ( $dataset->genelists )

 } #sub store_dataset

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub level1_adhore
# Submits the tasks to perform a level1 adhore.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub level1_adhore
 {
  my $central=$_[0];
  my (@task_list,@genelist_ids);

  print STDERR "Issuing level 2 multiplicon detection...";
  foreach my $list (values %{$central->genelist_ids} )
   {
    push( @genelist_ids, values %$list );
   } #foreach my $list (values %genomes)
  @genelist_ids= sort {$a <=> $b} ( @genelist_ids );
  
  foreach my $i ( 0 .. $#genelist_ids )
   {
    foreach my $j ( $i .. $#genelist_ids )
     {
      push( @task_list, [ $genelist_ids[$i], 'genelist', $genelist_ids[$j] ] );
     } #foreach my $j ( $i .. $#genelist_ids )
   } #foreach my $i ( 0 .. $#genelist_ids )
  
  $central->issue_tasks( @task_list );
  print STDERR "\nWaiting for results...";
  $central->get_results;
  print STDERR "\n";
  
 } #sub level1_adhore

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub profile_detection
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub profile_detection
 {
  my ($central,$dataset);
  $central=$_[0];
  $dataset=$central->dataset;
  
  print STDERR "\nHigher level detection\n";
  while (my $multiplicon = $dataset->next_multiplicon)
   {
    my (@new_multiplicons);
    print STDERR $dataset->multiplicons_left." multiplicons to evaluate - ";
    $central->current_level( $multiplicon->level );
    print STDERR "now at level ".$multiplicon->level."...";
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

    $central->save_multiplicon( $multiplicon );
    $multiplicon->mask_segments;
    $central->store_masking( $multiplicon );

    print STDERR " issuing search...";
    $settings->level_2_only || $central->profile_search;
    print STDERR "\n";
   } #while (my $multiplicon = $dataset->next_multiplicon)
  
 } #sub profile_detection

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub profile_search
# Submits the tasks to scan the dataset with the current profile.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub profile_search
 {
  my $central=$_[0];
  my (@genelist_ids,@multiplicons);
  
  foreach my $list (values %{$central->genelist_ids} )
   {
    push( @genelist_ids, values %$list );
   } #foreach my $list (values %genomes)
  @genelist_ids= sort {$a <=> $b} ( @genelist_ids );

  $central->issue_tasks( map { [ $central->parent_id, 'profile', $_ ] } (@genelist_ids) );
  
  $central->get_results;
  
 } #sub profile_search

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub save_multiplicon
# Stores an evaluated multiplicon in the database along with its
# baseclusters, anchorpoints, profile segments and elements. Also sets
# the parent and parent_id attributes of the central object.  These
# are used keep track of the parent multiplicon when doing higher
# level multiplicon detection.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub save_multiplicon
 {
  my ($central,$dataset,$multiplicon,$multiplicon_id,$order, $parent_id);
  ($central,$multiplicon)=@_;
  
  $central->store_multiplicon->execute(
   ($multiplicon->level == 2) ? ( $multiplicon->x_object->genome, $multiplicon->x_object->listname, undef ) :
   (undef, undef, $central->parent_id), $multiplicon->y_list->genome, $multiplicon->y_list->listname, $multiplicon->level,
   $multiplicon->number_of_anchorpoints, $multiplicon->profile->size, $multiplicon->begin_x, $multiplicon->end_x,
   $multiplicon->begin_y, $multiplicon->end_y );

  ($multiplicon_id)=$central->connection->selectrow_array("SELECT last_insert_id()");
  $central->current_level( $central->current_level + 1 );
  $central->parent_id( $multiplicon_id );
  $central->parent( $multiplicon );
  
  foreach my $basecluster ( $multiplicon->get_baseclusters )
   {
    $central->store_basecluster->execute( $multiplicon_id, $basecluster->number_of_anchorpoints, $basecluster->orientation, 
     $basecluster->was_twisted, $basecluster->random_probability );
    my ($basecluster_id)= $central->connection->selectrow_array("SELECT last_insert_id()");
    
    foreach my $anchorpoint ( $basecluster->get_anchorpoints )
     {
      $central->store_anchorpoint->execute( $multiplicon_id, $basecluster_id, $anchorpoint->gene_x->ID, $anchorpoint->gene_y->ID, 
       $anchorpoint->x, $anchorpoint->y, $anchorpoint->is_real_anchorpoint );
     } #foreach my $anchorpoint ( $basecluster->get_anchorpoints )

   } #foreach my $basecluster ( $multiplicon->get_baseclusters )
  
  $order = 0;
  foreach my $segment ( @{$multiplicon->profile->segments} )
   {
    my ($segment_id,$first,$last);
    ($segment_id) = $central->connection->selectrow_array("SELECT COUNT(*) FROM segments");
    $segment_id++;

    foreach my $element ( @{$segment->remapped} )
     {
      defined( $element ) || next;
      $first=$element->gene;
      last;
     } #foreach my $element ( @${$segment->remapped} )
    $last=$first;

    foreach my $i ( 0 .. $#{$segment->remapped} )
     {
      my $element = ${$segment->remapped}[$i];
      defined( $element ) || next;
      $central->store_list_element->execute( $segment_id, $element->gene->ID, $i, $element->orientation, $element->masked);
      ( $element->gene->coordinate < $first->coordinate ) && ( $first = $element->gene );
      ( $element->gene->coordinate > $last->coordinate ) && ( $last = $element->gene );
     } #foreach my $i ( 0 .. $#{$segment->remapped} )

    $central->store_segment->execute($multiplicon_id, $segment->genome, $segment->listname, $first->ID, $last->ID, $order);
    $order++;
   } #foreach my $segment ( @{$multiplicon->profile->segments} )
  
 } #sub save_multiplicon

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub store_masking
# Updates the masked attributes for newly masked list elements in the
# primary gene lists.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub store_masking
 {
  my ($central)=@_;
  
  foreach my $xy ( 'x', 'y' )
   {
    my ($object,$segment_id,$begin,$end);
    
    if ($xy eq 'x')
     {
      $object =$central->parent->x_object;
      $begin = $central->parent->begin_x;
      $end = $central->parent->end_x;
     } #if ($xy eq 'x')
    else 
     {
      $object = $central->parent->y_list;
      $begin = $central->parent->begin_y;
      $end = $central->parent->end_y;
     } #else 
    ( ref( $object ) eq 'iADHoRe::genelist' ) || next;
    
    $segment_id = ${$central->genelist_ids}{ $object->genome }{ $object->listname };
    
    $central->mask_elements->execute( $segment_id, $begin, $end );
    
   } #foreach my $xy ( 'x', 'y' )

 } #sub store_masking
 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub issue_tasks
# Issues a series of tasks. Each task is described as an x_object and
# a y_list. These two will be compared by an agent using the ADHoRe
# algorithm and any detected multiplicons can be retrieved using the
# get_results method.
# Takes as input a 2 dimensional array, where each row consists of 3
# columns:
# - The x_object id
# - The x_object type
# - The y_list id
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub issue_tasks
 {
  my  ($central, $add_task, $dbase);
  my @tasks;
  $central=shift( @_ );
  @tasks= @_;
  $dbase=$central->connection;
  
  $dbase->do("LOCK TABLES _tasks WRITE");
  $add_task=$dbase->prepare("INSERT INTO _tasks (x_object, x_object_type, y_list) VALUES (?, ?, ?)");
  foreach my $task ( @tasks )
   {
    $add_task->execute( @$task );
   } #foreach my $task ( @tasks )
  $add_task->finish; 
  $dbase->do("UNLOCK TABLES");
  
 } #sub issue_tasks

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_results
# Gets any detected multiplicons from issued tasks
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_results
 {
  my ($central, $dbase,$tasks_left,$multiplicon_table,$basecluster_table,$anchorpoint_table,$dataset);
  my (%multiplicons,%baseclusters,%anchorpoints);
  
  $central=$_[0];
  $dbase=$central->connection;
  $dataset=$central->dataset;
  
  #First, wait until all the tasks have been finished
  $tasks_left=1;
  until ( $tasks_left == 0 )
   {
    $dbase->do("LOCK TABLES _tasks READ");
    ($tasks_left)=$dbase->selectrow_array("SELECT COUNT(*) FROM _tasks");
    $dbase->do("UNLOCK TABLES");
    $tasks_left && sleep( 1 );
   } #until ( $tasks_finished )
  
  $multiplicon_table=$dbase->selectall_hashref("SELECT * FROM _new_multiplicons", "id");
  $basecluster_table=$dbase->selectall_hashref("SELECT * FROM _new_baseclusters", "id");
  $anchorpoint_table=$dbase->selectall_hashref("SELECT * FROM _new_anchorpoints", "id");
  
  foreach my $record ( values %$multiplicon_table )
   {
    my ($x_object,$y_list,$multiplicon);
    $x_object = ( $$record{'level'} == 2 )
     ? ${$central->genelists}{ $$record{'genome_x'} }{ $$record{'list_x'} }
     : $central->parent->profile;
    $y_list = ${$central->genelists}{ $$record{'genome_y'} }{ $$record{'list_y'} };
    $multiplicon=iADHoRe::multiplicon->new
     ( $x_object, $y_list, @{$record}{qw(begin_x end_x begin_y end_y number_of_anchorpoints level)}, [], undef, undef );
    $multiplicons{ $$record{'id'} } = $multiplicon;
   } #foreach my $record ( values %$multiplicon_table )
  
  foreach my $record ( values %$basecluster_table )
   {
    my ( $multiplicon, $basecluster );
    $multiplicon = $multiplicons{ $$record{'multiplicon'} };
    $basecluster = iADHoRe::basecluster->new
     ( $multiplicon->x_object, $multiplicon->y_list, @{$record}{qw(begin_x end_x begin_y end_y number_of_anchorpoints)}, 
       [], @{$record}{qw(random_probability orientation)}, $multiplicon, undef, undef, undef, undef, undef, undef, undef, 
       undef, undef, $$record{'was_twisted'} );
    push( @{$multiplicon->baseclusters}, $basecluster );
    $multiplicon->number_of_baseclusters( scalar( @{$multiplicon->baseclusters} ) );
    $baseclusters{ $$record{'id'} } = $basecluster;
   } #foreach my $record ( values %$basecluster_table )
  
  foreach my $record ( values %$anchorpoint_table )
   {
    my ( $multiplicon, $basecluster, $anchorpoint, $gene_hash );
    $multiplicon = $multiplicons{ $$record{'multiplicon'} };
    $basecluster = $baseclusters{ $$record{'basecluster'} };
    $gene_hash = $central->gene_hash;
    $anchorpoint=iADHoRe::anchorpoint->new
     ( $$gene_hash{ $$record{'gene_x'} }, $$gene_hash{ $$record{'gene_y'} }, 
       @{$record}{qw(coord_x coord_y is_real_anchorpoint)}, $basecluster, $multiplicon );
    push( @{$basecluster->anchorpoints}, $anchorpoint );   
   } #foreach my $record ( values %$anchorpoint_table )
  
  $dataset->add_multiplicons( sort by_multiplicon_size (values %multiplicons) );
  
  foreach my $table ( qw(_new_anchorpoints _new_baseclusters _new_multiplicons) )
   {
    $dbase->do("TRUNCATE TABLE $table");
   } #foreach my $table ( qw(_new_anchorpoints _new_baseclusters _new_multiplicons) )
  
  print STDERR " ".scalar( keys %multiplicons )." multiplicons found.";
   
 } #sub get_results

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub output
# When all is done, write the output files.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub output
 {
  my ($central,$dbase);
  $central=$_[0];
  $dbase=$central->connection;
  
  print STDERR "Creating output files";

  $central->substitute_coordinates;
  
  foreach my $table ( qw(anchorpoints baseclusters genes list_elements multiplicons segments) )
   {
    my ($fields, $sth);
    $fields=$dbase->selectcol_arrayref("SHOW COLUMNS FROM $table");
    $sth=$dbase->prepare("SELECT * FROM $table");
    $sth->execute;
    open (OUT, ">".$settings->output_path."$table.txt" );
    print OUT join("\t", @$fields )."\n";
    while (my @row = $sth->fetchrow_array)
     {
      print OUT join("\t", map { defined( $_ ) ? $_ : '' } (@row) )."\n";
     } #while (my @row = $sth->fetchrow_array)
    close OUT;
    print STDERR ".";
   } #foreach my $table ( qw(anchorpoints baseclusters genes list_elements multiplicons segments) )

  $dbase->do("DROP TABLE _agents");
  $dbase->do("DROP TABLE _tasks");
  $dbase->do("DROP TABLE _new_multiplicons");
  $dbase->do("DROP TABLE _new_baseclusters");
  $dbase->do("DROP TABLE _new_anchorpoints");
  $dbase->do("ALTER TABLE list_elements DROP COLUMN masked");

  print STDERR "\n";

 } #sub output

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub substitute_coordinates
# Replaces to remapped coordinates in the multiplicons table by the
# unmapped coordinates.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub substitute_coordinates
 {
  my ($central,$genes,$multiplicons,$update);
  my (%unmapped);
  $central=$_[0];
  $genes=$central->connection->selectall_hashref("SELECT * FROM genes", 'id');
  
  foreach my $record (values %$genes )
   {
    $unmapped{ $$record{'genome'} }{ $$record{'list'} }{ $$record{'remapped_coordinate'} } = $$record{'coordinate'};
   } #foreach my $record (values %$genes )
  
  $update=$central->connection->prepare
    ("UPDATE multiplicons SET begin_x = ?, end_x = ?, begin_y = ?, end_y = ? WHERE id = ?");

  $multiplicons=$central->connection->selectall_hashref("SELECT * FROM multiplicons", "id");
  
  foreach my $record (values %$multiplicons)  
   {
    my ($begin_x,$end_x,$begin_y,$end_y);
    if ( defined( $$record{'genome_x'} ) )
     {
      $begin_x= $unmapped{ $$record{'genome_x'} }{ $$record{'list_x'} }{ $$record{'begin_x'} };
      $end_x= $unmapped{ $$record{'genome_x'} }{ $$record{'list_x'} }{ $$record{'end_x'} };
     } #if ( defined( $$record{'genome_x'} ) )
    else 
     {
      $begin_x=$$record{'begin_x'};
      $end_x=$$record{'end_x'};
     } #else 

    $begin_y= $unmapped{ $$record{'genome_y'} }{ $$record{'list_y'} }{ $$record{'begin_y'} };
    $end_y= $unmapped{ $$record{'genome_y'} }{ $$record{'list_y'} }{ $$record{'end_y'} };

    $update->execute( $begin_x, $end_x, $begin_y, $end_y, $$record{'id'} )
   } #foreach my $record (values %$multiplicons)  
  
 } #sub substitute_coordinates

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub parent_id
# Accessor method for the parent_id attribute
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub parent_id
 {
  my ($central,$new_id,$level);
  ($central,$new_id)=@_;
  
  $level = $central->current_level;
  defined( $new_id ) && ( $$central{'parent_id'}{$level}= $new_id );
  
  return $$central{'parent_id'}{$level};
 } #sub parent_id

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub parent
# Accessor method for the parent attribute
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub parent
 {
  my ($central,$level,$new_parent);
  ($central,$new_parent)=@_;

  $level = $central->current_level;
  defined( $new_parent ) && ( $$central{'parent'}{$level}= $new_parent );
  
  return $$central{'parent'}{$level};
 } #sub parent

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub terminate_agents
# Terminates all agents
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub terminate_agents
 {
  my $central=$_[0];
  $central->connection->do("LOCK TABLES _agents WRITE");
  $central->connection->do("UPDATE _agents SET status = 'terminated'");
  $central->connection->do("UNLOCK TABLES");
 } #sub terminate_agents

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
# sub DESTROY
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub DESTROY
 {
  my $agent=$_[0];
  $agent->connection->disconnect;
 } #sub DESTROY
 
1;
