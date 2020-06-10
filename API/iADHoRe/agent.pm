package iADHoRe::agent;
#Created by Cedric Simillion on Wed May 18 18:23:06 CEST 2005

use strict;
use DBI;
use iADHoRe::genelist;
use iADHoRe::list_element;

our ($settings,$AUTOLOAD);

$settings=iADHoRe::settings->get_settings;

#===============================================================================
# Package methods
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub contact_central
# Contacts the central and returns an agent handle;
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub contact_central
 {
  my ($package,$dbase,$dbase_name,$dbase_host,$login,$password,$agent,$database_found,$agent_id,$dataset,$run_host);
  ($package,$dataset,$password)=@_;
  
  $dbase_host = $settings->database_host;
  $dbase_host || die "No database host specified!\n";
  $dbase_name = $settings->database_name;
  $dbase_name || die "No database name specified!\n";
  $login = $settings->database_login ? $settings->database_login : $ENV{'USER'};
  
  #Contact the database server
  $dbase=DBI->connect("dbi:mysql:host=$dbase_host",$login,$password)
   || die "Could not connect to database server $dbase_host as $login: $!\n";
  $dbase->{PrintError}=1;
  $dbase->{RaiseError}=1;

  #Wait until the database has been created by the central
  $database_found=0;
  do #until ($database_found)
   {
    my $db_list=$dbase->selectcol_arrayref("SHOW DATABASES");
    ( grep { $_ eq $dbase_name } ( @$db_list ) ) ? ( $database_found=-1 ) : sleep( 1 );
   } until ($database_found);

  #Connect to the database and register the agent
  $run_host = `hostname`;
  chomp $run_host;
  $dbase->do("USE $dbase_name");
  $dbase->do("LOCK TABLES _agents WRITE");
  $dbase->do("INSERT INTO _agents (host, status) VALUES ('$run_host', 'idle')");
  ($agent_id)=$dbase->selectrow_array("SELECT last_insert_id() FROM _agents");
  $dbase->do("UNLOCK TABLES");
  
  $agent = { 'connection' => $dbase, #Contains the connection to the database
             'dataset' => $dataset, #The iADHoRe dataset (genes + genelists);
             'id' => $agent_id, #The identififier of the agent.
	     'task_id' => undef, #The id of the current task.
	     'last_x_object_id' => -1, #The identifier of the last x_object
	     'last_x_object_type' => '', #Self-explenatory
	     'last_x_object' => undef, #Self-explenatory
	     'gene_hash' => {}, # A hash table to retrieve gene objects by their database identifiers.
	     'select_elements' => undef,    # -+
	     'select_multiplicon' => undef, #  +-> often-used SQL statement handles
	     'select_segment_ids' => undef, # -+
	     };
  
  bless $agent, $package;

  #Construct the gene hash table
  $agent->map_gene_ids;
  
  $agent->prepare_queries;

  return ($agent);

 } #sub contact_central

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub map_gene_ids
# Fills in the gene_hash table of an agent handle object. This table
# is used to link gene objects to their database identifiers.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub map_gene_ids
 {
  my ($agent,$dataset,$gene_hash); 
  $agent=$_[0];
  $dataset=$agent->dataset;
  $gene_hash=$agent->gene_hash;
  
  foreach my $genelist ( $dataset->genelists )
   {
    foreach my $list_element ( @{$genelist->elements} )
     {
      $$gene_hash{ $list_element->gene->ID }= $list_element->gene;
     } #foreach my $list_element ( @{$genelist->elements} )
   } #foreach my $genelist ( @{$dataset->genelists} )
   
 } #sub map_gene_ids

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub prepare_queries
# Prepares a few often-used SQL queries
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub prepare_queries
 {
  my ($agent,$dbase,$sql);
  $agent=$_[0];
  $dbase=$agent->connection;
  
  $sql="SELECT gene, position, orientation, masked FROM list_elements WHERE segment = ?";
  $agent->select_elements( $dbase->prepare( $sql ) );
  
  $sql="SELECT * FROM multiplicons WHERE id = ?";
  $agent->select_multiplicon( $dbase->prepare( $sql ) );
  
  $sql="SELECT id FROM segments WHERE multiplicon = ? ORDER BY `order` ASC";
  $agent->select_segment_ids( $dbase->prepare( $sql ) );
  
 } #sub prepare_queries

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_task
# Gets a new task from the central. a task consists of an x_object
# and a y_list. 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_task
 {
  my ($agent,$dbase,$status,$x_object,$y_list,$x_object_id,$y_list_id,$x_object_type,$task_id);
  $agent=$_[0];
  $dbase=$agent->connection;

  until ( defined $task_id )
   {
    #Abort if the status of the agent has been set to 'terminated'
    ( $agent->status eq 'terminated' ) && return;

    #Try to get a task description from the database
    $dbase->do("LOCK TABLES _tasks WRITE");
    ($task_id, $x_object_id, $x_object_type, $y_list_id)=
      $dbase->selectrow_array("SELECT id, x_object, x_object_type, y_list FROM _tasks WHERE status = 'queued'");

    #If a task was found, check it out.
    if ( defined $task_id )
     {
      $dbase->do("UPDATE _tasks SET status = 'running', agent = ".$agent->id." WHERE id = $task_id");
      $dbase->do("UNLOCK tables");
      $agent->status( 'busy' );
      $agent->task_id( $task_id );     
     } #if ( defined $task_id )
    else #If no task could be found, sleep for a second.
     {
      $dbase->do("UNLOCK tables");
      sleep( 1 );
     } #else 

   } #until ( defined $task_id )
  
  #Only retrieve the x_object if it's different then the previous one.
  unless ( ( $x_object_id == $agent->last_x_object_id ) && ( $x_object_type eq $agent->last_x_object_type ) )
   {
    $x_object= ($x_object_type eq 'genelist') ? $agent->get_list( $x_object_id ) : $agent->get_profile( $x_object_id );
    $agent->last_x_object_id( $x_object_id );
    $agent->last_x_object_type( $x_object_type );
    $agent->last_x_object( $x_object );
   } #unless ( ( $x_object_id == $agent->last_x_object_id ) && ( $x_object_type eq $agent->last_x_object_type ) )
  else
   {
    $x_object = $agent->last_x_object;
   } #else

  $y_list = $agent->get_list( $y_list_id );
  
  return ( $x_object, $y_list ); 
  
 } #sub get_task


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub status
# Retrieves the current status of an agent from the database.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub status
 {
  my ($agent,$dbase,$status,$new_status);
  ($agent,$new_status)=@_;
  $dbase=$agent->connection;

  if ( defined $new_status )
   {
    $dbase->do("LOCK TABLES _agents WRITE");
    $dbase->do("UPDATE _agents SET status = '$new_status' WHERE id = ".$agent->id);
    $status=$new_status;
    $dbase->do("UNLOCK TABLES");
   } #if ( defined $new_status )
  else
   {
    $dbase->do("LOCK TABLES _agents READ");
    ($status)=$dbase->selectrow_array("SELECT status FROM _agents WHERE id = ".$agent->id);
    $dbase->do("UNLOCK TABLES");
   } #else

  return ($status);
 } #sub status

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub switch_to_idle_if_busy
# Switches to status of an agent to 'idle' if it previously was set to
# 'busy'
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub switch_to_idle_if_busy
 {
  my ($agent,$dbase,$status);
  $agent=$_[0];
  $dbase=$agent->connection;
  
  $dbase->do("LOCK TABLES _agents WRITE");
  ($status)=$dbase->selectrow_array("SELECT status FROM _agents WHERE id = ".$agent->id);
  ($status eq 'busy') && $dbase->do("UPDATE _agents SET status = 'idle' WHERE id = ".$agent->id);
  $dbase->do("UNLOCK TABLES");
  
 } #sub switch_to_idle_if_busy

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_list
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_list
 {
  my ($agent,$list_id,$remapped,$genelist,$element_table,$gene_table,$genome,$listname);
  
  ($agent,$list_id)=@_;

  $gene_table=$agent->gene_hash;
  
  $agent->select_elements->execute( $list_id );
  $element_table=$agent->select_elements->fetchall_hashref("gene");
  
  foreach my $element (values %$element_table)
   {
    $$remapped[ $$element{'position'} ]=
     iADHoRe::list_element->new( $$gene_table{ $$element{'gene'} }, $$element{'orientation'}, $$element{'masked'} );
   } #foreach my $element (values %$element_table)
  
  foreach my $element ( @$remapped )
   {
    defined( $element ) || next;
    $genome=$element->gene->genome;
    $listname=$element->gene->list->listname;
    last;
   } #foreach my $element ( @$remapped )
   
  $genelist=iADHoRe::genelist->new( $listname, $genome, [], $remapped, undef );
  
  return $genelist;
 } #sub get_list

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub get_profile
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub get_profile
 {
  my ($agent,$profile_id, $profile, $multiplicon, $record, $segments, $segment_ids);
  ($agent,$profile_id)=@_;
  
  $agent->select_multiplicon->execute( $profile_id );
  $record=$agent->select_multiplicon->fetchrow_hashref;
  
  $multiplicon=iADHoRe::multiplicon->new
    ( undef, undef, @{$record}{qw(begin_x end_x begin_y end_y number_of_anchorpoints level )} );
  
  $agent->select_segment_ids->execute( $profile_id );
  @$segment_ids= map { $$_[0] } ( @{$agent->select_segment_ids->fetchall_arrayref} );
  
  foreach my $id ( @$segment_ids )
   {
    push( @$segments, $agent->get_list( $id ) );
   } #foreach my $id ( @$segment_ids )
  
  $profile=iADHoRe::profile->new( $multiplicon, -1, $segments );
  $multiplicon->profile( $profile );
  
  return( $profile );
  
 } #sub get_profile

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub submit_results
# Submits the output of a task to the database. The output consists
# of a set of multiplicons that are each linked to one or more
# baseclusters and a set of anchorpoints.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub submit_results
 {
  my ($agent,$dbase,$store_multiplicon,$store_basecluster,$store_anchorpoint);
  $agent=shift( @_ );
  my @multiplicons= @_;
  
  scalar( @multiplicons ) || return;
  
  $dbase=$agent->connection;
  
  $dbase->do("LOCK TABLES _new_multiplicons WRITE, _new_baseclusters WRITE, _new_anchorpoints WRITE");

  $store_multiplicon=$dbase->prepare
   ("INSERT INTO _new_multiplicons
       ( genome_x, list_x, parent, genome_y, list_y, level, number_of_anchorpoints, begin_x, end_x, begin_y, end_y )
      VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ? )");
  $store_basecluster=$dbase->prepare
   ("INSERT INTO _new_baseclusters
       (multiplicon, number_of_anchorpoints, orientation, was_twisted, random_probability, begin_x, end_x, begin_y, end_y )
      VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)");
  $store_anchorpoint=$dbase->prepare
   ("INSERT INTO _new_anchorpoints (multiplicon, basecluster, gene_x, gene_y, coord_x, coord_y, is_real_anchorpoint)
      VALUES (?, ?, ?, ?, ?, ?, ?)");
        
  foreach my $multiplicon ( @multiplicons ) 
   {
    my $multiplicon_id;
    my (@store_values);

    ($agent->last_x_object_type eq 'genelist')
      ? push( @store_values, $multiplicon->x_object->genome, $multiplicon->x_object->listname, undef )
      : push( @store_values, undef, undef, $agent->last_x_object_id );
    push( @store_values, $multiplicon->y_list->genome, $multiplicon->y_list->listname );
    push( @store_values, $multiplicon->level, $multiplicon->number_of_anchorpoints, $multiplicon->begin_x,
           $multiplicon->end_x, $multiplicon->begin_y, $multiplicon->end_y );
    
    $store_multiplicon->execute( @store_values );
    ($multiplicon_id)=$dbase->selectrow_array("SELECT last_insert_id() FROM _new_multiplicons");
    
    foreach my $basecluster ( $multiplicon->get_baseclusters )
     {
      my ($basecluster_id);
      $store_basecluster->execute( $multiplicon_id, $basecluster->number_of_anchorpoints, $basecluster->orientation,
        $basecluster->was_twisted, $basecluster->random_probability, $basecluster->begin_x, $basecluster->end_x,
	$basecluster->begin_y, $basecluster->end_y );
      ($basecluster_id)=$dbase->selectrow_array("SELECT last_insert_id() FROM _new_baseclusters");
      
      foreach my $anchorpoint ( $basecluster->get_anchorpoints )	
       {
        $store_anchorpoint->execute( $multiplicon_id, $basecluster_id, $anchorpoint->gene_x->ID, $anchorpoint->gene_y->ID,
	  $anchorpoint->x, $anchorpoint->y, $anchorpoint->is_real_anchorpoint );
       } #foreach my $anchorpoint ( $basecluster->get_anchorpoints )	

     } #foreach my $basecluster ( $multiplicon->get_baseclusters )
    
   } #foreach my $multiplicon ( @multiplicons ) 

  $dbase->do("UNLOCK TABLES");
  
 } #sub submit_results

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub finish_task
# Remove the task from the database and set the agent status to idle.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub finish_task
 {
  my $agent=$_[0];
  $agent->connection->do("DELETE FROM _tasks WHERE id = ".$agent->task_id);
  $agent->switch_to_idle_if_busy;
  $agent->task_id( undef );
 } #sub finish_task

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub AUTOLOAD
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub AUTOLOAD
 {
  my ($object, $setvalue, $attribute);

  $object=$_[0];
  $setvalue=$_[1];
  ($attribute)= ($AUTOLOAD=~ /::([^:]+)$/);
  
  ($attribute eq 'DESTROY') && return;
  
  exists( $$object{$attribute} ) || die "Attribute $attribute not defined for ".ref($object).".\n";

  defined($setvalue) && ( $$object{$attribute}=$setvalue );
  return ( $$object{$attribute} );
  
 } #sub AUTOLOAD

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub DESTROY
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub DESTROY
 {
  my $agent=$_[0];
  $agent->connection->do("LOCK TABLES _agents WRITE");
  $agent->connection->do("DELETE FROM _agents WHERE id = ".$agent->id);
  $agent->connection->do("UNLOCK TABLES");

  $agent->select_elements->finish;
  $agent->select_multiplicon->finish;
  $agent->select_segment_ids->finish;
  
  $agent->connection->disconnect;
 } #sub DESTROY

1;
