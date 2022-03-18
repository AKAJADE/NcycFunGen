#Import taxon information to database based on TaxID from GenBank
#/usr/bin/perl

#use strict;
use DBI;
use autodie;

#my $dir = "workdir";
#chdir $dir;

my $host="server_address";#server
my $dbname="database_name";#database name
my $dbuser="user_name";#database name
my $dbpass="password";#database user password
my $table="table_name"; # requrired table with taxnomy info in database

#my $taxfile = " ";
my $taxfile = "gi_taxid_ncbi_taxnonomy.txt";#table of GI and TaxID from GenBank
&connectdb(\$dbh,$host,$dbuser,$dbpass,$dbname);

open(IN,"$taxfile");
my %db_title;
my %taxid_gi;
my %taxid_tax;
my @added;
my $j=0;
while(<IN>){	
	if(/^(Gi?\S+)/){		
		chomp($_);
		@added = split("\t",$_);
		print "Querying titles from the DB database...\n";
		my $statement = "SELECT COLUMN_NAME,ORDINAL_POSITION FROM information_schema.COLUMNS WHERE TABLE_NAME= '$table' AND TABLE_SCHEMA = '$dbname'";
		my $sth = $dbh->prepare($statement) or die "Can't prepare $statement: $dbh->errstr\n";
		$sth-> execute();
		while(my @titles =  $sth->fetchrow_array()){
			my ($id,$pos) = @titles;
			#print "ID: $id\tPos: $pos\n";
			$db_title{$id} = $pos;
		}				
		my $i;
		for($i=1;$i<=$#added;$i++){
			unless(exists $db_title{$added[$i]}){
				print "Adding $added[$i] to the DB table $table\n";
				my $statement = "ALTER TABLE $table ADD `$added[$i]` text ";
				my $sth = $dbh->prepare($statement) or die "Can't prepare $statement: $dbh->errstr\n";
				$sth-> execute();				
			}
		}
	}elsif(/(\d+\S+)/){			
		chomp($_);
		my @temp = split("\t",$_);
		if(exists $taxid_tax{"$temp[1]"}){
			push @{$taxid_gi{$temp[1]}},$temp[0];
		}else{
			push @{$taxid_gi{$temp[1]}},$temp[0];
			shift(@temp);
			my $id1=shift(@temp);
			$taxid_tax{$id1} = join("\t",@temp);
		}
		
	}
}
close(IN);

for my $id (sort keys %taxid_tax){
	$j++;
	my $gi_range = join(",",@{$taxid_gi{$id}});
	print "Now is $id, $j in ".scalar(keys %taxid_tax)."...\n";
	my @taxinfo = split("\t",$taxid_tax{$id});	
	my $statement = "UPDATE $table SET `$added[1]`='$id',`$added[2]`='$taxinfo[0]',`$added[3]`='$taxinfo[1]',`$added[4]`='$taxinfo[2]',`$added[5]`='$taxinfo[3]',`$added[6]`='$taxinfo[4]',`$added[7]`='$taxinfo[5]',`$added[8]`='$taxinfo[6]',`$added[9]`='$taxinfo[7]' WHERE GI IN($gi_range)";
	my $sth = $dbh->prepare($statement) or die "Can't prepare $statement: $dbh->errstr\n";
	$sth-> execute();
}
#my @taxinfo = split("\t",$_);
#print "Updating $j records $taxinfo[0] into DB..\n";		
#my $statement = "UPDATE $table SET `$added[$i]`='$taxinfo[$i]' WHERE GI=$taxinfo[0]";
#my $statement = "UPDATE $table SET `$added[1]`='$taxinfo[1]',`$added[2]`='$taxinfo[2]',`$added[3]`='$taxinfo[3]',`$added[4]`='$taxinfo[4]',`$added[5]`='$taxinfo[5]',`$added[6]`='$taxinfo[6]',`$added[7]`='$taxinfo[7]',`$added[8]`='$taxinfo[8]',`$added[9]`='$taxinfo[9]' WHERE GI=$taxinfo[0]";
#my $statement = "insert into $table (GI,`$added[1]`, `$added[2]`, `$added[3]`,`$added[4]` , `$added[5]`, `$added[6]`, `$added[7]`,`$added[8]`, `$added[9]`) values ('$taxinfo[0]','$taxinfo[1]','$taxinfo[2]','$taxinfo[3]','$taxinfo[4]','$taxinfo[5]','$taxinfo[6]','$taxinfo[7]','$taxinfo[8]','$taxinfo[9]') ";
#my $sth = $dbh->prepare($statement) or die "Can't prepare $statement: $dbh->errstr\n";
#$sth-> execute();

&disconnectdb(\$dbh);

sub connectdb() {
	my ($dbhptr,$host,$dbuser,$dbpass,$objdbname)=@_;
	$$dbhptr = DBI->connect("DBI:mysql:$objdbname:$host","root",$dbpass);
	return 1;
}

sub disconnectdb() {
	my $dbhptr=shift(@_);
	$$dbhptr->disconnect;
	return 1;
}