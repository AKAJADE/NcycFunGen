#read HMMER out file and save to database
use DBI;
use FunctionGalaxy qw($HMMSEARCH_NCBI);
use autodie;
use Getopt::Long;

# my ($input,$emzyme,$input1,$input2,$input3,$input4,$input5,$database1,$output1);
# &GetOptions (
# "in=s" => \$input,
# "em=s" => \$emzyme,
# "i1=s" => \$input1,
# "i2=s" => \$input2,
# "i3=s" => \$input3,
# "i4=s" => \$input4,
# "i5=s" => \$input5,
# "s=s" => \$database1,
# "o1=s" => \$output1,
# );

# my $host = $input1;
# my $dbname = $input2;
# my $dbuser = $input3;
# my $dbpass = $input4;
# my $table = $input5;
# my $obj_hmms_gobal = $input; 
# my $gene_emzyme_file = $emzyme;

# my $fasta_db_file;
# if(!$database1 || $database1=~/none/i){
	# $fasta_db_file = $HMMSEARCH_NCBI;
# }else{
	# $fasta_db_file = $database1;
# }

my $host="server";#server address
my $dbname="database_name";#database name
my $dbuser="user_name";#database user name
my $dbpass="password";#database password
#my $table="seqs";
my $table = "table_name";#database table name
my $obj_hmms_gobal = "hmmsearch.out"; #HMMSEARCH out file
my $gene_emzyme_file = "gene_information.txt";#gene information
my $fasta_db_file = "combined_protein.fasta";#protein sequences from Genbank
my $output1="output_file.txt";#output file

#================================#
#read hmmsearch results and save to db
#================================#
my ($genes, $acc2gene, $acc2values, $acc2level) = &read_hmmsearch($obj_hmms_gobal, $gene_emzyme_file, $fasta_db_file);

#===============================#
&write_hmm_2db();


#print "Press <Enter> to continue...";
#<>;

#########################################


sub write_hmm_2db(){
	
	#my $test = "hmm_test.txt";
	#open(TEST,">$test") or die "!!!";
	my @gene_list = @$genes;
	my %acc2genes  = %$acc2gene;
	my %acc2hmm    = %$acc2values;
	my %acc2levels = %$acc2level;
	my $acc_num;
	
	
	my $table_seq = $table;#输出的表
	&connectdb(\$dbh,$host,$dbuser,$dbpass,$dbname);
	
	#print "Saving data to DB...\n";
	foreach my $gene(@gene_list){
		my $i=0;
		#print O"\n$gene sequences\n";
		foreach my $acc(keys %acc2genes){
			if($acc2genes{$acc} eq $gene and $acc ne ""){
				$i++;
				my @hmm_rep = @{$acc2hmm{$acc}};
				#$acc=~/\w+\.\d*/;
				#$gi=~/gi\|(\d+)\|/;
				#$acc_num = $1;
				#$acc_num = $acc;
				#print "$acc\t$acc_num\t$i\n";#test
				#print TEST scalar(keys %acc2genes)."\t$gene\t"."\t$acc\t$acc_num\n";#test
				my $statement="insert into $table_seq (gene_name, full_evalue, full_score, full_bias, best_evalue, best_score, best_bias, dom_exp, dom_N, prot_acc, annotation, Level) values (?,?,?,?,?,?,?,?,?,?,?,?) ";
				my $ssth = $dbh->prepare($statement) or die "Can't prepare $statement: $dbh->errstr\n";
				my $rv = $ssth->execute($gene, $hmm_rep[0], $hmm_rep[1], $hmm_rep[2], $hmm_rep[3], $hmm_rep[4], $hmm_rep[5], $hmm_rep[6], $hmm_rep[7], $hmm_rep[8], $hmm_rep[9], $acc2levels{$acc} ) or die "can't execute the query: $ssth->errstr";
				
			}
		}
	}
	&disconnectdb(\$dbh);		
	#close(TEST);
}



sub read_hmmsearch(){
	my ($in_file, $list_file, $db_file)=@_;
	
	my (%gene2key_words,%gene2neg_words);
	open(F,"<$list_file") or die "Couldn't open $list_file to read";#information files
	<F>;
	while(<F>){
		chomp();
		s/\"//g;
		my ($gene, $posit_info, $neg_info) = split(/\t/, $_);
		#my ($gene, $antibiotics, $temp, $emzyme, $anti_type) = split(/\t/, $_);
		##$gene2emzyme{$gene} = $emzyme;
		# my @anti = split(/,/, $antibiotics);
		# my @emzy = split(/\s+/, $emzyme);
		# push @anti, $anti_type if (defined($anti_type) and $anti_type ne "");
		# push @anti, @emzy;
		# @anti = map{lc} @anti;
		my @positive = split(/,/, $posit_info);
		my @negative = split(/,/, $neg_info);
		@positive = map{lc} @positive;
		@negative = map{lc} @negative;
		$gene2key_words{$gene} = \@positive;
		$gene2neg_words{$gene} = \@negative;
	}
	close(F);
	#print join(";", @{$gene2key_words{"mefa"}});
	#die;
	
	my (%acc2anno);
	open(F,"<$db_file") or die "Couldn't open $db_file to read";
	<F>;
	#print "Reading all seqs from Genbank DB file...\n";
	while(<F>){
		chomp();
		if(/^>(\S+) (.*)$/){
			$acc2anno{$1} = $2;
		}
	}
	close(F);
		
	my (@genes, %acc2gene, %acc2evalue, %acc2values, %acc2level);
	open(File_in,"<$in_file");
	
	my $c = 1e-5;
	
	#print "Reading HMMER results...\n";
	while(<File_in>){
		my $str_temp=$_;
		until($str_temp=~/^Query:\s+/ or eof(File_in)){
			$str_temp=<File_in>;
		}
		if(eof(File_in)){last;}
		$str_temp=~/^Query:\s+(\w+)\s+\[/;
		my $gene = $1;
		#die "$gene";
		push @genes, $gene;
		<File_in>;<File_in>;<File_in>;<File_in>;
		$str_temp=<File_in>;
		
		until($str_temp eq "\n"){
			chomp($str_temp);
			if($_!~/------/){
				$str_temp=~s/\n+?//g;
				my @values = split(/\s+/, $str_temp);
				#print join(",",@values);die;
				shift(@values);
				my $e_value = 0+$values[0];
				my $a = $#values;
				for(my $i=10;$i<=$a;$i++){
					pop(@values);
				}
				my $seq_acc = $values[8];
				$values[9] = $acc2anno{$seq_acc};
				my $annotation = lc($values[9]);
							
				#die "$annotation, $gene";
				if ($annotation and $e_value<=$c) {
					my @anno = split(/\s+|\-/, $annotation);
					my @neg_words = @{$gene2neg_words{$gene}};
					#if($#neg_words == -1){die "Couldn't find negtive words for $gene";}
					my %shared_neg;
					foreach $e(@neg_words){
						if($annotation=~/$e/i){
						$shared_neg{$e} = 1;
						}
					}		
					my @shared_neg_words = keys %shared_neg;
					#print "$annotation\n", join(",", @neg_words),"\n",join(",", @shared_neg_words),"\n";
					if($#shared_neg_words >= 0){
						if(!exists($acc2evalue{$seq_acc})){
							$acc2evalue{$seq_acc} = $e_value;
							$acc2gene{$seq_acc} = $gene;
							$acc2values{$seq_acc} = \@values;
							$acc2level{$seq_acc} = "Homologous";
							print "$seq_acc\n", join(",", @neg_words),"\n",join(",", @shared_neg_words),"\n";
						} elsif($acc2evalue{$seq_acc} > $e_value){
							$acc2evalue{$seq_acc} = $e_value;
							$acc2gene{$seq_acc} = $gene;
							$acc2values{$seq_acc} = \@values;
							$acc2level{$seq_acc} = "Homologous";
							}
					}
				}
				if($annotation=~/$gene\d*\s+/i and $e_value<=$c and $acc2level{$seq_acc} ne "Homologous"){
					$acc2evalue{$seq_acc} = $e_value;
					$acc2gene{$seq_acc} = $gene;
					$acc2values{$seq_acc} = \@values;
					$acc2level{$seq_acc} = "High_gene";
				}elsif($annotation=~/$gene\d*\s+/i and $e_value < 0.01 and $acc2level{$seq_acc} ne "Homologous"){
					$acc2evalue{$seq_acc} = $e_value;
					$acc2gene{$seq_acc} = $gene;
					$acc2values{$seq_acc} = \@values;
					$acc2level{$seq_acc} = "Low_gene";
				}elsif($e_value<=$c and $acc2level{$seq_acc} ne "Homologous"){
					#die "good";
					my @anno = split(/\s+|\-/, $annotation);
					my @key_words = @{$gene2key_words{$gene}};
					#if($#key_words == -1){die "Couldn't find key words for $gene";}
					#my @shared_words = grep( $anno{$_}, @key_words );
					my (%shared, %union);
					foreach $e(@key_words){
						if($annotation=~/$e/i){
							$shared{$e} = 1;
						}
					}
					
					my @shared_words = keys %shared;
					if($#shared_words >= 0){
						# if($annotation=~/hypothetical/i){
							# print "$annotation\n", join(",", @key_words),"\n",join(",", @shared_words),"\n";
							# exit;
						# }
						
						if(!exists($acc2evalue{$seq_acc})){
							$acc2evalue{$seq_acc} = $e_value;
							$acc2gene{$seq_acc} = $gene;
							$acc2values{$seq_acc} = \@values;
							$acc2level{$seq_acc} = "High";
						}elsif($acc2evalue{$seq_acc} > $e_value and $acc2level{$seq_acc} ne "High_gene" and $acc2level{$seq_acc} ne "Homologous"){
							$acc2evalue{$seq_acc} = $e_value;
							$acc2gene{$seq_acc} = $gene;
							$acc2values{$seq_acc} = \@values;
							$acc2level{$seq_acc} = "High";
						}
					}
					else{
						if(!exists($acc2evalue{$seq_acc})){
							$acc2evalue{$seq_acc} = $e_value;
							$acc2gene{$seq_acc} = $gene;
							$acc2values{$seq_acc} = \@values;
							$acc2level{$seq_acc} = "Middle";
						}elsif($acc2evalue{$seq_acc} > $e_value and $acc2level{$seq_acc} ne "High_gene" and $acc2level{$seq_acc} ne "Homologous"){
							$acc2evalue{$seq_acc} = $e_value;
							$acc2gene{$seq_acc} = $gene;
							$acc2values{$seq_acc} = \@values;
							$acc2level{$seq_acc} = "Middle";
						}
					}
				}elsif($e_value<0.01){
					my @anno = split(/\s+|\-/, $annotation);
					my @key_words = @{$gene2key_words{$gene}};
					#if($#key_words == -1){die "Could find key words for $gene";}
					#my @shared_words = grep( $anno{$_}, @key_words );
					my (%shared, %union);
					foreach $e(@key_words){
						if($annotation=~/$e/i){
							$shared{$e} = 1;
						}
					}
					my @shared_words = keys %shared;
					if($#shared_words >= 0){
						if(!exists($acc2evalue{$seq_acc})){
							$acc2evalue{$seq_acc} = $e_value;
							$acc2gene{$seq_acc} = $gene;
							$acc2values{$seq_acc} = \@values;
							$acc2level{$seq_acc} = "Low";
						}elsif($acc2evalue{$seq_acc} > $e_value and $acc2level{$seq_acc} ne "Low_gene" and $acc2level{$seq_acc} ne "Homologous"){
							$acc2evalue{$seq_acc} = $e_value;
							$acc2gene{$seq_acc} = $gene;
							$acc2values{$seq_acc} = \@values;
							$acc2level{$seq_acc} = "Low";
						}
					}
				}
				
			}
			
			$str_temp = <File_in>;
		}
	}
	close(File_in);
	
	#@genes, %acc2gene, %acc2values, %acc2level;
	open(O,">$output1");
	foreach my $gene(@genes){
		print O"\n$gene sequences\n";
		foreach my $acc(keys %acc2gene){
			if($acc2gene{$acc} eq $gene){
				print O"$acc\t",join("\t",@{$acc2values{$acc}}),"\t$acc2level{$acc}\n";
			}
		}
	}
	close(O);
	
	return(\@genes, \%acc2gene, \%acc2values, \%acc2level);
}

sub connectdb() {
	my ($dbhptr,$host,$dbuser,$dbpass,$objdbname)=@_;
	$$dbhptr = DBI->connect("DBI:mysql:$objdbname:$host",$dbuser,$dbpass);
	return 1;
}

sub disconnectdb() {
	my $dbhptr=shift(@_);
	$$dbhptr->disconnect;
	return 1;
}
