#for use: perl download_CDS_Genbank_2022.pl
#!/usr/local/bin/perl -w

use LWP::Simple;
use Getopt::Long;
use Bio::Factory::FTLocationFactory;
use Bio::DB::GenPept;
use Bio::DB::GenBank;
use Bio::Tools::CodonTable;
use DBI;

my $host="server_name";#server_name
my $dbname="NcycFunGen";#database name
my $dbuser="user_name";#user name of serbver
my $dbpass="password";#password of server
my $table="pmoA";#table name in database

my $dbh;
&connectdb(\$dbh,$host,$dbuser,$dbpass,$dbname);

my ($gene,$work_dir);
#my $result=&GetOptions("g=s"  => \$gene,"dir=s"  => \$work_dir);
$gene="pmoA";#gene name
$work_dir='work_dir';#work_dir

#my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
my $utils = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils";#GenBank link

chdir $work_dir;
my $file_log=$gene.".log";

#================================#
#main
#read the prot_acc in our db
my ($nuc_acc, $prot_seq)=&get_prot_acc();
my @a_nuc_acc=sort @{$nuc_acc};
#================================#
#find nucleotide acc number in NCBI#
my ($cds_acc)=&get_nuc_acc(@a_nuc_acc);
#================================#
my($cds_seq,$cds2nuc_acc)=&get_CDS_seq(\@a_nuc_acc,$prot_seq);

#================================#
#save all the original data, CDS_seq and nuc_acc to db
#&save_CDS();
#================================#


&disconnectdb(\$dbh);	


sub uniq {
    my %count;
    foreach my $e(@_){
		$count{$e}++;
	}
	my @diff;
	foreach my $e (keys %count) {
		if ($count{$e} == 1) 
		{
			push @diff, $e;
		}
	}
	return @diff;
}


sub get_nuc_acc(){
	open(File_log,">>$file_log");
	print File_log"==========================\n";
	print File_log"Now obtain the nucleic acc list in NBCI.\n";
	my @prot_acc=@_;
	#print "$#prot_acc\n";
	my (%cds_accs);
	for($i=0;$i<=$#prot_acc;$i+=200)#based on internet speed
	{
		my $prot_acc_200="";
		for(my $j=0;$j<200;$j++){
			if (exists($prot_acc[$i+$j])){
			   $prot_acc_200.="id=$prot_acc[$i+$j]&";
			}
		}
		$prot_acc_200=substr($prot_acc_200,0,length($prot_acc_200)-1);
		#print "$prot_acc_200\n";
		my $db  = "Protein";
    	my $cmd = 'neighbor';
    	my $elink = "$utils/elink.fcgi?dbfrom=protein&db=nuccore&cmd=$cmd&";
    	my $elink_result = get($elink . $prot_acc_200);
    	print File_log"$elink$prot_acc_200\n";
    	#print File_log"$elink_result\n";
		#die;
		my $pos_s=index($elink_result,"<LinkSet>");
		my $pos_e=index($elink_result,"</LinkSet>");
		my $str_elink=substr($elink_result,$pos_s,$pos_e-$pos_s);
		while($str_elink ne ""){
			if($str_elink =~ m|<IdList>.*<Id>(\w+\.\d*)</Id>.*</IdList>.*<Link>.*<Id>(\w+\.\d*)</Id>.*</Link>|s){
				$protein=$1;
				$nuc=$2;
				#$protein=substr($protein,7,length($protein)-14);
				#$nuc=substr($nuc,8,length($nuc)-16);
				$cds_accs{$protein}=$nuc;
				print "$protein, $nuc\n";
			}else{
				print "$str_elink";
			}
			$pos_s=index($elink_result,"<LinkSet>",$pos_s+1);
			$pos_e=index($elink_result,"</LinkSet>",$pos_e+1);
			$str_elink=substr($elink_result,$pos_s,$pos_e-$pos_s);
		}
	}
	my $i_aa_cds=keys %cds_accs;
	my $i_aa=scalar(@prot_acc);
	print File_log"Total $i_aa_cds proteins were found relative nucleotides acc numbers.\n";
	if($i_aa!=$i_aa_cds){
		print File_log"Some proteins could not find relative nt acc:";
		foreach my $aa_acc(@prot_acc){
			if(!exists($cds_accs{$aa_acc})){
				print File_log"$aa_acc ";
			}
		}
		print File_log"\n";
	}
	print File_log"Nucleic gi list was saved.\n";
	close(File_log);
	return (\%cds_accs);
}

sub get_CDS_seq(){
	my ($aa_acc,$prot_seq)=@_;
	open(LOG,">>$file_log");
	my %aa_seq=%{$prot_seq};
	my $gp = Bio::DB::GenPept->new;
    my $gb = Bio::DB::GenBank->new;
    # factory to turn strings into Bio::Location objects
    my $loc_factory = Bio::Factory::FTLocationFactory->new;
	my $myCodonTable  = Bio::Tools::CodonTable -> new ( -id => 11 );
	my (%cds_seqs,@false_acc,%cds2nuc_acc);
	my $i=0;
	my $i_total=scalar(@{$aa_acc});
	my $temp_file=$gene."_CDS.temp";
	if(-e $temp_file){
		 open(TEMP,"<$temp_file");
		 while(<TEMP>){
		 	chomp;
			my ($prot_acc,$prot_seq, $nuc_acc,$cds_seq)=split(/\t/,$_);
			$cds_seqs{$prot_acc}=$cds_seq;
			$cds2nuc_acc{$prot_acc}=$nuc_acc;
			#shift(@{$aa_acc});
			
			$i++;
		 }
		 close(TEMP);
		 my @aa_acc_list = @$aa_acc;
		 my @pre_acc = sort keys %cds_seqs;
		 my (%shared, %union);
		 foreach $e(@aa_acc_list, @pre_acc){
			$union{$e}++ && $shared{$e}++
		 }
		 my @left = keys %union;
		 $aa_acc = \@left;
	}
		
	open(TEMP,">>$temp_file");
    FOR1:foreach my $protein_acc(@{$aa_acc}){
		my $b_find=0;
		my $i_times=0;
        Label1:my $prot_obj; 
		if($prot_obj = $gp->get_Seq_by_id($protein_acc)){
    		if(!defined($prot_obj)){
    			$i_times++;
    			if($i_times>=3){
    				print LOG"We tried to 3 times in $protein_acc, but it doesn't work.\n";
    				next FOR1;
    			}else{
    				goto Label1;
    			}
    		}
		}else{
			$i_times++;
			if($i_times>=3){
				print LOG"We tried to 3 times in $protein_acc, but it doesn't work.\n";
				next FOR1; 
			}else{
				goto Label1;
			}
		}
        foreach my $feat ( $prot_obj->top_SeqFeatures ) {
		   if ($feat->primary_tag eq 'CDS' ) {
               my @coded_by = $feat->each_tag_value('coded_by');
			   if($coded_by[0]=~/^complement\((.*)\)/){
			   		$coded_by[0]=$1;
			   }
               my ($nuc_acc,$loc_str) = split /\:/, $coded_by[0];
			   #$nuc_acc=~s/\.\d*//;
			   if($nuc_acc=~/join\(/){
			   		$nuc_acc=~s/join\(//;
			   }
			   $cds2nuc_acc{$protein_acc}=$nuc_acc;
			   my ($cds_start,$cds_stop);
			   if($loc_str=~/\.\.>/){
			   		$loc_str=~/(\d+).*\.\.>(\d+)/;
			   		($cds_start,$cds_stop)=($1,$2);
			   }else{
			   		$loc_str=~/(\d+).*\.\.(\d+)/;
			   		($cds_start,$cds_stop)=($1,$2);
			   }
			   if(!defined($cds_start) or !defined($cds_start)){                                                                  
			   		print LOG"The CDS position of protein $protein_acc is wrong!\n";
					next FOR1;
			   }
               my $db  = "nuccore";
    		   my $type = "fasta";
    		   my $elink = "$utils/efetch.fcgi?db=$db&id=$nuc_acc&rettype=$type&seq_start=$cds_start&seq_stop=$cds_stop";
    		   my $elink_result = get($elink);			   
			   (my $cds_seq_all=$elink_result)=~s/>.*?\n//s;
			   $cds_seq_all=~s/\n//mg;
			   my $cds_seq;
			   if($cds_seq_all=~/(.*)>/){
			   		$cds_seq=$1;
			   }else{
			   		$cds_seq=$cds_seq_all;
			   }
			   my $len=length($cds_seq);
			   if(exists($aa_seq{$protein_acc})){
				   if($len>0){
						# my $cds_seq15=substr($cds_seq,3,15);
						# my $aa1 = $myCodonTable->translate($cds_seq15);
						# my $ori_aa_seq=substr($aa_seq{$protein_acc},1,5);
						# if($ori_aa_seq ne $aa1){
							# my $cds_seq2=reverse $cds_seq;
							# $cds_seq2=~tr/ATGCatgc/TACGtacg/;
							# my $cds_seq15=substr($cds_seq2,3,15);
							# my $aa2 = $myCodonTable->translate($cds_seq15);
							# if($ori_aa_seq eq $aa2){
								# $cds_seq=$cds_seq2;
							# }else{
								# my $cds_seq15=substr($cds_seq,4,15);
								# my $aa1 = $myCodonTable->translate($cds_seq15);
								# $cds_seq15=substr($cds_seq,2,15);
								# my $aa2 = $myCodonTable->translate($cds_seq15);
								# if($ori_aa_seq eq $aa1 or $ori_aa_seq eq $aa2){
								# }else{
									# print LOG"The CDS seq of protein $protein_acc is wrong!\n";
									# next FOR1;
								# }
							# }
						# }
						my $tot_shift=0;
						for(my $shift_nucl=0;$shift_nucl<=5;$shift_nucl++){
							my $cds_seq15=substr($cds_seq,$shift_nucl,15);
							my $aa1 = $myCodonTable->translate($cds_seq15);
							for(my $shift_prot=0;$shift_prot<=2;$shift_prot++){
								my $ori_aa_seq=substr($aa_seq{$protein_acc},$shift_prot,5);
								if($ori_aa_seq ne $aa1){
									$tot_shift++;
									#print "$protein_acc\t$tot_shift\t$ori_aa_seq\t$aa1\t$cds_seq15\n";
									my $cds_seq2=reverse $cds_seq;
									$cds_seq2=~tr/ATGCatgc/TACGtacg/;
									my $cds_seq15=substr($cds_seq2,$shift_nucl,15);
									my $aa2 = $myCodonTable->translate($cds_seq15);
									if($ori_aa_seq eq $aa2){
										$cds_seq=$cds_seq2;
										last;
									}else{
										$tot_shift++;
										#print "$protein_acc\t$tot_shift\t$ori_aa_seq\t$aa2\t$cds_seq15\n";
									}
								}else{
									last;
								}
							}
						}
						if($tot_shift eq "12"){
							print LOG"The CDS seq of protein $protein_acc is wrong!\n";
							next FOR1;
						}
						$i++;
						$b_find=1;
						$cds_seqs{$protein_acc}=$cds_seq;
						print TEMP"$protein_acc\t$aa_seq{$protein_acc}\t$nuc_acc\t$cds_seq\n";
						#print "Find $protein_acc in $nuc_acc. CDS length=$len.($i/$i_total)\t$cds_seq\n";	
						print LOG"Find $protein_acc in $nuc_acc. CDS length=$len.($i/$i_total)\n";
				   }else{
						print LOG"Wrong in $protein_acc in $nuc_acc.\n";
				   }
				}else{
					print LOG"No protein sequences for $protein_acc and couldn't assess CDS accuracy.\n";
				}
           }
        }
		if($b_find==0){push @false_acc,$protein_acc;}
	}
	close(TEMP);
	print LOG"Total $i CDS sequences were found.\n";
	if($#false_acc!=-1){
		print LOG"Protein GI number:",join(",",@false_acc),"can't find any cds sequences.\n";
	}
	print LOG"Finish read CDS sequences from NCBI.\n";
	close(LOG);
	return (\%cds_seqs,\%cds2nuc_acc);
}

sub save_CDS(){
	open(File_log,">>$file_log");
	print File_log "==========================\n";
	print File_log "Now save nucleic gi numbers and CDS sequences to db.\n";
	
	my $temp_file=$gene."_CDS.temp";
	if(-e $temp_file){
		 open(TEMP,"<$temp_file");
		 
		 my $statement="select ID, prot_acc from $table where gene_name=?";
		 my $ssth = $dbh->prepare($statement)  or die "Can't prepare $statement: $dbh->errstr\n";
		 my $rv = $ssth->execute($gene) or die "can't execute the query $statement: $ssth->errstr";
		 my %prot_acc2id;
		 while (my ($id, $prot_acc_0) = $ssth->fetchrow_array()) {
			$prot_acc2id{$prot_acc_0} = $id;
		 }
		 
		 while(<TEMP>){
		 	chomp;
			my ($prot_acc,$prot_seq, $nuc_acc,$cds_seq)=split(/\t/,$_);
			
			my $statement="update $table set Prot_seq=\"$prot_seq\", Nuc_acc=\"$nuc_acc\", CDS_seq =\"$cds_seq\" where ID=$prot_acc2id{$prot_acc}";
			my $ssth = $dbh->prepare($statement)  or die "Can't prepare $statement: $dbh->errstr\n";
			my $rv = $ssth->execute() or die "can't execute the query $statement: $ssth->errstr";
			
			$i++;
		 }
		 close(TEMP);
	}else{
		die "$gene.temp is not existed";
	}
	
	print File_log "Total $i sequences saved.\nFinish saving to db.\n";
	#my $temp_file=$gene."_CDS.temp";
	#unlink($temp_file);
	close(File_log);
}

sub get_prot_acc(){
	open(File_log,">>$file_log");
	print File_log"==========================\n";
	print File_log"Now obtain the protein acc_number list in our database.\n";
	
	my $statement="select prot_acc,gene_name from $table where gene_name=?";
	my $ssth = $dbh->prepare($statement)  or die "Can't prepare $statement: $dbh->errstr\n";
	my $rv = $ssth->execute($gene) or die "can't execute the query: $ssth->errstr";
	
	my @prot_acc;
	while (my ($id, $gene_name) = $ssth->fetchrow_array()) {
		if(defined($id) and $id ne ""){push @prot_acc, $id;}
	}
	
	$prot_acc_list=join(',',@prot_acc);
	print File_log $prot_acc_list."\n";
	print File_log "Total ".($#prot_acc+1)." prot_accs.\n";
	print File_log"Protein prot_acc list was saved.\n";
	
	#===obtain prot seq from Genbank=======#
	my (%cds_accs);
	for($i=0;$i<=$#prot_acc;$i+=200)#based on internet speed
	{
		my $prot_acc_200="";
		for(my $j=0;$j<200;$j++){
			if (exists($prot_acc[$i+$j])){
			   $prot_acc_200.="$prot_acc[$i+$j],";
			}
		}
		$prot_acc_200=substr($prot_acc_200,0,length($prot_acc_200)-1);
		my $elink = "$utils/efetch.fcgi?db=sequences&id=$prot_acc_200&rettype=fasta&retmode=text";
		my $elink_result = get($elink);
		print File_log"$elink\n";
		
		while($elink_result=~/>(\w+\.\d*)(.*?)\n(.*?)\n\n/sg){
			my $pro_seq = $3;
			my $nuc_acc_id = $1;
			$pro_seq=~s/\n//mg;
			$cds_accs{$nuc_acc_id} = $pro_seq;
			print "$nuc_acc_id,$pro_seq\n";
		}
	}
	
	my @found_seqs = sort keys %cds_accs;
	my @nofound = uniq(@found_seqs, @prot_acc);
	if(scalar(@nofound)>0){
		my $i_p = scalar @nofound;
		print File_log "Couldn't find $i_p protein sequences for ", join(",", @nofound), "\n";
	}else{
		print File_log "Found all protein sequences from Genebank\n";
	}
	close(File_log);
	#die;
	return (\@prot_acc,\%cds_accs);
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


