#similarity calculation between two aligned sequences
#use strict;
use Getopt::Long;

my ($input,$align,$output1,$output2);
&GetOptions (
"in=s" => \$input,
"a=s" => \$align,
"o1=s" => \$output1,
"o2=s" => \$output2,
);

if($align eq "muscle"){
	open (N,"$input") or die "!!"; 
	open (S,">muscleline.txt") or die "!!";
	my $mark = 0; 
	while(my $line = <N>){
		chomp($line);
		if($line =~ /^>/){
			if($mark eq 0){print S "$line\n";$mark = 10}else{print S "\n$line\n"}
		}else{print S "$line"}
	}
	close (N);
	close (S);

	open (N,"muscleline.txt") or die "!!";
}elsif($align eq "standard"){
	open (N,"input") or die "!!";
}

my @allseqname;
my $seqname;
my %h;#sequence names and sequences
while(my $line = <N>){
	chomp($line);
	if($line =~ /^>/){$seqname = $line}else{
		$h{$seqname} = $line;
		@allseqname = (@allseqname, $seqname);
	}
}
close (N);

#seuqence similarity calculation
open (O,">$output1") or die "!!";
open (P,">$output2") or die "!!";
my $i = 0;
my $seqnum = @allseqname;#sequence number
for(@allseqname){
#print "$_\n";
	my $seq1 = $h{$_};
	my $j = $i+1;
	while($j < $seqnum){
		my $seq2 = $h{$allseqname[$j]};
		my $simmain = &similarity($seq1,$seq2);
		print O "$allseqname[$i]\t$h{$allseqname[$i]}\t$allseqname[$j]\t$h{$allseqname[$j]}\t$simmain\n";
		print P "$allseqname[$i]\t$allseqname[$j]\t$simmain\n";
		$j++;
	}
	$i++
}
close (O);
close (P);



#similarity calculation between two aligned sequences
sub similarity{
	my @seq11 = split("",$_[0]);
	my @seq22 = split("",$_[1]);
	my $len = @seq11;
	my $i = 0;
		
	while($i < $len){
		if ($seq11[$i] eq "-" && $seq22[$i] eq "-"){
			undef $seq11[$i]; 
			undef $seq22[$i];
		}$i++
	}	
	my $i = 0;
	my @seq111;
	for(@seq11){if($_ ne ""){push(@seq111,$_)}}
	my @seq222;
	for(@seq22){if($_ ne ""){push(@seq222,$_)}}
	
	$len = @seq111;
	my $samebasenum = 0;
	$i = 0;
	for (@seq111){
		if ($_ eq $seq222[$i]){$samebasenum++}
		$i++
	}
	my $simsub = $samebasenum/$len;
}