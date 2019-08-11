use Bio::SeqIO;
$output_folder=$ARGV[2];
open(out1,'>',$output_folder . 'denovo_CDS_v1.fa');#CDS file
open(out2,'>',$output_folder . 'denovo_protein_v1.fa');#Protein file
open(out3,'>',$output_folder . 'denovo_intron_v1.fa');#Intron file
my $stream_Trans = Bio::SeqIO->new(-format => 'fasta',-file   => $ARGV[1]);#Original genome file
while(my $result_Trans=$stream_Trans->next_seq()){
	$idthis=$result_Trans->id();
	$seq{$idthis}=$result_Trans->seq();
}
local $/='# end gene';
open(infile,$ARGV[0]);#Gff3 file produced by augustus
@data=<infile>;
#print "$data[1]";

foreach$keys(@data){
	if($keys=~/start gene g(\d+)/){
		$count=$1;
	}
	@arr_introns=undef;
	@arr_introns_end=undef;
	@arr_introns_start=undef;
	if($keys=~/\tstart_codon\t/ and $keys=~/\tstop_codon\t/){
		$keys=~/(NODE.*?)\t/;
		$tag=$1;
		while($keys=~/\tintron\t(\d+)\t(\d+)\t.*?\t(.*?)\t/g){
			$start_intron=$1;
			$end_intron=$2;
			$frame_intron=$3;
			$length_intron=$end_intron+1-$start_intron;
			#if($length_intron<45){print "$length_intron\n";}
			$intron_seq=substr($seq{$tag},$start_intron-1,$length_intron);
			if($frame_intron eq '+'){
				$intron_seq=$intron_seq;
			}else{
				$intron_seq=&reverse_complement($intron_seq);
			}
			push @arr_introns,$intron_seq;
			push @arr_introns_start,$start_intron;
			push @arr_introns_end,$end_intron;
		}
		#print "\n";
		$ok=1;
		$keys=~/\tstart_codon\t(\d+)\t(\d+)\t\.\t(.*?)\t/;
		$num1=$1;
		$num2=$2;
		$frame=$3;
		if($frame eq '+'){
			$start_pos=$num1;
		}else{
			$end_pos=$num2;
		}
		$keys=~/\tstop_codon\t(\d+)\t(\d+)\t\.\t(.*?)\t/;
		$num1=$1;
		$num2=$2;
		$frame=$3;
		if($frame eq '+'){
			$end_pos=$num2;
		}else{
			$start_pos=$num1;
		}
		
		#if($start_pos>$end_pos){$int=$start_pos;$start_pos=$end_pos;$end_pos=$int;}
		#print "$start_pos\t$end_pos\t$frame\t$tag\n";
		if($frame eq '+'){
			$CDS=substr($seq{$tag},$start_pos-1,$end_pos-$start_pos+1);
		}else{
			$CDS=substr($seq{$tag},$start_pos-1,$end_pos-$start_pos+1);
			$CDS=&reverse_complement($CDS);
		}
	}else{
		$ok=0;
	}
	if($ok==1){
		foreach$keyx(keys@arr_introns){
			if($arr_introns[$keyx] ne undef){
				print out3 ">$tag|$arr_introns_start[$keyx]..$arr_introns_end[$keyx]\n$arr_introns[$keyx]\n";
				$seq_this_kill=$arr_introns[$keyx];
				$CDS=~s/$seq_this_kill//;
			}
		}
		#print "$keys";
		$qualified_contigs++;
		if($keys=~/protein sequence \= \[(.*?)\]/ms){
		$proteinx=$1;
		$proteinx=~s/#//g;
		$proteinx=~s/\n//g;
		$proteinx=~s/ //g;
		print out1 ">$tag|CDS|$start_pos..$end_pos|$frame\n$CDS\n";
		print out2 ">$tag|PROTEIN|$start_pos..$end_pos|$frame\n$proteinx\n";
		}
			
	}

}
local $/="\n";

sub reverse_complement{
	my($seqa)=@_;
	my(@seq_split,@new_seq);
	@seq_split=split('',$seqa);
	@seq_split=reverse(@seq_split);
	foreach $value(@seq_split){
		if ($value eq 'A'){
			push @new_seq,'T';
		}elsif($value eq 'T'){
			push @new_seq,'A';
		}elsif ($value eq 'C'){
			push @new_seq,'G';
		}elsif ($value eq 'G'){
			push @new_seq,'C';
		}else {
			push @new_seq,'N';
		}
	}
	$seqa=join('',@new_seq);
	return($seqa);
}
