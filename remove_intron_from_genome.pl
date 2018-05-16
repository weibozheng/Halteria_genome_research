use Bio::SeqIO;
use Bio::SearchIO;
open(out,'>C:\Users\wbzhe\Stentor\pool\pool_types\Stentor_genome_SPCA_v4_minor_p3_intronfree.fasta');
my $stream_Trans = Bio::SeqIO->new(-format => 'fasta',-file   => 'C:\Users\wbzhe\Stentor\pool\pool_types\Stentor_genome_SPCA_v4_minor_p3.fasta');#contigs
while(my $result_Trans=$stream_Trans->next_seq()){
	$idthis=$result_Trans->id();
	$idthis=~/(NODE.*)/;#title feature of fasta file
	$seq{$1}=$result_Trans->seq();
}
my $stream_Trans = Bio::SeqIO->new(-format => 'fasta',-file   => 'C:\Users\wbzhe\Stentor\intron\introns.fa');#intron file
while(my $result_Trans=$stream_Trans->next_seq()){
	$idthis=$result_Trans->id();
	$idthis=~/(.*)\|(\d+)\.\.(\d+)/;
	$seq2=$result_Trans->seq();
	$tag=$1;
	#$intron_start=$2;
	#$intron_end=$3;
	$seq{$tag}=~s/$seq2//;
	
}
foreach my $key (sort {&rule} keys %seq){
	print out ">$key\n$seq{$key}\n";
}

sub rule{
	my $count_a=length($seq{$b});
	my $count_b=length($seq{$a});
	
	if($count_a>$count_b){return 1;}else{return -1;}
	
}
