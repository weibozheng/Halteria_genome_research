#!/usr/bin/perl
use Bio::SeqIO;
use Bio::SearchIO;
open(out,'>C:\Users\wbzhe\Stentor\Stentor_genome_SPCA_v4_intronfree_CDS_INR.fasta');#predicted genes
open(out1,'>C:\Users\wbzhe\Stentor\Stentor_genome_SPCA_v4_intronfree_intron_INR.fasta');#potential introns which are not detected by rna-seq
my $stream_Trans = Bio::SeqIO->new(-format => 'fasta',-file   => 'C:\Users\wbzhe\Stentor\Stentor_genome_SPCA_v4_intronfree.fasta');#original contig file
while(my $result_Trans=$stream_Trans->next_seq()){
	$idthis=$result_Trans->id();
	$seq{$idthis}=$result_Trans->seq();
}
my $searchin = Bio::SearchIO->new(-format => 'blastxml', -file   => 'C:\Users\wbzhe\Stentor\Stentor_SPCA_v4_intronfree_blastx_swissprot_v2.xml');#blastxml file
while( my $query = $searchin->next_result) {
	$query_id=$query->query_description();
	$query_length=$query->query_length();
	$hsp_length=0;
	$identity=0;
	$new_hit=1;
	while( my $hit = $query->next_hit()){
		if(1){
		$new_hit=0;
		$hit_desc=$hit->name();
		$hit_length=$hit->length();
		$best_start=10000000;
		$best_end=0;
		@arr_Start=undef;
		@arr_End=undef;
			while(my $hsp=$hit->next_hsp()){
				$start=$hsp->start('hit');
				$end=$hsp->end('hit');
				$start_q=$hsp->start('query');
				$end_q=$hsp->end('query');
				$identity=$hsp->num_identical();
				$strand=$hsp->strand('query');
				$hit_true_length=$end-$start+1;
				if($start_q<=$best_start){$best_start=$start_q;}
				if($end_q>=$best_end){$best_end=$end_q;}
				push @arr_Start,$start_q;
				push @arr_End,$end_q;#recording all positions
			}
			@arr_Start=sort {$a<=>$b} @arr_Start;
			@arr_End=sort {$a<=>$b} @arr_End;
			$size=@arr_Start;
			#print "$size\n";
			@arr_introns=undef;
			@arr_introns_end=undef;
			@arr_introns_start=undef;
			foreach$i(keys@arr_Start){
				if($arr_Start[$i] != undef){
					if(($i+1)<$size){
						$inter_start=$arr_End[$i]+1;
						$inter_end=$arr_Start[$i+1]-1;
						$intron_ok=0;
						if($inter_end<$inter_start){
							if(($arr_End[$i+1]-$arr_Start[$i]+1)%3!=0){
								for($j_1=$inter_end;$j_1<=$inter_start;$j_1++){
									$local_A=substr($seq{$query_id},$j_1-1,2);
									if($local_A eq "GT"){
										for($k_1=$inter_start;$k_1>=$inter_end;$k_1--){
											$local_B=substr($seq{$query_id},$k_1-2,2);
											if($local_B eq "AG"){
												$distance_intron=$k_1-$j_1+1;
												$seq_intron=substr($seq{$query_id},$j_1-1,$distance_intron);
												if(($arr_End[$i+1]-$arr_Start[$i]+1-$distance_intron)%3==0 and $distance_intron>6){
													$intron_ok=1;
													last;
												}else{
													$intron_ok=0;
												}
											}
										}
										if($intron_ok==1){
											last;
										}
									}
									if($local_A eq "CT"){
										for($k_1=$inter_start;$k_1>=$inter_end;$k_1--){
											$local_B=substr($seq{$query_id},$k_1-2,2);
											if($local_B eq "AC"){
												$distance_intron=$k_1-$j_1+1;
												$seq_intron=substr($seq{$query_id},$j_1-1,$distance_intron);
												if(($arr_End[$i+1]-$arr_Start[$i]+1-$distance_intron)%3==0 and $distance_intron>6){
													$intron_ok=1;
													last;
												}else{
													$intron_ok=0;
												}
											}
										}
										if($intron_ok==1){
											last;
										}
									}
								}
								if($intron_ok==1){
									push @arr_introns,$seq_intron;
									push @arr_introns_start,$j_1;
									push @arr_introns_end,$k_1;
									#print "$query_id\t$j_1\t$k_1\t$seq_intron\t$arr_Start[$i]\t$arr_End[$i]\t$arr_Start[$i+1]\t$arr_End[$i+1]\t";
								}else{
									#print "break";
								}
								#<stdin>;
							}else{
								$intron_ok=1;
							}
						}else{
							if(($arr_End[$i+1]-$arr_Start[$i]+1)%3!=0){
								for($j_1=$inter_start;$j_1<=$inter_end;$j_1++){
									$local_A=substr($seq{$query_id},$j_1-1,2);
									if($local_A eq "GT"){
										for($k_1=$inter_end;$k_1>=$inter_start;$k_1--){
											$local_B=substr($seq{$query_id},$k_1-2,2);
											if($local_B eq "AG"){
												$distance_intron=$k_1-$j_1+1;
												$seq_intron=substr($seq{$query_id},$j_1-1,$distance_intron);
												if(($arr_End[$i+1]-$arr_Start[$i]+1-$distance_intron)%3==0 and $distance_intron>6){
													$intron_ok=1;
													last;
												}else{
													$intron_ok=0;
												}
											}
										}
										if($intron_ok==1){
											last;
										}
									}
									if($local_A eq "CT"){
										for($k_1=$inter_end;$k_1>=$inter_start;$k_1--){
											$local_B=substr($seq{$query_id},$k_1-2,2);
											if($local_B eq "AC"){
												$distance_intron=$k_1-$j_1+1;
												$seq_intron=substr($seq{$query_id},$j_1-1,$distance_intron);
												if(($arr_End[$i+1]-$arr_Start[$i]+1-$distance_intron)%3==0 and $distance_intron>6){
													$intron_ok=1;
													last;
												}else{
													$intron_ok=0;
												}
											}
										}
										if($intron_ok==1){
											last;
										}
									}
								}
								if($intron_ok==1){
									push @arr_introns,$seq_intron;
									push @arr_introns_start,$j_1;
									push @arr_introns_end,$k_1;
									#print "$query_id\t$j_1\t$k_1\t$seq_intron\t$arr_Start[$i]\t$arr_End[$i]\t$arr_Start[$i+1]\t$arr_End[$i+1]\t";
								}else{
									#print "break";
								}
							}else{
								$intron_ok=1;
							}
						}
					}
				}
			}
			if($intron_ok==1){
				
			}else{
				$longest=0;
				foreach$key_i(keys@arr_Start){
					$length_this=$arr_End[$key_i]-$arr_Start[$key_i]+1;
					if($longest<$length_this){
						$longest_start=$arr_Start[$key_i];
						$longest_end=$arr_End[$key_i];
					}
				}
				$best_start=$longest_start;
				$best_end=$longest_end;
			}
			
			
			if(1){#modify_ATG_TGA
					if($strand==-1){
						$start_codon=substr($seq{$query_id},$best_end-3,3);
						#print "$start_codon\n";
						if ($start_codon ne 'CAT'){
						for($i=3;$i<600;$i=$i+3){
							$new_end=$best_end+$i;
							$new_start_codon=substr($seq{$query_id},$new_end-3,3);
							if($new_start_codon eq 'CAT'){
								$dis=$new_end-$best_end;
								$best_end=$new_end;
								
								#print "$dis\n";
								last;
							}
						}
					}
					$stop_codon=substr($seq{$query_id},$best_start-1,3);
					#print "$stop_codon\n";
					if($best_start<600){
						$uplimit=$best_start;
					}else{
						$uplimit=600;
					}
					if ($stop_codon ne 'TCA'){
						for($i=3;$i<$uplimit;$i=$i+3){
							$new_start=$best_start-$i;
							$new_stop_codon=substr($seq{$query_id},$new_start-1,3);
							if($new_stop_codon eq 'TCA'){
								$dis=$best_start-$new_start;
								$best_start=$new_start;
									
								#print "$dis\n";
								last;
							}
						}
					}
					#print "##################################################################################";
					$seq_gene=substr($seq{$query_id},$best_start-1,$best_end-$best_start+1);
					if($intron_ok==1){
						foreach$key_introns(keys@arr_introns){
							if($arr_introns[$key_introns] ne undef){
								$seq_modify=$arr_introns[$key_introns];
								$seq_gene=~s/$seq_modify//;
								$seq_modify=&reverse_complement($seq_modify);
								print out1 ">$query_id\|intron\|$arr_introns_start[$key_introns]..$arr_introns_end[$key_introns]\n$seq_modify\n";
							}
						}
					}
					$seq_gene=&reverse_complement($seq_gene);
					print out ">$query_id\|$hit_desc\|$best_start..$best_end\|-\n$seq_gene\n";
				}else{
                                        if($best_start<600){
                                                $uplimit=$best_start;
                                        }else{
                                                $uplimit=600;
                                        }
					#print "$uplimit\n";
					$start_codon=substr($seq{$query_id},$best_start-1,3);
					if ($start_codon ne 'ATG'){
						for($i=3;$i<$uplimit;$i=$i+3){
							$new_start=$best_start-$i;
							$new_start_codon=substr($seq{$query_id},$new_start-1,3);
							if($new_start_codon eq 'ATG'){
								$dis=$best_start-$new_start;
								$best_start=$new_start;
								last;
							}
						}
					}
					$stop_codon=substr($seq{$query_id},$best_end-3,3);
					if ($stop_codon ne 'TGA'){
						for($i=3;$i<600;$i=$i+3){
							$new_end=$best_end+$i;
							$new_stop_codon=substr($seq{$query_id},$new_end-3,3);
							if($new_stop_codon eq 'TGA'){
								$dis=$new_end-$best_end;
								$best_end=$new_end;
								last;
							}
						}
					}
					$seq_gene=substr($seq{$query_id},$best_start-1,$best_end-$best_start+1);
					if($intron_ok==1){
						foreach$key_introns(keys@arr_introns){
							if($arr_introns[$key_introns] ne undef){
								$seq_modify=$arr_introns[$key_introns];
								$seq_gene=~s/$seq_modify//;
								print out1 ">$query_id\|intron\|$arr_introns_start[$key_introns]..$arr_introns_end[$key_introns]\n$seq_modify\n";
							}
						}
					}
					print out ">$query_id\|$hit_desc\|$best_start..$best_end\|+\n$seq_gene\n";
				}
			}
	
		}
	}
	
}
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
