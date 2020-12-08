#!/usr/bin/perl
# also removes those barcodes
use warnings;
use strict;

# perl GBS_fastq_Demultiplexer_vX.pl barcodes.txt R1.fastq R2.fastq thenameIwantStuckonEveryfile

unless (@ARGV == 4) {
	print "usage perl GBS_fastq_Demultiplexer_vX.pl barcodes.txt R1.fastq R2.fastq thenameIwantStuckonEveryfile\n";
	die;
}
my $bar = $ARGV[0];
my $fastq_1 = $ARGV[1];
my $fastq_2 = $ARGV[2];
my $out = $ARGV[3];
#make true to print R1 and R2 to seperate files
my $print_mates = 1;
#my $print_mates;
# get the barcodes
my %bar;
my %c_fail;
my %ph;

open IN, $bar;
while (<IN>){
	chomp;
	my @a = split (/\t/,$_);
	$bar{$a[1]}=$a[0]; 
	#this is to remove old ones	
	#print "\t$a[1]\t$a[0]";
	if ($print_mates){
		my $filenameR1 = "$out"."_$a[0]"."_R1.fastq";
		my $filenameR2 = "$out"."_$a[0]"."_R2.fastq";
		if (-e $filenameR1) {
			system ("rm $filenameR1");
		}
		if (-e $filenameR2) {
			system ("rm $filenameR2");
		}
	}else{
		if (-e "$out"."_$a[0].fastq") {
			system ("rm $out"."_$a[0].fastq");
		}
	}
}
system  ("rm $out"."_nobar*.fastq");
close IN;

open FAST1, $fastq_1;
open FAST2, $fastq_2;

my $c;
my $line_tracker;
my %read_info;
while (my $seq1 = <FAST1>){
	my $seq2 = (<FAST2>);
	chomp $seq1;
	chomp $seq2;
	$c++;
	$line_tracker++;
	if ($line_tracker == 1){
		$read_info{"R1"}{"header"} = $seq1;
		$read_info{"R2"}{"header"} = $seq2;
	} 
	if ($line_tracker == 2){
		$read_info{"R1"}{"seq"} = $seq1;
		$read_info{"R2"}{"seq"} = $seq2;
	} 
	if ($line_tracker == 4){
		$line_tracker = '';
		$read_info{"R1"}{"qual"} = $seq1;
		$read_info{"R2"}{"qual"} = $seq2;
	
		my $read = $read_info{"R1"}{"seq"};
		my $qual = $read_info{"R1"}{"qual"};	
            	my $read2 = $read_info{"R2"}{"seq"} ;
		my $qual2 = $read_info{"R2"}{"qual"};
		my $bar_number;
		my $bcseq;
		foreach my $i (4..9){
			my $tmp = $read;
			my $bc = substr($tmp,0,$i);
			if ($bar{$bc}){
				# check that the RE site is there
				my $re_site = substr($tmp,($i),5);
				if ($re_site eq "TGCAG"){
					$bar_number = $bar{$bc};
					$bcseq = $bc; 	
				}
			}
		}	

		if ($bar_number){
			$read =~ s/$bcseq//;
			my $l = length($bcseq);
			$qual = substr($read_info{"R1"}{"qual"},$l);
		}else{
			$bar_number = "nobar";
		}
		# clean the end of the first read
		# what might be there:
		# CTGCAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGA
		if ($read=~/CTGCA/){
			my $adapter = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG";
			my @rd = split "CTGCA",$read;
			my $end = $rd[1];
			if ($end){
                  		if(length($end)>3){
					if (length($end)>length($adapter)){
						#subset the end and see if it matches
						my $tmp = substr($end,0,length($adapter));
						if ($tmp eq $adapter){
							# only use the start of the read
							$read = $rd[0];
							#trim al	so the qual
							my $tmp_qual = substr ($qual,0,length($read));
							$qual = $tmp_qual;
							$c_fail{"R1"}++;
						}
					}else{
						#subset the adapter and see if matches
						my $tmp  = substr($adapter,0,length($end));
		       	        	      	if ($tmp eq $end){   
		       	                  		# only use the start of the read
		      	 	            		$read = $rd[0];
		       		                  	my $tmp_qual = substr ($qual,0,length($read));
		       	        	            	$qual = $tmp_qual;
							$c_fail{"R1"}++;
		       	                	}
					}
				}
			}
		}
		#highly redundant code, for the second read

		if (($read2=~/CTGCA/)and($bcseq)){
		#	print "befor:\t$read2\n";
			my $revBC = reverse $bcseq;
			$revBC =~ tr/A|T|G|C/T|A|C|G/;
	#		print "after\t$revBC\n";					
                        my $adapter = $revBC."AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
                        my @rd = split "CTGCA",$read2;
                        my $end = $rd[1];
		#	print "ad:\t\t$adapter\nend\t\t$end\n";
                        if ($end){
                                if(length($end)>3){
                                        if (length($end)>length($adapter)){
                                                #subset the end and see if it matches
                                                my $tmp = substr($end,0,length($adapter));
                                                if ($tmp eq $adapter){
                                                        # only use the start of the read
                                                	$read2 = $rd[0];
                                                	#trim also the qual
                                                      	my $tmp_qual = substr ($qual2,0,length($read2));
                                                      	$qual2 = $tmp_qual;
							$c_fail{"R2"}++;
                                                }
                                        }else{
                                                #subset the adapter and see if matches
                                                my $tmp  = substr($adapter,0,length($end));
                                                if ($tmp eq $end){
                                                      	# only use the start of the read
                                                      	$read2 = $rd[0];
                                                     	 my $tmp_qual = substr ($qual2,0,length($read2));
                                                	$qual2 = $tmp_qual;
							$c_fail{"R2"}++;
                                                }
                                        } 
                                }
                        }
                }

		my $outfile =  "$out"."_$bar_number";
		$read =~ s/\./N/g;
		$read2 =~ s/\./N/g;
		if((length($read)>29) && (length($read2))>29){
			my $tmp = $read_info{"R1"}{"header"}."\n$read\n+\n$qual\n";
			if ($print_mates){
				push(@{$ph{$outfile."_R1.fastq"}},$tmp);
			}else{
				push(@{$ph{$outfile.".fastq"}},$tmp);
			}
			$tmp = $read_info{"R2"}{"header"}."\n$read2\n+\n$qual2\n";
                        if ($print_mates){
                                push(@{$ph{$outfile."_R2.fastq"}},$tmp);
                        }else{
                                push(@{$ph{$outfile.".fastq"}},$tmp);
                        } 
		}# if it is smaller then you dont use it 
		%read_info = ();
	}
	#this waits for 1million lines to print out to save the HD some agony
	if ($c==1000000) {
		foreach my $file (keys %ph){
			open OUT, ">>$file";			
			foreach my $read (@{$ph{$file}}){			
				print OUT "$read";
			}
		}
		%ph = ();
		$c = 0;
	}
}

foreach my $file (keys %ph){
    open OUT, ">>$file";             
    foreach my $read (@{$ph{$file}}){
        print OUT "$read";
    }
}


