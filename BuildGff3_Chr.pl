# Plant genome assembly software
#author Bi Changwei NJFU
#version 1.0 
# Please direct questions and support issues to <511281975@qq.com>
#!/usr/bin/perl
use strict;
use Bio::SeqIO;
use  Bio::Seq;
my @scafname;
my @line;
my @scafflag;
my @scaflength;
my $scaflength;
my $scafflag;
my $scafname;
my $j=1;
system 'mkdir Wlchr';
open OUT,"> result1.gff3";
open OUT1,"> result1.fa";
for(my $i=1;$i<=19;$i++)
{
	system 'mkdir Wlchr/wl'.$i;
	my $path="chrname/chr$i.txt";
	open(IN,$path);
	while(my $line=<IN>)
	{
		chomp;
		@line=split (/\s+/,$line);
		$scafname=$line[0];
		push(@scafname,$scafname);
		$scaflength=$line[1];
		push(@scaflength,$scaflength);
		$scafflag=$line[2];
		push(@scafflag,$scafflag);
	}	

	#usecommand.pl

	foreach my $item(@scafname)
	{
		my $cmdStr="grep -E '" . $item . "\\>' " . "../final/Willow.gene.gff3 >Wlchr/wl".$i."/WL" . $j . "\n";
		$j=$j+1;
		system($cmdStr);
	}	       
		$j=1;

	#############
	#revGff3.pl
	#############
	#print "@scaflength\n";

	for(my $m=0;$m<@scaflength;$m++)
	{
		if($scafflag[$m] eq 1)
                {
             		my $path="Wlchr/wl".$i."/WL".($m+1);
                        my $path1="Wlchr/wl".$i."/WL".($m+1)."rev";
                        #print $path,"\n",$path1,"\n";
                        open(FILE,"<",$path)||die"cannot open the file: $!\n";
                        open(FILE1,">",$path1);

                        while (<FILE>)
                        {
                        	$_=~s/\s*$ //;
	                        my ($chr, $t1, $t2, $start, $end, $t3, $orientation, $t4, $geneid) = split ("\t", $_);

        	                if($orientation eq "+")
                	        {       $orientation = "-";}
                        	elsif($orientation eq "-")
	                        {       $orientation = "+";}

        	                my $newstart=$scaflength[$m]-$end+1;
          	                my $newend=$scaflength[$m]-$start+1;
                                my $record=$chr."\t".$t1."\t".$t2."\t".$newstart."\t".$newend."\t".$t3."\t".$orientation."\t".$t4."\t".$geneid;

                                print FILE1 $record;

                        }

                        close(FILE);
                        close(FILE1);
                 }
          }
        
	############
	#BuildGff3.pl
        # print "@scafflag\n";

	my $off=0;
	for(my $n=0;$n<@scaflength;$n++)
        {
                if($n>0)
                {$off=$off+100+$scaflength[($n-1)];}
                                #print $off."\n";
                if($scafflag[$n]==0)
                {
                	$path="Wlchr/wl".$i."/WL".($n+1);
                }
                else
                {
                        $path="Wlchr/wl".$i."/WL".($n+1)."rev";
                }
                open(FILE,"<",$path)||die"cannot open the file: $!\n";
                while (<FILE>){
                        $_=~s/\s*$//;
                        my ($chr, $t1, $t2, $start, $end, $t3, $orientation, $t4, $geneid) = split ("\t", $_);
			if($i<=9)
			{
                        	my $record="chr0$i"."\t".$t1."\t".$t2."\t".($start+$off)."\t".($end+$off)."\t".$t3."\t".$orientation."\t".$t4."\t".$geneid."\n";
                        	print OUT $record;
			}else{
				my $record="chr$i"."\t".$t1."\t".$t2."\t".($start+$off)."\t".($end+$off)."\t".$t3."\t".$orientation."\t".$t4."\t".$geneid."\n";
                                print OUT $record;
			}

                }

                close FILE;
       }

	my $ChrTotalLength=0;
	my $ChrSeq="R";
        my $SeqN="N" x 100;
	my $chrName;
	if($i<=9)
	{
        	 $chrName=">chr0$i";
	}else{
		 $chrName=">chr$i";
	}
        my $p=0;
        my $first=1;
	foreach my $line(@scafname)
	{
        	my $catchseq_seqio_obj = Bio::SeqIO->new(-file=>"../final/Willow.fa", -format=>'fasta');
                while(my $seq_obj = $catchseq_seqio_obj->next_seq)
                {
				my $display_name = $seq_obj->display_name;
				my $seq = $seq_obj->seq;                                   
                                my $seq_length = $seq_obj->length;
                                if ($display_name eq $line)
                                {
                                        $ChrTotalLength =$ChrTotalLength + $seq_length+100;
				        if($scafflag[$p] eq 1)
                                        {
                                                $seq = reverse $seq;
                                                $seq =~tr/ATGCatgc/TACGtacg/;
                                        }
                                        $p=$p+1;
                                        $ChrSeq= $ChrSeq.$seq.$SeqN;
                                        last;
                                }
                }
        }
#print $ChrTotalLength . "\n";
        print OUT1 $chrName . "\n";
	$ChrSeq=substr($ChrSeq,1,$ChrTotalLength-100);
	print OUT1 $ChrSeq."\n";
	undef @scafname;
        undef @scafflag ;
        undef @scaflength;
}
