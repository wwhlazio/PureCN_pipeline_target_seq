
#use cutadapt to cut adapters
#not finished need to check a reverse directorion
#filter reads with both reads matching and one reads quality ratio > 0.5

#perl reads_clear_hash_v5.pl "/sc/orga/projects/zhuj05a/Wenhui/HBV/data/HPV_raw_data/C1/CSZPE0131109934-42/131205_I297_FCC34KWACXX_L1_CSZPE0131109934-42_1.fq" "/sc/orga/projects/zhuj05a/Wenhui/HBV/data/HPV_raw_data/C1/CSZPE0131109934-42/131205_I297_FCC34KWACXX_L1_CSZPE0131109934-42_2.fq" "/sc/orga/projects/zhuj05a/Wenhui/HBV/data/HPV_raw_data/C1/CSZPE0131109934-42/1.adapter.list" "/sc/orga/projects/zhuj05a/Wenhui/HBV/data/HPV_raw_data/C1/CSZPE0131109934-42/2.adapter.list" "70" "file_not_passed_v5_1.fq" "file_not_passed_v5_2.fq" "file_passed_v5_1.fq" "file_passed_v5_2.fq"

use POSIX;

$fastq1=$ARGV[0];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/data/HPV_raw_data/C1/CSZPE0131109934-42/131205_I297_FCC34KWACXX_L1_CSZPE0131109934-42_1.fq";
$fastq2=$ARGV[1];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/data/HPV_raw_data/C1/CSZPE0131109934-42/131205_I297_FCC34KWACXX_L1_CSZPE0131109934-42_2.fq";

#$pl1=$ARGV[2];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/data/HPV_raw_data/C1/CSZPE0131109934-42/1.adapter.list";
#$pl2=$ARGV[3];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/data/HPV_raw_data/C1/CSZPE0131109934-42/2.adapter.list";

$thrq=$ARGV[2];#70; #64+5

$file_not_pass1=$ARGV[3];#"file_not_passed_v5_1.fq";
$file_not_pass2=$ARGV[4];#"file_not_passed_v5_2.fq";
$file_pass1=$ARGV[5];#"file_passed_v5_1.fq";
$file_pass2=$ARGV[6];#"file_passed_v5_2.fq";


#print $pl2,"\n";
#open(logfile1,$pl1) or die"I can't read the file!";
#open(logfile2,$pl2) or die"I can't read the file!";

#$npr=0;
#$u=0;
#for $line(<logfile1>){
#        @a=split("\t",$line);
#        if($u>0){
#                @b=split("\/",$a[0]);
#                $pr[$npr]=$b[0];
#                $npr=$npr+1;
#        }
#        $u=$u+1;
#}
#print "number of pelluted reads in file1:$npr\n";

#$u=0;
#for $line(<logfile2>){
#        @a=split("\t",$line);
#        if($u>0){
#                @b=split("\/",$a[0]);
#                $ind=0;
#                for $i(0..$npr-1){
#                        if($b[0] eq $pr[$i]){
#                                $ind=1;
#                                last;
#                        }
#                }
#                if($ind==0){
#                        $pr[$npr]=$b[0];
#                        $npr=$npr+1;
#                }
#        }
#        $u=$u+1;
#}
#print "number of pelluted reads after adding file2:$npr\n";
#close(logfile1);
#close(logfile2);



#my $uu = "0" x $npr;

#my @rid=split(//,$uu);


my @parax;

for $i(0..37){
	push(@parax,$i*128);
}

#for $i(0..$npr-1){
#	$pr[$i] =~ s/\\/\\\\/g;
#	my @psx=unpack("C*",$pr[$i]);
#	my $index = 0;
#	my @idd;
#	@idd = map{$psx[ $_] * $parax[$_]} 0..$#parax;
#	for $j(0..$#parax){
#		$index=$index+$idd[$j];
#		$u=$index % $npr;
#		$index=$u;
#	}
#	if($rid[$index] eq "0"){
#		$rid[$index]=$pr[$i];
#	}
#	else{
#		$rid[$index]=$rid[$index]."\n".$pr[$i];
#	}
	
#}

#$xx=0;
#for $i(0..$npr-1){
#	if($rid[$i] ne "0"){
#		@axx=split(/\n/,$rid[$i]);
#		if($#axx==0){
#			$xx=$xx+$#axx+1;
#		}else{
#			$xx=$xx+$#axx+1;
#		}		
#		
#	}	
#}
#print "xx:$xx\n";


#$file_not_pass1="file_not_passed_v5_1.fq";
#$file_not_pass2="file_not_passed_v5_2.fq";
open(file1,">",$file_not_pass1) or die"I can't write to the file!";
open(file2,">",$file_not_pass2) or die"I can't write to the file!";


open(logfile1,$fastq1) or die"I can't read the file!";
open(logfile2,$fastq2) or die"I can't read the file!";
print "start loading fastq1\n;";
my $nr1=0;
my $u=0;
my @r1;
my $xxb=0;
for $line1(<logfile1>){
	@a=split(/\n/,$line1);
	
	if($u==0){
		@b=split(/\//,$a[0]);
		@c=split(/\@/,$b[0]);
		$c[1] =~ s/\\/\\\\/g;
	
#		@psx=unpack("C*",$c[1]);
#		@idd = map{$psx[ $_] * $parax[$_]} 0..$#parax;
#        	$index=0;
#		$uu=0;
#		for $i(0..$#parax){
#                	$index=$index+$idd[$i];
#                	$uu=$index % $npr;
#                	$index=$uu;
#        	}
#		if($rid[$index] eq "0"){
			$r1[4*$nr1]=$a[0];
#			$xxb=1	
#		}
#		else{
#			@xx=split(/\n/,$rid[$index]);
#			$ll=0;
#			for $i(0..$#xx){
#				if($c[1] eq $xx[$i]){
#					$ll=1;
#					$xxb=0;
#					print file1 $a[0],"\n";
##					print "OOK:$i:$index:",$a[0],"\n";
#					last;
#				}
#			}
#			if($ll==0){
#				$r1[4*$nr1]=$a[0];
#				$xxb=1;
#			}
#		}
	}
	elsif($u==1){
#		if($xxb==1){
			$r1[4*$nr1+1]=$a[0];
#		}
#		elsif($xxb==0){
#			print file1 $a[0],"\n";
#		}
	}	
	elsif($u==2){
#		if($xxb==1){
			$r1[4*$nr1+2]=$a[0];
#		}
#		elsif($xxb==0){
#			print file1 $a[0],"\n";
#		}
	}
	elsif($u==3){
#		if($xxb==1){
			$r1[4*$nr1+3]=$a[0];
#		}
#		elsif($xxb==0){
#			print file1 $a[0],"\n";
#		}
	}
	$u=$u+1;
	if($u==4){
#		if($xxb == 1){
			$nr1=$nr1+1;
#		}
		$u=0;
#		$xxb=0;
	}	
}
print "finish loading fastq1.number of reads is $nr1\n";

print "start loading fastq2\n";
my $nr2=0;
my $u=0;
my @r2;
$xxb=0;
for $line1(<logfile2>){
        @a=split(/\n/,$line1);
        if($u==0){
		@b=split(/\//,$a[0]);
		@c=split(/\@/,$b[0]);
		$c[1] =~ s/\\/\\\\/g;
		
#		@psx=unpack("C*",$c[1]);
#		@idd = map{$psx[ $_] * $parax[$_]} 0..$#parax;
#		$index=0;
#                $uu=0;
#			
#		for $i(0..$#parax){
#                        $index=$index+$idd[$i];
#                        $uu=$index % $npr;
#                        $index=$uu;
#                }
#                if($rid[$index] eq "0"){
                        $r2[4*$nr2]=$a[0];
#                        $xxb=1
#                }
#                else{
#                        @xx=split(/\n/,$rid[$index]);
#                        $ll=0;
#                        for $i(0..$#xx){
#                                if($c[1] eq $xx[$i]){
#                                        $ll=1;
#                                        $xxb=0;
#					print file2 $a[0],"\n";
#                                        last;
#                                }
#                        }
#                        if($ll==0){
#                                $r2[4*$nr2]=$a[0];
#                                $xxb=1;
#                        }
#                }
		
	
        }
        elsif($u==1){
#		if($xxb==1){
                	$r2[4*$nr2+1]=$a[0];
#		}
#		elsif($xxb==0){
#			print file2 $a[0],"\n";
#		}
        }
        elsif($u==2){
#		if($xxb==1){
                	$r2[4*$nr2+2]=$a[0];
#		}
#		elsif($xxb==0){
#			print file2 $a[0],"\n";
#		}
        }
        elsif($u==3){
#		if($xxb==1){
                	$r2[4*$nr2+3]=$a[0];
#		}
#		elsif($xxb==0){
#			print file2 $a[0],"\n";
#		}
        }
        $u=$u+1;
        if($u==4){
#		if($xxb==1){	
                	$nr2=$nr2+1;
#		}
                $u=0;
#		$xxb=0;
        }
}
print "finish loading fastq2.number of reads is $nr2\n";

if($nr1!=$nr2){
	print "fatal error! input not paired\n";
}

print "start cleaing:\n";
#$file_not_pass1="file_not_passed_v5_1.fq";
#$file_not_pass2="file_not_passed_v5_2.fq";
#open(file1,">",$file_not_pass1) or die"I can't write to the file!";
#open(file2,">",$file_not_pass2) or die"I can't write to the file!";
#$file_pass1="file_passed_v5_1.fq";
#$file_pass2="file_passed_v5_2.fq";
open(file11,">",$file_pass1) or die"I can't write to the file!";
open(file21,">",$file_pass2) or die"I can't write to the file!";



$rsnr=0; #reserved number foreads

$rsnr_pr=ceil($nr1*2/3);
print "size:$nr1","\t","$rsnr_pr\n";

$vv = "0" x $rsnr_pr;
@unique_idx=split(//,$vv);


my @para1;

for $i(0..9){
	push(@para1,5**$i);	
}


for $i(0..$nr1-1){
	$r1[4*$i+3] =~ s/\\/\\\\/g;
	$r2[4*$i+3] =~ s/\\/\\\\/g;
	my @ps1=unpack("C*",$r1[4*$i+3]);
        my @ps2=unpack("C*",$r2[4*$i+3]);
	$r1[4*$i+3] =~ s/\\\\/\\/g;
        $r2[4*$i+3] =~ s/\\\\/\\/g;

        my @ps11=grep {$_<$thrq} @ps1;
        my @ps21=grep {$_<$thrq} @ps2;
        my $nbp1=scalar(@ps1);
	my $nbp2=scalar(@ps2);
        my $nbp1_1=scalar(@ps11);
	my $nbp1_2=scalar(@ps21);
	my $u=$nbp1_1/$nbp1;
	my $v=$nbp1_2/$nbp2;
        if($u<0.5 && $v<0.5){
		my @id1=split(//,$r1[4*$i+1]);
		my @id2=split(//,$r2[4*$i+1]);
		
#		@id1=unpack("C*",$r1[4*$i+1]);	
#		@id2=unpack("C*",$r2[4*$i+1]);
		push(@id1,@id2);
		if(scalar(@id1)!=302){
			print "error:read length error!\n";
		}
	
		s/A/0/ for @id1;
		s/C/1/ for @id1;
		s/G/2/ for @id1;
		s/T/3/ for @id1;
		s/N/4/ for @id1;
		my @id3;
		my $index=0;
		for $j(0..19){
			@id3 = map{$id1[10*$j + $_] * $para1[$_]} 0..$#para1;
			for $k(0..$#para1){
				$index=$index+$id3[$j];
			}
			$index = $index % $rsnr_pr;
		}
#		print $index,"\n";
#		print @id1,"\n";
#		my @id3 = map {$id1[$_] * $para1[$_]} 0..$#id1;
#		print @id3,"\n";
#		my $id = unpack  "%123d*", pack("d*",@id3);
#		my $id=0;
#		for $j(0..$#id1){
#			$id=$id+$id3[$j];
#		}
#		print $id,"\n";
		
#		$index=$id % $rsnr_pr;
#		print "index:$index\n";
	
		if($unique_idx[$index] eq "0"){
			$unique_idx[$index] = $r1[4*$i+1].$r2[4*$i+1];
#			print "**************************\n $unique_idx[$index]\n****************************\n";
			print file11 $r1[4*$i],"\n";
                        print file11 $r1[4*$i+1],"\n";
                        print file11 $r1[4*$i+2],"\n";
                        print file11 $r1[4*$i+3],"\n";

                        print file21 $r2[4*$i],"\n";
                        print file21 $r2[4*$i+1],"\n";
                        print file21 $r2[4*$i+2],"\n";
                        print file21 $r2[4*$i+3],"\n";
		}
		elsif($unique_idx[$index] ne "0"){
				
				@a=split(/\n/,$unique_idx[$index]);
				$l = 0;
				$v = $r1[4*$i+1].$r2[4*$i+1];
				
#				print $v,"\n";	
				for $j(0..$#a){
#					print "##################\n";
#					print "$v\n";
#					print "$a[$j]\n";
#					print "##################\n";
					if($v eq $a[$j]){
						$l=1;
						last;
					}
				}
#				print $l,"\n#######################\n";
				if($l == 0){
					$unique_idx[$index]=$unique_idx[$index]."\n".$v;
					print file11 $r1[4*$i],"\n";
                        		print file11 $r1[4*$i+1],"\n";
                        		print file11 $r1[4*$i+2],"\n";
                        		print file11 $r1[4*$i+3],"\n";

                        		print file21 $r2[4*$i],"\n";
                        		print file21 $r2[4*$i+1],"\n";
                        		print file21 $r2[4*$i+2],"\n";
                        		print file21 $r2[4*$i+3],"\n";
					
				
#				@set_pos=$index+1..$rsnr_pr-1;
#				@u=0..$index-1;
#				push(@set_pos, @u);
#				for $j(0..$#set_pos){
#					if($unique_reads1[$set_pos[$j]]==0){
#						$unique_reads1[$set_pos[$j]]=$r1[4*$i]."\n".$r1[4*$i+1]."\n".$r1[4*$i+2]."\n".$r1[4*$i+3];
#                				$unique_reads2[$set_pos[$j]]=$r2[4*$i]."\n".$r2[4*$i+1]."\n".$r2[4*$i+2]."\n".$r2[4*$i+3];
#						last;
#					}
#				}
				}
				if($l==1){
				 	print file1 $r1[4*$i],"\n";
                        	 	print file1 $r1[4*$i+1],"\n";
                        	 	print file1 $r1[4*$i+2],"\n";
                         	 	print file1 $r1[4*$i+3],"\n";
#					print $r1[4*$i],"\n";

                         	 	print file2 $r2[4*$i],"\n";
                        	 	print file2 $r2[4*$i+1],"\n";
                        	 	print file2 $r2[4*$i+2],"\n";
                        	 	print file2 $r2[4*$i+3],"\n";
#					print $r2[4*$i],"\n";
				}
		}	
#			else{
#				$unique_idx[$index]=$id;
#				@set_pos=$index+1..$rsnr_pr-1;
#                                @u=0..$index-1;
#                                push(@set_pos, @u);
#				for $j(0..$#set_pos-1){
#                                        if($unique_reads1[$set_pos[$j]]==0){
#                                                $unique_reads1[$set_pos[$j]]=$r1[4*$i]."\n".$r1[4*$i+1]."\n".$r1[4*$i+2]."\n".$r1[4*$i+3];
#                                                $unique_reads2[$set_pos[$j]]=$r2[4*$i]."\n".$r2[4*$i+1]."\n".$r2[4*$i+2]."\n".$r2[4*$i+3];
#                                                last;
#                                        }
#                                }
#			}
	}	
	else{
        	print file1 $r1[4*$i],"\n";
                print file1 $r1[4*$i+1],"\n";
                print file1 $r1[4*$i+2],"\n";
                print file1 $r1[4*$i+3],"\n";
#		print $r1[4*$i],"\n";

                print file2 $r2[4*$i],"\n";
                print file2 $r2[4*$i+1],"\n";
                print file2 $r2[4*$i+2],"\n";
                print file2 $r2[4*$i+3],"\n";
#		print $r2[4*$i],"\n";
	}

	if($i%100000==0){
		print "$i reads has been filtered\n";
	}
	
}
close(file1);
close(file2);
close(file11);
close(file21);

print "clean done!\n";

