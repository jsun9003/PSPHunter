#!/usr/bin/perl
my ($id,$email)=@ARGV;
#my $identi=0.3;

my $home_dir="/home/PSPHunter/";
my $dir="/var/www/psphunter.stemcellding.org/job/$id/";
my $progam_dir=$home_dir."scripts/";
my $sendmail_file="/var/www/psphunter.stemcellding.org/sendmail.pl";
my $sendmail_filetome="/var/www/psphunter.stemcellding.org/sendmailtome.pl";
my $emailsend_file_check ="$dir"."emailcheck.txt";

my $jsun_email="o.sj\@qq.com";
my $own_flag=0;
if (-e $emailsend_file_check)
{
    print "you will get an email latter;";
    $own_flag++;
}else{
    system "/home/PSPHunter/software/miniconda3/bin/perl $sendmail_filetome $jsun_email $id 2";
}

####Split fasta into different truncations
my $seq_fa=$dir."temp.fasta";

my $fasta_dir=$dir."fasta2/";
mkdir $fasta_dir;

my $seq="";
open(IN,$seq_fa) or die "$!";
while (<IN>)
{
    chomp;
    print $_,"\tRaw seq Each Line\n";
    if($_ =~/>/)
    {
       print $_,"\tHeader sequence Infor\n";
       next;
    }else{
       print $seq,"\tSeq\n";
       $_=~s/\n|\r|\s//g;
       $seq=$seq."".$_;
    }
}

close IN;
#my $seq=$lines[1];
$seq=~s/-|\*//g;

my $len=length $seq;
$seq=uc $seq;
print $seq,"\tfinal seq\n";
#print $len,"\t$seq\n";
my $frag="20";
my @list=();
for(my $i=0;$i<$len-$frag+1;$i++)
{
    my $new=substr($seq,0,$i).substr($seq,$i+$frag);
    my $temp="query"."_$i";
    push(@list,$temp);
    open(OUT,">$fasta_dir"."$temp.fasta") or die "$!";
    print OUT ">query"."_$i\n$new\n";
    close OUT;
}

my $ML_dir=$dir."ML2/";
mkdir $ML_dir;

open(OUT,">$ML_dir"."Intest.txt") or die "$!";
foreach my $human (@list)
{
    print OUT "$human\t0\n";
}
close OUT;

###word2vec hash
my %wv=();
open(IN, "$home_dir"."wordvec/uniprot_sprot70_size60.txt") or die "$!";
while (<IN>)
{
    chomp;
    next if /^Uniprot/;
    my @item=split/\t/;
    my @ref=@item[1..$#item];
    $wv{$item[0]}=\@ref;
}
close IN;

my $o_dir=$ML_dir."word2vec/";
mkdir $o_dir;

my $size=60;
open(OUT, ">$o_dir"."Intestinput.txt") or die "$!";
open(FA, "$ML_dir"."Intest.txt") or die "$!";
while (<FA>)
{
    chomp;
    my ($uni,$tag)=split/\s+/,$_;
    my @fea=();
    my @vec=();
    my $bacfasta="$fasta_dir"."$uni.fasta";

    if (!-s $bacfasta) {
        print $bacfasta,"\t","no bac protein\n";
        for(my $j=1;$j<=$size;$j++) #39
        {
            push(@vec,0);
        }
    }else{
        open(F,$bacfasta) or die "$!";
        chomp(my @bacteria_seq=<F>);
        close F;
        my @aa=split//,$bacteria_seq[1];
        for(my $i=0;$i<(scalar @aa)-2;$i++)
        {
            my $tri=$aa[$i].$aa[$i+1].$aa[$i+2];
            if (exists $wv{$tri})
            {
                my @tvec=@{$wv{$tri}};
                for(my $j=0;$j<@tvec;$j++)
                {
                    $vec[$j]+=$tvec[$j];
                }
            }
            
        }
    }
    push(@fea,@vec);
    print OUT "$tag\t",join "\t",@fea,"\n";
}
close OUT;
close FA;

###Model Test
my $python=$progam_dir."Intest-apply-test.py";
for(my $i=1;$i<=100;$i++)
{
    my $f1_dir=$home_dir."/train/$i/word2vec70_60/";
    my $f2_dir=$o_dir;
    system "/home/PSPHunter/software/miniconda3/bin/python $python $f1_dir $f2_dir $i";
}

#Avg Prob

open(IN,"$ML_dir"."Intest.txt") or die "$!";
my $n=0;
my %tag=();
while (<IN>)
{
    chomp;
    $n++;
    my ($uni,$tag)=split;
    $tag{$n}=$uni;
}

close IN;

open(OUT,">$o_dir"."InProbAVG.txt") or die "$!";
my %hash=();
for(my $i=1;$i<=100;$i++)
{
    open(IN,"$o_dir"."InProbphase$i.txt") or die "$!";
    my $k=0;
    while (<IN>)
    {
        chomp;
        $k++;
        if (exists $hash{$k})
        {
            $hash{$k}+=$_;
        }else{
            $hash{$k}=$_;
        }
    }
    close IN;
}

foreach my $key (sort {$hash{$b}<=>$hash{$a}} keys %hash)
{
    my $prob=sprintf "%0.3f",$hash{$key}/100;
    print OUT $key,"\t$prob\n";
}
close OUT;

####detect driving region

open(IN,"$o_dir"."InProbAVG.txt") or die "$!";
my @prob=();

%hash=();
while (<IN>)
{
   chomp;
    my ($index,$score)=split;
    push(@prob,$score);
    $hash{$index}=$score
}
close IN;

my $length=scalar @prob;

#my $bin=int($length/20);
#将蛋白质划分成bin，3，5，10，相较于平均值的偏离程度，(x-avg)/avg sort
#考虑如何合并
my $avg=&Avg(@prob);

#print $length,"\t",$avg,"\n";
my %rank=();
foreach my $key (sort {$a<=>$b} keys %hash)
{
   if ($key <= $length)
   {
      #my $prob=&lag($key);
      my $prob=$hash{$key};
      # my $bias=($avg-$prob)/$avg;
      my $bias=$avg-$prob;
      $rank{$key}=$bias;
   }
}

open(OUT,">$o_dir"."temp.bed") or die "$!";
$n=0;
foreach my $key (sort {$rank{$b}<=>$rank{$a}} keys %rank)
{
   #last if $rank{$key}<0;
    $n++;
    if($length >2000)
    {
      last if $n>=int ($length*0.01);# >20
      }elsif($length <=2000 && $length >1000)
    {
      last if $n>=int ($length*0.02);## 20-40
    }elsif($length > 500 && $length<=1000)
    {
      last if $n>=int ($length*0.04); ##20-40
    }else{
      last if $n>=int ($length*0.05);###25
    }
    
   print OUT "query\t",$key,"\t",$key+20,"\t",$rank{$key},"\n";
}
close OUT;

my $par1="sort -k1,1 -k2,2n ".$o_dir."temp.bed >".$o_dir."temp.sort.bed";
my $par2="/home/PSPHunter/software/miniconda3/bin/bedtools merge -i ".$o_dir."temp.sort.bed >".$o_dir."temp.drRegion";
system("$par1");
system("$par2");

%rank=();
open(IN,"$o_dir"."temp.drRegion") or die "$!";
while (<IN>)
{
   chomp;
   my ($uni,$s,$e)=split/\t/;
   my $sum=0;
   my $m=0;
   for(my $i=$s;$i<=$e;$i++)
   {
      if (exists $hash{$i})
      {
         $sum+=$hash{$i};
         $m++;
      }else{
         print $_,"\t$i\n";
      }
   }
   my $avg=$sum/$m;
   $rank{$_}=$avg;
}
close IN;

open(OUT,">$o_dir"."drRegion.rank") or die "$!";
#print OUT "#Total Length\t$length\n";
my $limit=3;
my %region=();
$n=0;
foreach my $key (sort {$rank{$a}<=>$rank{$b}} keys %rank)
{
    $n++;
    my ($query,$s,$e,$score)=split/\t/,$key;
    print OUT $key,"\n";
    for(my $i=$s;$i<=$e;$i++)
    {
        $region{$i}++;
    }
    last if $n>=$limit;
}
close OUT;

open(OUT,">$dir"."final.txt") or die "$!";
print OUT "#Residue in Purple denoted driving residues\n";
print OUT "Pos\tAA\tProb\tDRegion\n";

my @aa=split//,$seq;
for(my $i=0;$i<@aa;$i++)
{
    my $index=$i+1;
    next if $aa[$i] eq "-";
    if (exists $region{$index})
    {
        if (exists $hash{$i-9})
        {
            print OUT $index,"\t$aa[$i]\t$hash{$i-9}\t1\n";
        }else{
            print OUT $index,"\t$aa[$i]\t-\t1\n";
        }
    }else{
        if (exists $hash{$i-9})
        {
           if($aa[$i]){ print OUT $index,"\t$aa[$i]\t$hash{$i-9}\t0\n";}
        }else{
            if($aa[$i]){print OUT $index,"\t$aa[$i]\t-\t0\n";}
        } 
    }
}

close OUT;

sub Avg
{
   my @prob=@_;
   my $sum=0;
   my $n=0;
   foreach my $prob(@prob)
   {
      $sum+=$prob;
      $n++;
   }
   my $avg=$sum/$n;
   return $avg;
}


my $sendmail=$progam_dir."sendmail_myemail.py";

open(DO,">$dir"."done_predictDr.txt")||die"$!";
foreach my $key(keys %done){
  print DO $key,",",$done{$key},"\n";
}
close DO;

print $email,"\tEmail To\n";
if (-e $emailsend_file_check)
{
    print "you will get an email latter;"
}else{
    system "/home/PSPHunter/software/miniconda3/bin/perl $sendmail_file $email $id 2";
}

system "rm -rf $fasta_dir";
system "rm -rf $ML_dir";
#system "zip -j $dir/PSPHunter_results.zip $dir/*.out $dir/*.pdb";
