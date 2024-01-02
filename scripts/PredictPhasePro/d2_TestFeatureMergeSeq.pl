use warnings;
use strict;

my $abs_dir="/data1/juns/phase2/";
my $data_dir="/data1/juns/PathHost/data_raw/";
my $fasta_dir="/data1/juns/PathHost/data_raw/fasta_pathhost/";
my $hhm_dir="/data1/juns/PathHost/data_raw/hhm_pathhost/";
my $pssm_dir="/data1/juns/PathHost/data_raw/pssm_pathhost/";
my $choice="$ARGV[0]";#"scaffold";
my %class=(
            "C"=>"1","W"=>"1","H"=>"1","M"=>"1","Y"=>"1","F"=>"1","I"=>"1","V"=>"1","A"=>"1","L"=>"1",
            "G"=>"2","P"=>"2",
            "S"=>"3","Q"=>"3","T"=>"3","N"=>"3",
            "E"=>"4","R"=>"4","K"=>"4","D"=>"4"
           );
my $f_dir=$abs_dir."ML/PHApply/scaffold/$choice/test/";
my $o_dir=$f_dir."mergeSeq/";
mkdir $o_dir;

my %hash=();
open(FA, "$data_dir"."Pathhost_ASA_tag.txt") or die "$!";
while (<FA>)
{
    chomp;
    my ($uni,$pert1,$pert2,$pert3,$pert4)=split/\s+/;
    my @ref=($pert1,$pert2,$pert3,$pert4);
    $hash{$uni}=\@ref;
}
close FA;

open(FA, "$abs_dir"."Mut/NonPS/"."ASA_SS.txt") or die "$!";
while (<FA>)
{
    chomp;
    my ($uni,$pert1,$pert2,$pert3,$pert4)=split/\s+/;
    my @ref=($pert1,$pert2,$pert3,$pert4);
    $hash{$uni}=\@ref;
}
close FA;

open(FA, "$abs_dir"."Mut/PhaSePred/"."ASA_SS.txt") or die "$!";
while (<FA>)
{
	chomp;
	my ($uni,$pert1,$pert2,$pert3,$pert4)=split/\s+/;
	my @ref=($pert1,$pert2,$pert3,$pert4);
	$hash{$uni}=\@ref;
}
close FA;

open(FA, "$abs_dir"."Mut/BIB/"."ASA_SS.txt") or die "$!";
while (<FA>)
{
	chomp;
	my ($uni,$pert1,$pert2,$pert3,$pert4)=split/\s+/;
	my @ref=($pert1,$pert2,$pert3,$pert4);
	$hash{$uni}=\@ref;
}
close FA;

my %c=();
my $disorder_avg=0;
my $an=0;
my $sum=0;
open(FA, "/data1/juns/phase/IDR/SPINE-D-local.txt") or die "$!";
while (<FA>)
{
    chomp;
    next if /^uniprot/;
    my ($uni,$dis,$l,$count,$pert)=split/\s+/;
    next if exists $c{$uni};
    $c{$uni}=$pert;
    $sum+=$pert;
    $an++;
}
$disorder_avg=$sum/$an;

close FA;

open(FA, "$abs_dir"."Mut/NonPS/"."SPINE-D-local.txt") or die "$!";
while (<FA>)
{
    chomp;
    next if /^uniprot/;
    my ($uni,$dis,$l,$count,$pert)=split/\s+/;
    next if exists $c{$uni};
    $c{$uni}=$pert;
}
close FA;

open(FA, "$abs_dir"."Mut/PhaSePred/"."SPINE-D-local.txt") or die "$!";
while (<FA>)
{
	chomp;
	next if /^uniprot/;
	my ($uni,$dis,$l,$count,$pert)=split/\s+/;
	next if exists $c{$uni};
	$c{$uni}=$pert;
}
close FA;

open(FA, "$abs_dir"."Mut/BIB/"."SPINE-D-local.txt") or die "$!";
while (<FA>)
{
	chomp;
	next if /^uniprot/;
	my ($uni,$dis,$l,$count,$pert)=split/\s+/;
	next if exists $c{$uni};
	$c{$uni}=$pert;
}
close FA;

my %rna=();
my $rna_avg=0;
$an=0;
$sum=0;
open(FA, "/data1/juns/phase/NAR_binding/RNA_binding.txt") or die "$!";
while (<FA>)
{
    chomp;
    next if /^uniprot/;
    my ($uni,$dis,$l,$count,$pert)=split/\s+/;
    next if exists $rna{$uni};
    $rna{$uni}=$pert;
    $sum+=$pert;
    $an++;
}
$rna_avg=$sum/$an;
close FA;

my %dna=();
my $dna_avg=0;
$an=0;
$sum=0;
open(FA, "/data1/juns/phase/NAR_binding/DNA_binding.txt") or die "$!";
while (<FA>)
{
    chomp;
    next if /^uniprot/;
    my ($uni,$dis,$l,$count,$pert)=split/\s+/;
    next if exists $dna{$uni};
    $dna{$uni}=$pert;
    $sum+=$pert;
    $an++;
}
$dna_avg=$sum/$an;
close FA;

open(FA, "$abs_dir"."Mut/NonPS/"."NAR_binding/DNA_binding.txt") or die "$!";
while (<FA>)
{
    chomp;
    next if /^uniprot/;
    my ($uni,$l,$count,$pert)=split/\s+/;
    next if exists $dna{$uni};
    $dna{$uni}=$pert;
}
close FA;

open(FA, "$abs_dir"."Mut/NonPS/"."NAR_binding/RNA_binding.txt") or die "$!";
while (<FA>)
{
    chomp;
    next if /^uniprot/;
    my ($uni,$l,$count,$pert)=split/\s+/;
    next if exists $rna{$uni};
    $rna{$uni}=$pert;
}
close FA;

open(FA, "$abs_dir"."Mut/PhaSePred/"."NAR_binding/DNA_binding.txt") or die "$!";
while (<FA>)
{
	chomp;
	next if /^uniprot/;
	my ($uni,$l,$count,$pert)=split/\s+/;
	next if exists $dna{$uni};
	$dna{$uni}=$pert;
}
close FA;

open(FA, "$abs_dir"."Mut/PhaSePred/"."NAR_binding/RNA_binding.txt") or die "$!";
while (<FA>)
{
	chomp;
	next if /^uniprot/;
	my ($uni,$l,$count,$pert)=split/\s+/;
	next if exists $rna{$uni};
	$rna{$uni}=$pert;
}
close FA;

open(FA, "$abs_dir"."Mut/BIB/"."NAR_binding/DNA_binding.txt") or die "$!";
while (<FA>)
{
	chomp;
	next if /^uniprot/;
	my ($uni,$l,$count,$pert)=split/\s+/;
	next if exists $dna{$uni};
	$dna{$uni}=$pert;
}
close FA;

open(FA, "$abs_dir"."Mut/BIB/"."NAR_binding/RNA_binding.txt") or die "$!";
while (<FA>)
{
	chomp;
	next if /^uniprot/;
	my ($uni,$l,$count,$pert)=split/\s+/;
	next if exists $rna{$uni};
	$rna{$uni}=$pert;
}
close FA;

my %cs=();
open(FB, "$data_dir"."cscore_quantile.txt") or die "$!";
while (<FB>)
{
    chomp;
    my @item=split/\t/;
    my @ref=@item[1..$#item];
    $cs{$item[0]}=\@ref;
}
close FB;

open(FB, "$abs_dir"."Mut/NonPS/cscore_quantile.txt") or die "$!";
while (<FB>)
{
    chomp;
    my @item=split/\t/;
    my @ref=@item[1..$#item];
    $cs{$item[0]}=\@ref;
}
close FB;

open(FB, "$abs_dir"."Mut/PhaSePred/cscore_quantile.txt") or die "$!";
while (<FB>)
{
	chomp;
	my @item=split/\t/;
	my @ref=@item[1..$#item];
	$cs{$item[0]}=\@ref;
}
close FB;

open(FB, "$abs_dir"."Mut/BIB/cscore_quantile.txt") or die "$!";
while (<FB>)
{
	chomp;
	my @item=split/\t/;
	my @ref=@item[1..$#item];
	$cs{$item[0]}=\@ref;
}
close FB;


#Spscore
my %ps=();
my $ps_avg=0;
$an=0;
$sum=0;
open(IN,"/data1/juns/phase/PScore/pscore.txt") or die "$!";
while (<IN>)
{
    chomp;
    next if /^uniprot/;
    my @item=split/\t/;
    $ps{$item[0]}=$item[1];
    $sum+=$item[1];
    $an++;
}
close IN;
$ps_avg=$sum/$an;
#predicted mutation
#predicted mutation
my %mu=();
my $dis_avg=0;
my $neu_avg=0;
$an=0;
my $sum1=0;
my $sum2=0;
open(IN,"/data1/juns/phase/Disease/rhapsody_length.txt") or die "$!";
while (<IN>)
{
    chomp;
    next if /^uniprot/;
    my @item=split/\t/;
    $mu{$item[0]}=$item[1];
    $sum1+=$item[1];
    $an++;
}
close IN;
$dis_avg=$sum1/$an;

my %ptm=();
my $phos_avg=0;
my $meth_avg=0;
my $nitro_avg=0;
my $palm_avg=0;
$an=0;
$sum1=0;
$sum2=0;
my $sum3=0;
my $sum4=0;
open(FA, "/data1/juns/phase/PTM/predict/PTM_length.txt") or die "$!";
while (<FA>)
{
    chomp;
    next if /^Uniprot/;
    my ($uni,$type,$ace,$meth,$phos,$sumo,$nitro,$palm)=split/\s+/;
    next if exists $ptm{$uni};
    my @vec=($phos,$meth,$nitro,$palm);
    $ptm{$uni}=\@vec;
    $sum1+=$phos;
    $sum2+=$meth;
    $sum3+=$nitro;
    $sum4+=$palm;
    $an++;
}
$phos_avg=$sum1/$an;
$meth_avg=$sum2/$an;
$nitro_avg=$sum3/$an;
$palm_avg=$sum4/$an;

$an=0;
$sum=0;
my $le_avg=0;
my %le=();
open(IN,"/data1/juns/phase/proteinLength/proteinLength.txt") or die "$!";
while (<IN>)
{
    chomp;
    next if /^uniprot/;
    my @item=split/\t/;
    $le{$item[0]}=$item[1];
    $sum+=$item[1];
    $an++;
}
close IN;
$le_avg=$sum/$an;

open(FA, "$abs_dir"."Mut/NonPS/"."SPINE-D-local.txt") or die "$!";
while (<FA>)
{
    chomp;
    next if /^uniprot/;
    my ($uni,$dis,$l,$count,$pert)=split/\s+/;
    next if exists $c{$uni};
    $le{$uni}=$l;
}
close FA;

open(FA, "$abs_dir"."Mut/PhaSePred/"."SPINE-D-local.txt") or die "$!";
while (<FA>)
{
	chomp;
	next if /^uniprot/;
	my ($uni,$dis,$l,$count,$pert)=split/\s+/;
	next if exists $c{$uni};
	$le{$uni}=$l;
}
close FA;

open(FA, "$abs_dir"."Mut/BIB/"."SPINE-D-local.txt") or die "$!";
while (<FA>)
{
	chomp;
	next if /^uniprot/;
	my ($uni,$dis,$l,$count,$pert)=split/\s+/;
	next if exists $c{$uni};
	$le{$uni}=$l;
}
close FA;

my %wv=();
my $win="70";
my $size="60";
open(IN, "/data1/juns/word2v/uniprot_sprot"."$win"."_size$size.txt") or die "$!";
my $c=0;
while (<IN>)
{
    chomp;
    next if /^Uniprot/;
    my @item=split/\t/;
    my @ref=@item[1..$#item];
    $wv{$item[0]}=\@ref;
}
close IN;

open(OUT, ">$o_dir"."Intestinput.txt") or die "$!";
open(FA, "$f_dir"."Intest.txt") or die "$!";
while (<FA>)
{
    chomp;
    my ($uni,$tag)=split/\s+/,$_;
    my @vec=();
    my $bacfasta="$fasta_dir"."$uni.fasta";
    if (!-s $bacfasta) {
        print $bacfasta,"\t","no bac protein\n";
        for(my $j=1;$j<=20;$j++) #39
        {
            push(@vec,0);
        }
    }else{
        open(F,$bacfasta) or die "$!";
        chomp(my @bacteria_seq=<F>);
        close F;
        my $vector1=&sequence_feature1($bacteria_seq[1]);
        $vector1=~s/^\s+//;
        my @vec1=split/\t/,$vector1;
        my $vector2=&sequence_feature2($bacteria_seq[1]);
        $vector2=~s/^\s+//;
        my @vec2=split/\t/,$vector2;
        push(@vec,@vec1,@vec2);
    }
    
    my @vec1=();
    my $flag=0;
    my $bacpssm="$hhm_dir"."$uni.hhm";
    if (!-s $bacpssm) {
        print $bacpssm,"\tno bac pssm\n";
        for(my $j=1;$j<=4;$j++)
        {
            push(@vec1,0);
        }
    }else{
        @vec1=@{&hhm_feature($bacpssm)};
    }
    my @vec2=();
    $bacpssm="$pssm_dir"."$uni.pssm";
    if (!-s $bacpssm) {
        print $bacpssm,"\tno bac pssm\n";
        for(my $j=1;$j<=4;$j++)
        {
            push(@vec2,0);
        }
    }else{
        @vec2=@{&pssm_feature($bacpssm)};
    }
    my @vec3=();
    if (exists $cs{$uni})
    {
        @vec3=@{$cs{$uni}};
    }else{
        for(my $j=1;$j<=5;$j++)
        {
            push(@vec3,0);
        }
    }
    my @fea=();
    push(@fea,@vec,@vec1,@vec2,@vec3);
    
    if (exists $hash{$uni})
    {
       push (@fea,@{$hash{$uni}})
    }else{
       push (@fea,qw(0 0 0 0))
    }
    
    if (exists $c{$uni})
    {
        push(@fea,$c{$uni});
    }else{
        push(@fea,$disorder_avg);
    }
    
     if (exists $rna{$uni})
    {
        push(@fea,$rna{$uni});
    }else{
        push(@fea,$rna_avg);
    }
    
    if (exists $dna{$uni})
    {
        push(@fea,$dna{$uni});
    }else{
        push(@fea,$dna_avg);
    }
    if (exists $ps{$uni})
    {
        push(@fea,$ps{$uni});
    }else{
        push(@fea,$ps_avg);
    }
    
    my $vec1="";
    if (exists $mu{$uni})
    {
        $vec1=$mu{$uni};
    }else{
        $vec1=$dis_avg;
    }
    @vec1=();
    if (exists $ptm{$uni})
    {
        @vec1=@{$ptm{$uni}};
    }else{
        @vec1=($phos_avg,$meth_avg,$nitro_avg,$palm_avg);
    }
    push(@fea,$vec1,@vec1);
    
    if (exists $le{$uni})
    {
        push(@fea,$le{$uni});
    }else{
        push(@fea,$le_avg);
    }
    @vec=();
    $bacfasta="$fasta_dir"."$uni.fasta";
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

#feature
sub sequence_feature2
{
    my %freq=();
    for(my $i=1;$i<=4;$i++)
    {
        for(my $j=1;$j<=4;$j++)
        {
            my $ref=$i.$j;
            $freq{$ref}=0;
        }
    }
    my $sequnce=$_[0];
    my @elements=split//,$sequnce;
    my $new_seq="";
    for(my $i=0;$i<@elements;$i++)
    {
        if (exists $class{$elements[$i]})
        {
            $new_seq=$new_seq.$class{$elements[$i]};
        }
    }
    
    
    for(my $i=0;$i<(length $new_seq)-1;$i++)
    {
        my $triple=substr($new_seq,$i,2);
        $freq{$triple}++;
    }
    
    #my @vector=();
    my $max=0;
    my $min=0;
    for(my $i=1;$i<=4;$i++)
    {
        for(my $j=1;$j<=4;$j++)
        {
            my $ref=$i.$j;
            my $frequncy=$freq{$ref};
            if ($frequncy>$max)
            {
                $max=$frequncy;
            }
            if ($frequncy<$min)
            {
                $min=$frequncy;
            }
        }
    }
    #print $new_seq,"\t$min\t$max\n";
    my $result="";
    for(my $i=1;$i<=4;$i++)
    {
        for(my $j=1;$j<=4;$j++)
        {
            my $ref=$i.$j;
            my $frequncy=sprintf "%0.3f",($freq{$ref}-$min)/$max;
            $result=$result."\t".$frequncy;
        }
    }
    return $result;
}

sub sequence_feature1
{
    my %freq=();
    for(my $i=1;$i<=4;$i++)
    {
        $freq{$i}=0;
    }
    my $sequnce=$_[0];
    my @elements=split//,$sequnce;
    my $new_seq="";
    for(my $i=0;$i<@elements;$i++)
    {
        if (exists $class{$elements[$i]})
        {
            $new_seq=$new_seq.$class{$elements[$i]};
        }
    }
    for(my $i=0;$i<(length $new_seq);$i++)
    {
        my $triple=substr($new_seq,$i,1);
        $freq{$triple}++;
    }
    my $len=length $sequnce;
    my $result="";
    for(my $i=1;$i<=4;$i++)
    {
        my $ref=$i;
        my $frequncy=sprintf "%0.3f",$freq{$ref}/$len;
        $result=$result."\t".$frequncy;
    }
    return $result;
}

sub hhm_feature
{
    my $file=$_[0];
    #print $file,"\n";
    my @twenty=();
    my %mtx=();
    for(my $i=0;$i<20;$i++)
    {
        $twenty[$i]=0;
    }
    my $l=0;
    open(INS,"$file") or die "$!";
    chomp(my @hhm=<INS>);
	close INS;
    for(my $i=0;$i<@hhm;$i++)
    {
        last if $hhm[$i]=~/^\/\//;
        next if $i<=10;
        next if $hhm[$i]=~/^>/;
        my @item=split/\s+/,$hhm[$i];
        if ($item[0] && $item[1])
        {
            next if $item[0] eq "NULL";
            if ($item[0]=~/[A-Z]/ && $item[1]=~/[0-9]/)
            {
                my @hhm_line=@item[2..21];
                for(my $k=0;$k<=$#hhm_line;$k++)
                {
                    if($hhm_line[$k] ne "*")
                    {
                        $hhm_line[$k]=2**(-($hhm_line[$k]/1000));
                    }else
                    {
                        $hhm_line[$k]=0;
                    }
                }
                my @n_val=@{&zscore(@hhm_line)};
                $l++;
                for(my $j=0;$j<@n_val;$j++)
                {
                    $twenty[$j]+=$n_val[$j];
                }
                $i=$i+2;
            }
        }
    }
    close INS;
    
    for(my $j=0;$j<@twenty;$j++)
    {
        $twenty[$j]=sprintf "%0.3f",$twenty[$j]/$l;
    }
    my $h=$twenty[1]+$twenty[18]+$twenty[6]+$twenty[10]+$twenty[19]+$twenty[4]+$twenty[7]+$twenty[17]+$twenty[0]+$twenty[9];
    my $c=$twenty[3]+$twenty[14]+$twenty[8]+$twenty[2];
    my $p=$twenty[15]+$twenty[13]+$twenty[3]+$twenty[16]+$twenty[11];
    my $g=$twenty[5]+$twenty[12];
    my @res=($h,$c,$p,$g);
    return \@res;
}

sub pssm_feature
{
    my $file=$_[0];
    #print $file,"\n";
    my @twenty=();
    my %mtx=();
    for(my $i=0;$i<20;$i++)
    {
        $twenty[$i]=0;
    }
    my $l=0;
    open(INS,"$file") or die "$!";
    while (<INS>)
    {
        #chomp;
        my $raw=$_;
        $raw=~s/^\s+//g;
        my @raw_inf=split/ +/,$raw;
        my $num=scalar(@raw_inf);
        #print $num,"\n";
        if($num>42)
        {
            my @val=@raw_inf[2..21];
            my @n_val=@{&zscore(@val)};
            $l++;
            for(my $j=0;$j<@n_val;$j++)
            {
                $twenty[$j]+=$n_val[$j];
            }
        }
    }
    close INS;
    
    for(my $j=0;$j<@twenty;$j++)
    {
        $twenty[$j]=sprintf "%0.3f",$twenty[$j]/$l;
    }
    my $h=$twenty[4]+$twenty[17]+$twenty[8]+$twenty[12]+$twenty[18]+$twenty[13]+$twenty[9]+$twenty[19]+$twenty[0]+$twenty[10];
    my $c=$twenty[6]+$twenty[1]+$twenty[11]+$twenty[3];
    my $p=$twenty[15]+$twenty[5]+$twenty[16]+$twenty[2];
    my $g=$twenty[7]+$twenty[14];
    my @res=($h,$c,$p,$g);
    return \@res;
}

sub zscore
{
    my @tempsub=@_;
    my $total=0;
    foreach my $item (@tempsub)
    {
        $total+=$item;
    }
    my $avgsub=$total/scalar @tempsub;
    
    my $sdsub=0;
    foreach my $item (@tempsub){
        $sdsub+=($item-$avgsub)**2;
    }
    my $stdsub=sqrt($sdsub/($#tempsub));
   # print $avgsub,"\t$stdsub\n";
    my @norm=();
    foreach my $item (@tempsub)
    {
        my $norm=0;
        if ($stdsub ==0)
        {
            $norm=0;
        }else{
            $norm=($item-$avgsub)/$stdsub;
            #$norm=sprintf"%.4f",1/(1+exp(-$norm));
        }
        push(@norm,$norm);
    }
    return \@norm;
}



