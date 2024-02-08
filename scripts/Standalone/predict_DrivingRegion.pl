#!/usr/bin/env perl
use strict;
use 5.010;
use warnings;

use Getopt::Long;
use File::Basename;
use IPC::Cmd qw[can_run run];
use File::Temp qw/ tempfile tempdir /;
use FindBin qw($Bin);

use Cwd;

#usage information
my $usage = <<USAGE;
Usage: predict_DrivingRegion.pl -i in.fa  -o outfile
 Options:
  -i    input fasta file, default STDIN
  -o    output fasta file, default STDOUT
USAGE

my $in_fasta = '';
my $outfile  = '';
die $usage
  unless GetOptions(
    "i:s" => \$in_fasta,
    "o:s" => \$outfile,
  );

#Enviment checking
can_run('python')   or die 'python is not installed!';
can_run('sort')     or die 'sort is not installed!';
can_run('bedtools') or die 'bedtools is not installed!';
my $progam_dir  = $Bin . "../scripts/";
my $working_dir = getcwd;

my $wordvec_file = $Bin . "/../../datasets/wordvec/uniprot_sprot70_size60.txt";
die "Cannot dectect wordvec_file:$wordvec_file" unless -e $wordvec_file;

my $train_wvc = $Bin . "/../../datasets/train/";
die "Cannot dectect training data set:$train_wvc" unless -d $train_wvc;

#############################################
#input and output file initiation
#############################################

my $in_fasta_fh;
if ( $in_fasta && $in_fasta ne '-' ) {
    die "input fastafile does not exists\n" unless -e $in_fasta;
    open $in_fasta_fh, "<", $in_fasta or die "cannot open file $in_fasta:$!\n";
}
else {
    if (@ARGV) {
        if ( $ARGV[0] eq '-' ) {
            $in_fasta_fh = *STDIN;
        }
        elsif ( -e $ARGV[0] ) {
            open $in_fasta_fh, "<", $ARGV[0]
              or die "cannot open file $ARGV[0]:$!\n";
        }
        else {
            die "$ARGV[0] does not exists\n";
        }
    }
    else {
        $in_fasta_fh = *STDIN;
    }
}

my $outfile_fh;
if ( $outfile && $outfile ne '-' ) {
    open $outfile_fh, ">", $outfile or die "cannot open file $outfile:$!\n";
}
else {
    $outfile_fh = *STDOUT;
}

####################################################
#Read fasta files and store to hash %seqs
####################################################
# my %seqs = ();
my %seqs;
warn "Reading fasta file...\n";
while (<$in_fasta_fh>) {
    s/\r?\n//;
    my $id = "";
    if (/^(>.*?)\s/) {
        $id = $1;
        warn "\tReading $id...\n";
        if ( exists $seqs{$id} ) {
            warn "duplicated squence id:$id";
        }
        $seqs{$id} = '';
    }
    elsif (/^\s*$/) {

        #preprocess sequence
        s/\n|\r|\s//g;
        s/-|\*//g;
        $seqs{$id} .= uc $_;
    }
}

####################################################
#Interate Each Sequence
####################################################
foreach my $id ( keys %seqs ) {
    warn "Processing sequencing $id\n";

    my $seq = $seqs{$id};
    my $len = length $seq;

    ####Split fasta into different truncations

    my $fasta_dir = "tmp.$id.fasta/";
    mkdir $fasta_dir;

    my $frag = "20";
    my @list = ();
    for ( my $i = 0 ; $i < $len - $frag + 1 ; $i++ ) {
        my $new  = substr( $seq, 0, $i ) . substr( $seq, $i + $frag );
        my $temp = "query" . "_$i";
        push( @list, $temp );
        open my $fasta_fh, ">$fasta_dir" . "$temp.fasta" or die "$!";
        print $fasta_fh ">query" . "_$i\n$new\n";
        close $fasta_fh;
    }

    my $ML_dir = $working_dir . "ML2/";
    mkdir $ML_dir;

    open my $intest_fh, ">$ML_dir" . "Intest.txt" or die "$!";
    foreach my $human (@list) {
        print $intest_fh "$human\t0\n";
    }
    close $intest_fh;

    ###word2vec hash
    my %wv = ();
    open my $wv_fh, "<", $wordvec_file or die "$!";
    while (<$wv_fh>) {
        chomp;
        next if /^Uniprot/;
        my @item = split /\t/;
        my @ref  = @item[ 1 .. $#item ];
        $wv{ $item[0] } = \@ref;
    }
    close $wv_fh;

    my $o_dir = $ML_dir . "word2vec/";
    mkdir $o_dir;

    my $size = 60;
    open my $Intestinput_fh, ">$o_dir" . "Intestinput.txt" or die "$!";
    open my $fa_fh, "<", "$ML_dir" . "Intest.txt" or die "$!";
    while (<$fa_fh>) {
        chomp;
        my ( $uni, $tag ) = split /\s+/;
        my @fea      = ();
        my @vec      = ();
        my $bacfasta = "$fasta_dir" . "$uni.fasta";

        if ( !-s $bacfasta ) {
            print $bacfasta, "\t", "no bac protein\n";
            for ( my $j = 1 ; $j <= $size ; $j++ )    #39
            {
                push( @vec, 0 );
            }
        }
        else {
            open my $tmp_fh, $bacfasta or die "$!";
            chomp( my @bacteria_seq = <$tmp_fh> );
            close $tmp_fh;
            my @aa = split //, $bacteria_seq[1];
            for ( my $i = 0 ; $i < ( scalar @aa ) - 2 ; $i++ ) {
                my $tri = $aa[$i] . $aa[ $i + 1 ] . $aa[ $i + 2 ];
                if ( exists $wv{$tri} ) {
                    my @tvec = @{ $wv{$tri} };
                    for ( my $j = 0 ; $j < @tvec ; $j++ ) {
                        $vec[$j] += $tvec[$j];
                    }
                }

            }
        }
        push( @fea, @vec );
        print $Intestinput_fh "$tag\t", join "\t", @fea, "\n";
    }
    close $Intestinput_fh;
    close $fa_fh;

    ###Model Test
    for ( my $i = 1 ; $i <= 100 ; $i++ ) {
        my $f1_dir = $train_wvc . "$i/word2vec70_60/";
        my $f2_dir = $o_dir;
        system "python $Bin/Intest-apply-test.py $f1_dir $f2_dir $i";
    }

    #Avg Prob
    open my $ml_out_fh, "<", "$ML_dir" . "Intest.txt" or die "$!";
    my $n   = 0;
    my %tag = ();
    while (<$ml_out_fh>) {
        chomp;
        $n++;
        my ( $uni, $tag ) = split;
        $tag{$n} = $uni;
    }

    close $ml_out_fh;

    open my $inprob_avg_fh, ">", "$o_dir" . "InProbAVG.txt" or die "$!";
    my %hash = ();
    for ( my $i = 1 ; $i <= 100 ; $i++ ) {
        open my $tmp_fh, "$o_dir" . "InProbphase$i.txt" or die "$!";
        my $k = 0;
        while (<IN>) {
            chomp;
            $k++;
            if ( exists $hash{$k} ) {
                $hash{$k} += $_;
            }
            else {
                $hash{$k} = $_;
            }
        }
        close $tmp_fh;
    }

    foreach my $key ( sort { $hash{$b} <=> $hash{$a} } keys %hash ) {
        my $prob = sprintf "%0.3f", $hash{$key} / 100;
        print $inprob_avg_fh $key, "\t$prob\n";
    }
    close $inprob_avg_fh;

    ####detect driving region

    open $inprob_avg_fh, "<", "$o_dir" . "InProbAVG.txt" or die "$!";
    my @prob = ();

    %hash = ();
    while (<IN>) {
        chomp;
        my ( $index, $score ) = split;
        push( @prob, $score );
        $hash{$index} = $score;
    }
    close $inprob_avg_fh;

    my $length = scalar @prob;

#my $bin=int($length/20);
#将蛋白质划分成bin，3，5，10，相较于平均值的偏离程度，(x-avg)/avg sort
#考虑如何合并
    my $avg = &Avg(@prob);

    #print $length,"\t",$avg,"\n";
    my %rank = ();
    foreach my $key ( sort { $a <=> $b } keys %hash ) {
        if ( $key <= $length ) {

            #my $prob=&lag($key);
            my $prob = $hash{$key};

            # my $bias=($avg-$prob)/$avg;
            my $bias = $avg - $prob;
            $rank{$key} = $bias;
        }
    }

    open my $bed_tmp, ">", "$o_dir" . "temp.bed" or die "$!";
    $n = 0;
    foreach my $key ( sort { $rank{$b} <=> $rank{$a} } keys %rank ) {

        #last if $rank{$key}<0;
        $n++;
        if ( $length > 2000 ) {
            last if $n >= int( $length * 0.01 );    # >20
        }
        elsif ( $length <= 2000 && $length > 1000 ) {
            last if $n >= int( $length * 0.02 );    ## 20-40
        }
        elsif ( $length > 500 && $length <= 1000 ) {
            last if $n >= int( $length * 0.04 );    ##20-40
        }
        else {
            last if $n >= int( $length * 0.05 );    ###25
        }

        print $bed_tmp "query\t", $key, "\t", $key + 20, "\t", $rank{$key},
          "\n";
    }
    close $bed_tmp;

    my $par1 =
      "sort -k1,1 -k2,2n " . $o_dir . "temp.bed >" . $o_dir . "temp.sort.bed";
    my $par2 =
        "bedtools merge -i "
      . $o_dir
      . "temp.sort.bed >"
      . $o_dir
      . "temp.drRegion";
    system("$par1");
    system("$par2");

    %rank = ();
    open my $drRegion_fh, "$o_dir" . "temp.drRegion" or die "$!";
    while (<$drRegion_fh>) {
        chomp;
        my ( $uni, $s, $e ) = split /\t/;
        my $sum = 0;
        my $m   = 0;
        for ( my $i = $s ; $i <= $e ; $i++ ) {
            if ( exists $hash{$i} ) {
                $sum += $hash{$i};
                $m++;
            }
            else {
                print $_, "\t$i\n";
            }
        }
        my $avg = $sum / $m;
        $rank{$_} = $avg;
    }
    close $drRegion_fh;

    open my $drRegion_rank_fh, ">$o_dir" . "drRegion.rank" or die "$!";

    #print $drRegion_rank_fh "#Total Length\t$length\n";
    my $limit  = 3;
    my %region = ();
    $n = 0;
    foreach my $key ( sort { $rank{$a} <=> $rank{$b} } keys %rank ) {
        $n++;
        my ( $query, $s, $e, $score ) = split /\t/, $key;
        print $drRegion_rank_fh $key, "\n";
        for ( my $i = $s ; $i <= $e ; $i++ ) {
            $region{$i}++;
        }
        last if $n >= $limit;
    }
    close $drRegion_rank_fh;

    print $outfile_fh "#Sequecing ID:$id\n";
    print $outfile_fh "#Residue in Purple denoted driving residues\n";
    print $outfile_fh "Pos\tAA\tProb\tDRegion\n";

    my @aa = split //, $seq;
    for ( my $i = 0 ; $i < @aa ; $i++ ) {
        my $index = $i + 1;
        next if $aa[$i] eq "-";
        if ( exists $region{$index} ) {
            if ( exists $hash{ $i - 9 } ) {
                print $outfile_fh $index, "\t$aa[$i]\t$hash{$i-9}\t1\n";
            }
            else {
                print $outfile_fh $index, "\t$aa[$i]\t-\t1\n";
            }
        }
        else {
            if ( exists $hash{ $i - 9 } ) {
                if ( $aa[$i] ) {
                    print $outfile_fh $index, "\t$aa[$i]\t$hash{$i-9}\t0\n";
                }
            }
            else {
                if ( $aa[$i] ) {
                    print $outfile_fh $index, "\t$aa[$i]\t-\t0\n";
                }
            }
        }
    }

    close $outfile_fh;

    system "rm -rf $fasta_dir";
    system "rm -rf $ML_dir";

}

sub Avg {
    my @prob = @_;
    my $sum  = 0;
    my $n    = 0;
    foreach my $prob (@prob) {
        $sum += $prob;
        $n++;
    }
    my $avg = $sum / $n;
    return $avg;
}
