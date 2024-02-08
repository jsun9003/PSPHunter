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
Usage: predict_MutationEffect.pl -i in.fa  -o outfile
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
can_run('python') or die 'python is not installed!';
my $progam_dir  = $Bin . "../scripts/";
my $working_dir = getcwd;

my $wordvec_file = $Bin . "/../../datasets/wordvec/uniprot_sprot70_size60.txt";
die "Cannot dectect wordvec_file:$wordvec_file" unless -e $wordvec_file;

my $train_wvc = $Bin . "/../../datasets/train/";
die "Cannot dectect training data set:$train_wvc" unless -d $train_wvc;

my @aa = qw(C W H M Y F I V A L G P S Q T N E R K D);

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

    my @list  = ();
    my @seq   = split //, $seq;
    my %index = ();
    for ( my $i = 0 ; $i < @seq ; $i++ ) {
        if ( $seq[$i] ~~ @aa ) {
            $index{ $i + 1 } = $seq[$i];
            for ( my $j = 0 ; $j < @aa ; $j++ ) {
                my $raw = $seq;
                substr( $raw, $i, 1, $aa[$j] );
                my $temp = "query" . "_$i" . "_$aa[$j]";
                push( @list, $temp );
                open my $tmp_fa_out_fh, ">$fasta_dir" . "$temp.fasta"
                  or die "$!";
                print $tmp_fa_out_fh ">query" . "_$i" . "_$aa[$j]", "\n$raw\\n";
                close $tmp_fa_out_fh;
            }
        }
    }

    my $ML_dir = "tmp.$id.ML3/";
    mkdir $ML_dir;
    my $ml_intest_file = $ML_dir . "Intest.txt";

    open my $ml_out_fh, ">", $ml_intest_file or die "$!";
    foreach my $human (@list) {
        print $ml_out_fh "$human\t0\n";
    }
    close $ml_out_fh;

    ###word2vec hash
    my %wv = ();
    open my $wv_in_fh, "<", $wordvec_file
      or die "$!";
    while (<$wv_in_fh>) {
        chomp;
        next if /^Uniprot/;
        my @item = split /\t/;
        my @ref  = @item[ 1 .. $#item ];
        $wv{ $item[0] } = \@ref;
    }
    close $wv_in_fh;

    my $o_dir = $ML_dir . "word2vec/";
    mkdir $o_dir;

    my $size = 60;
    open my $intest_input_fh, ">", $o_dir . "Intestinput.txt"
      or die "$!";
    open my $ml_intest_fh, $ml_intest_file or die "$!";
    while (<$ml_intest_fh>) {
        chomp;
        my ( $uni, $tag ) = split /\s+/, $_;
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
            open my $bac_fh, "<", $bacfasta or die "$!";
            chomp( my @bacteria_seq = <$bac_fh> );
            close $bac_fh;
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
        print $intest_input_fh "$tag\t", join "\t", @fea, "\n";
    }
    close $intest_input_fh;
    close $ml_intest_fh;

    ###Model Test
    for ( my $i = 1 ; $i <= 100 ; $i++ ) {
        my $f1_dir = $train_wvc . "$i/word2vec70_60/";
        my $f2_dir = $o_dir;
        system "python $Bin/Intest-apply-test.py $f1_dir $f2_dir $i";
    }

    #Avg Prob

    open $ml_out_fh, "<", $ml_intest_file or die "$!";
    my $n   = 0;
    my %tag = ();
    while (<$ml_out_fh>) {
        chomp;
        $n++;
        my ( $uni, $tag ) = split;
        $tag{$n} = $uni;
    }

    close $ml_out_fh;

    open my $inprob_avg_fh, ">", "$working_dir/$id/InProbAVG.txt" or die "$!";
    my %hash = ();
    for ( my $i = 1 ; $i <= 100 ; $i++ ) {
        open my $in_prob_phase_fh, "<", "$o_dir" . "InProbphase$i.txt"
          or die "$!";
        my $k = 0;
        while (<$in_prob_phase_fh>) {
            chomp;
            $k++;
            if ( exists $hash{$k} ) {
                $hash{$k} += $_;
            }
            else {
                $hash{$k} = $_;
            }
        }
        close $in_prob_phase_fh;
    }

    my %prob = ();
    foreach my $key ( sort { $hash{$b} <=> $hash{$a} } keys %hash ) {
        my $prob = sprintf "%0.3f", $hash{$key} / 100;
        my ( $query, $index, $aa ) = split /_/, $tag{$key};
        $prob{ $index + 1 }{$aa} = $prob;
        print $outfile_fh $key, "\t$prob\n";
    }
    close $outfile_fh;

    print $inprob_avg_fh "SeqIndex\tAA\t";
    foreach my $aa ( sort @aa ) {
        print $inprob_avg_fh $aa, "\t";
    }
    print $inprob_avg_fh "\n";

    foreach my $key1 ( sort { $a <=> $b } keys %prob ) {
        print $inprob_avg_fh $key1, "\t$index{$key1}\t";
        foreach my $key2 ( sort keys %{ $prob{$key1} } ) {
            print $inprob_avg_fh $prob{$key1}{$key2}, "\t";
        }
        print $inprob_avg_fh "\n";
    }
    close $inprob_avg_fh;

    system "rm -rf $fasta_dir";
    system "rm -rf $ML_dir";

}
