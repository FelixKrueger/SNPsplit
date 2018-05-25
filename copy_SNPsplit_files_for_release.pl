#!/usr/bin/env perl
use warnings;
use strict;
use File::Copy "cp";

my $dir = shift@ARGV;

die "Please provide a directory to copy files to!\n\n" unless ($dir);

unless (-d $dir){
  warn "Specified directory '$dir' doesn't exist. Creating it for you...\n\n";
  mkdir $dir or die "Failed to create directory: $!\n\n";
}

my @files = ('CHANGELOG.md','SNPsplit','tag2sort','SNPsplit_User_Guide.pdf','license.txt','SNPsplit_genome_preparation');

foreach my $file (@files){
  copy_and_warn($file);
}

sub copy_and_warn{
    my $file = shift;
    warn "Now copying '$file' to $dir\n";
    cp($file,"$dir/") or die "Copy failed: $!";
    
}

@files = ('SNPsplit','tag2sort','SNPsplit_genome_preparation');

foreach my $file (@files){
    set_permissions($file);
}

sub set_permissions{
    my $file = shift;
    warn "Setting permissions for ${dir}/$file\n";
    chmod 0755, "${dir}/$file";
}


### Taring up the folder
$dir =~ s/\/$//;
warn "Tar command:\ntar czvf ${dir}.tar.gz $dir\n\n";
system ("tar czvf ${dir}.tar.gz $dir/");
