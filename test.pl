#! /usr/bin/perl
#########################################################################
#        This Perl script is Copyright (c) 2002, Peter J Billam         #
#               c/o P J B Computing, www.pjb.com.au                     #
#                                                                       #
#     This script is free software; you can redistribute it and/or      #
#            modify it under the same terms as Perl itself.             #
#########################################################################

require './WalshTransform.pm'; import Math::WalshTransform;

my $taken  = 0;
my $failed = 0;
my $eps = .000000001;

for ($n=0; $n<=1023; $n++) { push @f, rand 9.9; }

my @H = &fht(@f); @h = &fhtinv(@H);
print STDERR "Hadamard transform and inverse: "; &unequal(\@f,\@h);

my @W = &fwt(@f); @w = &fwtinv(@W);
print STDERR "Walsh transform and inverse: "; &unequal(\@f,\@w);

my @HW = &hadamard2walsh(@H);
print STDERR "Hadamard to Walsh:   "; &unequal(\@W,\@HW);

my @WH = &walsh2hadamard(@W);
print STDERR "Walsh to Hadamard:   "; &unequal(\@H,\@WH);

my @f1; for ($n=0; $n<=511; $n++) { push @f1, rand 9.9; }
my @f2; for ($n=0; $n<=511; $n++) { push @f2, rand 9.9; }
my @lc1 = &Math::WalshTransform::old_logical_convolution(\@f1, \@f2);
my @lc2 = &logical_convolution(\@f1, \@f2);
print STDERR "Logical Convolution: "; &unequal(\@lc1,\@lc2);

if ($failed) {
	warn "failed $failed tests out of $taken\n"; exit 1;
} else {
	warn "all $taken tests succeeded\n"; exit 0;
}
# --------------------------- infrastructure ----------------
sub unequal { my ($xref, $yref) = @_;
	my @x = @$xref;
	my @y = @$yref;
	$taken++;
	if (scalar @x != scalar @y) {
		print STDERR "test $taken: unequal sized arrays ".scalar @x." and ",
		scalar @y," elements\n";
		$failed ++;
		return 0;
	}
	my $i; for ($i=$[; $i<=$#x; $i++) {
		if (abs($x[$i]-$y[$i]) > $eps) {
			print STDERR (sprintf "unequal array element %d, %g versus %g\n",
			$i,$x[$i],$y[$i]);
			$failed ++;
			return 0;
		}
	}
	print STDERR "ok\n"; return 1;
}

__END__

=pod

=head1 NAME

test.pl - Perl script to test Math::WalshTransform.pm

=head1 SYNOPSIS

 perl test.pl

=head1 DESCRIPTION

This script tests Math::WalshTransform.pm

=head1 AUTHOR

Peter J Billam  http://www.pjb.com.au/comp/contact.html

=head1 SEE ALSO

Math::WalshTransform.pm , http://www.pjb.com.au/ , perl(1).

=cut

