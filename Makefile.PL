use ExtUtils::MakeMaker;

WriteMakefile(
	NAME      => 'Math::WalshTransform',
	VERSION_FROM => 'WalshTransform.pm',
	EXE_FILES => [ glob "bin/*" ],
	AUTHOR    => 'PJB, Peter Billam, www.pjb.com.au/comp/contact.html',
	ABSTRACT  => 'Fast Hadamard and Walsh Transforms',
	dist      => {COMPRESS=>'gzip -9f', SUFFIX => 'gz'},
);
