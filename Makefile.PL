# -*- perl -*-

# try and be helpful if PDL is not available
eval "use PDL::Core::Dev;";
if ( $@ ) {
    die <<'EOD';

ERROR

  The Astro::Cosmology module requires the PDL package to be installed
  before it will build. PDL can be obtained from CPAN or from the PDL 
  web site at http://pdl.perl.org/

EOD
}


use ExtUtils::MakeMaker;
use Cwd;

$package = [ 'cosmo.pd', Cosmology, Astro::Cosmology ];

%hash = pdlpp_stdargs( $package );

$hash{VERSION_FROM} = 'cosmo.pd';

# what code do we want compiled and linked in?
@wanted = qw( romberg utils );
for $file ( @wanted ) {
    $hash{OBJECT} .= " $file\$(OBJ_EXT)";
}

push @{ $hash{LIBS} }, "-lm";

if ( $ExtUtils::MakeMaker::Version >= 5.4301 ) {
    $hash{ABSTRACT} = 'calculate cosmological distances, volumes and times';
    $hash{AUTHOR}   = 'Doug Burke (djburke@cpan.org)';
}

# given that this file uses part of the PDL distribution
# I'm not sure this is really any use
$hash{PREREQ_PM} = { PDL => 0 };

WriteMakefile( %hash );

sub MY::postamble {
    pdlpp_postamble($package);
}
