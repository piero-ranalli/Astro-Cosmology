# -*-perl-*-
#
# thanks to Brad Holden for providing an independent check
# of the Astro::Cosmology code. Since the algorithm is
# essentially taken from the same source it's not completely
# independent, but was written by someone else and
# uses a different integration scheme for the calculations.
#
# Note:
#   . only checks the distance calculations
#
# To do:
#   . check time and volume calculations
#

use strict;
$|++;

use Test;

plan tests => 4;

use PDL;
use PDL::Math;
use Astro::Cosmology;

# time for some actual comments
# This "module" is written so Doug has something to test against
# for his Astro::Cosmology package

# this routine determines whether or not you have flat universe or
# if you have to use the correction when going between the comoving
# distance and the interesting distances (lum, ang, etc.)

my $_tol = 1e-5;

sub _not_flat {

  my $matter = shift;
  my $lambda = shift;

  my $tot = $matter+$lambda;

  if (abs(1-$tot) < $_tol) {
    return(0);
  } elsif ($tot > 1) {
    return(-1); # negative curvature
  } else {
    return(1); # positive curvature
  }

}

=pod

=head1 B<d_c>

Inputs

=over 4

=item $z

Input redshift

=item $matter

Omega matter - defaults to 0.3

=item $lambda

Omega lambda - defaults to 0.7

=item $nsteps

Number of integration steps - defaults to 10000

=back

Outputs

=over 4

=item comoving distance in units of the horizon distance

=back

=cut

sub d_c {

  my $z = shift;
  my $matter = 0.3;
  my $lambda = 0.7;
  my $nsteps = 10000;

  $matter = shift if(@_);
  $lambda = shift if(@_);
  $nsteps = shift if(@_);

  my $curve = 1.0 - ($lambda + $matter);

  # here I do the actual integration, slowly, in PDL

  my $zs = xvals $nsteps;
  $zs *= $z/$nsteps;
  $zs += 0.5*$z/$nsteps;

  # below is the E(z) term from David Hogg's writeup

  my $ezs = $matter*(1+$zs)**3;
  $ezs += $curve*(1+$zs)**2;
  $ezs += $lambda;

  # below is integrand for the comoving distance

  $ezs = sqrt($ezs);
  $ezs = 1.0/$ezs;
  $ezs *= $z/$nsteps;

  return(intover($ezs));
}

=pod

=head1 B<d_m>

Inputs

=over 4

=item $z

Input redshift

=item $matter

Omega matter - defaults to 0.3

=item $lambda

Omega lambda - defaults to 0.7

=item $nsteps

Number of integration steps - defaults to 10000

=back

Outputs

=over 4

=item transverse comoving distance in units of the horizon distance

=back

=cut

sub d_m {

  my $z = shift;
  my $matter = 0.3;
  my $lambda = 0.7;
  my $nsteps = 10000;

  $matter = shift if(@_);
  $lambda = shift if(@_);
  $nsteps = shift if(@_);

  # here I have to "correct" the comoving distance depending on the
  # the curvature of the universe.  So, I calculate the curvature,
  # calculate the comoving distance along the line of sight and
  # only then do compute whether the curvature is positive, negatuve
  # or "zero", which means less than the tolerance variable $_tol
  #
  # _not_flat is the routine that checks the curvature, is uses the
  # input lambda and matter, not the computed curvature

  my $curve = 1.0- ($lambda+$matter);
  my $dc = d_c($z,$matter,$lambda,$nsteps);

  if (!_not_flat($matter,$lambda)) {
    return($dc);
  } elsif (_not_flat($matter,$lambda) < 0) {
    # negative curvature
    $curve = abs($curve);
    return(sin($dc*sqrt($curve))/sqrt($curve));
  } else {
    return(sinh($dc*sqrt($curve))/sqrt($curve));
  }

}

=pod

=head1 B<lum_z>

Inputs

=over 4

=item $z

Input redshift



=item $matter

Omega matter - defaults to 0.3



=item $lambda

Omega lambda - defaults to 0.7



=item $nsteps

Number of integration steps - defaults to 10000

=back

Outputs

=over 4

=item luminosity distance in units of the horizon distance

=back

=cut

# wrapper with correct 1+z

sub lum_z {

  my $z = shift;
  my $matter = 0.3;
  my $lambda = 0.7;
  my $nsteps = 10000;

  $matter = shift if(@_);
  $lambda = shift if(@_);
  $nsteps = shift if(@_);

#  $nsteps = floor($nsteps);

  my $dm = d_m($z,$matter,$lambda,$nsteps);
  return((1+$z)*$dm);

}

=pod

=head1 B<ang_z>

Inputs

=over 4

=item $z

Input redshift



=item $matter

Omega matter - defaults to 0.3



=item $lambda

Omega lambda - defaults to 0.7



=item $nsteps

Number of integration steps - defaults to 10000

=back

Outputs

=over 4

=item angular diameter distance in units of the horizon distance

=back

=cut

# wrapper with correct 1+z

sub ang_z {

  my $z = shift;
  my $matter = 0.3;
  my $lambda = 0.7;
  my $nsteps = 10000;

  $matter = shift if(@_);
  $lambda = shift if(@_);
  $nsteps = shift if(@_);

#  $nsteps = floor($nsteps);

  my $dm = d_m($z,$matter,$lambda,$nsteps);
  return($dm/(1+$z));

}


# Test Doug's code
# first the right universe according to the Nova PBS series
my $diff;
my $abstol = 1.0e-5; # (this is the default tolerance for the romberg integration in Astro::Cosmology)

my $sn = Astro::Cosmology->new({Matter=>0.3, Lambda=>0.7, H0=>0});
my $z = 1.235;
print $sn,"\n";
my $m_ld = $sn->lum_dist($z);
my $l_ld = lum_z($z,0.3,0.7,1000000);
print "Module: ",$m_ld,"\n";
print "Mine: ",$l_ld,"\n";
$diff = abs($l_ld-$m_ld);
print "Difference: ",sprintf("%.3g",$diff),"\n";
ok( $diff < $abstol );

# an open universe, coincidentally the best from my thesis....

my $opn = Astro::Cosmology->new({Matter=>0.32, Lambda=>0.0, H0=>0});
print $opn,"\n";
#$z = 0.9464;
$m_ld = $opn->lum_dist($z);
$l_ld = lum_z($z,0.32,0.0,1000000);
print "Module: ",$m_ld,"\n";
print "Mine: ",$l_ld,"\n";
$diff = abs($l_ld-$m_ld);
print "Difference: ",sprintf("%.3g",$diff),"\n";
ok( $diff < $abstol );

# a wacky closed universe

my $closed = Astro::Cosmology->new({Matter=>1.3, Lambda=>0.0, H0=>0});
print $closed,"\n";
$m_ld = $closed->lum_dist($z);
$l_ld = lum_z($z,1.3,0.0,1000000);
print "Module: ",$m_ld,"\n";
print "Mine: ",$l_ld,"\n";
$diff = abs($l_ld-$m_ld);
print "Difference: ",sprintf("%.3g",$diff),"\n";
ok( $diff < $abstol );

# a really wacky open universe

my $weird = Astro::Cosmology->new({Matter=>1.3, Lambda=>-1.0, H0=>0});
print $weird,"\n";
$m_ld = $weird->lum_dist($z);
$l_ld = lum_z($z,1.3,-1.0,1000000);
print "Module: ",$m_ld,"\n";
print "Mine: ",$l_ld,"\n";
$diff = abs($l_ld-$m_ld);
print "Difference: ",sprintf("%.3g",$diff),"\n";
ok( $diff < $abstol );

