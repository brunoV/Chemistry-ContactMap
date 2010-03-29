package Chemistry::ContactMap::Distance;

use warnings;
use strict;
use PDL::NiceSlice;
use PDL::LiteF;
use Carp qw/croak/;
use base 'Chemistry::ContactMap';
use List::Util qw//;

=head1 NAME

Chemistry::ContactMap::Distance - Calculate contacts between macromolecules
using the closest distance method.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

    use Chemistry::MacroMol;
    use Chemistry::ContactMap::Distance;

    $cmap = Chemistry::ContactMap::Distance->new;
    $cmap->calculate($mol1, $mol2)

    # Are residues 4 of $mol1 and 200 of $mol2 making contact?
    print "Yes!" if ( $cmap(3, 199) == 1 );

    $cmap->print;
    $cmap->save('contact_map.cmp');

    $cmap2 = Chemistry::ContactMap::Distance->read('contact_map2.cmp');

    $cmap_average = ( $cmap + $cmap2 ) / 2;

This module implements the closest atom algorithm for determining
non-covalent contacts between macromolecules. See ContactMap for
a more general introduction.

Contacts are calculated using the closest atom method, which defines
determines them if the closest atom-atom distance between two residues is
below an arbitrary threshold. Typically, this distance cutoff is set at 6
Angstroms (for a comparative review, see Fischer, T. et. al: Assessing
methods for identifying pair-wise atomic contacts across binding interfaces,
J. Struc.  Biol., vol 153, p. 103-112, 2006), but it could be set at any
arbitrary value.

This package is meant to be used to analyze protein interfaces, whether
between chains of an oligomer, or between two polypeptide chains of a
protein-protein complex.

=head1 METHODS

=head2 $cmap->radius

Sets or gets the value for the 'radius' attribute. Default is 6 Angstroms.
The radius defines a sphere around each atom; every other atom of
the interacting partner with a center mass within this sphere is
considered to be in contact with the first atom. The default value
gives a compromise between finding all the actual contacts while
minimizing false positives. See the bibliographic reference in the
description for more details. If the value is set to '-1', the
distance between every residue-residue pair is returned when calling
the 'calculate' method.

=cut

sub radius {
   my ( $self, $radius ) = @_;
   if ($radius) {
      unless (
         $radius =~ /^[0-9]+  # Starts with digits
	(\.[0-9]+)*$          # Optionally, a decimal part with one point
	                      # And one or more digits behind.
	|                     # Alternatively,
	^-1$                  # The number -1.
	/x
          )
      {
         croak "Distance threshold is not a number\n";
      }

      $self->{radius} = $radius;
      return 1;
   } else {
      return $self->{radius};
   }
}

=head2 $cmap->calculate($mol1, $mol2)

Calculates the contacts (putative non-covalent interactions) between
$mol1 and $mol2, two Chemistry::MacroMol objects using the closest
atom approach.

=cut

sub calculate {

   # Calculate the actual contact values.
   # Args: two MacroMol objects and a radius.
   # If not defined, it tries to get the arguments
   # from the object's properties. When called with arguments,
   # they are stored in the object's proper attributes, overriding
   # them if they were set.

   my $self = shift;
   my (%args) = @_;

   my ( $mol1, $mol2 ) = ( $args{mol1}, $args{mol2} );
   my $radius = $args{radius};

   if ( $mol1 && $mol2 ) {
      $self->structures( $mol1, $mol2 );
   } elsif ( $self->structures ) {
      ( $mol1, $mol2 ) = $self->structures;
   } else {
      croak "Need two Chemistry::MacroMol objects!\n";
   }

   if ($radius) {
      $self->radius($radius);
   } elsif ( $self->radius ) {
      $radius = $self->radius;
   } else {
      croak "Radius attribute has not been set.\n
	 You should never have gotten here.\n";
   }

   if ( $radius == -1 ) {

      # Return distance matrix.
      $self->_set_PDL( _distance_matrix( $mol1, $mol2 ) );
   } else {

      # Return contact map.
      $self->_set_PDL( _distance_matrix( $mol1, $mol2 ) < $radius );

   }

   return 1;
}

=head2 $cmap->contacts

Returns a list of Chemistry::Bond elements that represent the calculated
contacts. This are automatically created after the "calculate" method
is called.

=cut

sub contacts {
   my $self = shift;
   if ( !$self->SUPER::contacts ) {
      $self->_set_bonds_from_pdl;
   }
   $self->SUPER::contacts;
}

# Private methods.

sub _distance_matrix {

   # Calculate a distance matrix between two Chemistry::MacroMol
   # objects. Distances are the minimum atom-atom distance for
   # each residue.

   # TODO split into smaller matrixes dynamically
   # to reduce memory consuption. Important to reduce
   # the chances of swapping.

   my ( $self, $other ) = @_;

   my $x = _get_coords($self);
   my $y = _get_coords($other);

   # Create dummy dimensions to allow
   # proper threading in the substraction
   my $xx = $x->dummy(1)->dummy(2);

   # Calculate all-versus-all square difference
   # between atom coordinates (memory hungry).
   # Operations are performed inplace to avoid
   # duplication of large datasets.
   my $D = ( $xx - $y )->inplace->power( 2, 0 );

   # Return the the lowest atom-atom distance for each
   # residue-residue pair.
   return sumover($D)->xchg( 3, 0 )->minimum->minimum->inplace->sqrt;
}

sub _get_coords {

   # Return a 3D piddle with coordinates.
   # Takes a Chemistry::MacroMol as an argument.
   # 1st D: xyz values
   # 2nd D: residues
   # 3rd D: atoms of each residue.

   my $self = shift;
   my @res  = $self->domains;

   # Create the coordinate matrix filled with
   # "bad values" (missing values).

   my $nres
       = List::Util::max( map { $_->attr('pdb/sequence_number') } @res );
   my $matrix = zeroes( 3, $nres + 1, 14 );
   $matrix->inplace->setvaltobad(0);

   #   for ( my $j = 1; $j < max(@res); $j++ ) {
   foreach my $res (@res) {
      my @atoms = $res->atoms;
      my $j     = $res->attr('pdb/sequence_number');
      for ( my $k = 0; $k < 14; $k++ ) {
         if ( $atoms[$k] ) {
            $matrix ( 0, $j, $k ) .= $atoms[$k]->x3;
            $matrix ( 1, $j, $k ) .= $atoms[$k]->y3;
            $matrix ( 2, $j, $k ) .= $atoms[$k]->z3;
         }
      }
   }
   return $matrix;
}

sub _set_bonds_from_pdl {

   # create Chemistry::Bond objects from contact information.
   my $self = shift;
   my ( $mol1, $mol2 ) = $self->structures;

   my ( @residues1, @residues2 );

   map { $residues1[ $_->attr('pdb/sequence_number') ] = $_ }
       $mol1->domains;
   map { $residues2[ $_->attr('pdb/sequence_number') ] = $_ }
       $mol2->domains;

   my $idx = whichND( $self > 0 );

   foreach my $c ( 0 .. $idx->dim(1) - 1 ) {

      my ( $distance, $atom1, $atom2 )
          = $residues1[ $idx->at( 0, $c ) ]
          ->distance( $residues2[ $idx->at( 1, $c ) ] );

      my $bond = Chemistry::Bond->new(
         type  => 'non-covalent',
         order => 0,
         atoms => [ $atom1, $atom2 ],
      );

      # Apply modifications to the user's molecule.
      $mol1->add_bond($bond);

      # Add the bond to the 'contact' attribute
      # of the ContactMap object
      $self->_add_contacts($bond);
   }

   # Appy modifications.
   $self->structures( $mol1, $mol2 );
}

=head1 AUTHOR

Bruno Vecchi, C<< <vecchi at cpan.org> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-contactmap at rt.cpan.org>,
or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=ContactMap>.  I will be
notified, and then you'll automatically be notified of progress on your bug as
I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Chemistry::ContactMap::Distance


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=ContactMap>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/ContactMap>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/ContactMap>

=item * Search CPAN

L<http://search.cpan.org/dist/ContactMap/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 COPYRIGHT & LICENSE

Copyright 2008 Bruno Vecchi, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.


=cut

1;    # End of Chemistry::ContactMap::Distance
