package Chemistry::ContactMap;

# ABSTRACT: An interface for ContactMap classes

use Mouse::Role;
use Chemistry::ContactMap::Types 'MolPair';

requires '_build_contacts';

has structures => (
    is  => 'ro',
    isa => MolPair,
    required => 1,
    coerce   => 1
);

has contacts => ( is => 'ro', lazy_build => 1 );

1;

__END__

=head1 SYNOPSYS

    {
        package My::Class;
        use Mouse;

        with 'Chemistry::ContactMap';

        ...
    }

=head1 DESCRIPTION

L<Chemistry::ContactMap> is a Mouse role that implements a minimal
interface for ContactMap classes. It currently adds simple C<structures>
and C<contacts> attributes, and enforces the implementation of
C<contacts>' lazy builder.

=attr structures

Holds two read-only and required C<Chemistry::MacroMol> objects whose
contacts will be computed. It coerces from opened filehandles to PDB
files, and also from strings with the contents of said PDB files.

Thus, these will all work:

    My::Class->new( structures => \@fhs     );
    My::Class->new( structures => \@strings );
    My::Class->new( structures => \@mols    );

=contacts

Read only, lazily-built attribute without a type. You are encouraged to
add one in the consuming class. Also, you must implement a
C<_build_contacts> method that will populate this attribute when its
getter method is called.
