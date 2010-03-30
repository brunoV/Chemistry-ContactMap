package Chemistry::ContactMap::Types;

use MouseX::Types::Mouse qw(ArrayRef);
use MouseX::Types -declare => [qw( MolPair )];
# ABSTRACT: Specific types for Chemistry::ContactMap.
use List::MoreUtils 'all';

subtype MolPair, as ArrayRef, where {
    @{$_} == 2 and all { $_->isa('Chemistry::MacroMol') } @$_;
};

1;

__END__

=head1 DESCRIPTION

This is a type library for L<Chemistry::ContactMap> related modules. You
don't have to use it directly.
