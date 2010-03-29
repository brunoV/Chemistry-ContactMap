package Chemistry::ContactMap::Types;

use MouseX::Types::Mouse qw(ArrayRef);
use MouseX::Types -declare => [qw( MolPair )];
use List::MoreUtils 'all';

subtype MolPair, as ArrayRef, where {
    @{$_} == 2 and all { $_->isa('Chemistry::MacroMol') } @$_;
};

1;
