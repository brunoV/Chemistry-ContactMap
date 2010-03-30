package Chemistry::ContactMap::Types;

# ABSTRACT: Specific types for Chemistry::ContactMap.

use MouseX::Types::Mouse qw(ArrayRef Str FileHandle);
use MouseX::Types -declare =>
  [qw( MolPair ArrayRefofFileHandles ArrayRefofStrs )];
use List::MoreUtils 'all';
use Chemistry::File::PDB;

subtype ArrayRefofStrs,        as ArrayRef[Str];
subtype ArrayRefofFileHandles, as ArrayRef[FileHandle];
subtype MolPair, as ArrayRef, where {
    @{$_} == 2 and all { $_->isa('Chemistry::MacroMol') } @$_;
};

coerce MolPair,
    from ArrayRefofFileHandles, via { _fhs_to_molpair($_)  },
    from ArrayRefofStrs,        via { _strs_to_molpair($_) };

sub _str_to_mol {

    my ($mol) = Chemistry::MacroMol->read( @_ );

    return $mol;
}

sub _strs_to_molpair {
    my @mols = map { _str_to_mol($_) } @$_;

    return \@mols;
}

sub _fh_to_molpair {

    my $mol = Chemistry::File::PDB->read_mol($_[0]);

    return $mol;
}

sub _fhs_to_molpair {
    my @mols = map { _fh_to_molpair($_) } @$_;

    return \@mols;
}

1;

__END__

=head1 DESCRIPTION

This is a type library for L<Chemistry::ContactMap> related modules. You
don't have to use it directly.
