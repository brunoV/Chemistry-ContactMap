use Test::More;
use Test::Exception;
use Chemistry::MacroMol;

{
    package My::ContactMap;
    use Mouse;

    with 'Chemistry::ContactMap';

    sub _build_contacts {
        return "I have contacts!";
    }


}

my @mols = map { Chemistry::MacroMol->new } (0..1);

my $cmap = My::ContactMap->new( structures => \@mols );
isa_ok $cmap, 'My::ContactMap';

dies_ok { My::ContactMap->new } 'Dies with no structures';

dies_ok { My::ContactMap->new( structures => [ 'foo', 'bar' ] ) }
'Dies with bad arguments';

dies_ok { My::ContactMap->new( structures => [ $mols[0] ] ) }
'Dies with bad arguments';

is $cmap->contacts, "I have contacts!";

done_testing();
