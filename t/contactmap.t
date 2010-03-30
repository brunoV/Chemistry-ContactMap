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

my $strings = [ 't/data/A.pdb', 't/data/B.pdb'          ];
my $mols    = [ map { Chemistry::MacroMol->new } (0..1) ];

my $fhs;
for (0 .. $#$strings) {
    open( $fhs->[$_], '<', $strings->[$_] );
}

foreach my $input_type ($strings, $fhs, $mols) {
    $cmap = My::ContactMap->new( structures => $input_type );
    isa_ok $cmap, 'My::ContactMap';
}

dies_ok { My::ContactMap->new } 'Dies with no structures';

dies_ok { My::ContactMap->new( structures => [ 'foo', 'bar' ] ) }
'Dies with bad arguments';

dies_ok { My::ContactMap->new( structures => [ $mols[0] ] ) }
'Dies with bad arguments';

is $cmap->contacts, "I have contacts!";

done_testing();
