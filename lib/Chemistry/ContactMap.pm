package Chemistry::ContactMap;

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
