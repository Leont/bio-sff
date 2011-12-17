package Bio::SFF::Index;

use Moo;

has manifest => (
	is => 'ro',
	required => 1,
);

has _offsets => (
	is => 'ro',
	isa => sub { ref($_[0]) eq 'HASH' },
	init_arg => 'offsets',
	required => 1,
);

sub offset_of {
	my ($self, $name) = @_;
	return $self->_offsets->{$name};
}

1;

#ABSTRACT: SFF index object

__END__

=begin Pod::Coverage

offset_of

=end Pod::Coverage

